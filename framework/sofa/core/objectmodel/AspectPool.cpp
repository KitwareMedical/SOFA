/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                              SOFA :: Framework                              *
*                                                                             *
* Authors: The SOFA Team (see Authors.txt)                                    *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/

#include "AspectPool.h"
#include <sofa/helper/system/thread/CircularQueue.inl>
#include <iostream>

namespace sofa
{

namespace core
{

namespace objectmodel
{

/**
 * Constructor.
 * It is private to only allow AspectPool to use it.
 */
Aspect::Aspect(AspectPool& pool, int id)
    : pool(pool), id(id), counter(0)
{
}

/**
 * Destructor
 */
Aspect::~Aspect()
{
}


void Aspect::add_ref()
{
    counter.inc();
}

void Aspect::release()
{
    if(counter.dec_and_test_null())
    {
        // release the aspect from the pool.
        pool.release(id);
    }
}

void intrusive_ptr_add_ref(Aspect* a)
{
    a->add_ref();
}
void intrusive_ptr_release(Aspect* a)
{
    a->release();
}


/**
 * Constructor: creates a new aspect pool.
 */
AspectPool::AspectPool()
{
    // Create all aspects and fill the list of free aspects.
    std::cout << "AspectPool"<<this<<": filling " << SOFA_DATA_MAX_ASPECTS << " aspects" << std::endl;
    aspects.resize(SOFA_DATA_MAX_ASPECTS);
    for(int i = 0; i < SOFA_DATA_MAX_ASPECTS; ++i)
    {
        aspects[i] = new Aspect(*this, i);
        AtomicInt aspectID(i);
        freeAspects.push(aspectID);
    }
    //std::cout << "AspectPool"<<this<<": " << freeAspects.size() << " aspects available" << std::endl;
}

/**
 * Destructor.
 */
AspectPool::~AspectPool()
{
    for(unsigned int i = 0; i < aspects.size(); ++i)
    {
        delete aspects[i];
    }
}

void AspectPool::setReleaseCallback(const boost::function<void (int)>& callback)
{
    releaseCallback = callback;
}

/**
 * Request a new aspect.
 * The returned object should stay alive as long as the aspect is in use.
 * It it possible to duplicate the AspectRef if several threads/algorithm use
 * the same aspect.
 * If no aspect remains available, null pointer is returned.
 */
AspectRef AspectPool::allocate()
{
    AspectRef ref;
    AtomicInt aspectID;
    std::cout << "AspectPool"<<this<<": allocate" << std::endl;
    std::cout << "AspectPool"<<this<<": " << freeAspects.size() << " aspects available" << std::endl;
    if(freeAspects.pop(aspectID))
    {
        std::cout << "AspectPool"<<this<<": aspect " << aspectID << " allocated" << std::endl;
        ref = aspects[aspectID];
    }
    else
        std::cout << "AspectPool"<<this<<": no aspect available" << std::endl;
    return ref;
}

/**
 * Release the aspect having the specified number.
 * It makes the number immediately available to satisfy subsequent AspectPool::allocate
 * requests.
 */
void AspectPool::release(int id)
{
    std::cout << "AspectPool"<<this<<": release aspect " << id << std::endl;
    if(releaseCallback != 0)
    {
        releaseCallback(id);
    }
    AtomicInt aspectID(id);
    freeAspects.push(aspectID);
}

AspectRef AspectPool::getAspect(int id)
{
    if ((unsigned)id >= aspects.size())
    {
        return AspectRef();
    }
    else
    {
        return aspects[id];
    }
}


AspectBuffer::AspectBuffer(AspectPool& pool)
    : pool(pool), latestID(-1), availableID(-1)
{
}

AspectBuffer::~AspectBuffer()
{
    // we can't call clear() here as the pool may have been destroyed before us
    // we rely on the code using this buffer to call clear at the right place.
    //clear();
}

void AspectBuffer::clear()
{
    int id = latestID.exchange(-1);
    if (id != -1)
        pool.getAspect(id)->release(); // remove the ref token that was implicitly held by latestID
    id = availableID.exchange(-1);
    if (id != -1)
        pool.getAspect(id)->release(); // remove the ref token that was implicitly held by availableID
}

AspectRef AspectBuffer::allocate()
{
    int id = availableID.exchange(-1);
    //std::cerr << "availableID: " << id << " -> " << -1 << std::endl;
    if (id == -1)
        return pool.allocate();
    else
    {
        AspectRef a = pool.getAspect(id);
        a->release(); // remove the ref token that was implicitly hold by availableID
        return a;
    }
}

/// Send a new version, overriding the latest if it was not already received (in which case it can be "recycled" using allocate)
void AspectBuffer::push(AspectRef id)
{
    int newID = -1;
    if (id)
    {
        newID = id->aspectID();
        id->add_ref(); // add a ref token that will be implicitly held by latestID
    }
    int prevID = latestID.exchange(newID);
    //std::cerr << "latestID: " << prevID << " -> " << newID << std::endl;
    if (prevID != -1)
    {
        // store prevID as the next available ID
        int freeID = availableID.exchange(prevID);
        //std::cerr << "availableID: " << freeID << " -> " << prevID << std::endl;
        if (freeID != -1)
            pool.getAspect(freeID)->release(); // remove the ref token that was implicitly held by availableID
    }
}

/// Receive the latest version, return true if one is available, or false otherwise (in which case id is unchanged)
bool AspectBuffer::pop(AspectRef& id)
{
    int newID = latestID.exchange(-1);
    //std::cerr << "latestID: " << newID << " -> " << -1 << std::endl;
    if (newID == -1) return false;
    if (id)
    {
        int prevID = id->aspectID();
        // store prevID as the next available ID
        id->add_ref(); // add a ref token that will be implicitly held by availableID
        int freeID = availableID.exchange(prevID);
        //std::cerr << "availableID: " << freeID << " -> " << prevID << std::endl;
        if (freeID != -1)
            pool.getAspect(freeID)->release(); // remove the ref token that was implicitly held by availableID
    }
    // replace id by newID
    id = pool.getAspect(newID);
    id->release(); // remove the ref token that was implicitly held by latestID
    return true;
}

} // namespace objectmodel

} // namespace core

} // namespace sofa
