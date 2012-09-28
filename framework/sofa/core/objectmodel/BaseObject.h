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
#ifndef SOFA_CORE_OBJECTMODEL_BASEOBJECT_H
#define SOFA_CORE_OBJECTMODEL_BASEOBJECT_H

#include <sofa/core/objectmodel/Base.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/objectmodel/BaseObjectDescription.h>
#include <sofa/core/objectmodel/Link.h>
#ifdef SOFA_SMP
#include <sofa/defaulttype/SharedTypes.h>
#include <sofa/core/objectmodel/Context.h>
#include <sofa/core/objectmodel/BaseObjectTasks.h>
#include <sofa/helper/set.h>
#endif
#ifdef SOFA_SMP_NUMA
#include <numa.h>
#endif


namespace sofa
{
using helper::vector;

namespace core
{

// forward declaration of referenced classes
namespace topology
{
class Topology;
} // namespace topology

namespace visual
{
class VisualParams;
class DisplayFlags;
}


namespace objectmodel
{

class Event;
class BaseNode;

/**
 *  \brief Base class for simulation objects.
 *
 *  An object defines a part of the functionnality in the simulation
 *  (stores state data, specify topology, compute forces, etc).
 *  Each simulation object is related to a context, which gives access to all available external data.
 *  It is able to process events, if listening enabled (default is false).
 *
 */
class SOFA_CORE_API BaseObject : public virtual Base
#ifdef SOFA_SMP
    , public BaseObjectTasks
#endif
{
public:
    SOFA_CLASS(BaseObject, Base);
protected:
    BaseObject();

    virtual ~BaseObject();
public:

    /// @name control
    ///   Basic control
    /// @{

    /// Pre-construction check method called by ObjectFactory.
    template<class T>
    static bool canCreate(T* /*obj*/, BaseContext* /*context*/, BaseObjectDescription* /*arg*/)
    {
        return true;
    }

    /// Construction method called by ObjectFactory.
    template<class T>
    static typename T::SPtr create(T*, BaseContext* context, BaseObjectDescription* arg)
    {
        typename T::SPtr obj = sofa::core::objectmodel::New<T>();
        if (context) context->addObject(obj);
        if (arg) obj->parse(arg);
        return obj;
    }

    /// Parse the given description to assign values to this object's fields and potentially other parameters
    virtual void parse ( BaseObjectDescription* arg );

    /// Initialization method called at graph creation and modification, during top-down traversal.
    virtual void init();

    /// Initialization method called at graph creation and modification, during bottom-up traversal.
    virtual void bwdInit();

    /// Update method called when variables used in precomputation are modified.
    virtual void reinit();

    /// Save the initial state for later uses in reset()
    virtual void storeResetState();

    /// Reset to initial state
    virtual void reset();

    /// Called just before deleting this object
    /// Any object in the tree bellow this object that are to be removed will be removed only after this call,
    /// so any references this object holds should still be valid.
    virtual void cleanup();

    /// @}

    /// Render internal data of this object, for debugging purposes.
    virtual void draw(const core::visual::VisualParams*)
    {
#ifndef SOFA_DEPRECATE_OLD_API
        draw();
#endif
    }
    ///@}
#ifndef SOFA_DEPRECATE_OLD_API
    virtual void draw() {}
#endif

    /// @name Context accessors
    /// @{

    //void setContext(BaseContext* n);

    const BaseContext* getContext() const;

    BaseContext* getContext();

    const BaseObject* getMaster() const;

    BaseObject* getMaster();


    typedef MultiLink<BaseObject, BaseObject, BaseLink::FLAG_DOUBLELINK|BaseLink::FLAG_STRONGLINK> LinkSlaves;
    typedef LinkSlaves::Container VecSlaves;

    const VecSlaves& getSlaves() const;

    BaseObject* getSlave(const std::string& name) const;

    virtual void addSlave(BaseObject::SPtr s);

    virtual void removeSlave(BaseObject::SPtr s);

    virtual void copyAspect(int destAspect, int srcAspect);

    virtual void releaseAspect(int aspect);

    /// @}

    /// @name Component accessors
    /// @{

    /// Local search of an object of the given type
    template<class T>
    typename T::SPtr searchLocal() const { const BaseContext* context = getContext(); return typename T::SPtr(context->get<T>(BaseContext::Local)); }
    /// Upward search of an object of the given type, starting from the local context
    template<class T>
    typename T::SPtr searchUp() const { const BaseContext* context = getContext(); return typename T::SPtr(context->get<T>(BaseContext::SearchUp)); }
    /// Downward search of an object of the given type, starting from the local context
    template<class T>
    typename T::SPtr searchDown() const { const BaseContext* context = getContext(); return typename T::SPtr(context->get<T>(BaseContext::SearchDown)); }
    /// Search of an object of the given type, starting from the root
    /// @todo or only in the root ?
    template<class T>
    typename T::SPtr searchFromRoot() const { const BaseContext* context = getContext(); return typename T::SPtr(context->get<T>(BaseContext::SearchRoot)); }
    /// Search of an object of the given type, in the parents of the local context
    /// @todo is this an upward search starting from each parent, or only a local search in each parent ?
    template<class T>
    typename T::SPtr searchInParents() const { const BaseContext* context = getContext(); return typename T::SPtr(context->get<T>(BaseContext::SearchParents)); }

    /// Local search of all objects of the given type
    template<class T>
    vector<typename T::SPtr> searchAllLocal() const
    {
        vector<T*> v;
        const BaseContext* context = getContext();
        context->get<T>(&v,BaseContext::Local);
        vector<typename T::SPtr> vp;
        for( unsigned i=0; i<v.size(); i++ )
        {
            vp.push_back(typename T::SPtr(v[i]));
        }
        return vp;
    }
    /// Upward search of all objects of the given type, starting from the local context
    template<class T>
    vector<typename T::SPtr> searchAllUp() const
    {
        vector<T*> v;
        const BaseContext* context = getContext();
        context->get<T>(&v,BaseContext::SearchUp);
        vector<typename T::SPtr> vp;
        for( unsigned i=0; i<v.size(); i++ )
        {
            vp.push_back(typename T::SPtr(v[i]));
        }
        return vp;
    }
    /// Downward search of all objects of the given type, starting from the local context
    template<class T>
    vector<typename T::SPtr> searchAllDown() const
    {
        vector<T*> v;
        const BaseContext* context = getContext();
        context->get<T>(&v,BaseContext::SearchDown);
        vector<typename T::SPtr> vp;
        for( unsigned i=0; i<v.size(); i++ )
        {
            vp.push_back(typename T::SPtr(v[i]));
        }
        return vp;
    }
    /// Search of all objects of the given type, starting from the root
    template<class T>
    vector<typename T::SPtr> searchAllFromRoot() const
    {
        vector<T*> v;
        const BaseContext* context = getContext();
        context->get<T>(&v,BaseContext::SearchRoot);
        vector<typename T::SPtr> vp;
        for( unsigned i=0; i<v.size(); i++ )
        {
            vp.push_back(typename T::SPtr(v[i]));
        }
        return vp;
    }
    /// Search of all objects of the given type, in the parents of the local context
    /// @todo is this an upward search starting from each parent, or only a local search in each parent ?
    template<class T>
    vector<typename T::SPtr> searchAllInParents() const
    {
        vector<T*> v;
        const BaseContext* context = getContext();
        context->get<T>(&v,BaseContext::SearchParents);
        vector<typename T::SPtr> vp;
        for( unsigned i=0; i<v.size(); i++ )
        {
            vp.push_back(typename T::SPtr(v[i]));
        }
        return vp;
    }



    /// Local search of all objects of the given type with a given Tag
    template<class T>
    vector<typename T::SPtr> searchAllLocal(const Tag& t) const
    {
        vector<T*> v;
        const BaseContext* context = getContext();
        context->get<T>(&v,t,BaseContext::Local);
        vector<typename T::SPtr> vp;
        for( unsigned i=0; i<v.size(); i++ )
        {
            vp.push_back(typename T::SPtr(v[i]));
        }
        return vp;
    }
    /// Upward search of all objects of the given type with a given Tag, starting from the local context
    template<class T>
    vector<typename T::SPtr> searchAllUp(const Tag& t) const
    {
        vector<T*> v;
        const BaseContext* context = getContext();
        context->get<T>(&v,t,BaseContext::SearchUp);
        vector<typename T::SPtr> vp;
        for( unsigned i=0; i<v.size(); i++ )
        {
            vp.push_back(typename T::SPtr(v[i]));
        }
        return vp;
    }
    /// Downward search of all objects of the given typee with a given Tag, starting from the local context
    template<class T>
    vector<typename T::SPtr> searchAllDown(const Tag& t) const
    {
        vector<T*> v;
        const BaseContext* context = getContext();
        context->get<T>(&v,t,BaseContext::SearchDown);
        vector<typename T::SPtr> vp;
        for( unsigned i=0; i<v.size(); i++ )
        {
            vp.push_back(typename T::SPtr(v[i]));
        }
        return vp;
    }
    /// Search of all objects of the given typee with a given Tag, starting from the root
    template<class T>
    vector<typename T::SPtr> searchAllFromRoot(const Tag& t) const
    {
        vector<T*> v;
        const BaseContext* context = getContext();
        context->get<T>(&v,t,BaseContext::SearchRoot);
        vector<typename T::SPtr> vp;
        for( unsigned i=0; i<v.size(); i++ )
        {
            vp.push_back(typename T::SPtr(v[i]));
        }
        return vp;
    }
    /// Search of all objects of the given typee with a given Tag, in the parents of the local context
    /// @todo is this an upward search starting from each parent, or only a local search in each parent ?
    template<class T>
    vector<typename T::SPtr> searchAllInParents(const Tag& t) const
    {
        vector<T*> v;
        const BaseContext* context = getContext();
        context->get<T>(&v,t,BaseContext::SearchParents);
        vector<typename T::SPtr> vp;
        for( unsigned i=0; i<v.size(); i++ )
        {
            vp.push_back(typename T::SPtr(v[i]));
        }
        return vp;
    }




    /// Local search of all objects of the given type with a given TagSet
    template<class T>
    vector<typename T::SPtr> searchAllLocal(const TagSet& t) const
    {
        vector<T*> v;
        const BaseContext* context = getContext();
        context->get<T>(&v,t,BaseContext::Local);
        vector<typename T::SPtr> vp;
        for( unsigned i=0; i<v.size(); i++ )
        {
            vp.push_back(typename T::SPtr(v[i]));
        }
        return vp;
    }
    /// Upward search of all objects of the given type with a given TagSet, starting from the local context
    template<class T>
    vector<typename T::SPtr> searchAllUp(const TagSet& t) const
    {
        vector<T*> v;
        const BaseContext* context = getContext();
        context->get<T>(&v,t,BaseContext::SearchUp);
        vector<typename T::SPtr> vp;
        for( unsigned i=0; i<v.size(); i++ )
        {
            vp.push_back(typename T::SPtr(v[i]));
        }
        return vp;
    }
    /// Downward search of all objects of the given typee with a given TagSet, starting from the local context
    template<class T>
    vector<typename T::SPtr> searchAllDown(const TagSet& t) const
    {
        vector<T*> v;
        const BaseContext* context = getContext();
        context->get<T>(&v,t,BaseContext::SearchDown);
        vector<typename T::SPtr> vp;
        for( unsigned i=0; i<v.size(); i++ )
        {
            vp.push_back(typename T::SPtr(v[i]));
        }
        return vp;
    }
    /// Search of all objects of the given typee with a given TagSet, starting from the root
    template<class T>
    vector<typename T::SPtr> searchAllFromRoot(const TagSet& t) const
    {
        vector<T*> v;
        const BaseContext* context = getContext();
        context->get<T>(&v,t,BaseContext::SearchRoot);
        vector<typename T::SPtr> vp;
        for( unsigned i=0; i<v.size(); i++ )
        {
            vp.push_back(typename T::SPtr(v[i]));
        }
        return vp;
    }
    /// Search of all objects of the given typee with a given TagSet, in the parents of the local context
    /// @todo is this an upward search starting from each parent, or only a local search in each parent ?
    template<class T>
    vector<typename T::SPtr> searchAllInParents(const TagSet& t) const
    {
        vector<T*> v;
        const BaseContext* context = getContext();
        context->get<T>(&v,t,BaseContext::SearchParents);
        vector<typename T::SPtr> vp;
        for( unsigned i=0; i<v.size(); i++ )
        {
            vp.push_back(typename T::SPtr(v[i]));
        }
        return vp;
    }

    /// @}


    /// @name data access
    ///   Access to external data
    /// @{

    /// Current time
    double getTime() const;

    /// @}

    /// @name events
    ///   Methods related to Event processing
    /// @{

    Data<bool> f_listening;

    /// Handle an event
    virtual void handleEvent( Event* );

    /// Handle topological Changes
    /// @deprecated topological changes now rely on TopologyEngine
    virtual void handleTopologyChange() {}

    /// Handle topological Changes from a given Topology
    /// @deprecated topological changes now rely on TopologyEngine
    virtual void handleTopologyChange(core::topology::Topology* t);

    ///@}

    /// Bounding Box computation method.
    /// Default to empty method.
    virtual void computeBBox(const core::ExecParams* /* params */) {};

    /// Sets a source Object and parses it to collect dependent Data
    void setSrc(const std::string &v, std::vector< std::string > *attributeList=0);

    /// Sets a source Object and parses it to collect dependent Data
    /// Use it before scene graph insertion
    void setSrc(const std::string &v, const BaseObject *loader, std::vector< std::string > *attributeList=0);

    void* findLinkDestClass(const BaseClass* destType, const std::string& path, const BaseLink* link);

#ifdef SOFA_SMP
    void setPartition(Iterative::IterativePartition* p);
    Iterative::IterativePartition*  getPartition();
    Iterative::IterativePartition*  prepareTask();
#endif

protected:

    SingleLink<BaseObject, BaseContext, BaseLink::FLAG_DOUBLELINK> l_context;
    LinkSlaves l_slaves;
    SingleLink<BaseObject, BaseObject, BaseLink::FLAG_DOUBLELINK> l_master;

    // This method insures that context is never NULL (using BaseContext::getDefault() instead)
    // and that all slaves of an object share its context
    void changeContextLink(BaseContext* before, BaseContext*& after);

    /// This method insures that slaves objects have master and context links set correctly
    void changeSlavesLink(BaseObject::SPtr ptr, unsigned int /*index*/, bool add);

    // BaseNode can set the context of its own objects
    friend class BaseNode;

#ifdef SOFA_SMP
    Iterative::IterativePartition *partition_;
#endif
};

} // namespace objectmodel

} // namespace core

} // namespace sofa

#endif
