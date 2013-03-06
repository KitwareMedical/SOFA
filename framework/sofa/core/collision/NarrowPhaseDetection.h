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
#ifndef SOFA_COMPONENT_COLLISION_NARROWPHASEDETECTION_H
#define SOFA_COMPONENT_COLLISION_NARROWPHASEDETECTION_H

#include <sofa/core/collision/Detection.h>
#include <vector>
#include <map>
#include <algorithm>

namespace sofa
{

namespace core
{

namespace collision
{

/**
* @brief Given a set of potentially colliding pairs of models, compute set of contact points
*/

class NarrowPhaseDetection : virtual public Detection
{
public:
    SOFA_ABSTRACT_CLASS(NarrowPhaseDetection, Detection);

    typedef std::map< std::pair< core::CollisionModel*, core::CollisionModel* >, DetectionOutputVector* > DetectionOutputMap;
    typedef std::vector< std::pair< std::pair< core::CollisionModel*, core::CollisionModel* >, DetectionOutputVector** > > DetectionOutputVectors;
protected:
    /// Destructor
    virtual ~NarrowPhaseDetection() { }
public:
    /// Clear all the potentially colliding pairs detected in the previous simulation step
    virtual void beginNarrowPhase()
    {
        for (DetectionOutputVectors::iterator it = m_outputsVec.begin(); it != m_outputsVec.end(); it++)
        {
            DetectionOutputVector *do_vec = *(it->second);

            if (do_vec != 0)
                do_vec->clear();
        }
    }

    /// Add a new potentially colliding pairs of models
    virtual void addCollisionPair (const std::pair<core::CollisionModel*, core::CollisionModel*>& cmPair) = 0;

    /// Add a new list of potentially colliding pairs of models
    virtual void addCollisionPairs(const sofa::helper::vector< std::pair<core::CollisionModel*, core::CollisionModel*> >& v)
    {
        for (sofa::helper::vector< std::pair<core::CollisionModel*, core::CollisionModel*> >::const_iterator it = v.begin(); it!=v.end(); it++)
            addCollisionPair(*it);
    }

    virtual void endNarrowPhase()
    {
        DetectionOutputVectors::iterator it = m_outputsVec.begin();
        
        while (it != m_outputsVec.end())
        {
            DetectionOutputVector *do_vec = *(it->second);

            if (!do_vec || do_vec->empty())
            {
                /// @TODO Optimization
                DetectionOutputMap::iterator it_map = m_outputsMap.find(it->first);

                if (it_map != m_outputsMap.end())
                {
                        m_outputsMap.erase(it_map);
                }

                it = m_outputsVec.erase(it);

                if (do_vec) do_vec->release();
            }
            else
            {
                ++it;
            }
        }
    }

    //sofa::helper::vector<std::pair<core::CollisionElementIterator, core::CollisionElementIterator> >& getCollisionElementPairs() { return elemPairs; }

    const DetectionOutputMap& getDetectionOutputs()
    {
        return m_outputsMap;
    }

    const DetectionOutputVectors& getDetectionOutputsVector()
    {
        return m_outputsVec;
    }

    DetectionOutputVector*& getDetectionOutputs(CollisionModel *cm1, CollisionModel *cm2)
    {
        std::pair< CollisionModel*, CollisionModel* > cm_pair = std::make_pair(cm1, cm2);

        DetectionOutputMap::iterator it = m_outputsMap.find(cm_pair);

        if (it == m_outputsMap.end())
        {
            // new contact
            it = m_outputsMap.insert( std::make_pair(cm_pair, static_cast< DetectionOutputVector * >(0)) ).first;

            m_outputsVec.push_back( std::make_pair(cm_pair, &(it->second)) );
        }

        return it->second;
    }

protected:
    std::map<Instance, DetectionOutputMap> m_storedOutputsMap;

    DetectionOutputVectors m_outputsVec;
    std::map<Instance, DetectionOutputVectors> m_storedOutputsVec;

    virtual void changeInstanceNP(Instance inst)
    {
        m_storedOutputsMap[instance].swap(m_outputsMap);
        m_outputsMap.swap(m_storedOutputsMap[inst]);

        m_storedOutputsVec[instance].swap(m_outputsVec);
        m_outputsVec.swap(m_storedOutputsVec[inst]);
    }

private:
    DetectionOutputMap m_outputsMap;
};

} // namespace collision

} // namespace core

} // namespace sofa

#endif
