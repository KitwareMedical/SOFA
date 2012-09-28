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
#ifndef SOFA_CORE_COLLISION_PIPELINE_H
#define SOFA_CORE_COLLISION_PIPELINE_H

#include <sofa/core/CollisionModel.h>
#include <sofa/core/CollisionElement.h>
#include <sofa/core/collision/Intersection.h>
#include <sofa/core/collision/BroadPhaseDetection.h>
#include <sofa/core/collision/NarrowPhaseDetection.h>
#include <sofa/core/collision/DetectionOutput.h>
#include <sofa/core/collision/ContactManager.h>
#include <sofa/core/collision/CollisionGroupManager.h>

#include <sofa/helper/set.h>
#include <sofa/helper/vector.h>


namespace sofa
{

namespace core
{

namespace collision
{

/**
 * @brief Pipeline component gather list of collision models and control the sequence of computations
*/

class SOFA_CORE_API Pipeline : public virtual sofa::core::objectmodel::BaseObject
{
public:
    SOFA_ABSTRACT_CLASS(Pipeline, sofa::core::objectmodel::BaseObject);

protected:

    //sofa::helper::vector<DetectionOutput*> detectionOutputs;

    sofa::helper::vector<Intersection*> intersectionMethods;
    sofa::helper::vector<BroadPhaseDetection*> broadPhaseDetections;
    sofa::helper::vector<NarrowPhaseDetection*> narrowPhaseDetections;
    sofa::helper::vector<ContactManager*> contactManagers;
    sofa::helper::vector<CollisionGroupManager*> groupManagers;

    Intersection* intersectionMethod;
    BroadPhaseDetection* broadPhaseDetection;
    NarrowPhaseDetection* narrowPhaseDetection;
    ContactManager* contactManager;
    CollisionGroupManager* groupManager;

public:
    typedef NarrowPhaseDetection::DetectionOutputMap DetectionOutputMap;
protected:
    Pipeline();

    virtual ~Pipeline();
public:
    virtual void reset()=0;

    /// Remove collision response from last step
    virtual void computeCollisionReset()=0;
    /// Detect new collisions. Note that this step must not modify the simulation graph
    virtual void computeCollisionDetection()=0;
    /// Add collision response in the simulation graph
    virtual void computeCollisionResponse()=0;

    void computeCollisions()
    {
        computeCollisionReset();
        computeCollisionDetection();
        computeCollisionResponse();
    }

    //sofa::helper::vector<DetectionOutput*>& getDetectionOutputs() { return detectionOutputs; }

    /// Broad phase collision detection method accessor.
    const BroadPhaseDetection *getBroadPhaseDetection() const;

    /// Narrow phase collision detection method accessor.
    const NarrowPhaseDetection *getNarrowPhaseDetection() const;

    /// get the set of response available with the current collision pipeline
    virtual helper::set< std::string > getResponseList() const=0;
protected:
    /// Remove collision response from last step
    virtual void doCollisionReset() = 0;
    /// Detect new collisions. Note that this step must not modify the simulation graph
    virtual void doCollisionDetection(const sofa::helper::vector<core::CollisionModel*>& collisionModels) = 0;
    /// Add collision response in the simulation graph
    virtual void doCollisionResponse() = 0;
};

} // namespace collision

} // namespace core

} // namespace sofa

#endif
