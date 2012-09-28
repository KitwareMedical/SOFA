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
#include <sofa/core/collision/Pipeline.h>
//#include <sofa/component/collision/DiscreteIntersection.h>
//#include <sofa/simulation/tree/GNode.h>

namespace sofa
{

namespace core
{

namespace collision
{
//using namespace core::objectmodel;
//using namespace core::behavior;

Pipeline::Pipeline()
    : intersectionMethod(NULL),
      broadPhaseDetection(NULL),
      narrowPhaseDetection(NULL),
      contactManager(NULL),
      groupManager(NULL)
{
}

Pipeline::~Pipeline()
{
}


const BroadPhaseDetection *Pipeline::getBroadPhaseDetection() const
{
    return broadPhaseDetection;
}


const NarrowPhaseDetection *Pipeline::getNarrowPhaseDetection() const
{
    return narrowPhaseDetection;
}

#if 0

void Pipeline::init()
{
    simulation::tree::GNode* root = dynamic_cast<simulation::tree::GNode*>(getContext());
    if(root == NULL) return;
    intersectionMethods.clear();
    root->getTreeObjects<Intersection>(&intersectionMethods);
    intersectionMethod = (intersectionMethods.empty() ? NULL : intersectionMethods[0]);
    broadPhaseDetections.clear();
    root->getTreeObjects<BroadPhaseDetection>(&broadPhaseDetections);
    broadPhaseDetection = (broadPhaseDetections.empty() ? NULL : broadPhaseDetections[0]);
    narrowPhaseDetections.clear();
    root->getTreeObjects<NarrowPhaseDetection>(&narrowPhaseDetections);
    narrowPhaseDetection = (narrowPhaseDetections.empty() ? NULL : narrowPhaseDetections[0]);
    contactManagers.clear();
    root->getTreeObjects<ContactManager>(&contactManagers);
    contactManager = (contactManagers.empty() ? NULL : contactManagers[0]);
    groupManagers.clear();
    root->getTreeObjects<CollisionGroupManager>(&groupManagers);
    groupManager = (groupManagers.empty() ? NULL : groupManagers[0]);

    if (intersectionMethod==NULL)
        intersectionMethod = new sofa::component::collision::DiscreteIntersection;
}

void Pipeline::reset()
{
    computeCollisionReset();
}

void Pipeline::computeCollisionReset()
{
    simulation::tree::GNode* root = dynamic_cast<simulation::tree::GNode*>(getContext());
    if(root == NULL) return;
    if (broadPhaseDetection!=NULL && broadPhaseDetection->getIntersectionMethod()!=intersectionMethod)
        broadPhaseDetection->setIntersectionMethod(intersectionMethod);
    if (narrowPhaseDetection!=NULL && narrowPhaseDetection->getIntersectionMethod()!=intersectionMethod)
        narrowPhaseDetection->setIntersectionMethod(intersectionMethod);
    if (contactManager!=NULL && contactManager->getIntersectionMethod()!=intersectionMethod)
        contactManager->setIntersectionMethod(intersectionMethod);
    doCollisionReset();
}

void Pipeline::computeCollisionDetection()
{
    simulation::tree::GNode* root = dynamic_cast<simulation::tree::GNode*>(getContext());
    if(root == NULL) return;
    sofa::helper::vector<CollisionModel*> collisionModels;
    root->getTreeObjects<CollisionModel>(&collisionModels);
    doCollisionDetection(collisionModels);
}

void Pipeline::computeCollisionResponse()
{
    simulation::tree::GNode* root = dynamic_cast<simulation::tree::GNode*>(getContext());
    if(root == NULL) return;
    doCollisionResponse();
}

#endif

} // namespace collision

} // namespace core

} // namespace sofa

