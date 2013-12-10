/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 MGH, INRIA, USTL, UJF, CNRS                    *
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
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#include <sofa/helper/system/config.h>
#include <sofa/helper/FnDispatcher.inl>
#include <sofa/component/collision/DiscreteIntersection.inl>
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/collision/Intersection.inl>
#include <sofa/helper/proximity.h>
#include <iostream>
#include <algorithm>

namespace sofa
{

namespace component
{

namespace collision
{

using namespace sofa::defaulttype;
using namespace sofa::core::collision;

SOFA_DECL_CLASS(DiscreteIntersection)

int DiscreteIntersectionClass = core::RegisterObject("TODO-DiscreteIntersectionClass")
        .add< DiscreteIntersection >()
        ;


DiscreteIntersection::DiscreteIntersection()
{
    intersectors.add<CubeModel,       CubeModel,         DiscreteIntersection> (this);

    intersectors.add<SphereModel,     SphereModel,       DiscreteIntersection> (this);

    intersectors.add<CapsuleModel,CapsuleModel, DiscreteIntersection> (this);
    intersectors.add<CapsuleModel,SphereModel, DiscreteIntersection> (this);

    intersectors.add<OBBModel,OBBModel,DiscreteIntersection>(this);
    intersectors.add<SphereModel,OBBModel, DiscreteIntersection> (this);
    intersectors.add<CapsuleModel,OBBModel,DiscreteIntersection>(this);

    intersectors.add<RigidSphereModel,RigidSphereModel,DiscreteIntersection>(this);
    intersectors.add<SphereModel,RigidSphereModel, DiscreteIntersection> (this);
    intersectors.add<CapsuleModel,RigidSphereModel,DiscreteIntersection>(this);
    intersectors.add<RigidSphereModel,OBBModel,DiscreteIntersection>(this);

    //IntersectorFactory::getInstance()->addIntersectors(this);
}

/// Return the intersector class handling the given pair of collision models, or NULL if not supported.
ElementIntersector* DiscreteIntersection::findIntersector(core::CollisionModel* object1, core::CollisionModel* object2, bool& swapModels)
{
    return intersectors.get(object1, object2, swapModels);
}

bool DiscreteIntersection::testIntersection(Cube& cube1, Cube& cube2)
{    
    return BaseIntTool::testIntersection(cube1,cube2,getAlarmDistance());
}

int DiscreteIntersection::computeIntersection(Cube&, Cube&, OutputVector*)
{
    return 0; /// \todo
}



bool DiscreteIntersection::testIntersection(Capsule&, Capsule&){
    //TO DO
    return false;
}


int DiscreteIntersection::computeIntersection(Capsule & e1,Capsule & e2,OutputVector * contacts){
    return CapsuleIntTool::computeIntersection(e1,e2,getAlarmDistance(),getContactDistance(),contacts);
}




bool DiscreteIntersection::testIntersection(OBB &,OBB &){
    return false;
}

bool DiscreteIntersection::testIntersection(Capsule &,OBB &){
    return false;
}




int DiscreteIntersection::computeIntersection(OBB & box0, OBB & box1,OutputVector* contacts){
    return OBBIntTool::computeIntersection(box0,box1,getAlarmDistance(), getContactDistance(),contacts);
}

int DiscreteIntersection::computeIntersection(Capsule & cap, OBB & box,OutputVector* contacts){
    return CapsuleIntTool::computeIntersection(cap,box,getAlarmDistance(),getContactDistance(),contacts);
}



} // namespace collision

} // namespace component

namespace core
{
namespace collision
{
template class SOFA_BASE_COLLISION_API IntersectorFactory<component::collision::DiscreteIntersection>;
}
}

} // namespace sofa

