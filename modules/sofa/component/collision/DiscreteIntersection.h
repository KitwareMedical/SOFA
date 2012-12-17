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
#ifndef SOFA_COMPONENT_COLLISION_DISCRETEINTERSECTION_H
#define SOFA_COMPONENT_COLLISION_DISCRETEINTERSECTION_H

#include <sofa/core/collision/Intersection.h>
#include <sofa/core/collision/IntersectorFactory.h>
#include <sofa/helper/FnDispatcher.h>
#include <sofa/component/collision/SphereModel.h>
#include <sofa/component/collision/CubeModel.h>
#include <sofa/component/collision/CapsuleModel.h>
#include <sofa/component/collision/CapsuleIntTool.h>

namespace sofa
{

namespace component
{

namespace collision
{
class SOFA_BASE_COLLISION_API DiscreteIntersection : public core::collision::Intersection, public core::collision::BaseIntersector
{
public:
    SOFA_CLASS(DiscreteIntersection,sofa::core::collision::Intersection);
protected:
    DiscreteIntersection();
public:
    /// Return the intersector class handling the given pair of collision models, or NULL if not supported.
    /// @param swapModel output value set to true if the collision models must be swapped before calling the intersector.
    virtual core::collision::ElementIntersector* findIntersector(core::CollisionModel* object1, core::CollisionModel* object2, bool& swapModels);

    core::collision::IntersectorMap intersectors;
    typedef core::collision::IntersectorFactory<DiscreteIntersection> IntersectorFactory;

    bool testIntersection(Cube&, Cube&);

    template <class Sphere>
    bool testIntersection(Sphere&, Sphere&);
    template <class Sphere>
    bool testIntersection(Sphere&, Cube&);
    bool testIntersection(Capsule&, Capsule&);
    bool testIntersection(Capsule&, Sphere&);

    int computeIntersection(Cube&, Cube&, OutputVector*);    
    template <class Sphere>
    int computeIntersection(Sphere&, Sphere&, OutputVector*);    
    template <class Sphere>
    int computeIntersection(Sphere&, Cube&, OutputVector*);
    int computeIntersection(Capsule&, Capsule&,OutputVector* contacts);
    int computeIntersection(Capsule&, Sphere&,OutputVector* contacts);
};

} // namespace collision

} // namespace component

namespace core
{
namespace collision
{
#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_BUILD_BASE_COLLISION)
extern template class SOFA_BASE_COLLISION_API IntersectorFactory<component::collision::DiscreteIntersection>;
#endif
}
}

} // namespace sofa

#endif
