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
#define SOFA_COMPONENT_PROJECTIVECONSTRAINTSET_PARTIALLINEARMOVEMENTCONSTRAINT_CPP
#include <sofa/component/projectiveconstraintset/PartialLinearMovementConstraint.inl>
#include <sofa/core/behavior/ProjectiveConstraintSet.inl>
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/Vec3Types.h>
#include <sofa/defaulttype/RigidTypes.h>

#include <sofa/simulation/common/Node.h>

namespace sofa
{

namespace component
{

namespace projectiveconstraintset
{

//declaration of the class, for the factory
SOFA_DECL_CLASS(PartialLinearMovementConstraint)


int PartialLinearMovementConstraintClass = core::RegisterObject("translate given particles")
#ifndef SOFA_FLOAT
        .add< PartialLinearMovementConstraint<Vec3dTypes> >()
        .add< PartialLinearMovementConstraint<Vec2dTypes> >()
        .add< PartialLinearMovementConstraint<Vec1dTypes> >()
        .add< PartialLinearMovementConstraint<Vec6dTypes> >()
        .add< PartialLinearMovementConstraint<Rigid3dTypes> >()
#endif
#ifndef SOFA_DOUBLE
        .add< PartialLinearMovementConstraint<Vec3fTypes> >()
        .add< PartialLinearMovementConstraint<Vec2fTypes> >()
        .add< PartialLinearMovementConstraint<Vec1fTypes> >()
        .add< PartialLinearMovementConstraint<Vec6fTypes> >()
        .add< PartialLinearMovementConstraint<Rigid3fTypes> >()
#endif
        ;

#ifndef SOFA_FLOAT
template class SOFA_BOUNDARY_CONDITION_API PartialLinearMovementConstraint<Vec3dTypes>;
template class SOFA_BOUNDARY_CONDITION_API PartialLinearMovementConstraint<Vec2dTypes>;
template class SOFA_BOUNDARY_CONDITION_API PartialLinearMovementConstraint<Vec1dTypes>;
template class SOFA_BOUNDARY_CONDITION_API PartialLinearMovementConstraint<Vec6dTypes>;
template class SOFA_BOUNDARY_CONDITION_API PartialLinearMovementConstraint<Rigid3dTypes>;
#endif
#ifndef SOFA_DOUBLE
template class SOFA_BOUNDARY_CONDITION_API PartialLinearMovementConstraint<Vec3fTypes>;
template class SOFA_BOUNDARY_CONDITION_API PartialLinearMovementConstraint<Vec2fTypes>;
template class SOFA_BOUNDARY_CONDITION_API PartialLinearMovementConstraint<Vec1fTypes>;
template class SOFA_BOUNDARY_CONDITION_API PartialLinearMovementConstraint<Vec6fTypes>;
template class SOFA_BOUNDARY_CONDITION_API PartialLinearMovementConstraint<Rigid3fTypes>;
#endif

} // namespace projectiveconstraintset

} // namespace component

} // namespace sofa
