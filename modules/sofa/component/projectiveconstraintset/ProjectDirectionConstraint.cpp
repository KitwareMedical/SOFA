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
#define SOFA_COMPONENT_PROJECTIVECONSTRAINTSET_ProjectDirectionConstraint_CPP
#include <sofa/component/projectiveconstraintset/ProjectDirectionConstraint.inl>
#include <sofa/core/behavior/ProjectiveConstraintSet.inl>
#include <sofa/core/ObjectFactory.h>

#include <sofa/simulation/common/Node.h>

namespace sofa
{

namespace component
{

namespace projectiveconstraintset
{

using namespace sofa::defaulttype;
using namespace sofa::helper;


SOFA_DECL_CLASS(ProjectDirectionConstraint)

int ProjectDirectionConstraintClass = core::RegisterObject("Attach given particles to their initial positions")
#ifndef SOFA_FLOAT
        .add< ProjectDirectionConstraint<Vec3dTypes> >()
        .add< ProjectDirectionConstraint<Vec2dTypes> >()
//.add< ProjectDirectionConstraint<Vec1dTypes> >()
//.add< ProjectDirectionConstraint<Vec6dTypes> >()
//.add< ProjectDirectionConstraint<Rigid3dTypes> >()
//.add< ProjectDirectionConstraint<Rigid2dTypes> >()
#endif
#ifndef SOFA_DOUBLE
        .add< ProjectDirectionConstraint<Vec3fTypes> >()
        .add< ProjectDirectionConstraint<Vec2fTypes> >()
//.add< ProjectDirectionConstraint<Vec1fTypes> >()
//.add< ProjectDirectionConstraint<Vec6fTypes> >()
//.add< ProjectDirectionConstraint<Rigid3fTypes> >()
//.add< ProjectDirectionConstraint<Rigid2fTypes> >()
#endif
        ;

#ifndef SOFA_FLOAT
template class SOFA_BOUNDARY_CONDITION_API ProjectDirectionConstraint<Vec3dTypes>;
template class SOFA_BOUNDARY_CONDITION_API ProjectDirectionConstraint<Vec2dTypes>;
//template class SOFA_BOUNDARY_CONDITION_API ProjectDirectionConstraint<Vec1dTypes>;
//template class SOFA_BOUNDARY_CONDITION_API ProjectDirectionConstraint<Vec6dTypes>;
//template class SOFA_BOUNDARY_CONDITION_API ProjectDirectionConstraint<Rigid3dTypes>;
//template class SOFA_BOUNDARY_CONDITION_API ProjectDirectionConstraint<Rigid2dTypes>;
#endif
#ifndef SOFA_DOUBLE
template class SOFA_BOUNDARY_CONDITION_API ProjectDirectionConstraint<Vec3fTypes>;
template class SOFA_BOUNDARY_CONDITION_API ProjectDirectionConstraint<Vec2fTypes>;
//template class SOFA_BOUNDARY_CONDITION_API ProjectDirectionConstraint<Vec1fTypes>;
//template class SOFA_BOUNDARY_CONDITION_API ProjectDirectionConstraint<Vec6fTypes>;
//template class SOFA_BOUNDARY_CONDITION_API ProjectDirectionConstraint<Rigid3fTypes>;
//template class SOFA_BOUNDARY_CONDITION_API ProjectDirectionConstraint<Rigid2fTypes>;
#endif



} // namespace projectiveconstraintset

} // namespace component

} // namespace sofa

