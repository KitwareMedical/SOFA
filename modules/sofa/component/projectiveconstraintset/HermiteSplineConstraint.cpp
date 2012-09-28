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
#include <sofa/component/projectiveconstraintset/HermiteSplineConstraint.inl>
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/Vec3Types.h>
#include <sofa/defaulttype/RigidTypes.h>


namespace sofa
{

namespace component
{

namespace projectiveconstraintset
{

SOFA_DECL_CLASS(HermiteSplineConstraint)


int HermiteSplineConstraintClass = core::RegisterObject("Apply a hermite cubic spline trajectory to given points")
#ifndef SOFA_FLOAT
        .add< HermiteSplineConstraint<Vec3dTypes> >()
        .add< HermiteSplineConstraint<Rigid3dTypes> >()
#endif
#ifndef SOFA_DOUBLE
        .add< HermiteSplineConstraint<Vec3fTypes> >()
        .add< HermiteSplineConstraint<Rigid3fTypes> >()
#endif
        ;

#ifndef SOFA_FLOAT
template class HermiteSplineConstraint<Rigid3dTypes>;
template class HermiteSplineConstraint<Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
template class HermiteSplineConstraint<Rigid3fTypes>;
template class HermiteSplineConstraint<Vec3fTypes>;
#endif



} // namespace projectiveconstraintset

} // namespace component

} // namespace sofa

