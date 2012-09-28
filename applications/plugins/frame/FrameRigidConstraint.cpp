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
*                               SOFA :: Plugins                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#define FRAME_FrameRigidConstraint_CPP

#include "QuadraticTypes.h"
#include "AffineTypes.h"
#include "FrameRigidConstraint.inl"
#include <sofa/component/projectiveconstraintset/FixedConstraint.inl>
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace projectiveconstraintset
{

using namespace sofa::defaulttype;

SOFA_DECL_CLASS(FrameRigidConstraint);

int FrameRigidConstraintClass = core::RegisterObject("Rigidify a deformable frame")
#ifndef SOFA_FLOAT
        .add< FrameRigidConstraint<Affine3dTypes> >()
        .add< FrameRigidConstraint<Quadratic3dTypes> >()
#endif
#ifndef SOFA_DOUBLE
        .add< FrameRigidConstraint<Affine3fTypes> >()
        .add< FrameRigidConstraint<Quadratic3fTypes> >()
#endif
        ;

#ifndef SOFA_FLOAT
template class SOFA_FRAME_API FrameRigidConstraint<Affine3dTypes>;

template class SOFA_FRAME_API FrameRigidConstraint<Affine3fTypes>;
#endif
#ifndef SOFA_DOUBLE

template class SOFA_FRAME_API FrameRigidConstraint<Quadratic3dTypes>;

template class SOFA_FRAME_API FrameRigidConstraint<Quadratic3fTypes>;
#endif

} // namespace projectiveconstraintset

} // namespace component

} // namespace sofa
