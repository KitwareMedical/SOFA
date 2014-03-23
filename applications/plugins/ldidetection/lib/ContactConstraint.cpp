/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 3      *
*                (c) 2006-2008 MGH, INRIA, USTL, UJF, CNRS                    *
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
#include "ContactConstraint.inl"
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/Vec3Types.h>

namespace sofa
{

namespace component
{

namespace constraint
{

using namespace sofa::defaulttype;
using namespace sofa::helper;


SOFA_DECL_CLASS(ContactConstraint)

int ContactConstraintClass = core::RegisterObject("Maintain the length of the edges of an object within a given contact rate")
#ifndef SOFA_FLOAT
.add< ContactConstraint<Vec3dTypes> >()
// .add< ContactConstraint<Vec2dTypes> >()
// .add< ContactConstraint<Vec1dTypes> >()
// .add< ContactConstraint<Vec6dTypes> >()
#endif
#ifndef SOFA_DOUBLE
.add< ContactConstraint<Vec3fTypes> >()
// .add< ContactConstraint<Vec2fTypes> >()
// .add< ContactConstraint<Vec1fTypes> >()
// .add< ContactConstraint<Vec6fTypes> >()
#endif
;

#ifndef SOFA_FLOAT
template class ContactConstraint<Vec3dTypes>;
// template class ContactConstraint<Vec2dTypes>;
// template class ContactConstraint<Vec1dTypes>;
// template class ContactConstraint<Vec6dTypes>;
#endif
#ifndef SOFA_DOUBLE
template class ContactConstraint<Vec3fTypes>;
// template class ContactConstraint<Vec2fTypes>;
// template class ContactConstraint<Vec1fTypes>;
// template class ContactConstraint<Vec6fTypes>;
#endif 




} // namespace constraint

} // namespace component

} // namespace sofa

