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
#include <sofa/core/behavior/LMConstraint.inl>
#include <sofa/defaulttype/Vec3Types.h>
#include <sofa/defaulttype/RigidTypes.h>

namespace sofa
{

namespace core
{

namespace behavior
{

using namespace sofa::defaulttype;

template class SOFA_CORE_API LMConstraint<Vec3dTypes,Vec3dTypes>;
template class SOFA_CORE_API LMConstraint<Vec3dTypes,Vec2dTypes>;
template class SOFA_CORE_API LMConstraint<Vec3dTypes,Vec1dTypes>;
template class SOFA_CORE_API LMConstraint<Vec3dTypes,Rigid3dTypes>;
template class SOFA_CORE_API LMConstraint<Vec3dTypes,Rigid2dTypes>;

template class SOFA_CORE_API LMConstraint<Vec2dTypes,Vec3dTypes>;
template class SOFA_CORE_API LMConstraint<Vec2dTypes,Vec2dTypes>;
template class SOFA_CORE_API LMConstraint<Vec2dTypes,Vec1dTypes>;
template class SOFA_CORE_API LMConstraint<Vec2dTypes,Rigid3dTypes>;
template class SOFA_CORE_API LMConstraint<Vec2dTypes,Rigid2dTypes>;

template class SOFA_CORE_API LMConstraint<Vec1dTypes,Vec3dTypes>;
template class SOFA_CORE_API LMConstraint<Vec1dTypes,Vec2dTypes>;
template class SOFA_CORE_API LMConstraint<Vec1dTypes,Vec1dTypes>;
template class SOFA_CORE_API LMConstraint<Vec1dTypes,Rigid3dTypes>;
template class SOFA_CORE_API LMConstraint<Vec1dTypes,Rigid2dTypes>;

template class SOFA_CORE_API LMConstraint<Rigid3dTypes,Vec3dTypes>;
template class SOFA_CORE_API LMConstraint<Rigid3dTypes,Vec2dTypes>;
template class SOFA_CORE_API LMConstraint<Rigid3dTypes,Vec1dTypes>;
template class SOFA_CORE_API LMConstraint<Rigid3dTypes,Rigid3dTypes>;
template class SOFA_CORE_API LMConstraint<Rigid3dTypes,Rigid2dTypes>;

template class SOFA_CORE_API LMConstraint<Rigid2dTypes,Vec3dTypes>;
template class SOFA_CORE_API LMConstraint<Rigid2dTypes,Vec2dTypes>;
template class SOFA_CORE_API LMConstraint<Rigid2dTypes,Vec1dTypes>;
template class SOFA_CORE_API LMConstraint<Rigid2dTypes,Rigid3dTypes>;
template class SOFA_CORE_API LMConstraint<Rigid2dTypes,Rigid2dTypes>;

template class SOFA_CORE_API LMConstraint<Vec3fTypes,Vec3fTypes>;
template class SOFA_CORE_API LMConstraint<Vec3fTypes,Vec2fTypes>;
template class SOFA_CORE_API LMConstraint<Vec3fTypes,Vec1fTypes>;
template class SOFA_CORE_API LMConstraint<Vec3fTypes,Rigid3fTypes>;
template class SOFA_CORE_API LMConstraint<Vec3fTypes,Rigid2fTypes>;

template class SOFA_CORE_API LMConstraint<Vec2fTypes,Vec3fTypes>;
template class SOFA_CORE_API LMConstraint<Vec2fTypes,Vec2fTypes>;
template class SOFA_CORE_API LMConstraint<Vec2fTypes,Vec1fTypes>;
template class SOFA_CORE_API LMConstraint<Vec2fTypes,Rigid3fTypes>;
template class SOFA_CORE_API LMConstraint<Vec2fTypes,Rigid2fTypes>;

template class SOFA_CORE_API LMConstraint<Vec1fTypes,Vec3fTypes>;
template class SOFA_CORE_API LMConstraint<Vec1fTypes,Vec2fTypes>;
template class SOFA_CORE_API LMConstraint<Vec1fTypes,Vec1fTypes>;
template class SOFA_CORE_API LMConstraint<Vec1fTypes,Rigid3fTypes>;
template class SOFA_CORE_API LMConstraint<Vec1fTypes,Rigid2fTypes>;

template class SOFA_CORE_API LMConstraint<Rigid3fTypes,Vec3fTypes>;
template class SOFA_CORE_API LMConstraint<Rigid3fTypes,Vec2fTypes>;
template class SOFA_CORE_API LMConstraint<Rigid3fTypes,Vec1fTypes>;
template class SOFA_CORE_API LMConstraint<Rigid3fTypes,Rigid3fTypes>;
template class SOFA_CORE_API LMConstraint<Rigid3fTypes,Rigid2fTypes>;

template class SOFA_CORE_API LMConstraint<Rigid2fTypes,Vec3fTypes>;
template class SOFA_CORE_API LMConstraint<Rigid2fTypes,Vec2fTypes>;
template class SOFA_CORE_API LMConstraint<Rigid2fTypes,Vec1fTypes>;
template class SOFA_CORE_API LMConstraint<Rigid2fTypes,Rigid3fTypes>;
template class SOFA_CORE_API LMConstraint<Rigid2fTypes,Rigid2fTypes>;


//Need the combinations

} // namespace behavior

} // namespace core

} // namespace sofa
