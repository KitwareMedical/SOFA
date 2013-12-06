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
#include "Mapping.inl"
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>

namespace sofa
{

namespace core
{

using namespace sofa::defaulttype;
using namespace core;

#ifndef SOFA_FLOAT
template class SOFA_CORE_API Mapping< sofa::defaulttype::Vec1dTypes, sofa::defaulttype::Vec1dTypes >;
template class SOFA_CORE_API Mapping< sofa::defaulttype::Vec2dTypes, sofa::defaulttype::Vec2dTypes >;
template class SOFA_CORE_API Mapping< sofa::defaulttype::Vec3dTypes, sofa::defaulttype::Vec3dTypes >;
template class SOFA_CORE_API Mapping< sofa::defaulttype::Vec3dTypes, sofa::defaulttype::Vec1dTypes >;
template class SOFA_CORE_API Mapping< sofa::defaulttype::Rigid3dTypes, sofa::defaulttype::Vec3dTypes >;
template class SOFA_CORE_API Mapping< sofa::defaulttype::Rigid3dTypes, sofa::defaulttype::Rigid3dTypes >;
template class SOFA_CORE_API Mapping< sofa::defaulttype::Rigid2dTypes, sofa::defaulttype::Vec2dTypes >;
template class SOFA_CORE_API Mapping< sofa::defaulttype::Rigid2dTypes, sofa::defaulttype::Rigid2dTypes >;
template class SOFA_CORE_API Mapping< sofa::defaulttype::Vec3dTypes, sofa::defaulttype::ExtVec3fTypes >;  // this one is special; ExtVec3fTypes are used for output outside of Sofa.
#endif

#ifndef SOFA_DOUBLE
template class SOFA_CORE_API Mapping< sofa::defaulttype::Vec1fTypes, sofa::defaulttype::Vec1fTypes >;
template class SOFA_CORE_API Mapping< sofa::defaulttype::Vec2fTypes, sofa::defaulttype::Vec2fTypes >;
template class SOFA_CORE_API Mapping< sofa::defaulttype::Vec3fTypes, sofa::defaulttype::Vec3fTypes >;
template class SOFA_CORE_API Mapping< sofa::defaulttype::Vec3fTypes, sofa::defaulttype::Vec1fTypes >;
template class SOFA_CORE_API Mapping< sofa::defaulttype::Vec3fTypes, sofa::defaulttype::ExtVec3fTypes >;
template class SOFA_CORE_API Mapping< sofa::defaulttype::Rigid3fTypes, sofa::defaulttype::Vec3fTypes >;
template class SOFA_CORE_API Mapping< sofa::defaulttype::Rigid3fTypes, sofa::defaulttype::Rigid3fTypes >;
template class SOFA_CORE_API Mapping< sofa::defaulttype::Rigid2fTypes, sofa::defaulttype::Vec2fTypes >;
template class SOFA_CORE_API Mapping< sofa::defaulttype::Rigid2fTypes, sofa::defaulttype::Rigid2fTypes >;
#endif

#ifndef SOFA_FLOAT
#ifndef SOFA_DOUBLE
template class SOFA_CORE_API Mapping< sofa::defaulttype::Vec1dTypes, sofa::defaulttype::Vec1fTypes >;
template class SOFA_CORE_API Mapping< sofa::defaulttype::Vec1fTypes, sofa::defaulttype::Vec1dTypes > ;
template class SOFA_CORE_API Mapping< sofa::defaulttype::Vec2dTypes, sofa::defaulttype::Vec2fTypes >;
template class SOFA_CORE_API Mapping< sofa::defaulttype::Vec2fTypes, sofa::defaulttype::Vec2dTypes > ;
template class SOFA_CORE_API Mapping< sofa::defaulttype::Vec3dTypes, sofa::defaulttype::Vec3fTypes >;
template class SOFA_CORE_API Mapping< sofa::defaulttype::Vec3fTypes, sofa::defaulttype::Vec3dTypes > ;
template class SOFA_CORE_API Mapping< sofa::defaulttype::Rigid3dTypes, sofa::defaulttype::Vec3fTypes >;
template class SOFA_CORE_API Mapping< sofa::defaulttype::Rigid3fTypes, sofa::defaulttype::Vec3dTypes >;
template class SOFA_CORE_API Mapping< sofa::defaulttype::Rigid3dTypes, sofa::defaulttype::Rigid3fTypes >;
template class SOFA_CORE_API Mapping< sofa::defaulttype::Rigid3fTypes, sofa::defaulttype::Rigid3dTypes >;
template class SOFA_CORE_API Mapping< sofa::defaulttype::Rigid2dTypes, sofa::defaulttype::Vec2fTypes >;
template class SOFA_CORE_API Mapping< sofa::defaulttype::Rigid2fTypes, sofa::defaulttype::Vec2dTypes >;
template class SOFA_CORE_API Mapping< sofa::defaulttype::Rigid2dTypes, sofa::defaulttype::Rigid2fTypes >;
template class SOFA_CORE_API Mapping< sofa::defaulttype::Rigid2fTypes, sofa::defaulttype::Rigid2dTypes >;
#endif
#endif

} // namespace core

} // namespace sofa

