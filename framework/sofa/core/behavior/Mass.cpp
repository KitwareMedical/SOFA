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
#include <sofa/core/behavior/Mass.inl>

namespace sofa
{

namespace core
{

namespace behavior
{

using namespace sofa::defaulttype;

template class SOFA_CORE_API Mass<Vec3dTypes>;
template class SOFA_CORE_API Mass<Vec2dTypes>;
template class SOFA_CORE_API Mass<Vec1dTypes>;
template class SOFA_CORE_API Mass<Vec6dTypes>;
template class SOFA_CORE_API Mass<Rigid3dTypes>;
template class SOFA_CORE_API Mass<Rigid2dTypes>;

template class SOFA_CORE_API Mass<Vec3fTypes>;
template class SOFA_CORE_API Mass<Vec2fTypes>;
template class SOFA_CORE_API Mass<Vec1fTypes>;
template class SOFA_CORE_API Mass<Vec6fTypes>;
template class SOFA_CORE_API Mass<Rigid3fTypes>;
template class SOFA_CORE_API Mass<Rigid2fTypes>;

} // namespace behavior

} // namespace core

} // namespace sofa
