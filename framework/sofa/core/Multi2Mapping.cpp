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

#include <sofa/core/Multi2Mapping.inl>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>


namespace sofa
{

namespace core
{

using namespace sofa::defaulttype;
using namespace core::behavior;


template class Multi2Mapping< Vec3dTypes, Vec3fTypes, Vec3dTypes >;
template class Multi2Mapping< Vec3dTypes, Vec3fTypes, Vec3fTypes >;

template class Multi2Mapping< Vec3dTypes, Rigid3dTypes, Vec3dTypes >;
template class Multi2Mapping< Vec3dTypes, Rigid3fTypes, Vec3dTypes >;

template class Multi2Mapping< Vec3dTypes, Rigid3dTypes, Vec3fTypes >;
template class Multi2Mapping< Vec3dTypes, Rigid3fTypes, Vec3fTypes >;

template class Multi2Mapping< Vec3fTypes, Rigid3dTypes, Vec3dTypes >;
template class Multi2Mapping< Vec3fTypes, Rigid3fTypes, Vec3dTypes >;

template class Multi2Mapping< Vec3fTypes, Rigid3dTypes, Vec3fTypes >;
template class Multi2Mapping< Vec3fTypes, Rigid3fTypes, Vec3fTypes >;

template class Multi2Mapping< Vec3dTypes, Vec3fTypes, Rigid3dTypes >;
template class Multi2Mapping< Vec3dTypes, Vec3fTypes, Rigid3fTypes >;

template class Multi2Mapping< Vec3dTypes, Rigid3dTypes, Rigid3dTypes >;
template class Multi2Mapping< Vec3dTypes, Rigid3fTypes, Rigid3dTypes >;

template class Multi2Mapping< Vec3dTypes, Rigid3dTypes, Rigid3fTypes >;
template class Multi2Mapping< Vec3dTypes, Rigid3fTypes, Rigid3fTypes >;

template class Multi2Mapping< Vec3fTypes, Rigid3dTypes, Rigid3dTypes >;
template class Multi2Mapping< Vec3fTypes, Rigid3fTypes, Rigid3dTypes >;

template class Multi2Mapping< Vec3fTypes, Rigid3dTypes, Rigid3fTypes >;
template class Multi2Mapping< Vec3fTypes, Rigid3fTypes, Rigid3fTypes >;

template class Multi2Mapping< Vec1dTypes, Rigid3dTypes, Rigid3dTypes >;
template class Multi2Mapping< Vec1fTypes, Rigid3fTypes, Rigid3fTypes >;

template class Multi2Mapping< Vec1fTypes, Rigid3fTypes, Rigid3dTypes >;
template class Multi2Mapping< Vec1fTypes, Rigid3dTypes, Rigid3dTypes >;
template class Multi2Mapping< Vec1dTypes, Rigid3fTypes, Rigid3dTypes >;
template class Multi2Mapping< Vec1fTypes, Rigid3dTypes, Rigid3fTypes >;
template class Multi2Mapping< Vec1dTypes, Rigid3fTypes, Rigid3fTypes >;
template class Multi2Mapping< Vec1dTypes, Rigid3dTypes, Rigid3fTypes >;
}

}
