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

#include "MyMappingPendulumInPlane.inl"
#include <sofa/core/Mapping.inl>
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/Vec3Types.h>
#include <sofa/defaulttype/VecTypes.h>


namespace sofa
{

namespace component
{

namespace mapping
{

using namespace sofa::defaulttype;
using sofa::defaulttype::Vec3dTypes;
using sofa::defaulttype::Vec3fTypes;
using sofa::defaulttype::Vec2dTypes;
using sofa::defaulttype::Vec1dTypes;
using sofa::defaulttype::Vec1fTypes;



SOFA_DECL_CLASS(MyMappingPendulumInPlane)


int MyMappingPendulumInPlaneClass = core::RegisterObject("Mapping from an angle to a point in 2D")
#ifndef SOFA_FLOAT
        .add< MyMappingPendulumInPlane<Vec1dTypes,Vec3dTypes> >()
        .add< MyMappingPendulumInPlane<Vec1dTypes,Vec2dTypes> >()
#endif
#ifndef SOFA_DOUBLE
        .add< MyMappingPendulumInPlane<Vec1fTypes,Vec3fTypes> >()
#endif
        ;

#ifndef SOFA_FLOAT
template class MyMappingPendulumInPlane<Vec1dTypes,Vec3dTypes>;
template class MyMappingPendulumInPlane<Vec1dTypes,Vec2dTypes>;
#endif
#ifndef SOFA_DOUBLE
template class MyMappingPendulumInPlane<Vec1fTypes,Vec3fTypes>;
#endif



}	//mapping

}	//component

}	//sofa

