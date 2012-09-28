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
#define SOFA_COMPONENT_ENGINE_TRANSFORMPOSITION_CPP
#include <sofa/component/engine/TransformPosition.inl>
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/Vec3Types.h>
#include <sofa/defaulttype/RigidTypes.h>

namespace sofa
{

namespace component
{

namespace engine
{

SOFA_DECL_CLASS(TransformPosition)

int TransformPositionClass = core::RegisterObject("Transform position of 3d points")
#ifdef SOFA_FLOAT
        .add< TransformPosition<Vec3fTypes> >(true)
#else
        .add< TransformPosition<Vec3dTypes> >(true)
#ifndef SOFA_DOUBLE
        .add< TransformPosition<Vec3fTypes> >()
#endif //SOFA_DOUBLE
#endif
        ;

#ifndef SOFA_FLOAT
template class SOFA_ENGINE_API TransformPosition<Vec3dTypes>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
template class SOFA_ENGINE_API TransformPosition<Vec3fTypes>;
#endif //SOFA_DOUBLE


} // namespace constraint

} // namespace component

} // namespace sofa

