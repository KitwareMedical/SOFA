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
#define SOFA_COMPONENT_MAPPING_RIGIDMAPPING_CPP
#include <sofa/component/mapping/RigidMapping.inl>
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace mapping
{

SOFA_DECL_CLASS(RigidMapping)

using namespace defaulttype;

// Register in the Factory
int RigidMappingClass = core::RegisterObject("Set the positions and velocities of points attached to a rigid parent")
#ifndef SOFA_FLOAT
        .add< RigidMapping< Rigid3dTypes, Vec3dTypes > >()
        .add< RigidMapping< Rigid2dTypes, Vec2dTypes > >()
        .add< RigidMapping< Rigid3dTypes, ExtVec3fTypes > >()
#endif
#ifndef SOFA_DOUBLE
        .add< RigidMapping< Rigid3fTypes, Vec3fTypes > >()
        .add< RigidMapping< Rigid2fTypes, Vec2fTypes > >()
        .add< RigidMapping< Rigid3fTypes, ExtVec3fTypes > >()
#endif

#ifndef SOFA_FLOAT
#ifndef SOFA_DOUBLE
        .add< RigidMapping< Rigid3dTypes, Vec3fTypes > >()
        .add< RigidMapping< Rigid3fTypes, Vec3dTypes > >()
        .add< RigidMapping< Rigid2dTypes, Vec2fTypes > >()
        .add< RigidMapping< Rigid2fTypes, Vec2dTypes > >()
#endif
#endif
        ;

#ifndef SOFA_FLOAT
template class SOFA_RIGID_API RigidMapping< Rigid3dTypes, Vec3dTypes >;
template class SOFA_RIGID_API RigidMapping< Rigid2dTypes, Vec2dTypes >;
template class SOFA_RIGID_API RigidMapping< Rigid3dTypes, ExtVec3fTypes >;
#endif

#ifndef SOFA_DOUBLE
template class SOFA_RIGID_API RigidMapping< Rigid3fTypes, Vec3fTypes >;
template class SOFA_RIGID_API RigidMapping< Rigid2fTypes, Vec2fTypes >;
template class SOFA_RIGID_API RigidMapping< Rigid3fTypes, ExtVec3fTypes >;
#endif

#ifndef SOFA_FLOAT
#ifndef SOFA_DOUBLE
template class SOFA_RIGID_API RigidMapping< Rigid3dTypes, Vec3fTypes >;
template class SOFA_RIGID_API RigidMapping< Rigid3fTypes, Vec3dTypes >;
template class SOFA_RIGID_API RigidMapping< Rigid2dTypes, Vec2fTypes >;
template class SOFA_RIGID_API RigidMapping< Rigid2fTypes, Vec2dTypes >;
#endif
#endif

/// used to put a breakpoint, because gdb not always succeed in breaking within template methods
void rigidMappingDummyFunction()
{
    std::cerr << "rigidMappingDummyFunction()" << std::endl;
}

} // namespace mapping

} // namespace component

} // namespace sofa

