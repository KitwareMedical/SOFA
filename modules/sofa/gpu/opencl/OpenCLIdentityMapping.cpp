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
#include "OpenCLTypes.h"
#include "OpenCLIdentityMapping.inl"
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/Mapping.inl>

namespace sofa
{

namespace component
{

namespace mapping
{

using namespace sofa::defaulttype;
using namespace sofa::core;
using namespace sofa::core::behavior;
using namespace sofa::gpu::opencl;

template class  IdentityMapping< OpenCLVec3fTypes, OpenCLVec3fTypes>;
template class  IdentityMapping< OpenCLVec3fTypes, Vec3fTypes>;
template class  IdentityMapping< OpenCLVec3fTypes, Vec3dTypes>;
template class  IdentityMapping< Vec3fTypes, OpenCLVec3fTypes>;
template class  IdentityMapping< Vec3dTypes, OpenCLVec3fTypes>;

template class  IdentityMapping< OpenCLVec3fTypes, OpenCLVec3dTypes>;
template class  IdentityMapping< OpenCLVec3dTypes, OpenCLVec3fTypes>;
template class  IdentityMapping< OpenCLVec3dTypes, OpenCLVec3dTypes>;
template class  IdentityMapping< OpenCLVec3dTypes, Vec3fTypes>;
template class  IdentityMapping< OpenCLVec3dTypes, Vec3dTypes>;
template class  IdentityMapping< Vec3fTypes, OpenCLVec3dTypes>;
template class  IdentityMapping< Vec3dTypes, OpenCLVec3dTypes>;

template class  IdentityMapping< OpenCLVec3d1Types, ExtVec3fTypes>;
template class  IdentityMapping< OpenCLVec3dTypes, ExtVec3fTypes>;

template class  IdentityMapping< OpenCLVec3fTypes, ExtVec3fTypes>;
template class  IdentityMapping< OpenCLVec3f1Types, OpenCLVec3f1Types>;
template class  IdentityMapping< OpenCLVec3f1Types, Vec3dTypes>;
template class  IdentityMapping< OpenCLVec3f1Types, Vec3fTypes>;
template class  IdentityMapping< Vec3dTypes, OpenCLVec3f1Types>;
template class  IdentityMapping< Vec3fTypes, OpenCLVec3f1Types>;
template class  IdentityMapping< OpenCLVec3f1Types, ExtVec3fTypes>;
template class  IdentityMapping< OpenCLVec3f1Types, OpenCLVec3fTypes>;
template class  IdentityMapping< OpenCLVec3fTypes, OpenCLVec3f1Types>;

} // namespace mapping

} // namespace component

namespace gpu
{

namespace opencl
{

using namespace sofa::defaulttype;
using namespace sofa::core;
using namespace sofa::core::behavior;
using namespace sofa::component::mapping;

SOFA_DECL_CLASS(OpenCLIdentityMapping)

int IdentityMappingOpenCLClass = core::RegisterObject("Supports GPU-side computations using OPENCL")
        .add< IdentityMapping< OpenCLVec3fTypes, OpenCLVec3fTypes> >()
        .add< IdentityMapping< OpenCLVec3fTypes, Vec3fTypes> >()
        .add< IdentityMapping< OpenCLVec3fTypes, Vec3dTypes> >()
        .add< IdentityMapping< Vec3fTypes, OpenCLVec3fTypes> >()
        .add< IdentityMapping< Vec3dTypes, OpenCLVec3fTypes> >()

        .add< IdentityMapping< OpenCLVec3fTypes, OpenCLVec3dTypes> >()
        .add< IdentityMapping< OpenCLVec3dTypes, OpenCLVec3fTypes> >()
        .add< IdentityMapping< OpenCLVec3dTypes, OpenCLVec3dTypes> >()
        .add< IdentityMapping< OpenCLVec3dTypes, Vec3fTypes> >()
        .add< IdentityMapping< OpenCLVec3dTypes, Vec3dTypes> >()
        .add< IdentityMapping< Vec3fTypes, OpenCLVec3dTypes> >()
        .add< IdentityMapping< Vec3dTypes, OpenCLVec3dTypes> >()

        .add< IdentityMapping< OpenCLVec3d1Types, ExtVec3fTypes> >()
        .add< IdentityMapping< OpenCLVec3dTypes, ExtVec3fTypes> >()

        .add< IdentityMapping< OpenCLVec3fTypes, ExtVec3fTypes> >()
        .add< IdentityMapping< OpenCLVec3f1Types, OpenCLVec3f1Types> >()
        .add< IdentityMapping< OpenCLVec3f1Types, Vec3dTypes> >()
        .add< IdentityMapping< OpenCLVec3f1Types, Vec3fTypes> >()
        .add< IdentityMapping< Vec3dTypes, OpenCLVec3f1Types> >()
        .add< IdentityMapping< Vec3fTypes, OpenCLVec3f1Types> >()
        .add< IdentityMapping< OpenCLVec3f1Types, ExtVec3fTypes> >()
        .add< IdentityMapping< OpenCLVec3f1Types, OpenCLVec3fTypes> >()
        .add< IdentityMapping< OpenCLVec3fTypes, OpenCLVec3f1Types> >()
        ;

} // namespace opencl

} // namespace gpu

} // namespace sofa
