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
#include "CudaTypes.h"
#include "CudaBarycentricMapping.inl"
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>

namespace sofa
{

namespace component
{

namespace mapping
{

using namespace sofa::defaulttype;
using namespace sofa::core;
using namespace sofa::core::behavior;
using namespace sofa::gpu::cuda;

template class BarycentricMapping< CudaVec3fTypes, CudaVec3fTypes>;
template class BarycentricMapping< CudaVec3fTypes, ExtVec3fTypes>;
template class BarycentricMapping< CudaVec3f1Types, CudaVec3f1Types>;
template class BarycentricMapping< CudaVec3f1Types, CudaVec3fTypes>;
template class BarycentricMapping< CudaVec3fTypes, CudaVec3f1Types>;
template class BarycentricMapping< CudaVec3f1Types, ExtVec3fTypes>;


#ifndef SOFA_FLOAT
template class BarycentricMapping< Vec3dTypes, CudaVec3fTypes>;
template class BarycentricMapping< CudaVec3fTypes, Vec3dTypes>;
template class BarycentricMapping< Vec3dTypes, CudaVec3f1Types>;
template class BarycentricMapping< CudaVec3f1Types, Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
template class BarycentricMapping< Vec3fTypes, CudaVec3fTypes>;
template class BarycentricMapping< CudaVec3fTypes, Vec3fTypes>;
template class BarycentricMapping< Vec3fTypes, CudaVec3f1Types>;
template class BarycentricMapping< CudaVec3f1Types, Vec3fTypes>;
#endif

} // namespace mapping

} // namespace component

namespace gpu
{

namespace cuda
{

using namespace sofa::defaulttype;
using namespace sofa::core;
using namespace sofa::core::behavior;
using namespace sofa::component::mapping;

SOFA_DECL_CLASS(CudaBarycentricMapping)

int BarycentricMappingCudaClass = core::RegisterObject("Supports GPU-side computations using CUDA")
        .add< BarycentricMapping< CudaVec3fTypes, CudaVec3fTypes> >()
        .add< BarycentricMapping< CudaVec3fTypes, ExtVec3fTypes> >()
        .add< BarycentricMapping< CudaVec3f1Types, CudaVec3f1Types> >()
        .add< BarycentricMapping< CudaVec3f1Types, CudaVec3fTypes> >()
        .add< BarycentricMapping< CudaVec3fTypes, CudaVec3f1Types> >()
        .add< BarycentricMapping< CudaVec3f1Types, ExtVec3fTypes> >()
#ifndef SOFA_FLOAT
        .add< BarycentricMapping< Vec3dTypes, CudaVec3fTypes> >()
        .add< BarycentricMapping< CudaVec3fTypes, Vec3dTypes> >()
        .add< BarycentricMapping< Vec3dTypes, CudaVec3f1Types> >()
        .add< BarycentricMapping< CudaVec3f1Types, Vec3dTypes> >()
#endif
#ifndef SOFA_DOUBLE
        .add< BarycentricMapping< Vec3fTypes, CudaVec3fTypes> >()
        .add< BarycentricMapping< CudaVec3fTypes, Vec3fTypes> >()
        .add< BarycentricMapping< Vec3fTypes, CudaVec3f1Types> >()
        .add< BarycentricMapping< CudaVec3f1Types, Vec3fTypes> >()
#endif
#ifdef SOFA_GPU_CUDA_DOUBLE
        .add< BarycentricMapping< CudaVec3fTypes, CudaVec3dTypes> >()
        .add< BarycentricMapping< CudaVec3dTypes, CudaVec3fTypes> >()
        .add< BarycentricMapping< CudaVec3dTypes, CudaVec3dTypes> >()
        .add< BarycentricMapping< CudaVec3dTypes, Vec3fTypes> >()
        .add< BarycentricMapping< CudaVec3dTypes, Vec3dTypes> >()
        .add< BarycentricMapping< Vec3fTypes, CudaVec3dTypes> >()
        .add< BarycentricMapping< Vec3dTypes, CudaVec3dTypes> >()
//.add< BarycentricMapping< CudaVec3d1Types, ExtVec3fTypes> >()
//.add< BarycentricMapping< CudaVec3dTypes, ExtVec3fTypes> >()
#endif
        ;

} // namespace cuda

} // namespace gpu

} // namespace sofa
