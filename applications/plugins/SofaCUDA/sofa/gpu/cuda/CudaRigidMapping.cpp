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
#include "CudaRigidMapping.inl"
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
using namespace sofa::gpu::cuda;

template class RigidMapping< CudaRigid3fTypes, CudaVec3fTypes>;
template class RigidMapping< Rigid3fTypes, CudaVec3fTypes>;
template class RigidMapping< Rigid3dTypes, CudaVec3fTypes>;
//template class RigidMapping< CudaRigid3fTypes, Vec3dTypes>;
//template class RigidMapping< CudaRigid3fTypes, Vec3fTypes>;
template class RigidMapping< CudaRigid3fTypes, CudaVec3f1Types>;
template class RigidMapping< Rigid3fTypes, CudaVec3f1Types>;
template class RigidMapping< Rigid3dTypes, CudaVec3f1Types>;

#ifdef SOFA_GPU_CUDA_DOUBLE
template class RigidMapping< CudaRigid3fTypes, CudaVec3dTypes>;
template class RigidMapping< Rigid3fTypes, CudaVec3dTypes>;
template class RigidMapping< Rigid3dTypes, CudaVec3dTypes>;
//template class RigidMapping< CudaRigid3fTypes, Vec3dTypes>;
//template class RigidMapping< CudaRigid3fTypes, Vec3fTypes>;
template class RigidMapping< CudaRigid3fTypes, CudaVec3d1Types>;
template class RigidMapping< Rigid3fTypes, CudaVec3d1Types>;
template class RigidMapping< Rigid3dTypes, CudaVec3d1Types>;
#endif // SOFA_GPU_CUDA_DOUBLE
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

SOFA_DECL_CLASS(CudaRigidMapping)

int RigidMappingCudaClass = core::RegisterObject("Supports GPU-side computations using CUDA")
        .add< RigidMapping< CudaRigid3fTypes, CudaVec3fTypes> >()
        .add< RigidMapping< Rigid3fTypes, CudaVec3fTypes> >()
        .add< RigidMapping< Rigid3dTypes, CudaVec3fTypes> >()
//.add< RigidMapping< CudaRigid3fTypes, Vec3dTypes> >()
//.add< RigidMapping< CudaRigid3fTypes, Vec3fTypes> >()
        .add< RigidMapping< CudaRigid3fTypes, CudaVec3f1Types> >()
        .add< RigidMapping< Rigid3fTypes, CudaVec3f1Types> >()
        .add< RigidMapping< Rigid3dTypes, CudaVec3f1Types> >()
#ifdef SOFA_GPU_CUDA_DOUBLE
        .add< RigidMapping< CudaRigid3fTypes, CudaVec3dTypes> >()
        .add< RigidMapping< Rigid3fTypes, CudaVec3dTypes> >()
        .add< RigidMapping< Rigid3dTypes, CudaVec3dTypes> >()
//.add< RigidMapping< CudaRigid3fTypes, Vec3dTypes> >()
//.add< RigidMapping< CudaRigid3fTypes, Vec3fTypes> >()
        .add< RigidMapping< CudaRigid3fTypes, CudaVec3d1Types> >()
        .add< RigidMapping< Rigid3fTypes, CudaVec3d1Types> >()
        .add< RigidMapping< Rigid3dTypes, CudaVec3d1Types> >()
#endif // SOFA_GPU_CUDA_DOUBLE
        ;

} // namespace cuda

} // namespace gpu

} // namespace sofa
