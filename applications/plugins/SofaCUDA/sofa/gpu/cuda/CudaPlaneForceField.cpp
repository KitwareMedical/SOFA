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
#include "CudaPlaneForceField.inl"
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace gpu
{

namespace cuda
{

SOFA_DECL_CLASS(CudaPlaneForceField)

int PlaneForceFieldCudaClass = core::RegisterObject("Supports GPU-side computations using CUDA")
        .add< component::forcefield::PlaneForceField<CudaVec3fTypes> >()
        .add< component::forcefield::PlaneForceField<CudaVec3f1Types> >()
#ifdef SOFA_GPU_CUDA_DOUBLE
        .add< component::forcefield::PlaneForceField<CudaVec3dTypes> >()
        .add< component::forcefield::PlaneForceField<CudaVec3d1Types> >()
#endif // SOFA_GPU_CUDA_DOUBLE
        ;

} // namespace cuda

} // namespace gpu

} // namespace sofa
