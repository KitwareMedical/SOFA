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
#include "CudaFixedConstraint.inl"
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/Vec3Types.h>
#include <sofa/defaulttype/RigidTypes.h>

namespace sofa
{

namespace component
{

namespace projectiveconstraintset
{


template class FixedConstraint<gpu::cuda::CudaVec3fTypes>;
template class FixedConstraint<gpu::cuda::CudaVec3f1Types>;
template class FixedConstraint<gpu::cuda::CudaVec6fTypes>;
template class FixedConstraint<gpu::cuda::CudaRigid3fTypes>;
#ifdef SOFA_GPU_CUDA_DOUBLE
template class FixedConstraint<gpu::cuda::CudaVec3dTypes>;
template class FixedConstraint<gpu::cuda::CudaVec3d1Types>;
template class FixedConstraint<gpu::cuda::CudaVec6dTypes>;
template class FixedConstraint<gpu::cuda::CudaRigid3dTypes>;
#endif // SOFA_GPU_CUDA_DOUBLE

} // namespace projectiveconstraintset

} // namespace component

namespace gpu
{

namespace cuda
{


SOFA_DECL_CLASS(CudaFixedConstraint)

int FixedConstraintCudaClass = core::RegisterObject("Supports GPU-side computations using CUDA")
        .add< component::projectiveconstraintset::FixedConstraint<CudaVec3fTypes> >()
        .add< component::projectiveconstraintset::FixedConstraint<CudaVec3f1Types> >()
        .add< component::projectiveconstraintset::FixedConstraint<CudaVec6fTypes> >()
        .add< component::projectiveconstraintset::FixedConstraint<CudaRigid3fTypes> >()
#ifdef SOFA_GPU_CUDA_DOUBLE
        .add< component::projectiveconstraintset::FixedConstraint<CudaVec3dTypes> >()
        .add< component::projectiveconstraintset::FixedConstraint<CudaVec3d1Types> >()
        .add< component::projectiveconstraintset::FixedConstraint<CudaVec6dTypes> >()
        .add< component::projectiveconstraintset::FixedConstraint<CudaRigid3dTypes> >()
#endif // SOFA_GPU_CUDA_DOUBLE
        ;

} // namespace cuda

} // namespace gpu

} // namespace sofa
