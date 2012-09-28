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
#include "CudaLinearMovementConstraint.inl"

#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/ObjectFactory.h>
namespace sofa
{

// namespace component
// {
//
// namespace projectiveconstraintset
// {
//
// template class LinearMovementConstraint<gpu::cuda::CudaRigid3fTypes>;
//
// }// namespace projectiveconstraintset
//
// }// namespace component

namespace gpu
{

namespace cuda
{


SOFA_DECL_CLASS(CudaLinearMovementConstraint)

int LinearMovementConstraintCudaClass = core::RegisterObject("Supports GPU-side computations using CUDA")
// .add< component::projectiveconstraintset::LinearMovementConstraint<CudaVec3fTypes> >()
// .add< component::projectiveconstraintset::LinearMovementConstraint<CudaVec3f1Types> >()
        .add< component::projectiveconstraintset::LinearMovementConstraint<CudaVec6fTypes> >()
        .add< component::projectiveconstraintset::LinearMovementConstraint<CudaRigid3fTypes> >()
#ifdef SOFA_GPU_CUDA_DOUBLE
// .add< component::projectiveconstraintset::LinearMovementConstraint<CudaVec3dTypes> >()
// .add< component::projectiveconstraintset::LinearMovementConstraint<CudaVec3d1Types> >()
        .add< component::projectiveconstraintset::LinearMovementConstraint<CudaVec6dTypes> >()
        .add< component::projectiveconstraintset::LinearMovementConstraint<CudaRigid3dTypes> >()
#endif // SOFA_GPU_CUDA_DOUBLE
        ;

} // namespace cuda

} // namespace gpu

} // namespace sofa
