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
#ifndef SOFA_GPU_CUDA_CUDATRIANGULARFEMFORCEFIELDOPTIM_INL
#define SOFA_GPU_CUDA_CUDATRIANGULARFEMFORCEFIELDOPTIM_INL

#include "CudaTriangularFEMForceFieldOptim.h"
#include <sofa/component/forcefield/TriangularFEMForceFieldOptim.inl>

namespace sofa
{

namespace gpu
{

namespace cuda
{

extern "C"
{
void TriangularFEMForceFieldOptimCuda3f_addForce(unsigned int size, void* f, const void* x, const void* v,
    void* triangleState, const void* triangleInfo,
    unsigned int nbTriangles,
    const void* gpuTriangleInfo,
    float gamma, float mu);

void TriangularFEMForceFieldOptimCuda3f_addDForce(unsigned int size, void* f, const void* dx, float kFactor,
    const void* triangleState, const void* triangleInfo,
    unsigned int nbTriangles,
    const void* gpuTriangleInfo,
    float gamma, float mu); //, const void* dfdx);
}

} // namespace cuda

} // namespace gpu

namespace component
{

namespace forcefield
{

using namespace gpu::cuda;

template <>
void TriangularFEMForceFieldOptim<gpu::cuda::CudaVec3fTypes>::addForce(const core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecDeriv& d_f, const DataVecCoord& d_x, const DataVecDeriv& d_v)
{
    VecDeriv& f = *d_f.beginEdit();
    const VecCoord& x = d_x.getValue();
    const VecDeriv& v = d_v.getValue();

    VecTriangleState& triState = *triangleState.beginEdit();
    const VecTriangleInfo& triInfo = triangleInfo.getValue();
    const unsigned int nbTriangles = _topology->getNbTriangles();
    const InternalData::VecGPUTriangleInfo& gpuTriangleInfo = data.gpuTriangleInfo;
    const Real gamma = this->gamma;
    const Real mu = this->mu;

    f.resize(x.size());

    TriangularFEMForceFieldOptimCuda3f_addForce(x.size(), f.deviceWrite(), x.deviceRead(), v.deviceRead(),
        triState.deviceWrite(), 
        triInfo.deviceRead(), 
        nbTriangles,
        gpuTriangleInfo.deviceRead(),
        gamma, mu);
    
    triangleState.endEdit();
    d_f.endEdit();
}

template <>
void TriangularFEMForceFieldOptim<gpu::cuda::CudaVec3fTypes>::addDForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& d_df, const DataVecDeriv& d_dx)
{
    VecDeriv& df = *d_df.beginEdit();
    const VecDeriv& dx = d_dx.getValue();
    const Real kFactor = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());

    const VecTriangleState& triState = triangleState.getValue();
    const VecTriangleInfo& triInfo = triangleInfo.getValue();
    const unsigned int nbTriangles = _topology->getNbTriangles();
    const InternalData::VecGPUTriangleInfo& gpuTriangleInfo = data.gpuTriangleInfo;
    const Real gamma = this->gamma;
    const Real mu = this->mu;

    df.resize(dx.size());

    TriangularFEMForceFieldOptimCuda3f_addDForce(dx.size(), df.deviceWrite(), dx.deviceRead(), kFactor,
        triState.deviceRead(), 
        triInfo.deviceRead(), 
        nbTriangles,
        gpuTriangleInfo.deviceRead(),
        gamma, mu);
    
    d_df.endEdit();
}

} // namespace forcefield

} // namespace component

} // namespace sofa

#endif
