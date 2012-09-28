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
#ifndef SOFA_GPU_CUDA_CUDAPARTICLESREPULSIONFORCEFIELD_INL
#define SOFA_GPU_CUDA_CUDAPARTICLESREPULSIONFORCEFIELD_INL

#include "CudaParticlesRepulsionForceField.h"
#include <sofa/component/forcefield/ParticlesRepulsionForceField.inl>
//#include <sofa/gpu/cuda/CudaSpatialGridContainer.inl>

namespace sofa
{

namespace gpu
{

namespace cuda
{

extern "C"
{

    void ParticlesRepulsionForceFieldCuda3f_addForce (unsigned int size, const void* cells, const void* cellGhost, GPURepulsion3f* repulsion, void* f, const void* x, const void* v );
    void ParticlesRepulsionForceFieldCuda3f_addDForce(unsigned int size, const void* cells, const void* cellGhost, GPURepulsion3f* repulsion, void* f, const void* x, const void* dx);

#ifdef SOFA_GPU_CUDA_DOUBLE

    void ParticlesRepulsionForceFieldCuda3d_addForce (unsigned int size, const void* cells, const void* cellGhost, GPURepulsion3d* repulsion, void* f, const void* x, const void* v );
    void ParticlesRepulsionForceFieldCuda3d_addDForce(unsigned int size, const void* cells, const void* cellGhost, GPURepulsion3d* repulsion, void* f, const void* x, const void* dx);

#endif // SOFA_GPU_CUDA_DOUBLE
}

} // namespace cuda

} // namespace gpu

namespace component
{

namespace forcefield
{

using namespace gpu::cuda;


template <>
void ParticlesRepulsionForceField<gpu::cuda::CudaVec3fTypes>::addForce(const core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecDeriv& d_f, const DataVecCoord& d_x, const DataVecDeriv& d_v)
{
    if (grid == NULL) return;

    VecDeriv& f = *d_f.beginEdit();
    const VecCoord& x = d_x.getValue();
    const VecDeriv& v = d_v.getValue();

    grid->updateGrid(x);
    GPURepulsion3f repulsion;
    repulsion.d = distance.getValue();
    repulsion.d2 = repulsion.d*repulsion.d;
    repulsion.stiffness = stiffness.getValue();
    repulsion.damping = damping.getValue();
    f.resize(x.size());
    Grid::Grid* g = grid->getGrid();
    ParticlesRepulsionForceFieldCuda3f_addForce(
        g->getNbCells(), g->getCellsVector().deviceRead(), g->getCellGhostVector().deviceRead(),
        &repulsion, f.deviceWrite(), x.deviceRead(), v.deviceRead());

    d_f.endEdit();
}

template <>
void ParticlesRepulsionForceField<gpu::cuda::CudaVec3fTypes>::addDForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& d_df, const DataVecDeriv& d_dx)
{
    if (grid == NULL) return;

    VecDeriv& df = *d_df.beginEdit();
    const VecDeriv& dx = d_dx.getValue();
    double kFactor = mparams->kFactor();
    double bFactor = mparams->bFactor();

    const VecCoord& x = *this->mstate->getX();
    GPURepulsion3f repulsion;
    repulsion.d = distance.getValue();
    repulsion.d2 = repulsion.d*repulsion.d;
    repulsion.stiffness = (float)(stiffness.getValue()*kFactor);
    repulsion.damping = (float)(damping.getValue()*bFactor);
    df.resize(dx.size());
    Grid::Grid* g = grid->getGrid();
    ParticlesRepulsionForceFieldCuda3f_addDForce(
        g->getNbCells(), g->getCellsVector().deviceRead(), g->getCellGhostVector().deviceRead(),
        &repulsion, df.deviceWrite(), x.deviceRead(), dx.deviceRead());

    d_df.endEdit();
}


#ifdef SOFA_GPU_CUDA_DOUBLE

template <>
void ParticlesRepulsionForceField<gpu::cuda::CudaVec3dTypes>::addForce(const core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecDeriv& d_f, const DataVecCoord& d_x, const DataVecDeriv& d_v)
{
    if (grid == NULL) return;

    VecDeriv& f = *d_f.beginEdit();
    const VecCoord& x = d_x.getValue();
    const VecDeriv& v = d_v.getValue();

    grid->updateGrid(x);
    GPURepulsion3d repulsion;
    repulsion.d = distance.getValue();
    repulsion.d2 = repulsion.d*repulsion.d;
    repulsion.stiffness = stiffness.getValue();
    repulsion.damping = damping.getValue();
    f.resize(x.size());
    Grid::Grid* g = grid->getGrid();
    ParticlesRepulsionForceFieldCuda3d_addForce(
        g->getNbCells(), g->getCellsVector().deviceRead(), g->getCellGhostVector().deviceRead(),
        &repulsion, f.deviceWrite(), x.deviceRead(), v.deviceRead());

    d_f.endEdit();
}

template <>
void ParticlesRepulsionForceField<gpu::cuda::CudaVec3dTypes>::addDForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& d_df, const DataVecDeriv& d_dx)
{
    if (grid == NULL) return;

    VecDeriv& df = *d_df.beginEdit();
    const VecDeriv& dx = d_dx.getValue();
    double kFactor = mparams->kFactor();
    double bFactor = mparams->bFactor();

    const VecCoord& x = *this->mstate->getX();
    GPURepulsion3d repulsion;
    repulsion.d = distance.getValue();
    repulsion.d2 = repulsion.d*repulsion.d;
    repulsion.stiffness = stiffness.getValue()*kFactor;
    repulsion.damping = damping.getValue()*bFactor;
    df.resize(dx.size());
    Grid::Grid* g = grid->getGrid();
    ParticlesRepulsionForceFieldCuda3d_addDForce(
        g->getNbCells(), g->getCellsVector().deviceRead(), g->getCellGhostVector().deviceRead(),
        &repulsion, df.deviceWrite(), x.deviceRead(), dx.deviceRead());

    d_df.endEdit();
}

#endif // SOFA_GPU_CUDA_DOUBLE


} // namespace forcefield

} // namespace component

} // namespace sofa

#endif
