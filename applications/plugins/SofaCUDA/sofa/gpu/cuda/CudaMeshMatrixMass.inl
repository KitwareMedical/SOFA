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
#ifndef SOFA_COMPONENT_MASS_CUDAMESHMATRIXMASS_INL
#define SOFA_COMPONENT_MASS_CUDAMESHMATRIXMASS_INL

#include "CudaMeshMatrixMass.h"
#include <sofa/component/mass/MeshMatrixMass.inl>
#include <sofa/gpu/cuda/CudaTypes.h>

namespace sofa
{
namespace gpu
{
namespace cuda
{

using namespace sofa::gpu::cuda;

extern "C"
{
    void MeshMatrixMassCuda_addMDx2f(unsigned int size, float factor, float massLumpingCoeff,const void * vertexMass, const void* dx, void* res);
    void MeshMatrixMassCuda_addForce2f(int dim, void * f, const void * vertexMass, const double * gravity, float massLumpingCoeff);
    void MeshMatrixMassCuda_accFromF2f(int dim, void * acc, const void * f,  const void * vertexMass, float massLumpingCoeff);
}
}// cuda
}// gpu

namespace component
{

namespace mass
{

using namespace sofa::gpu::cuda;


template<>
void MeshMatrixMass<CudaVec2fTypes, float>::copyVertexMass()
{
    vector<MassType>& vertexInf = *(vertexMassInfo.beginEdit());
    data.vMass.resize(_topology->getNbPoints());

    for (int i=0; i<this->_topology->getNbPoints(); ++i)
        data.vMass[i] = (float) vertexInf[i];

    vertexMassInfo.endEdit();
}


template<>
void MeshMatrixMass<CudaVec2fTypes, float>::addMDx(const core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecDeriv& d_f, const DataVecDeriv& d_dx, double d_factor)
{
    VecDeriv& f = *d_f.beginEdit();
    const VecDeriv& dx = d_dx.getValue();
    const CudaVector<float>& vertexMass = data.vMass;

    MeshMatrixMassCuda_addMDx2f(dx.size(),(float) d_factor, (float) massLumpingCoeff, vertexMass.deviceRead() , dx.deviceRead(), f.deviceWrite());
    d_f.endEdit();
}


template<>
void MeshMatrixMass<CudaVec2fTypes, float>::addForce(const core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecDeriv& d_f, const DataVecCoord& /* */, const DataVecDeriv& /* */)
{
    VecDeriv& f = *d_f.beginEdit();
    const CudaVector<float>& vertexMass = data.vMass;
    Vec3d g ( this->getContext()->getGravity() );

    MeshMatrixMassCuda_addForce2f( vertexMass.size(), f.deviceWrite(), vertexMass.deviceRead(), g.ptr(), (float) massLumpingCoeff);
    d_f.endEdit();
}


template<>
void MeshMatrixMass<CudaVec2fTypes, float>::accFromF(const core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecDeriv& a, const DataVecDeriv& f)
{
    VecDeriv& _acc = *a.beginEdit();
    const VecDeriv& _f = f.getValue();
    const CudaVector<float>& vertexMass = data.vMass;

    MeshMatrixMassCuda_accFromF2f( vertexMass.size(), _acc.deviceWrite(), _f.deviceRead(), vertexMass.deviceRead(), (float) massLumpingCoeff);
    a.endEdit();
}


} // namespace mass
} // namespace component
} // namespace sofa

#endif // SOFA_COMPONENT_MASS_CUDAMESHMATRIXMASS_INL
