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
#include <sofa/component/forcefield/LinearForceField.h>
#include <sofa/core/behavior/ForceField.inl>

namespace sofa
{

namespace gpu
{

namespace cuda
{

template<class DataTypes>
class CudaKernelsLinearForceField;

}// namespace cuda

}// namespace gpu

namespace component
{

namespace forcefield
{


template<int N, class real>
class LinearForceFieldInternalData< gpu::cuda::CudaRigidTypes<N, real> >
{
public:
    typedef LinearForceFieldInternalData< gpu::cuda::CudaRigidTypes<N, real> > Data;
    typedef gpu::cuda::CudaRigidTypes<N, real> DataTypes;
    typedef LinearForceField<DataTypes> Main;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
    typedef typename Main::SetIndex SetIndex;
    typedef typename Main::SetIndexArray SetIndexArray;

    typedef gpu::cuda::CudaKernelsLinearForceField<DataTypes> Kernels;

    // vector indices of concerned dofs
    gpu::cuda::CudaVector< int > indices;

    int size;

    static void init(Main* m);

    static void addForce(Main* m, VecDeriv& f);

};// LinearForceFieldInternalData< CudaRigidTypes >


template<>
void LinearForceField< gpu::cuda::CudaRigid3fTypes >::init();

template<>
void LinearForceField< gpu::cuda::CudaRigid3fTypes >::addForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v);

template<>
double LinearForceField< gpu::cuda::CudaRigid3fTypes >::getPotentialEnergy(const core::MechanicalParams* /* PARAMS FIRST */, const DataVecCoord& ) const;

#ifdef SOFA_GPU_CUDA_DOUBLE
template<>
void LinearForceField< gpu::cuda::CudaRigid3dTypes >::init();

template<>
void LinearForceField< gpu::cuda::CudaRigid3dTypes >::addForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v);

template<>
double LinearForceField< gpu::cuda::CudaRigid3dTypes >::getPotentialEnergy(const core::MechanicalParams* /* PARAMS FIRST */, const DataVecCoord& ) const;
#endif // SOFA_GPU_CUDA_DOUBLE

} // namespace forcefield

} // namespace component

} // namespace sofa
