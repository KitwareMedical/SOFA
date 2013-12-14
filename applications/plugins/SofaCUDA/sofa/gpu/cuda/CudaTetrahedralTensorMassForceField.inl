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
#ifndef SOFA_GPU_CUDA_CUDATETRAHEDRALTENSORMASSFORCEFIELD_INL
#define SOFA_GPU_CUDA_CUDATETRAHEDRALTENSORMASSFORCEFIELD_INL

#include <sofa/gpu/cuda/CudaTetrahedralTensorMassForceField.h>
#include <sofa/component/forcefield/TetrahedralTensorMassForceField.inl>

namespace sofa
{
namespace gpu
{
namespace cuda
{
    extern "C"
    {
        void TetrahedralTensorMassForceFieldCuda3f_addForce(int nbPoints, int nbMaxEdgesPerNode, const void* neighbourhoodPoints, void* contribEdge, int nbEdges, void* f, const void* x, const void* initialPoints, const void* edgeInfo );
        void TetrahedralTensorMassForceFieldCuda3f_addDForce(int nbPoints, int nbMaxEdgesPerNode, const void* neighbourhoodPoints, void* contribEdge, int nbEdges, void* df, const void* dx, const void* edgeInfo, float kFactor );
#ifdef SOFA_GPU_CUDA_DOUBLE
        void TetrahedralTensorMassForceFieldCuda3d_addForce(int nbPoints, int nbMaxEdgesPerNode, const void* neighbourhoodPoints, void* contribEdge, int nbEdges, void* f, const void* x, const void* initialPoints, const void* edgeInfo );
        void TetrahedralTensorMassForceFieldCuda3d_addDForce(int nbPoints, int nbMaxEdgesPerNode, const void* neighbourhoodPoints, void* contribEdge, int nbEdges, void* df, const void* dx, const void* edgeInfo, double kFactor );
#endif
    }
} // namespace cuda
} // namespace gpu

namespace component
{

namespace forcefield
{

using namespace gpu::cuda;

    template <>
    void TetrahedralTensorMassForceField<gpu::cuda::CudaVec3fTypes>::addForce(const core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecDeriv& d_f, const DataVecCoord& d_x, const DataVecDeriv& /*d_v*/)
    {
		sofa::helper::AdvancedTimer::stepBegin("addForceTetraTensorMass");

        VecDeriv& f = *d_f.beginEdit();
        const VecCoord& x = d_x.getValue();

        int nbEdges=_topology->getNbEdges();
        int nbPoints=_topology->getNbPoints();

        edgeRestInfoVector& edgeInf = *(edgeInfo.beginEdit());

        TetrahedralTensorMassForceField_contribEdge.resize(6*nbEdges);
        TetrahedralTensorMassForceFieldCuda3f_addForce(nbPoints, TetrahedralTensorMassForceField_nbMaxEdgesPerNode, TetrahedralTensorMassForceField_neighbourhoodPoints.deviceRead(), TetrahedralTensorMassForceField_contribEdge.deviceWrite(), nbEdges,  f.deviceWrite(), x.deviceRead(), _initialPoints.deviceRead(), edgeInf.deviceRead());

        edgeInfo.endEdit();
        d_f.endEdit();
		sofa::helper::AdvancedTimer::stepEnd("addForceTetraTensorMass");
    }

    template <>
    void TetrahedralTensorMassForceField<gpu::cuda::CudaVec3fTypes>::addDForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& d_df, const DataVecDeriv& d_dx)
    {
		sofa::helper::AdvancedTimer::stepBegin("addDForceTetraTensorMass");

        VecDeriv& df = *d_df.beginEdit();
        const VecDeriv& dx = d_dx.getValue();
        Real kFactor = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());

        int nbEdges=_topology->getNbEdges();
        int nbPoints=_topology->getNbPoints();
        edgeRestInfoVector& edgeInf = *(edgeInfo.beginEdit());

        TetrahedralTensorMassForceField_contribEdge.resize(6*nbEdges);
        TetrahedralTensorMassForceFieldCuda3f_addDForce(nbPoints, TetrahedralTensorMassForceField_nbMaxEdgesPerNode, TetrahedralTensorMassForceField_neighbourhoodPoints.deviceRead(), TetrahedralTensorMassForceField_contribEdge.deviceWrite(), nbEdges,  df.deviceWrite(), dx.deviceRead(), edgeInf.deviceRead(), (float)kFactor);

        edgeInfo.endEdit();
        d_df.endEdit();

        sofa::helper::AdvancedTimer::stepEnd("addDForceTetraTensorMass");
    }


    template<>
    void TetrahedralTensorMassForceField<CudaVec3fTypes>::initNeighbourhoodPoints()
    {
        std::cout<<"(TetrahedralTensorMassForceField) GPU-GEMS activated"<<std::endl;

        /// Initialize the number max of edges per node
        TetrahedralTensorMassForceField_nbMaxEdgesPerNode = 0;

        /// Compute it
        for(int i=0;i<_topology->getNbPoints();++i)
        {
            if((int)_topology->getEdgesAroundVertex(i).size()>TetrahedralTensorMassForceField_nbMaxEdgesPerNode)
                TetrahedralTensorMassForceField_nbMaxEdgesPerNode = _topology->getEdgesAroundVertex(i).size();
        }

        /// Initialize the vector neighbourhoodPoints
        TetrahedralTensorMassForceField_neighbourhoodPoints.resize((_topology->getNbPoints())*TetrahedralTensorMassForceField_nbMaxEdgesPerNode);

        unsigned int edgeID;

        for (int i=0;i<_topology->getNbPoints();++i)
        {
            for(int j=0;j<TetrahedralTensorMassForceField_nbMaxEdgesPerNode;++j)
            {
                if(j>(int)_topology->getEdgesAroundVertex(i).size()-1)
                    TetrahedralTensorMassForceField_neighbourhoodPoints[i*TetrahedralTensorMassForceField_nbMaxEdgesPerNode+j] = -1;
                else
                {
                    edgeID = _topology->getEdgesAroundVertex(i)[j];
                    if(i == (int)_topology->getEdge(edgeID)[0])
                        TetrahedralTensorMassForceField_neighbourhoodPoints[i*TetrahedralTensorMassForceField_nbMaxEdgesPerNode+j] = 2*edgeID;   //v0
                    else
                        TetrahedralTensorMassForceField_neighbourhoodPoints[i*TetrahedralTensorMassForceField_nbMaxEdgesPerNode+j] = 2*edgeID+1; //v1
                }
            }
        }


    }


#ifdef SOFA_GPU_CUDA_DOUBLE
    template <>
    void TetrahedralTensorMassForceField<gpu::cuda::CudaVec3dTypes>::addForce(const core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecDeriv& d_f, const DataVecCoord& d_x, const DataVecDeriv& /*d_v*/)
    {
        VecDeriv& f = *d_f.beginEdit();
        const VecCoord& x = d_x.getValue();

        int nbEdges=_topology->getNbEdges();
        int nbPoints=_topology->getNbPoints();

        edgeRestInfoVector& edgeInf = *(edgeInfo.beginEdit());

        TetrahedralTensorMassForceField_contribEdge.resize(6*nbEdges);
        TetrahedralTensorMassForceFieldCuda3d_addForce(nbPoints, TetrahedralTensorMassForceField_nbMaxEdgesPerNode, TetrahedralTensorMassForceField_neighbourhoodPoints.deviceRead(), TetrahedralTensorMassForceField_contribEdge.deviceWrite(), nbEdges,  f.deviceWrite(), x.deviceRead(), _initialPoints.deviceRead(), edgeInf.deviceRead());

        edgeInfo.endEdit();
        d_f.endEdit();


    }

    template <>
    void TetrahedralTensorMassForceField<gpu::cuda::CudaVec3dTypes>::addDForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& d_df, const DataVecDeriv& d_dx)
    {
        VecDeriv& df = *d_df.beginEdit();
        const VecDeriv& dx = d_dx.getValue();
        Real kFactor = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());

        int nbEdges=_topology->getNbEdges();
        int nbPoints=_topology->getNbPoints();
        edgeRestInfoVector& edgeInf = *(edgeInfo.beginEdit());

        TetrahedralTensorMassForceField_contribEdge.resize(6*nbEdges);
        TetrahedralTensorMassForceFieldCuda3d_addDForce(nbPoints, TetrahedralTensorMassForceField_nbMaxEdgesPerNode, TetrahedralTensorMassForceField_neighbourhoodPoints.deviceRead(), TetrahedralTensorMassForceField_contribEdge.deviceWrite(), nbEdges,  df.deviceWrite(), dx.deviceRead(), edgeInf.deviceRead(), kFactor);

        edgeInfo.endEdit();
        d_df.endEdit();
    }

	template<>
	void TetrahedralTensorMassForceField<CudaVec3dTypes>::initNeighbourhoodPoints()
	{
		std::cout<<"(TetrahedralTensorMassForceField) GPU-GEMS activated"<<std::endl;

		/// Initialize the number max of edges per node
		TetrahedralTensorMassForceField_nbMaxEdgesPerNode = 0;

		/// Compute it
		for(int i=0;i<_topology->getNbPoints();++i)
		{
			if((int)_topology->getEdgesAroundVertex(i).size()>TetrahedralTensorMassForceField_nbMaxEdgesPerNode)
				TetrahedralTensorMassForceField_nbMaxEdgesPerNode = _topology->getEdgesAroundVertex(i).size();
		}

		/// Initialize the vector neighbourhoodPoints
		TetrahedralTensorMassForceField_neighbourhoodPoints.resize((_topology->getNbPoints())*TetrahedralTensorMassForceField_nbMaxEdgesPerNode);

		unsigned int edgeID;

        for (int i=0;i<_topology->getNbPoints();++i)
		{
			for(int j=0;j<TetrahedralTensorMassForceField_nbMaxEdgesPerNode;++j)
			{
				if(j>(int)_topology->getEdgesAroundVertex(i).size()-1)
					TetrahedralTensorMassForceField_neighbourhoodPoints[i*TetrahedralTensorMassForceField_nbMaxEdgesPerNode+j] = -1;
				else
				{
					edgeID = _topology->getEdgesAroundVertex(i)[j];
                    if((unsigned) i == _topology->getEdge(edgeID)[0])
						TetrahedralTensorMassForceField_neighbourhoodPoints[i*TetrahedralTensorMassForceField_nbMaxEdgesPerNode+j] = 2*edgeID;   //v0
					else
						TetrahedralTensorMassForceField_neighbourhoodPoints[i*TetrahedralTensorMassForceField_nbMaxEdgesPerNode+j] = 2*edgeID+1; //v1
				}
			}
		}


	}

#endif


} // namespace forcefield
} // namespace component
} // namespace sofa
#endif //SOFA_GPU_CUDA_CUDATETRAHEDRALTENSORMASSFORCEFIELD_INL
