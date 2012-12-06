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
#ifndef SOFA_COMPONENT_FORCEFIELD_TETRAHEDRALTENSORMASSFORCEFIELD_INL
#define SOFA_COMPONENT_FORCEFIELD_TETRAHEDRALTENSORMASSFORCEFIELD_INL

#include <sofa/component/forcefield/TetrahedralTensorMassForceField.h>
#include <sofa/core/visual/VisualParams.h>
#include <fstream> // for reading the file
#include <iostream> //for debugging
#include <sofa/helper/gl/template.h>
#include <sofa/component/topology/TopologyData.inl>
#include <sofa/core/behavior/ForceField.inl>

namespace sofa
{

namespace component
{

namespace forcefield
{

using namespace sofa::defaulttype;
using namespace	sofa::component::topology;
using namespace core::topology;

using core::topology::BaseMeshTopology;

typedef BaseMeshTopology::Tetra				Tetra;
typedef BaseMeshTopology::EdgesInTetrahedron		EdgesInTetrahedron;

typedef Tetra			Tetrahedron;
typedef EdgesInTetrahedron		EdgesInTetrahedron;


template< class DataTypes>
void TetrahedralTensorMassForceField<DataTypes>::TetrahedralTMEdgeHandler::applyCreateFunction(unsigned int, EdgeRestInformation &ei, const Edge &edge, const sofa::helper::vector<unsigned int> &, const sofa::helper::vector<double> &)
{
    if (ff)
    {

        unsigned int u,v;
        /// set to zero the stiffness matrix
        for (u=0; u<3; ++u)
        {
            for (v=0; v<3; ++v)
            {
                ei.DfDx[u][v]=0;
            }
        }

    }
    ei.vertices[0] =(float) edge[0];
    ei.vertices[1] =(float) edge[1];
}

template< class DataTypes>
void TetrahedralTensorMassForceField<DataTypes>::TetrahedralTMEdgeHandler::applyTetrahedronCreation(const sofa::helper::vector<unsigned int> &tetrahedronAdded,
        const sofa::helper::vector<Tetrahedron> &,
        const sofa::helper::vector<sofa::helper::vector<unsigned int> > &,
        const sofa::helper::vector<sofa::helper::vector<double> > &)
{
    if (ff)
    {

        unsigned int i,j,k,l,u,v;

        typename DataTypes::Real val1,volume;
        typename DataTypes::Real lambda=ff->getLambda();
        typename DataTypes::Real mu=ff->getMu();
        typename DataTypes::Real lambdastar, mustar;
        typename DataTypes::Coord point[4],shapeVector[4];

        const typename DataTypes::VecCoord *restPosition=ff->mstate->getX0();

        edgeRestInfoVector& edgeData = *(ff->edgeInfo.beginEdit());

        for (i=0; i<tetrahedronAdded.size(); ++i)
        {

            /// get a reference on the edge set of the ith added tetrahedron
            const EdgesInTetrahedron &te= ff->_topology->getEdgesInTetrahedron(tetrahedronAdded[i]);
            ///get a reference on the vertex set of the ith added tetrahedron
            const Tetrahedron &t= ff->_topology->getTetrahedron(tetrahedronAdded[i]);
            // store points
            for(j=0; j<4; ++j)
                point[j]=(*restPosition)[t[j]];
            /// compute 6 times the rest volume
            volume=dot(cross(point[1]-point[0],point[2]-point[0]),point[0]-point[3]);
            // store shape vectors
            for(j=0; j<4; ++j)
            {
                if ((j%2)==0)
                    shapeVector[j]=cross(point[(j+2)%4] - point[(j+1)%4],point[(j+3)%4] - point[(j+1)%4])/volume;
                else
                    shapeVector[j]= -cross(point[(j+2)%4] - point[(j+1)%4],point[(j+3)%4] - point[(j+1)%4])/volume;
            }

            lambdastar=lambda*fabs(volume)/6;
            mustar=mu*fabs(volume)/6;


            for(j=0; j<6; ++j)
            {
                /// local indices of the edge
                k = ff->_topology->getLocalEdgesInTetrahedron(j)[0];
                l = ff->_topology->getLocalEdgesInTetrahedron(j)[1];

                Mat3 &m=edgeData[te[j]].DfDx;

                val1= dot(shapeVector[k],shapeVector[l])*mustar;
                // print if obtuse tetrahedron along that edge
                if (ff->f_printLog.getValue())
                {
                    if (val1<0)
                        ff->serr<<"negative cotangent["<<tetrahedronAdded[i]<<"]["<<j<<"]"<<ff->sendl;
                }

                if (ff->_topology->getEdge(te[j])[0]!=t[l])
                {
                    for (u=0; u<3; ++u)
                    {
                        for (v=0; v<3; ++v)
                        {
                            m[u][v]+= lambdastar*shapeVector[l][u]*shapeVector[k][v]+mustar*shapeVector[k][u]*shapeVector[l][v];
                        }
                        m[u][u]+=val1;
                    }
                }
                else
                {
                    for (u=0; u<3; ++u)
                    {
                        for (v=0; v<3; ++v)
                        {
                            m[v][u]+= lambdastar*shapeVector[l][u]*shapeVector[k][v]+mustar*shapeVector[k][u]*shapeVector[l][v];
                        }
                        m[u][u]+=val1;
                    }
                }


            }
        }
        ff->edgeInfo.endEdit();
    }
}

template< class DataTypes>
void TetrahedralTensorMassForceField<DataTypes>::TetrahedralTMEdgeHandler::applyTetrahedronDestruction(const sofa::helper::vector<unsigned int> &tetrahedronRemoved)
{
    if (ff)
    {

        unsigned int i,j,k,l,u,v;

        typename DataTypes::Real val1,volume;
        typename DataTypes::Real lambda=ff->getLambda();
        typename DataTypes::Real mu=ff->getMu();
        typename DataTypes::Real lambdastar, mustar;
        typename DataTypes::Coord point[4],shapeVector[4];

        const typename DataTypes::VecCoord *restPosition=ff->mstate->getX0();

        edgeRestInfoVector& edgeData = *(ff->edgeInfo.beginEdit());

        for (i=0; i<tetrahedronRemoved.size(); ++i)
        {

            /// get a reference on the edge set of the ith added tetrahedron
            const EdgesInTetrahedron &te= ff->_topology->getEdgesInTetrahedron(tetrahedronRemoved[i]);
            ///get a reference on the vertex set of the ith added tetrahedron
            const Tetrahedron &t= ff->_topology->getTetrahedron(tetrahedronRemoved[i]);
            // store points
            for(j=0; j<4; ++j)
                point[j]=(*restPosition)[t[j]];
            /// compute 6 times the rest volume
            volume=dot(cross(point[1]-point[0],point[2]-point[0]),point[0]-point[3]);
            // store shape vectors
            for(j=0; j<4; ++j)
            {
                if ((j%2)==0)
                    shapeVector[j]=cross(point[(j+2)%4] - point[(j+1)%4],point[(j+3)%4] - point[(j+1)%4])/volume;
                else
                    shapeVector[j]= -cross(point[(j+2)%4] - point[(j+1)%4],point[(j+3)%4] - point[(j+1)%4])/volume;
            }

            lambdastar=lambda*fabs(volume)/6;
            mustar=mu*fabs(volume)/6;


            for(j=0; j<6; ++j)
            {
                /// local indices of the edge
                k = ff->_topology->getLocalEdgesInTetrahedron(j)[0];
                l = ff->_topology->getLocalEdgesInTetrahedron(j)[1];

                Mat3 &m=edgeData[te[j]].DfDx;

                val1= dot(shapeVector[k],shapeVector[l])*mustar;
                // print if obtuse tetrahedron along that edge
                if (ff->f_printLog.getValue())
                {
                    if (val1<0)
                        ff->serr<<"negative cotangent["<<tetrahedronRemoved[i]<<"]["<<j<<"]"<<ff->sendl;
                }

                if (ff->_topology->getEdge(te[j])[0]!=t[l])
                {
                    for (u=0; u<3; ++u)
                    {
                        for (v=0; v<3; ++v)
                        {
                            m[u][v]-= lambdastar*shapeVector[l][u]*shapeVector[k][v]+mustar*shapeVector[k][u]*shapeVector[l][v];
                        }
                        m[u][u]-=val1;
                    }
                }
                else
                {
                    for (u=0; u<3; ++u)
                    {
                        for (v=0; v<3; ++v)
                        {
                            m[v][u]-= lambdastar*shapeVector[l][u]*shapeVector[k][v]+mustar*shapeVector[k][u]*shapeVector[l][v];
                        }
                        m[u][u]-=val1;
                    }
                }


            }

        }
        ff->edgeInfo.endEdit();
    }
}


template <class DataTypes> TetrahedralTensorMassForceField<DataTypes>::TetrahedralTensorMassForceField()
    : _initialPoints(0)
    , updateMatrix(true)
    , f_poissonRatio(initData(&f_poissonRatio,(Real)0.3,"poissonRatio","Poisson ratio in Hooke's law"))
    , f_youngModulus(initData(&f_youngModulus,(Real)1000.,"youngModulus","Young modulus in Hooke's law"))
    , lambda(0)
    , mu(0)
	, atomicGPU(initData(&atomicGPU, (bool) false, "atomicGPU","True if the GPU can handle atomic operations (CUDA version > 2.0)"))
    , edgeInfo(initData(&edgeInfo, "edgeInfo", "Internal edge data"))
{
    edgeHandler = new TetrahedralTMEdgeHandler(this, &edgeInfo);
}

template <class DataTypes> TetrahedralTensorMassForceField<DataTypes>::~TetrahedralTensorMassForceField()
{
    if(edgeHandler) delete edgeHandler;
}

template <class DataTypes> void TetrahedralTensorMassForceField<DataTypes>::init()
{
    this->Inherited::init();

    _topology = this->getContext()->getMeshTopology();

    edgeInfo.createTopologicalEngine(_topology,edgeHandler);
    edgeInfo.linkToTetrahedronDataArray();
    edgeInfo.registerTopologicalData();


    if (_topology->getNbTetrahedra()==0)
    {
        serr << "ERROR(TetrahedralTensorMassForceField): object must have a Tetrahedral Set Topology."<<sendl;
        return;
    }
    updateLameCoefficients();

    edgeRestInfoVector& edgeInf = *(edgeInfo.beginEdit());


    /// prepare to store info in the edge array
    edgeInf.resize(_topology->getNbEdges());


    if (_initialPoints.size() == 0)
    {
        // get restPosition
        const VecCoord& p = *this->mstate->getX0();
        _initialPoints=p;
    }

    int i;
    // set edge tensor to 0
    for (i=0; i<_topology->getNbEdges(); ++i)
    {
        edgeHandler->applyCreateFunction(i, edgeInf[i],
                _topology->getEdge(i),  (const sofa::helper::vector< unsigned int > )0,
                (const sofa::helper::vector< double >)0);
    }
    // create edge tensor by calling the tetrahedron creation function
    sofa::helper::vector<unsigned int> tetrahedronAdded;
    for (i=0; i<_topology->getNbTetrahedra(); ++i)
        tetrahedronAdded.push_back(i);

    edgeHandler->applyTetrahedronCreation(tetrahedronAdded,
            (const sofa::helper::vector<Tetrahedron>)0,
            (const sofa::helper::vector<sofa::helper::vector<unsigned int> >)0,
            (const sofa::helper::vector<sofa::helper::vector<double> >)0);

    edgeInfo.endEdit();

    /// FOR CUDA
    /// Save the neighbourhood for points (in case of CudaTypes, non atomic)
    if(!atomicGPU.getValue())
    {
        this->initNeighbourhoodPoints();
    }
}

template <class DataTypes>
void TetrahedralTensorMassForceField<DataTypes>::initNeighbourhoodPoints() {}

template <class DataTypes>
void TetrahedralTensorMassForceField<DataTypes>::addForce(const core::MechanicalParams* /* mparams */ /* PARAMS FIRST */, DataVecDeriv& d_f, const DataVecCoord& d_x, const DataVecDeriv& /* d_v */)
{
    VecDeriv& f = *d_f.beginEdit();
    const VecCoord& x = d_x.getValue();

    unsigned int v0,v1;
    int nbEdges=_topology->getNbEdges();

    EdgeRestInformation *einfo;

    edgeRestInfoVector& edgeInf = *(edgeInfo.beginEdit());
    edgeRestInfoVector cudaEdgeInfo;
    cudaEdgeInfo.resize(edgeInf.size());


    for(unsigned int i=0; i<edgeInf.size(); i++)
        cudaEdgeInfo[i] = edgeInf[i];

    Deriv force;
    Coord dp0,dp1,dp;

    for(int i=0; i<nbEdges; i++ )
    {
        einfo=&edgeInf[i];
        v0=_topology->getEdge(i)[0];
        v1=_topology->getEdge(i)[1];
        dp0=x[v0]-_initialPoints[v0];
        dp1=x[v1]-_initialPoints[v1];
        dp = dp1-dp0;

        f[v1]+=einfo->DfDx*dp;
        f[v0]-=einfo->DfDx.transposeMultiply(dp);
    }

    edgeInfo.endEdit();
    d_f.endEdit();
}


template <class DataTypes>
void TetrahedralTensorMassForceField<DataTypes>::addDForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& d_df, const DataVecDeriv& d_dx)
{
    VecDeriv& df = *d_df.beginEdit();
    const VecDeriv& dx = d_dx.getValue();
    double kFactor = mparams->kFactor();

    unsigned int v0,v1;
    int nbEdges=_topology->getNbEdges();

    EdgeRestInformation *einfo;

    edgeRestInfoVector& edgeInf = *(edgeInfo.beginEdit());

    Deriv force;
    Coord dp0,dp1,dp;

    for(int i=0; i<nbEdges; i++ )
    {
        einfo=&edgeInf[i];
        v0=_topology->getEdge(i)[0];
        v1=_topology->getEdge(i)[1];
        dp0=dx[v0];
        dp1=dx[v1];
        dp = dp1-dp0;

        df[v1]+= (einfo->DfDx*dp) * kFactor;
        df[v0]-= (einfo->DfDx.transposeMultiply(dp)) * kFactor;
    }
    edgeInfo.endEdit();

    d_df.endEdit();
}


template<class DataTypes>
void TetrahedralTensorMassForceField<DataTypes>::updateLameCoefficients()
{
    lambda= f_youngModulus.getValue()*f_poissonRatio.getValue()/((1-2*f_poissonRatio.getValue())*(1+f_poissonRatio.getValue()));
    mu = f_youngModulus.getValue()/(2*(1+f_poissonRatio.getValue()));
//	serr << "initialized Lame coef : lambda=" <<lambda<< " mu="<<mu<<sendl;
}


template<class DataTypes>
void TetrahedralTensorMassForceField<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    if (!vparams->displayFlags().getShowForceFields()) return;
    if (!this->mstate) return;

    if (vparams->displayFlags().getShowWireFrame())
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

// 	VecCoord& x = *this->mstate->getX();
// 	int nbTriangles=_topology->getNbTriangles();

    /*
    	glDisable(GL_LIGHTING);

    	glBegin(GL_TRIANGLES);
    	for(i=0;i<nbTriangles; ++i)
    	{
    		int a = _topology->getTriangle(i)[0];
    		int b = _topology->getTriangle(i)[1];
    		int c = _topology->getTriangle(i)[2];

    		glColor4f(0,1,0,1);
    		helper::gl::glVertexT(x[a]);
    		glColor4f(0,0.5,0.5,1);
    		helper::gl::glVertexT(x[b]);
    		glColor4f(0,0,1,1);
    		helper::gl::glVertexT(x[c]);
    	}
    	glEnd();

    */
    if (vparams->displayFlags().getShowWireFrame())
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FORCEFIELD_TETRAHEDRALTENSORMASSFORCEFIELD_INL
