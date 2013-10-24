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
#ifndef SOFA_COMPONENT_FORCEFIELD_FASTTETRAHEDRALCOROTATIONALFORCEFIELD_INL
#define SOFA_COMPONENT_FORCEFIELD_FASTTETRAHEDRALCOROTATIONALFORCEFIELD_INL

#include <sofa/component/forcefield/FastTetrahedralCorotationalForceField.h>
#include <sofa/core/visual/VisualParams.h>
#include <fstream> // for reading the file
#include <iostream> //for debugging
#include <sofa/helper/gl/template.h>
#include <sofa/component/topology/TopologyData.inl>
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/helper/decompose.h>

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

typedef Tetra			        Tetrahedron;
typedef EdgesInTetrahedron		EdgesInTetrahedron;


template< class DataTypes>
void FastTetrahedralCorotationalForceField<DataTypes>::FTCFTetrahedronHandler::applyCreateFunction(unsigned int tetrahedronIndex,
        TetrahedronRestInformation &my_tinfo,
        const Tetrahedron &,
        const sofa::helper::vector<unsigned int> &,
        const sofa::helper::vector<double> &)
{
    if (ff)
    {
        const std::vector< Tetrahedron > &tetrahedronArray=ff->_topology->getTetrahedra() ;
        //		const std::vector< Edge> &edgeArray=ff->_topology->getEdges() ;
        unsigned int j,k,l,m,n;
        typename DataTypes::Real lambda=ff->getLambda();
        typename DataTypes::Real mu=ff->getMu();
        typename DataTypes::Real volume,val;
        typename DataTypes::Coord point[4]; //shapeVector[4];
        const typename DataTypes::VecCoord *restPosition=ff->mstate->getX0();

        ///describe the indices of the 4 tetrahedron vertices
        const Tetrahedron &t= tetrahedronArray[tetrahedronIndex];
//    BaseMeshTopology::EdgesInTetrahedron te=ff->_topology->getEdgesInTetrahedron(tetrahedronIndex);


        // store the point position
        for(j=0; j<4; ++j)
            point[j]=(*restPosition)[t[j]];
        /// compute 6 times the rest volume
        volume=dot(cross(point[1]-point[0],point[2]-point[0]),point[0]-point[3]);
        /// store the rest volume
        my_tinfo.restVolume=volume/6;
        mu*=fabs(volume)/6;
        lambda*=fabs(volume)/6;

        // store shape vectors at the rest configuration
        for(j=0; j<4; ++j)
        {
            if ((j%2)==0)
                my_tinfo.shapeVector[j]=cross(point[(j+2)%4] - point[(j+1)%4],point[(j+3)%4] - point[(j+1)%4])/volume;
            else
                my_tinfo.shapeVector[j]= -cross(point[(j+2)%4] - point[(j+1)%4],point[(j+3)%4] - point[(j+1)%4])/volume;
        }

        /// compute the edge stiffness of the linear elastic material
        for(j=0; j<6; ++j)
        {
            Edge e=ff->_topology->getLocalEdgesInTetrahedron(j);
            k=e[0];
            l=e[1];

            // store the rest edge vector
            my_tinfo.restEdgeVector[j]=point[l]-point[k];

            // the linear stiffness matrix using shape vectors and Lam� coefficients
            val=mu*dot(my_tinfo.shapeVector[l],my_tinfo.shapeVector[k]);
            for(m=0; m<3; ++m)
            {
                for(n=0; n<3; ++n)
                {
                    my_tinfo.linearDfDx[j][m][n]=lambda*my_tinfo.shapeVector[k][n]*my_tinfo.shapeVector[l][m]+
                            mu*my_tinfo.shapeVector[l][n]*my_tinfo.shapeVector[k][m];

                    if (m==n)
                    {
                        my_tinfo.linearDfDx[j][m][m]+=(Real)val;
                    }
                }
            }
            my_tinfo.transposedLinearDfDx[j]=my_tinfo.linearDfDx[j].transposed();
        }
        // compute the rotation matrix of the initial tetrahedron for the QR decomposition
        computeQRRotation(my_tinfo.restRotation,my_tinfo.restEdgeVector);
    }
}

template <class DataTypes> FastTetrahedralCorotationalForceField<DataTypes>::FastTetrahedralCorotationalForceField()
    : edgeInfo(initData(&edgeInfo, "edgeInfo", "Internal edge data"))
    , tetrahedronInfo(initData(&tetrahedronInfo, "tetrahedronInfo", "Internal tetrahedron data"))
    , _initialPoints(0)
    , updateMatrix(true)
    , f_method(initData(&f_method,std::string("large"),"method","\"large\" (by QR) or \"polar\" displacements"))
    , f_poissonRatio(initData(&f_poissonRatio,(Real)0.3,"poissonRatio","Poisson ratio in Hooke's law"))
    , f_youngModulus(initData(&f_youngModulus,(Real)1000.,"youngModulus","Young modulus in Hooke's law"))
    , lambda(0)
    , mu(0)
    , tetrahedronHandler(NULL)
{
    tetrahedronHandler = new FTCFTetrahedronHandler(this,&tetrahedronInfo);
}

template <class DataTypes> FastTetrahedralCorotationalForceField<DataTypes>::~FastTetrahedralCorotationalForceField()
{
    if (tetrahedronHandler) delete tetrahedronHandler;
}

template <class DataTypes> void FastTetrahedralCorotationalForceField<DataTypes>::init()
{
    //	serr << "initializing FastTetrahedralCorotationalForceField" << sendl;
    this->Inherited::init();

    _topology = this->getContext()->getMeshTopology();

    if (_topology->getNbTetrahedra()==0)
    {
        serr << "ERROR(FastTetrahedralCorotationalForceField): object must have a Tetrahedral Set Topology."<<sendl;
        return;
    }
    updateLameCoefficients();


    if (f_method.getValue() == "polar")
        decompositionMethod= POLAR_DECOMPOSITION;
    else if (f_method.getValue() == "large")
        decompositionMethod= QR_DECOMPOSITION;
    else
    {
        serr << "cannot recognize method "<< f_method.getValue() << ". Must be either large or polar" << sendl;
    }


    helper::vector<TetrahedronRestInformation>& tetrahedronInf = *(tetrahedronInfo.beginEdit());
    tetrahedronInf.resize(_topology->getNbTetrahedra());


    helper::vector<EdgeRestInformation>& edgeInf = *(edgeInfo.beginEdit());

    /// prepare to store info in the edge array
    edgeInf.resize(_topology->getNbEdges());
    edgeInfo.createTopologicalEngine(_topology);
    edgeInfo.registerTopologicalData();
    edgeInfo.endEdit();

    if (_initialPoints.size() == 0)
    {
        // get restPosition
        const VecCoord& p = *this->mstate->getX0();
        _initialPoints=p;
    }

    int i;

    /// initialize the data structure associated with each tetrahedron
    for (i=0; i<_topology->getNbTetrahedra(); ++i)
    {
        tetrahedronHandler->applyCreateFunction(i,tetrahedronInf[i],_topology->getTetrahedron(i),
                (const helper::vector< unsigned int > )0,
                (const helper::vector< double >)0);
    }
    /// set the call back function upon creation of a tetrahedron
    tetrahedronInfo.createTopologicalEngine(_topology,tetrahedronHandler);
    tetrahedronInfo.registerTopologicalData();
    tetrahedronInfo.endEdit();

    updateTopologyInfo=true;
}


template <class DataTypes>
void FastTetrahedralCorotationalForceField<DataTypes>::updateTopologyInformation()
{
    int i;
    unsigned int j;

    int nbEdges=_topology->getNbEdges();
    int nbTetrahedra=_topology->getNbTetrahedra();

    EdgeRestInformation *einfo;
    TetrahedronRestInformation *tetinfo;

    helper::vector<typename FastTetrahedralCorotationalForceField<DataTypes>::TetrahedronRestInformation>& tetrahedronInf = *(tetrahedronInfo.beginEdit());
    helper::vector<typename FastTetrahedralCorotationalForceField<DataTypes>::EdgeRestInformation>& edgeInf = *(edgeInfo.beginEdit());


    for(i=0; i<nbEdges; i++ )
    {
        einfo=&edgeInf[i];
        einfo->v[0]=_topology->getEdge(i)[0];
        einfo->v[1]=_topology->getEdge(i)[1];
    }

    for(i=0; i<nbTetrahedra; i++ )
    {
        tetinfo=&tetrahedronInf[i];
        /// describe the jth edge index of triangle no i
        const EdgesInTetrahedron &tea= _topology->getEdgesInTetrahedron(i);
        /// describe the jth vertex index of triangle no i
        const Tetrahedron &ta= _topology->getTetrahedron(i);

        for (j=0; j<4; ++j)
        {
            tetinfo->v[j]=ta[j];
        }
        for (j=0; j<6; ++j)
        {
            /// store the pointer to the local edge information
            tetinfo->edgeInfo[j]=& edgeInf[tea[j]];
            /// store the information about the orientation of the edge : 1 if the edge orientation matches the orientation in getLocalEdgesInTetrahedron
            /// ie edgesInTetrahedronArray[6][2] = {{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}};
            if (ta[ _topology->getLocalEdgesInTetrahedron(j)[0]]== _topology->getEdge(tea[j])[0])
                tetinfo->edgeOrientation[j]=1;
            else
                tetinfo->edgeOrientation[j]= -1;
        }

    }
    updateTopologyInfo=false;
    edgeInfo.endEdit();
    tetrahedronInfo.endEdit();
}
template<class DataTypes>
void FastTetrahedralCorotationalForceField<DataTypes>::computeQRRotation( Mat3x3 &r, const Coord *dp)
{
    // first vector on first edge
    // second vector in the plane of the two first edges
    // third vector orthogonal to first and second

    Coord edgex = dp[0];
    edgex.normalize();

    Coord edgey = dp[1];

    Coord edgez = cross( edgex, edgey );
    edgez.normalize();

    edgey = cross( edgez, edgex );
    edgey.normalize();

    r[0][0] = edgex[0];
    r[0][1] = edgex[1];
    r[0][2] = edgex[2];
    r[1][0] = edgey[0];
    r[1][1] = edgey[1];
    r[1][2] = edgey[2];
    r[2][0] = edgez[0];
    r[2][1] = edgez[1];
    r[2][2] = edgez[2];
}

template <class DataTypes>
void FastTetrahedralCorotationalForceField<DataTypes>::addForce(const sofa::core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecDeriv &  dataF, const DataVecCoord &  dataX , const DataVecDeriv & /*dataV*/ )
{

    VecDeriv& f        = *(dataF.beginEdit());
    const VecCoord& x  =   dataX.getValue()  ;


    unsigned int j,k,l;
    int nbTetrahedra=_topology->getNbTetrahedra();
    int i;
    const unsigned int edgesInTetrahedronArray[6][2] = {{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}};


    if (updateTopologyInfo)
    {
        updateTopologyInformation();
    }
    helper::vector<TetrahedronRestInformation>& tetrahedronInf = *(tetrahedronInfo.beginEdit());
    TetrahedronRestInformation *tetinfo;

    Coord dp[6],sv;
    Mat3x3 deformationGradient,S,R;

    for(i=0; i<nbTetrahedra; i++ )
    {
        tetinfo=&tetrahedronInf[i];

        for (j=0; j<6; ++j)
        {
            dp[j]=x[tetinfo->v[edgesInTetrahedronArray[j][1]]]-x[tetinfo->v[edgesInTetrahedronArray[j][0]]];
        }

        if (decompositionMethod==POLAR_DECOMPOSITION)
        {
            // compute the deformation gradient
            // deformation gradient = sum of tensor product between vertex position and shape vector
            // optimize by using displacement with first vertex
            sv=tetinfo->shapeVector[1];
            for (k=0; k<3; ++k)
            {
                for (l=0; l<3; ++l)
                {
                    deformationGradient[k][l]=dp[0][k]*sv[l];
                }
            }
            for (j=1; j<3; ++j)
            {
                sv=tetinfo->shapeVector[j+1];
                for (k=0; k<3; ++k)
                {
                    for (l=0; l<3; ++l)
                    {
                        deformationGradient[k][l]+=dp[j][k]*sv[l];
                    }
                }
            }
            // polar decomposition of the transformation
            helper::Decompose<Real>::polarDecomposition(deformationGradient,R);
        }
        else
        {

            /// perform QR decomposition
            computeQRRotation(S,dp);
            R=S.transposed()*tetinfo->restRotation;

        }

        // store transpose of rotation
        tetinfo->rotation=R.transposed();
        Coord force[4];


        for (j=0; j<6; ++j)
        {

            // displacement in the rest configuration
            dp[j]=tetinfo->rotation*dp[j]-tetinfo->restEdgeVector[j];

            // force on first vertex in the rest configuration
            force[edgesInTetrahedronArray[j][1]]+=tetinfo->linearDfDx[j]*dp[j];
            // force on second vertex in the rest configuration
            force[edgesInTetrahedronArray[j][0]]-=tetinfo->transposedLinearDfDx[j]*dp[j];
        }
        for (j=0; j<4; ++j)
        {
            f[tetinfo->v[j]]+=R*force[j];
        }


    }

    updateMatrix=true; // next time assemble the matrix
    tetrahedronInfo.endEdit();

    dataF.endEdit();

}


template <class DataTypes>
void FastTetrahedralCorotationalForceField<DataTypes>::addDForce(const sofa::core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv&   datadF , const DataVecDeriv&   datadX )
{
    VecDeriv& df       = *(datadF.beginEdit());
    const VecCoord& dx =   datadX.getValue()  ;
    Real kFactor = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());

    unsigned int j;
    int i;
    //	const std::vector< Edge> &edgeArray=_topology->getEdges() ;
    int nbEdges=_topology->getNbEdges();

    if (updateMatrix==true)
    {
        // the matrix must be stored in edges
        helper::vector<TetrahedronRestInformation>& tetrahedronInf = *(tetrahedronInfo.beginEdit());
        helper::vector<EdgeRestInformation>& edgeInf = *(edgeInfo.beginEdit());

        TetrahedronRestInformation *tetinfo;
        EdgeRestInformation *einfo;
        int nbTetrahedra=_topology->getNbTetrahedra();
        Mat3x3 tmp;

        updateMatrix=false;
        // reset all edge matrices
        for(einfo=&edgeInf[0],i=0; i<nbEdges; i++,einfo++ )
        {
            einfo->DfDx.clear();
            einfo->reverseDfDx.clear();
        }

        for(i=0; i<nbTetrahedra; i++ )
        {
            tetinfo=&tetrahedronInf[i];


            for (j=0; j<6; ++j)
            {

                einfo=tetinfo->edgeInfo[j];

                // test if the tetrahedron edge has the same orientation as the global edge
                tmp=tetinfo->linearDfDx[j]*tetinfo->rotation;

                if (tetinfo->edgeOrientation[j]==1)
                {
                    // store the two edge matrices since the stiffness matrix is not symmetric
                    einfo->DfDx+=tetinfo->rotation.transposed()*tmp;
                    einfo->reverseDfDx+= tmp.transposed()*tetinfo->rotation;

                }
                else
                {
                    einfo->reverseDfDx+=tetinfo->rotation.transposed()*tmp;
                    einfo->DfDx+= tmp.transposed()*tetinfo->rotation;

                }
            }
        }

        tetrahedronInfo.endEdit();
        edgeInfo.endEdit();
    }

    const helper::vector<EdgeRestInformation> &edgeInf= edgeInfo.getValue();
    unsigned int v0,v1;
    const EdgeRestInformation *einfo;
    Coord deltax;
    // use the already stored matrix
    for(i=0; i<nbEdges; i++ )
    {
        einfo=&edgeInf[i];

        v0=einfo->v[0];
        v1=einfo->v[1];

        deltax= dx[v1] -dx[v0];
        // einfo->reverseDfDx is not the transposed of einfo->DfDx
        df[v1]+= (einfo->DfDx*deltax) * kFactor;
        df[v0]-= (einfo->reverseDfDx*deltax) * kFactor;
    }

    datadF.endEdit();
}


template<class DataTypes>
void FastTetrahedralCorotationalForceField<DataTypes>::updateLameCoefficients()
{
    lambda= f_youngModulus.getValue()*f_poissonRatio.getValue()/((1-2*f_poissonRatio.getValue())*(1+f_poissonRatio.getValue()));
    mu = f_youngModulus.getValue()/(2*(1+f_poissonRatio.getValue()));

}


template<class DataTypes>
void FastTetrahedralCorotationalForceField<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
#ifndef SOFA_NO_OPENGL
    if (!vparams->displayFlags().getShowForceFields()) return;
    if (!this->mstate) return;

    if (vparams->displayFlags().getShowWireFrame())
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);


    if (vparams->displayFlags().getShowWireFrame())
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
#endif /* SOFA_NO_OPENGL */
}

} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FORCEFIELD_FASTTETRAHEDRALCOROTATIONALFORCEFIELD_INL
