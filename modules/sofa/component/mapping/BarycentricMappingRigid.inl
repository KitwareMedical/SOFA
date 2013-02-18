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
#ifndef SOFA_COMPONENT_MAPPING_BARYCENTRICMAPPINGRIGID_INL
#define SOFA_COMPONENT_MAPPING_BARYCENTRICMAPPINGRIGID_INL

#include <sofa/component/mapping/BarycentricMappingRigid.h>
#include <sofa/component/mapping/BarycentricMapping.inl>

#include <sofa/helper/decompose.inl>

#ifdef SOFA_IP_TRACES
#include <sofa/component/mapping/_IP_MapTraceMacros.h>
#endif

namespace sofa
{

namespace component
{

namespace mapping
{


template <class In, class Out>
void BarycentricMapperTetrahedronSetTopologyRigid<In,Out>::clear ( int reserve )
{
    helper::vector<MappingData>& vectorData = *(map.beginEdit());
    vectorData.clear(); if ( reserve>0 ) vectorData.reserve ( reserve );
    map.endEdit();
}

template <class In, class Out>
int BarycentricMapperTetrahedronSetTopologyRigid<In,Out>::addPointInTetra ( const int tetraIndex, const SReal* baryCoords )
{
    helper::vector<MappingData>& vectorData = *(map.beginEdit());
    vectorData.resize ( map.getValue().size() +1 );
    MappingData& data = *vectorData.rbegin();
    map.endEdit();
    data.in_index = tetraIndex;
    data.baryCoords[0] = ( Real ) baryCoords[0];
    data.baryCoords[1] = ( Real ) baryCoords[1];
    data.baryCoords[2] = ( Real ) baryCoords[2];
    return map.getValue().size()-1;
}

template<class In, class Out>
int BarycentricMapperTetrahedronSetTopologyRigid<In,Out>::addPointOrientationInTetra( const int tetraIndex, const Matrix3 baryCoorsOrient )
{
    //storing the frame in 3 maps: one direction vector in one map  (3 coor. elements inside a map)
    // IPTR_BARCPP_ADDOR("addPointOrientation BEGIN" << endl);
    helper::vector<MappingOrientData>& vectorData = *(mapOrient.beginEdit());
    vectorData.resize ( vectorData.size() +1 );
    MappingOrientData& data = *vectorData.rbegin();
    for (unsigned int dir = 0; dir < 3; dir++)
    {
        data[dir].in_index = tetraIndex;
        // IPTR_BARCPP_ADDOR("baryCoords of vector["<<dir<<"]: ");
        for (unsigned int coor = 0; coor < 3; coor++)
        {
            data[dir].baryCoords[coor] = ( Real ) baryCoorsOrient[coor][dir];
            //IPNTR_BARCPP_ADDOR(data[dir].baryCoords[coor] << " ");
        }
        //IPNTR_BARCPP_ADDOR(endl);

    }
    mapOrient.endEdit();
    // IPTR_BARCPP_ADDOR("addPointOrientation END" << endl);
    return map.getValue().size()-1;
}


template<class In, class Out>
void BarycentricMapperTetrahedronSetTopologyRigid<In,Out>::init(const typename Out::VecCoord& out, const typename In::VecCoord& in)
{
#ifdef SOFA_IP_TRACES
    IPTR_BARCPP_INIT("BarycentricMapperTetrahedronSetTopology::init BEGIN " << endl);
    IPTR_BARCPP_INIT("out size = " << out.size() << endl);
#endif

    _fromContainer->getContext()->get ( _fromGeomAlgo );

    int outside = 0;
    const sofa::helper::vector<topology::Tetrahedron>& tetrahedra = this->fromTopology->getTetrahedra();

    sofa::helper::vector<Matrix3> bases;
    sofa::helper::vector<Vector3> centers;

    clear ( out.size() );
    bases.resize ( tetrahedra.size() );
    centers.resize ( tetrahedra.size() );
    for ( unsigned int t = 0; t < tetrahedra.size(); t++ )
    {
        Mat3x3d m,mt;
        m[0] = in[tetrahedra[t][1]]-in[tetrahedra[t][0]];
        m[1] = in[tetrahedra[t][2]]-in[tetrahedra[t][0]];
        m[2] = in[tetrahedra[t][3]]-in[tetrahedra[t][0]];
        mt.transpose ( m );
        bases[t].invert ( mt );
        centers[t] = ( in[tetrahedra[t][0]]+in[tetrahedra[t][1]]+in[tetrahedra[t][2]]+in[tetrahedra[t][3]] ) *0.25;
    }

    //find the closest tetra for given points of beam
    for ( unsigned int i=0; i<out.size(); i++ )
    {
        Vector3 pos = out[i].getCenter(); // Rigid3dTypes::GetCPos(out[i]); // pos = defaulttype::Rigid3dType::getCPos(out[i]);

        //associate the point with the tetrahedron, point in Barycentric coors wrt. the closest tetra, store to an associated structure
        Vector3 coefs;
        int index = -1;
        double distance = 1e10;
        for ( unsigned int t = 0; t < tetrahedra.size(); t++ )
        {
            Vec3d v = bases[t] * ( pos - in[tetrahedra[t][0]] );
            double d = std::max ( std::max ( -v[0],-v[1] ),std::max ( -v[2],v[0]+v[1]+v[2]-1 ) );
            if ( d>0 ) d = ( pos-centers[t] ).norm2();
            if ( d<distance )
            {
                coefs = v; distance = d; index = t;
            }
        }
        if ( distance>0 ) ++outside;

        //convert the orientation to basis given by closest tetrahedron
        Quat quatA = out[i].getOrientation();
        //initRigidOrientation[i]=quatA;
        Matrix3 orientationMatrix, orientationMatrixBary;
        quatA.toMatrix(orientationMatrix);    //direction vectors stored as columns
        orientationMatrix.transpose(); //to simplify the multiplication below
        for (unsigned c=0; c < 3; c++)
        {
            orientationMatrixBary[c]=bases[index]*orientationMatrix[c];
            orientationMatrixBary[c].normalize();
        }
        orientationMatrixBary.transpose();  //to get the directions as columns

        //store the point and orientation in local coordinates of tetrahedra frame
        addPointInTetra ( index, coefs.ptr() );
        addPointOrientationInTetra(index, orientationMatrixBary);
    }
#ifdef SOFA_IP_TRACES
    IPTR_BARCPP_INIT("BarycentricMapperTetrahedronSetTopology::init END" << endl);
#endif
}


template<class In, class Out>
void BarycentricMapperTetrahedronSetTopologyRigid<In,Out>::apply( typename Out::VecCoord& out, const typename In::VecCoord& in )
{
#ifdef SOFA_IP_TRACES
    IPTR_BARCPP_APPLY( "BarycentricMapperTetrahedronSetTopology<SPEC>::apply BEGIN" << endl);
#endif

    actualTetraPosition=in;
    //get number of point being mapped
    unsigned int nbPoints = map.getValue().size();
    out.resize (nbPoints);
    const sofa::helper::vector<topology::Tetrahedron>& tetrahedra = this->fromTopology->getTetrahedra();
    const sofa::helper::vector<MappingData >& map = this->map.getValue();
    const sofa::helper::vector<MappingOrientData >& mapOrient = this->mapOrient.getValue();

    typename In::VecCoord inCopy = in;

    for ( unsigned int i=0; i<map.size(); i++ )
    {
        //get barycentric coors for given point (node of a beam in this case)
        const Real fx = map[i].baryCoords[0];
        const Real fy = map[i].baryCoords[1];
        const Real fz = map[i].baryCoords[2];
        int index = map[i].in_index;
        const topology::Tetrahedron& tetra = tetrahedra[index];

        Vector3 rotatedPosition= in[tetra[0]] * ( 1-fx-fy-fz ) + in[tetra[1]] * fx + in[tetra[2]] * fy + in[tetra[3]] * fz ;
        Out::setCPos(out[i] , rotatedPosition); // glPointPositions[i] );
    }

    //sofa::helper::vector<Vector3> vectors
    sofa::helper::vector< sofa::defaulttype::Mat<12,3> > rotJ;
    rotJ.resize(map.size());
    //point running over each DoF (assoc. with frame) in the out vector; get it from the mapOrient[0]
    for (unsigned int point = 0; point < mapOrient.size(); point++)
    {
        int index = mapOrient[point][0].in_index;
        const topology::Tetrahedron& tetra = tetrahedra[index];

        //compute the rotation of the rigid point using the "basis" approach
        Matrix3 orientationMatrix, polarMatrixQ; // orthogMatrix
        Matrix3 m,basis;
        m[0] = in[tetra[1]]-in[tetra[0]];
        m[1] = in[tetra[2]]-in[tetra[0]];
        m[2] = in[tetra[3]]-in[tetra[0]];
        basis.transpose ( m );

        for (unsigned int dir = 0; dir < 3; dir++)   //go through the three maps
        {
            Vector3 inGlobal;
            inGlobal[0] = mapOrient[point][dir].baryCoords[0];
            inGlobal[1] = mapOrient[point][dir].baryCoords[1];
            inGlobal[2] = mapOrient[point][dir].baryCoords[2];

            orientationMatrix[dir]= basis*inGlobal;
        }

        orientationMatrix.transpose();
        helper::Decompose<Matrix3::Real>::polarDecomposition(orientationMatrix, polarMatrixQ);
        Quat quatA;
        quatA.fromMatrix(polarMatrixQ);
        Out::setCRot(out[point], quatA);
    }
} //apply


template<class In, class Out>
void BarycentricMapperTetrahedronSetTopologyRigid<In,Out>::applyJT( typename In::VecDeriv& out, const typename Out::VecDeriv& in )
{
    const sofa::helper::vector<topology::Tetrahedron>& tetrahedra = this->fromTopology->getTetrahedra();
    const sofa::helper::vector<MappingData >& map = this->map.getValue();
    typename core::behavior::MechanicalState<Out>* mechanicalObject;
    this->getContext()->get(mechanicalObject);

    const typename  Out::VecCoord& pX = *mechanicalObject->getX();

    // TODO: use mapOrient
    //const sofa::helper::vector<MappingOrientData >& mapOrient = this->mapOrient.getValue();
#ifdef SOFA_IP_TRACES
    IPTR_BARCPP_APPLYJT( "BarycentricMapperTetrahedronSetTopology::applyJT BEGIN " << endl);
    //unsigned int nbTetra = this->fromTopology->getTetrahedra().size();
    //IPTR_BARCPP_APPLYJT( "  number of tetrahedra = " << nbTetra );
#endif

    actualPos.clear();
    actualPos.resize(map.size());
    if ((!maskTo)||(maskTo&& !(maskTo->isInUse())) )
    {
        maskFrom->setInUse(false);
        for ( unsigned int i=0; i<map.size(); i++ )
        {
            //get velocity of the DoF
            const defaulttype::Rigid3dTypes::DPos v = defaulttype::Rigid3dTypes::getDPos(in[i]);

            //get its coordinated wrt to the associated tetra with given index
            const OutReal fx = ( OutReal ) map[i].baryCoords[0];
            const OutReal fy = ( OutReal ) map[i].baryCoords[1];
            const OutReal fz = ( OutReal ) map[i].baryCoords[2];
            int index = map[i].in_index;
            const topology::Tetrahedron& tetra = tetrahedra[index];

            out[tetra[0]] += v * ( 1-fx-fy-fz );
            out[tetra[1]] += v * fx;
            out[tetra[2]] += v * fy;
            out[tetra[3]] += v * fz;

            actualPos[i] = pX[i];


            //compute the linear forces for each vertex from the torque, inspired by rigid mapping
            Vector3 torque = getVOrientation(in[i]);

            for (unsigned int ti = 0; ti<4; ti++) {
                Vector3 lever;
                for (size_t dim = 0; dim < 3; dim++)
                    lever[dim] = actualTetraPosition[tetra[ti]][dim]-actualPos[i][dim];
                out[tetra[ti]] -= cross(lever,torque);
                //std::cout << "Force[" << tetra[ti] << "]: " << out[tetra[ti]] << std::endl;
            }


#ifdef SOFA_IP_TRACES
            //IPTR_BARCPP_APPLYJT("torque = " << torque  << endl); //  " |"<<torqueMagnitude<<"| "<< endl); // in point = " << glPointPositions[i] << endl);
            //IPTR_BARCPP_APPLYJT("   change of force due to the torque = " << cross(actualTetraPosition[tetra[0]],torque) << endl);*/
#endif
        }
    }
    else
    {
        typedef helper::ParticleMask ParticleMask;
        const ParticleMask::InternalStorage &indices=maskTo->getEntries();


        ParticleMask::InternalStorage::const_iterator it;
        for (it=indices.begin(); it!=indices.end(); it++)
        {
            const int i=(int)(*it);
            const defaulttype::Rigid3dTypes::DPos v = defaulttype::Rigid3dTypes::getDPos(in[i]);
            const OutReal fx = ( OutReal ) map[i].baryCoords[0];
            const OutReal fy = ( OutReal ) map[i].baryCoords[1];
            const OutReal fz = ( OutReal ) map[i].baryCoords[2];
            int index = map[i].in_index;
            const topology::Tetrahedron& tetra = tetrahedra[index];
            out[tetra[0]] += v * ( 1-fx-fy-fz );
            out[tetra[1]] += v * fx;
            out[tetra[2]] += v * fy;
            out[tetra[3]] += v * fz;

            //compute the linear forces for each vertex from the torque, inspired by rigid mapping
            Vector3 torque = getVOrientation(in[i]);
            //if (torque.norm() > 10e-6) {
            for (unsigned int ti = 0; ti<4; ti++)
                out[tetra[ti]] -= cross(actualTetraPosition[tetra[ti]],torque);
            //}
            maskFrom->insertEntry(tetra[0]);
            maskFrom->insertEntry(tetra[1]);
            maskFrom->insertEntry(tetra[2]);
            maskFrom->insertEntry(tetra[3]);
        }
    }

    actualOut.clear();
    actualOut.resize(out.size());
    for (size_t i = 0; i < out.size(); i++)
        for (size_t j = 0; j < out[i].size(); j++)
            actualOut[i][j] = 0.1*out[i][j];

}  //applyJT


template<class In, class Out>
void BarycentricMapperTetrahedronSetTopologyRigid<In,Out>::applyJ( typename Out::VecDeriv& out, const typename In::VecDeriv& in )
{
#ifdef SOFA_IP_TRACES
    IPTR_BARCPP_APPLYJ( "BarycentricMapperTetrahedronSetTopology::applyJ BEGIN " << endl);
    IPTR_BARCPP_APPLYJ( "Velocities: " << in << endl);
#endif
    const sofa::helper::vector<topology::Tetrahedron>& tetrahedra = this->fromTopology->getTetrahedra();
    const sofa::helper::vector<MappingData >& map = this->map.getValue();
    // TODO: use mapOrient
    //const sofa::helper::vector<MappingOrientData >& mapOrient = this->mapOrient.getValue();
    out.resize ( map.size() );


    if ((!maskTo)||(maskTo&& !(maskTo->isInUse())) )
    {
        for ( unsigned int i=0; i<map.size(); i++ )
        {
            const Real fx = map[i].baryCoords[0];
            const Real fy = map[i].baryCoords[1];
            const Real fz = map[i].baryCoords[2];
            int index = map[i].in_index;
            const topology::Tetrahedron& tetra = tetrahedra[index];
            Vector3 actualDPos = in[tetra[0]] * ( 1-fx-fy-fz )
                    + in[tetra[1]] * fx
                    + in[tetra[2]] * fy
                    + in[tetra[3]] * fz;
            Out::setDPos(out[i], actualDPos);

            Vector3 actualDRot(0,0,0);
            for (unsigned int vert = 0; vert < 4; vert++)
            {
                //if (in[tetra[vert]].norm() > 10e-6)
                actualDRot += cross(actualTetraPosition[tetra[vert]], in[tetra[vert]]);
            }

            Out::setDRot(out[i], actualDRot);


        }
    }
    else
    {
        typedef helper::ParticleMask ParticleMask;
        const ParticleMask::InternalStorage &indices=maskTo->getEntries();


        ParticleMask::InternalStorage::const_iterator it;
        for (it=indices.begin(); it!=indices.end(); it++)
        {
            const int i=(int)(*it);
            const Real fx = map[i].baryCoords[0];
            const Real fy = map[i].baryCoords[1];
            const Real fz = map[i].baryCoords[2];
            int index = map[i].in_index;
            const topology::Tetrahedron& tetra = tetrahedra[index];
            Out::setDPos(out[i] , in[tetra[0]] * ( 1-fx-fy-fz )
                    + in[tetra[1]] * fx
                    + in[tetra[2]] * fy
                    + in[tetra[3]] * fz );

            Vector3 actualDRot(0,0,0);
            for (unsigned int vert = 0; vert < 4; vert++)
            {
                //if (in[tetra[vert]].norm() > 10e-6)
                actualDRot += cross(actualTetraPosition[tetra[vert]], in[tetra[vert]]);
            }

            Out::setDRot(out[i], actualDRot);
        }
    }
#ifdef SOFA_IP_TRACES
    IPTR_BARCPP_APPLYJ( "BarycentricMapperTetrahedronSetTopology::applyJ END " << endl);
#endif
} //applyJ


template<class In, class Out>
const sofa::defaulttype::BaseMatrix* BarycentricMapperTetrahedronSetTopologyRigid<In,Out>::getJ(int outSize, int inSize)
{    
    //if (matrixJ && !updateJ && matrixJ->rowBSize() == (unsigned)outSize && matrixJ->colBSize() == (unsigned)inSize)
    //    return matrixJ;
    if (outSize > 0 && map.getValue().size() == 0)
    {
        std::cout << "Maps not created yet" << std::endl;
        return NULL; // error: maps not yet created ?
    }
    if (!matrixJ)
    {        
        matrixJ = new MatrixType;
    }

    if (matrixJ->rowBSize() != (unsigned)outSize || matrixJ->colBSize() != (unsigned)inSize)
    {
        //std::cout << "Resizing to " << outSize*NOut  << " X " << inSize*NIn << std::endl;
        matrixJ->resize(outSize*NOut, inSize*NIn);
    }
    else
        matrixJ->clear();

    const sofa::helper::vector<topology::Tetrahedron>& tetrahedra = this->fromTopology->getTetrahedra();
    const sofa::helper::vector<MappingData >& map = this->map.getValue();
    // TODO: use mapOrient
    //const sofa::helper::vector<MappingOrientData >& mapOrient = this->mapOrient.getValue();


    for (unsigned int beamNode = 0; beamNode < map.size(); beamNode++)
    {
        //linear forces
        const OutReal fx = ( OutReal ) map[beamNode].baryCoords[0];
        const OutReal fy = ( OutReal ) map[beamNode].baryCoords[1];
        const OutReal fz = ( OutReal ) map[beamNode].baryCoords[2];


        int index = map[beamNode].in_index;
        const topology::Tetrahedron& tetra = tetrahedra[index];

        for (int dim = 0; dim < 3; dim++)
        {
            matrixJ->add(beamNode*6+dim, 3*tetra[0]+dim, 1-fx-fy-fz);
            matrixJ->add(beamNode*6+dim, 3*tetra[1]+dim, fx);
            matrixJ->add(beamNode*6+dim, 3*tetra[2]+dim, fy);
            matrixJ->add(beamNode*6+dim, 3*tetra[3]+dim, fz);
        }

        for (int vert = 0; vert < 4; vert++)
        {
            Vector3 v;
            for (size_t dim = 0; dim < 3; dim++)
                v[dim] = actualTetraPosition[tetra[vert]][dim] - actualPos[beamNode][dim];
            matrixJ->add(beamNode*6+3, 3*tetra[vert]+1, -v[2]);
            matrixJ->add(beamNode*6+3, 3*tetra[vert]+2, +v[1]);
            matrixJ->add(beamNode*6+4, 3*tetra[vert]+0, +v[2]);
            matrixJ->add(beamNode*6+4, 3*tetra[vert]+2, -v[0]);
            matrixJ->add(beamNode*6+5, 3*tetra[vert]+0, -v[1]);
            matrixJ->add(beamNode*6+5, 3*tetra[vert]+1, +v[0]);            
        }
    }

    matrixJ->compress();
    updateJ = false;

    return matrixJ;
} // getJ


template <class In, class Out>
void BarycentricMapperTetrahedronSetTopologyRigid<In,Out>::applyJT ( typename In::MatrixDeriv& out, const typename Out::MatrixDeriv& in )
{
    typename Out::MatrixDeriv::RowConstIterator rowItEnd = in.end();
    const sofa::helper::vector<topology::Tetrahedron>& tetrahedra = this->fromTopology->getTetrahedra();
    const sofa::helper::vector<MappingData >& map = this->map.getValue();
    // TODO: use mapOrient
    //const sofa::helper::vector<MappingOrientData >& mapOrient = this->mapOrient.getValue();

    for (typename Out::MatrixDeriv::RowConstIterator rowIt = in.begin(); rowIt != rowItEnd; ++rowIt)
    {
        typename Out::MatrixDeriv::ColConstIterator colItEnd = rowIt.end();
        typename Out::MatrixDeriv::ColConstIterator colIt = rowIt.begin();

        if (colIt != colItEnd)
        {
            typename In::MatrixDeriv::RowIterator o = out.writeLine(rowIt.index());

            for ( ; colIt != colItEnd; ++colIt)
            {
                unsigned indexIn = colIt.index();
                InDeriv data = (InDeriv) Out::getDPos(colIt.val());

                const OutReal fx = ( OutReal ) map[indexIn].baryCoords[0];
                const OutReal fy = ( OutReal ) map[indexIn].baryCoords[1];
                const OutReal fz = ( OutReal ) map[indexIn].baryCoords[2];
                int index = map[indexIn].in_index;
                const topology::Tetrahedron& tetra = tetrahedra[index];

                o.addCol (tetra[0], data * ( 1-fx-fy-fz ) );
                o.addCol (tetra[1], data * fx );
                o.addCol (tetra[2], data * fy );
                o.addCol (tetra[3], data * fz );
            }
        }
    }
}

template <class In, class Out>
void BarycentricMapperTetrahedronSetTopologyRigid<In,Out>::draw  (const core::visual::VisualParams* vparams,const typename Out::VecCoord& out, const typename In::VecCoord& in )
{
    const sofa::helper::vector<topology::Tetrahedron>& tetrahedra = this->fromTopology->getTetrahedra();
    const sofa::helper::vector<MappingData >& map = this->map.getValue();
    // TODO: use mapOrient
    //const sofa::helper::vector<MappingOrientData >& mapOrient = this->mapOrient.getValue();

    std::vector< Vector3 > points;
    {
        for ( unsigned int i=0; i<map.size(); i++ )
        {
            const Real fx = map[i].baryCoords[0];
            const Real fy = map[i].baryCoords[1];
            const Real fz = map[i].baryCoords[2];
            int index = map[i].in_index;
            const topology::Tetrahedron& tetra = tetrahedra[index];
            Real f[4];
            f[0] = ( 1-fx-fy-fz );
            f[1] = fx;
            f[2] = fy;
            f[3] = fz;
            for ( int j=0; j<4; j++ )
            {
                if ( f[j]<=-0.0001 || f[j]>=0.0001 )
                {
                    //                     glColor3f((float)f[j],1,(float)f[j]);
                    points.push_back ( Out::getCPos(out[i]) );
                    points.push_back ( in[tetra[j]] );
                }
            }
        }
    }
    vparams->drawTool()->drawLines ( points, 1, Vec<4,float> ( 0,1,0,1 ) );
    //std::cout << "Drawing" << std::endl;

    for ( unsigned int i=0; i<map.size(); i++ )
    {
        //get velocity of the DoF
        //const defaulttype::Rigid3dTypes::DPos v = defaulttype::Rigid3dTypes::getDPos(in[i]);

        //get its coordinated wrt to the associated tetra with given index
        //const OutReal fx = ( OutReal ) map[i].baryCoords[0];
        //const OutReal fy = ( OutReal ) map[i].baryCoords[1];
        //const OutReal fz = ( OutReal ) map[i].baryCoords[2];
        int index = map[i].in_index;
        const topology::Tetrahedron& tetra = tetrahedra[index];

        //out[tetra[0]] += v * ( 1-fx-fy-fz );
        //out[tetra[1]] += v * fx;
        //out[tetra[2]] += v * fy;
        //out[tetra[3]] += v * fz;

        //compute the linear forces for each vertex from the torque, inspired by rigid mapping
        //Vector3 torque = getVOrientation(in[i]);
        //if (torque.norm() > 10e-6) {
        //for (unsigned int ti = 0; ti<4; ti++)
        //    out[tetra[ti]] -= cross(actualTetraPosition[tetra[ti]],torque);
        //}

        for (size_t i = 0; i < actualPos.size(); i++) {
            glPointSize(10);
            glColor3d(1.0,0,0.0);
            glBegin(GL_POINTS);
            helper::gl::glVertexT(actualPos[i]);
            //std::cout << "DRW " << actualPos[i] << std::endl;
            glEnd();

        }


        for (unsigned int ti = 0; ti<4; ti++) {
            glPointSize(10);
            glColor3d(1.0,0,1.0);
            glBegin(GL_POINTS);
            helper::gl::glVertexT(actualTetraPosition[tetra[ti]]);
            glEnd();


            if (tetra[ti] < actualOut.size()) {
                //std::cout << "Drawing the linear force in " << tetra[ti] << ": " << actualOut[tetra[ti]] << std::endl;
                glLineWidth(3.0);
                glBegin(GL_LINES);
                helper::gl::glVertexT(actualTetraPosition[tetra[ti]]);
                helper::gl::glVertexT(actualTetraPosition[tetra[ti]]+actualOut[tetra[ti]]);
                glEnd();
            }

        }
        glEnd();
    }
}

} // namespace mapping

} // namespace component

} // namespace sofa

#endif
