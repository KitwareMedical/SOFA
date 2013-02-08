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
#ifndef SOFA_COMPONENT_MAPPING_TriangleStrainAverageMapping_INL
#define SOFA_COMPONENT_MAPPING_TriangleStrainAverageMapping_INL

#include "../deformationMapping/TriangleStrainAverageMapping.h"
#include <sofa/core/visual/VisualParams.h>
#include <map>
#include <iostream>
using std::cerr;
using std::endl;

namespace sofa
{

namespace component
{

namespace mapping
{

using namespace sofa::defaulttype;

template <class TIn, class TOut>
TriangleStrainAverageMapping<TIn, TOut>::TriangleStrainAverageMapping()
    : Inherit()
    , f_triangleIndices(initData(&f_triangleIndices, "triangleIndices", "For each node, index of the adjacent triangles."))
    , f_endIndices(initData(&f_endIndices,  "endIndices", "For each node, end index of its triangle list"))
    , f_weights(initData(&f_weights, "weights", "For each node, weight of each triangle in the average"))
{

}

template <class TIn, class TOut>
TriangleStrainAverageMapping<TIn, TOut>::~TriangleStrainAverageMapping()
{
}


template <class TIn, class TOut>
void TriangleStrainAverageMapping<TIn, TOut>::init()
{
    triangleContainer = this->core::objectmodel::BaseObject::searchUp<topology::TriangleSetTopologyContainer>();
    if( !triangleContainer ) serr<<"No TriangleSetTopologyContainer found ! "<<sendl;

    SeqTriangles triangles = triangleContainer->getTriangles();

    typedef std::map< unsigned, std::set<unsigned> > MapNodeToTriangles;
    MapNodeToTriangles nodeToTriangles;
    for(unsigned i=0; i<triangles.size(); i++ )
    {
        for(unsigned j=0; j<3; j++)
            nodeToTriangles[triangles[i][j]].insert(i);
    }

    helper::WriteAccessor< Data<vector<unsigned> > > triangleIndices(f_triangleIndices);
    helper::WriteAccessor< Data<vector<unsigned> > > endIndices(f_endIndices);
    helper::WriteAccessor< Data<vector<Real> > > weights(f_weights);
    unsigned endIndex=0;
    for( MapNodeToTriangles::const_iterator i=nodeToTriangles.begin(), iend=nodeToTriangles.end(); i!=iend; i++ ) // for each node
    {
        unsigned n = (*i).second.size();
        endIndex += n;
        endIndices.push_back( endIndex );
        for( std::set<unsigned>::const_iterator j=(*i).second.begin(), jend=(*i).second.end(); j!=jend; j++ ) // for each triangle connected to the node
        {
            triangleIndices.push_back(*j);
            weights.push_back( 1./n );  // todo: used triangle areas
        }
    }

    this->toModel->resize(endIndices.size());
    diagMat.resize(endIndices.size());

    jacobian.resizeBlocks(endIndices.size(),triangles.size());
    unsigned startIndex=0;
    for(unsigned i=0; i<endIndices.size(); i++ )
    {
        diagMat[i] = sqrt( (endIndices[i]-startIndex)/3.0 );  // sqrt of the area associated with the node. Todo: use triangle areas rather than assuming unit areas.

        jacobian.beginBlockRow(i);
        for( unsigned j=startIndex; j<endIndices[i]; j++)
        {
            jacobian.createBlock( triangleIndices[j], Block::Identity() * weights[j] * diagMat[i] );
        }
        jacobian.endBlockRow();

        startIndex = endIndices[i];
    }

//    cerr<<"TriangleStrainAverageMapping<TIn, TOut>::init, indices = "<< triangleIndices << endl;
//    cerr<<"TriangleStrainAverageMapping<TIn, TOut>::init, bounds = "<< endIndices << endl;
//    cerr<<"TriangleStrainAverageMapping<TIn, TOut>::init, weights = "<< weights << endl;
//    cerr<<"TriangleStrainAverageMapping<TIn, TOut>::init, diagMatrix = "<< diagMat << endl;

    baseMatrices.resize( 1 );
    baseMatrices[0] = &jacobian;

    this->Inherit::init();  // applies the mapping, so after the Data init
}



template <class TIn, class TOut>
void TriangleStrainAverageMapping<TIn, TOut>::mult( Data<OutVecCoord>& dOut, const Data<InVecCoord>& dIn )
{
    helper::WriteAccessor< Data<OutVecCoord> >  nodeValues = dOut;
    helper::ReadAccessor< Data<InVecCoord> >  triangleValues = dIn;
    helper::ReadAccessor< Data<vector<unsigned> > > triangleIndices(f_triangleIndices);
    helper::ReadAccessor< Data<vector<unsigned> > > endIndices(f_endIndices);
    helper::ReadAccessor< Data<vector<Real> > > weights(f_weights);

//    cerr<<"begin TriangleStrainAverageMapping<TIn, TOut>::mult, parent: " << triangleValues << endl;
    unsigned startIndex=0;
    for(unsigned i=0; i<endIndices.size(); i++ )
    {
//        cerr<<"TriangleStrainAverageMapping<TIn, TOut>::mult, node " << i << endl;
        nodeValues[i] = OutCoord();
        for( unsigned j=startIndex; j<endIndices[i]; j++)  // product with J
        {
            nodeValues[i] += triangleValues[triangleIndices[j]] * weights[j];
//            cerr<<"  TriangleStrainAverageMapping<TIn, TOut>::mult, add contribution from triangle " << triangleIndices[j] << " with weight " << weights[j] << ": " << triangleValues[triangleIndices[j]] * weights[j] << endl;
        }
        nodeValues[i] *= diagMat[i];                       // product with the diagonal matrix
        startIndex = endIndices[i];
    }
//    cerr<<"==== end TriangleStrainAverageMapping<TIn, TOut>::mult, final out: " << nodeValues << endl;
}


template <class TIn, class TOut>
void TriangleStrainAverageMapping<TIn, TOut>::apply(const core::MechanicalParams * , Data<OutVecCoord>& dOut, const Data<InVecCoord>& dIn)
{
    mult(dOut,dIn);
}

template <class TIn, class TOut>
void TriangleStrainAverageMapping<TIn, TOut>::applyJ(const core::MechanicalParams * , Data<OutVecDeriv>& dOut, const Data<InVecDeriv>& dIn)
{
    mult(dOut,dIn);
}

template <class TIn, class TOut>
void TriangleStrainAverageMapping<TIn, TOut>::applyJT(const core::MechanicalParams *, Data<InVecDeriv>& dIn, const Data<OutVecDeriv>& dOut)
{
    helper::ReadAccessor< Data<OutVecDeriv> >  nodeValues = dOut;
    helper::WriteAccessor< Data<InVecDeriv> >  triangleValues = dIn;
    helper::ReadAccessor< Data<vector<unsigned> > > triangleIndices(f_triangleIndices);
    helper::ReadAccessor< Data<vector<unsigned> > > endIndices(f_endIndices);
    helper::ReadAccessor< Data<vector<Real> > > weights(f_weights);

//    cerr<<"begin TriangleStrainAverageMapping<TIn, TOut>::applyJT, child vector : " << nodeValues << endl;
    unsigned startIndex=0;
    for(unsigned i=0; i<endIndices.size(); i++ )
    {
//        cerr<<"TriangleStrainAverageMapping<TIn, TOut>::applyJT, node " << i << endl;
        OutCoord val = nodeValues[i] * diagMat[i];          // product with the diagonal matrix
        for( unsigned j=startIndex; j<endIndices[i]; j++)   // product with J^T
        {
            triangleValues[triangleIndices[j]] +=  val * weights[j];
//           cerr<<"  TriangleStrainAverageMapping<TIn, TOut>::applyJT, add contribution to triangle " << triangleIndices[j] << " with weight " << weights[j] << ": " << val * weights[j] << endl;
        }
        startIndex = endIndices[i];
    }
//    cerr<<"==== end TriangleStrainAverageMapping<TIn, TOut>::applyJT, final in: " << triangleValues << endl;

}


template <class TIn, class TOut>
void TriangleStrainAverageMapping<TIn, TOut>::applyJT(const core::ConstraintParams*, Data<InMatrixDeriv>& , const Data<OutMatrixDeriv>& )
{
    //    cerr<<"TriangleStrainAverageMapping<TIn, TOut>::applyJT(const core::ConstraintParams*, Data<InMatrixDeriv>& , const Data<OutMatrixDeriv>& ) does nothing " << endl;
}


template <class TIn, class TOut>
const sofa::defaulttype::BaseMatrix* TriangleStrainAverageMapping<TIn, TOut>::getJ()
{
    return &jacobian;
}

template <class TIn, class TOut>
const vector<sofa::defaulttype::BaseMatrix*>* TriangleStrainAverageMapping<TIn, TOut>::getJs()
{
    return &baseMatrices;
}



} // namespace mapping

} // namespace component

} // namespace sofa

#endif
