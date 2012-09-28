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
#ifndef SOFA_COMPONENT_MAPPING_LogRigidMapping_INL
#define SOFA_COMPONENT_MAPPING_LogRigidMapping_INL

#include "../deformationMapping/LogRigidMapping.h"
#include <sofa/core/visual/VisualParams.h>
#include <iostream>

#include <sofa/core/Mapping.inl>

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
LogRigidMapping<TIn, TOut>::LogRigidMapping()
    : Inherit()

{
}

template <class TIn, class TOut>
LogRigidMapping<TIn, TOut>::~LogRigidMapping()
{
}


template <class TIn, class TOut>
void LogRigidMapping<TIn, TOut>::init()
{
    // TODO wrap this somehow ?
    baseMatrices.resize( 1 );
    baseMatrices[0] = &jacobian;

    this->Inherit::init();  // applies the mapping, so after the Data init
}


template <class TIn, class TOut>
void LogRigidMapping<TIn, TOut>::apply(const core::MechanicalParams * /*mparams*/ , Data<OutVecCoord>& dOut, const Data<InVecCoord>& dIn)
{
    helper::WriteAccessor< Data<OutVecCoord> >  out = dOut;
    helper::ReadAccessor< Data<InVecCoord> >  in = dIn;

    assert( in.size() <= out.size() );
    assert(  Nout == Nin  );

    // TODO is this needed on each apply ?
    this->jacobian.resizeBlocks(out.size(), in.size());

    se3_type se3;
    typedef unsigned int index_type;

    for(index_type d = 0; d < in.size(); ++d)
    {

        out[ d ] = se3.product_log( in[ d ] ).getVAll();

        typename se3_type::mat66 chunk = se3.product_dlog( in[d] );

        for(index_type i = 0; i < Nout; ++i)
        {
            index_type row = d * Nout + i;
            jacobian.beginRow( row );

            for(index_type j = 0; j < Nin; ++j)
            {
                index_type col = d * Nin + j;

                jacobian.insertBack( row, col, chunk(i, j) );

            }
        }
    }

    jacobian.compress();
}

template <class TIn, class TOut>
void LogRigidMapping<TIn, TOut>::applyJ(const core::MechanicalParams * /*mparams*/ , Data<OutVecDeriv>& dOut, const Data<InVecDeriv>& dIn)
{
    if( jacobian.rowSize() > 0 )
        jacobian.mult(dOut,dIn);
}

template <class TIn, class TOut>
void LogRigidMapping<TIn, TOut>::applyJT(const core::MechanicalParams * /*mparams*/ , Data<InVecDeriv>& dIn, const Data<OutVecDeriv>& dOut)
{
    if( jacobian.rowSize() > 0 )
        jacobian.addMultTranspose(dIn,dOut);
}

template <class TIn, class TOut>
void LogRigidMapping<TIn, TOut>::applyJT(const core::ConstraintParams*, Data<InMatrixDeriv>& , const Data<OutMatrixDeriv>& )
{
    //    cerr<<"LogRigidMapping<TIn, TOut>::applyJT does nothing " << endl;
}


template <class TIn, class TOut>
const sofa::defaulttype::BaseMatrix* LogRigidMapping<TIn, TOut>::getJ()
{
    return &jacobian;
}

template <class TIn, class TOut>
const vector<sofa::defaulttype::BaseMatrix*>* LogRigidMapping<TIn, TOut>::getJs()
{
    return &baseMatrices;
}

template <class TIn, class TOut>
void LogRigidMapping<TIn, TOut>::draw(const core::visual::VisualParams* /*vparams*/)
{
    // typename core::behavior::MechanicalState<In>::ReadVecCoord pos = this->getFromModel()->readPositions();
    // SeqEdges links = edgeContainer->getEdges();

    // vector< Vec3d > points;

    // for(unsigned i=0; i<links.size(); i++ ){
    //     points.push_back(pos[links[i][0]]);
    //     points.push_back(pos[links[i][1]]);
    // }
    // vparams->drawTool()->drawLines ( points, 1, Vec<4,float> ( 1,1,0,1 ) );
}


} // namespace mapping

} // namespace component

} // namespace sofa

#endif
