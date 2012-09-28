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
#ifndef SOFA_COMPONENT_MAPPING_JointRigidMapping_INL
#define SOFA_COMPONENT_MAPPING_JointRigidMapping_INL

#include "../deformationMapping/JointRigidMapping.h"
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
JointRigidMapping<TIn, TOut>::JointRigidMapping()
    : Inherit(),
      source(initData(&source, "source", "input dof and rigid offset for each output dof"))
{
}

template <class TIn, class TOut>
JointRigidMapping<TIn, TOut>::~JointRigidMapping()
{
}


template <class TIn, class TOut>
void JointRigidMapping<TIn, TOut>::init()
{
    // TODO wrap this somehow ?
    baseMatrices.resize( 1 );
    baseMatrices[0] = &jacobian;

    this->Inherit::init();  // applies the mapping, so after the Data init
}


// namespace impl {

//   // delta computations for different dof types
//   template<class U>
//   defaulttype::RigidCoord<3, U> delta(const defaulttype::RigidCoord<3, U>& lhs,
// 					const defaulttype::RigidCoord<3, U>& rhs) {
//     SE3<U> se3;

//     return se3.prod(lhs, se3.inv(rhs));
//   }

//   template<int N, class U>
//   defaulttype::Vec<N, U> delta(const defaulttype::Vec<N, U>& lhs,
// 				 const defaulttype::Vec<N, U>& rhs) {
//     return lhs - rhs;
//   }



//   template<class U>
//   typename SE3<U>::mat66 d_delta(const defaulttype::RigidCoord<3, U>& lhs,
// 				   const defaulttype::RigidCoord<3, U>& rhs) {
//     SE3<U> se3;

//     return se3.dR(se3.inv(rhs), lhs );
//   }


//   template<int N, class U>
//   Eigen::Matrix<U, N, N> d_delta(const defaulttype::Vec<N, U>&,
// 				   const defaulttype::Vec<N, U>& ) {
//     return Eigen::Matrix<U, N, N>::Identity();
//   }


// }


template <class TIn, class TOut>
void JointRigidMapping<TIn, TOut>::apply(const core::MechanicalParams * /*mparams*/ , Data<OutVecCoord>& dOut, const Data<InVecCoord>& dIn)
{
    helper::WriteAccessor< Data<OutVecCoord> >  out = dOut;
    helper::ReadAccessor< Data<InVecCoord> >  in = dIn;

    assert( in.size() <= out.size() );
    assert(  Nout == Nin  );
    assert( out.size() >= source.getValue().size() );

    SE3<Real> se3;

    // TODO is this needed on each apply ?
    this->jacobian.resizeBlocks(out.size(), in.size());

    typedef unsigned int uint;

    for(uint d = 0; d < out.size(); ++d)
    {

        const source_type& s = source.getValue()[d];

        out[ d ] = se3.prod( in[ s.first ], s.second );

        Eigen::Matrix<Real, Nout, Nin> chunk = se3.dR(s.second, in[ s.first ] );

        for(uint i = 0; i < Nout; ++i)
        {
            uint row = d * Nout + i;
            jacobian.beginRow( row );

            for(uint j = 0; j < Nin; ++j)
            {
                uint col = s.first * Nin + j;
                jacobian.insertBack(row, col, chunk(i, j));
            }
        }
    }

    jacobian.compress();
}

template <class TIn, class TOut>
void JointRigidMapping<TIn, TOut>::applyJ(const core::MechanicalParams * /*mparams*/ , Data<OutVecDeriv>& dOut, const Data<InVecDeriv>& dIn)
{
    if( jacobian.rowSize() > 0 )
        jacobian.mult(dOut,dIn);
}

template <class TIn, class TOut>
void JointRigidMapping<TIn, TOut>::applyJT(const core::MechanicalParams * /*mparams*/ , Data<InVecDeriv>& dIn, const Data<OutVecDeriv>& dOut)
{
    if( jacobian.rowSize() > 0 )
        jacobian.addMultTranspose(dIn,dOut);
}

template <class TIn, class TOut>
void JointRigidMapping<TIn, TOut>::applyJT(const core::ConstraintParams*, Data<InMatrixDeriv>& , const Data<OutMatrixDeriv>& )
{
    //    cerr<<"JointRigidMapping<TIn, TOut>::applyJT does nothing " << endl;
}


template <class TIn, class TOut>
const sofa::defaulttype::BaseMatrix* JointRigidMapping<TIn, TOut>::getJ()
{
    return &jacobian;
}

template <class TIn, class TOut>
const vector<sofa::defaulttype::BaseMatrix*>* JointRigidMapping<TIn, TOut>::getJs()
{
    return &baseMatrices;
}

template <class TIn, class TOut>
void JointRigidMapping<TIn, TOut>::draw(const core::visual::VisualParams* )
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
