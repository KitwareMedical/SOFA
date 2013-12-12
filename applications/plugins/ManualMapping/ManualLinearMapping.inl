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
#ifndef SOFA_COMPONENT_MAPPING_MANUALLINEARMAPPING_INL
#define SOFA_COMPONENT_MAPPING_MANUALLINEARMAPPING_INL

#include "ManualLinearMapping.h"
#include <sofa/defaulttype/VecTypes.h>

#include <sofa/core/Mapping.inl>


namespace sofa
{

namespace component
{

namespace mapping
{




// _matrixJ must have been filled before calling init()
template<class TIn, class TOut>
void ManualLinearMapping<TIn, TOut>::init()
{
    Inherit::init();

    // resize output
    this->toModel->resize( _matrixJ.rowSize() / NOut );
}



template <class TIn, class TOut>
void ManualLinearMapping<TIn, TOut>::apply(const core::MechanicalParams * /*mparams*/ /* PARAMS FIRST */, Data<VecCoord>& dOut, const Data<InVecCoord>& dIn)
{
//    if( !_matrixJ.colSize() && !_matrixJ.rowSize() ) return;

    _matrixJ.mult( *reinterpret_cast<Data<VecDeriv>*>(&dOut), *reinterpret_cast<const Data<InVecDeriv>*>(&dIn) );
}

template <class TIn, class TOut>
void ManualLinearMapping<TIn, TOut>::applyJ(const core::MechanicalParams * /*mparams*/ /* PARAMS FIRST */, Data<VecDeriv>& dOut, const Data<InVecDeriv>& dIn)
{
//    if( !_matrixJ.colSize() && !_matrixJ.rowSize() ) return;

    _matrixJ.mult( dOut, dIn );
}

template<class TIn, class TOut>
void ManualLinearMapping<TIn, TOut>::applyJT(const core::MechanicalParams * /*mparams*/ /* PARAMS FIRST */, Data<InVecDeriv>& dOut, const Data<VecDeriv>& dIn)
{
//    if( !_matrixJ.colSize() && !_matrixJ.rowSize() ) return;

    _matrixJ.addMultTranspose( dOut, dIn );
}

template <class TIn, class TOut>
void ManualLinearMapping<TIn, TOut>::applyJT(const core::ConstraintParams * /*cparams*/ /* PARAMS FIRST */, Data<InMatrixDeriv>& /*dOut*/, const Data<MatrixDeriv>& /*dIn*/)
{
//    serr<<SOFA_CLASS_METHOD<<"not yet implemented"<<sendl;
}


template <class TIn, class TOut>
const sofa::defaulttype::BaseMatrix* ManualLinearMapping<TIn, TOut>::getJ()
{
    return &_matrixJ;
}


template <class TIn, class TOut>
const typename ManualLinearMapping<TIn, TOut>::js_type* ManualLinearMapping<TIn, TOut>::getJs()
{
	return &js;
}


} // namespace mapping

} // namespace component

} // namespace sofa

#endif
