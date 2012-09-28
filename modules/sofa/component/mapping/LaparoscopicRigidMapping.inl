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
#ifndef SOFA_COMPONENT_MAPPING_LAPAROSCOPICRIGIDMAPPING_INL
#define SOFA_COMPONENT_MAPPING_LAPAROSCOPICRIGIDMAPPING_INL

#include <sofa/component/mapping/LaparoscopicRigidMapping.h>
#include <sofa/core/visual/VisualParams.h>

#include <sofa/core/Mapping.inl>

#include <string>


namespace sofa
{

namespace component
{

namespace mapping
{

template <class TIn, class TOut>
void LaparoscopicRigidMapping<TIn, TOut>::init()
{
    Inherit::init();
}

template <class TIn, class TOut>
void LaparoscopicRigidMapping<TIn, TOut>::apply(const core::MechanicalParams * /*mparams*/ /* PARAMS FIRST */, Data<OutVecCoord>& dOut, const Data<InVecCoord>& dIn)
{
    helper::WriteAccessor< Data<OutVecCoord> > out = dOut;
    helper::ReadAccessor< Data<InVecCoord> > in = dIn;

    out.resize(1);
    out[0].getOrientation() = in[0].getOrientation(); // * rotation.getValue();
    out[0].getCenter() = pivot.getValue() + in[0].getOrientation().rotate(sofa::defaulttype::Vector3(0,0,in[0].getTranslation()));
    currentRotation = in[0].getOrientation();
}

template <class TIn, class TOut>
void LaparoscopicRigidMapping<TIn, TOut>::applyJ(const core::MechanicalParams * /*mparams*/ /* PARAMS FIRST */, Data<OutVecDeriv>& dOut, const Data<InVecDeriv>& dIn)
{
    helper::WriteAccessor< Data<OutVecDeriv> > out = dOut;
    helper::ReadAccessor< Data<InVecDeriv> > in = dIn;

    out.resize(1);
    getVOrientation(out[0]) =in[0].getVOrientation(); //rotation * in[0].getVOrientation();
    getVCenter(out[0]) = currentRotation.rotate(sofa::defaulttype::Vector3(0,0,in[0].getVTranslation()));
}

template <class TIn, class TOut>
void LaparoscopicRigidMapping<TIn, TOut>::applyJT(const core::MechanicalParams * /*mparams*/ /* PARAMS FIRST */, Data<InVecDeriv>& dOut, const Data<OutVecDeriv>& dIn)
{
    helper::WriteAccessor< Data<InVecDeriv> > out = dOut;
    helper::ReadAccessor< Data<OutVecDeriv> > in = dIn;

    out[0].getVOrientation() += getVOrientation(in[0]); //rotation * in[0].getVOrientation();
    out[0].getVTranslation() += dot(currentRotation.rotate(sofa::defaulttype::Vector3(0,0,1)), getVCenter(in[0]));
}

template <class TIn, class TOut>
void LaparoscopicRigidMapping<TIn, TOut>::draw(const core::visual::VisualParams* vparams)
{
    if (!vparams->displayFlags().getShowMappings())
        return;
}

} // namespace mapping

} // namespace component

} // namespace sofa

#endif
