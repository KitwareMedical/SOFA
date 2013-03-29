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
*                               SOFA :: Plugins                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_COMPONENT_MAPPING_PrincipalStretchesMAPPING_H
#define SOFA_COMPONENT_MAPPING_PrincipalStretchesMAPPING_H

#include "../initFlexible.h"
#include "../strainMapping/BaseStrainMapping.h"
#include "../strainMapping/PrincipalStretchesJacobianBlock.h"

#include <sofa/helper/OptionsGroup.h>

namespace sofa
{
namespace component
{
namespace mapping
{

using helper::vector;


/** Deformation Gradient to Principal Stretches mapping.
*/



//////////////////////////////////////////////////////////////////////////////////
////  F -> U, F -> D
//////////////////////////////////////////////////////////////////////////////////

template <class TIn, class TOut>
class PrincipalStretchesMapping : public BaseStrainMappingT<defaulttype::PrincipalStretchesJacobianBlock<TIn,TOut> >
{
public:
    typedef defaulttype::PrincipalStretchesJacobianBlock<TIn,TOut> BlockType;
    typedef BaseStrainMappingT<BlockType > Inherit;

    SOFA_CLASS(SOFA_TEMPLATE2(PrincipalStretchesMapping,TIn,TOut), SOFA_TEMPLATE(BaseStrainMappingT,BlockType ));

    Data<bool> asStrain;
    Data<SReal> threshold;
//    Data<bool> f_PSDStabilization;


    virtual void reinit()
    {
        for(unsigned int i=0; i<this->jacobian.size(); i++)
        {
            this->jacobian[i].init( asStrain.getValue(), threshold.getValue()/*, f_PSDStabilization.getValue()*/ );
        }
        Inherit::reinit();
    }


protected:

    PrincipalStretchesMapping (core::State<TIn>* from = NULL, core::State<TOut>* to= NULL)
        : Inherit ( from, to )
        , asStrain(initData(&asStrain,false,"asStrain","compute principal stretches - 1"))
        , threshold(initData(&threshold,-std::numeric_limits<SReal>::max(),"threshold","threshold the principal stretches to ensure detF=J=U1*U2*U3 is not too close or < 0"))
//        , f_PSDStabilization(initData(&f_PSDStabilization,false,"PSDStabilization","project geometric stiffness sub-matrices to their nearest symmetric, positive semi-definite matrices"))
    {
    }

    virtual ~PrincipalStretchesMapping() { }

};



} // namespace mapping
} // namespace component
} // namespace sofa

#endif // SOFA_COMPONENT_MAPPING_PrincipalStretchesMAPPING_H
