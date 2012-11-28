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
#ifndef SOFA_COMPONENT_MAPPING_InvariantMAPPING_H
#define SOFA_COMPONENT_MAPPING_InvariantMAPPING_H

#include "../initFlexible.h"
#include "../strainMapping/BaseStrainMapping.h"
#include "../strainMapping/InvariantJacobianBlock.inl"


namespace sofa
{
namespace component
{
namespace mapping
{

using helper::vector;

/** Map deformation gradients to the square root of (deviatoric) invariants of the right Cauchy Green deformation tensor: sqrt(I1),sqrt(I2) and J
*/

template <class TIn, class TOut>
class SOFA_Flexible_API InvariantMapping : public BaseStrainMappingT<defaulttype::InvariantJacobianBlock<TIn,TOut> >
{
public:
    typedef defaulttype::InvariantJacobianBlock<TIn,TOut> BlockType;
    typedef BaseStrainMappingT<BlockType > Inherit;

    SOFA_CLASS(SOFA_TEMPLATE2(InvariantMapping,TIn,TOut), SOFA_TEMPLATE(BaseStrainMappingT,BlockType ));

    Data<bool> deviatoric;

    virtual void reinit()
    {
        for(unsigned int i=0; i<this->jacobian.size(); i++) this->jacobian[i].deviatoric=deviatoric.getValue();
        Inherit::reinit();
    }

protected:
    InvariantMapping (core::State<TIn>* from = NULL, core::State<TOut>* to= NULL)
        : Inherit ( from, to )
        , deviatoric(initData(&deviatoric,true,"deviatoric","use deviatoric tensor: C.J^(-2/3)"))
    {
    }

    virtual ~InvariantMapping()     { }

};


} // namespace mapping
} // namespace component
} // namespace sofa

#endif
