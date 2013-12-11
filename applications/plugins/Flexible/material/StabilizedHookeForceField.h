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
#ifndef SOFA_StabilizedHookeFORCEFIELD_H
#define SOFA_StabilizedHookeFORCEFIELD_H

#include "../initFlexible.h"
#include "../material/BaseMaterialForceField.h"
#include "../material/StabilizedHookeMaterialBlock.h"

#include <sofa/core/objectmodel/Event.h>
#include <sofa/simulation/common/AnimateEndEvent.h>

namespace sofa
{
namespace component
{
namespace forcefield
{

using helper::vector;

/** Apply Hooke's Law for isotropic homogeneous incompressible materials from principal stretches.
  * This is the stabilized formulation from "Energetically Consistent Invertible Elasticity", SCA'12
  *
  * The energy is : mu.sum_i((Ui-1)^2) + lambda/2(J-1)^2
*/

template <class _DataTypes>
class StabilizedHookeForceField : public BaseMaterialForceFieldT<defaulttype::StabilizedHookeMaterialBlock<_DataTypes> >
{
public:
    typedef defaulttype::StabilizedHookeMaterialBlock<_DataTypes> BlockType;
    typedef BaseMaterialForceFieldT<BlockType> Inherit;

    SOFA_CLASS(SOFA_TEMPLATE(StabilizedHookeForceField,_DataTypes),SOFA_TEMPLATE(BaseMaterialForceFieldT, BlockType));

    typedef typename Inherit::Real Real;

    /** @name  Material parameters */
    //@{
    Data<vector<Real> > _youngModulus;
    Data<vector<Real> > _poissonRatio;
    //@}

    virtual void reinit()
    {
        Real youngModulus=0,poissonRatio=0;
        for(unsigned int i=0; i<this->material.size(); i++)
        {
            if(i<_youngModulus.getValue().size()) youngModulus=_youngModulus.getValue()[i]; else if(_youngModulus.getValue().size()) youngModulus=_youngModulus.getValue()[0];
            if(i<_poissonRatio.getValue().size()) poissonRatio=_poissonRatio.getValue()[i]; else if(_poissonRatio.getValue().size()) poissonRatio=_poissonRatio.getValue()[0];

            assert( helper::isClamped<Real>( poissonRatio, -1+std::numeric_limits<Real>::epsilon(), 0.5-std::numeric_limits<Real>::epsilon() ) );

            this->material[i].init( youngModulus, poissonRatio );
        }
        Inherit::reinit();
    }

    void handleEvent(sofa::core::objectmodel::Event *event)
    {
        if ( dynamic_cast<simulation::AnimateEndEvent*>(event))
        {
            if(_youngModulus.isDirty() || _poissonRatio.isDirty()) reinit();
        }
    }

protected:
    StabilizedHookeForceField(core::behavior::MechanicalState<_DataTypes> *mm = NULL)
        : Inherit(mm)
        , _youngModulus(initData(&_youngModulus,vector<Real>((int)1,(Real)5000),"youngModulus","Young Modulus"))
        , _poissonRatio(initData(&_poissonRatio,vector<Real>((int)1,(Real)0),"poissonRatio","Poisson Ratio ]-1,0.5["))
//        , _viscosity(initData(&_viscosity,(Real)0,"viscosity","Viscosity (stress/strainRate)"))
    {
        this->f_listening.setValue(true);
    }

    virtual ~StabilizedHookeForceField()     {    }

};


}
}
}

#endif
