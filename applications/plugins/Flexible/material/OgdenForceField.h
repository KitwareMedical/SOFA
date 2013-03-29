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
#ifndef SOFA_OgdenFORCEFIELD_H
#define SOFA_OgdenFORCEFIELD_H

#include "../initFlexible.h"
#include "../material/BaseMaterialForceField.h"
#include "../material/OgdenMaterialBlock.h"

#include <sofa/core/objectmodel/Event.h>
#include <sofa/simulation/common/AnimateEndEvent.h>

namespace sofa
{
namespace component
{
namespace forcefield
{

using helper::vector;


/** Ogden compressible energy with N=3 (note that a more generic implementation for any N is possible)
    W = sum(1<=i=<N) mui/alphai (~U1^alphai+~U2^alphai+~U3^alphai-3) + sum(1<=i=<N) 1/di(J-1)^{2i} with J = U1*U2*U3 and ~Ui=J^{-1/3}Ui deviatoric principal stretches
  */
template <class _DataTypes>
class OgdenForceField : public BaseMaterialForceFieldT<defaulttype::OgdenMaterialBlock<_DataTypes> >
{
public:
    typedef defaulttype::OgdenMaterialBlock<_DataTypes> BlockType;
    typedef BaseMaterialForceFieldT<BlockType> Inherit;

    SOFA_CLASS(SOFA_TEMPLATE(OgdenForceField,_DataTypes),SOFA_TEMPLATE(BaseMaterialForceFieldT, BlockType));

    typedef typename Inherit::Real Real;

    /** @name  Material parameters */
    //@{
    Data<vector<Real> > f_mu1;
    Data<vector<Real> > f_mu2;
    Data<vector<Real> > f_mu3;
    Data<vector<Real> > f_alpha1;
    Data<vector<Real> > f_alpha2;
    Data<vector<Real> > f_alpha3;
    Data<vector<Real> > f_d1;
    Data<vector<Real> > f_d2;
    Data<vector<Real> > f_d3;
    Data<bool > f_PSDStabilization;
    //@}

    virtual void reinit()
    {
        Real mu1=0,mu2=0,mu3=0,alpha1=0,alpha2=0,alpha3=0,d1=0,d2=0,d3=0;
        for(unsigned int i=0; i<this->material.size(); i++)
        {
            if(i<f_mu1.getValue().size()) mu1=f_mu1.getValue()[i]; else if(f_mu1.getValue().size()) mu1=f_mu1.getValue()[0];
            if(i<f_mu2.getValue().size()) mu2=f_mu2.getValue()[i]; else if(f_mu2.getValue().size()) mu2=f_mu2.getValue()[0];
            if(i<f_mu3.getValue().size()) mu3=f_mu3.getValue()[i]; else if(f_mu3.getValue().size()) mu3=f_mu3.getValue()[0];
            if(i<f_alpha1.getValue().size()) alpha1=f_alpha1.getValue()[i]; else if(f_alpha1.getValue().size()) alpha1=f_alpha1.getValue()[0];
            if(i<f_alpha2.getValue().size()) alpha2=f_alpha2.getValue()[i]; else if(f_alpha2.getValue().size()) alpha2=f_alpha2.getValue()[0];
            if(i<f_alpha3.getValue().size()) alpha3=f_alpha3.getValue()[i]; else if(f_alpha3.getValue().size()) alpha3=f_alpha3.getValue()[0];
            if(i<f_d1.getValue().size()) d1=f_d1.getValue()[i]; else if(f_d1.getValue().size()) d1=f_d1.getValue()[0];
            if(i<f_d2.getValue().size()) d2=f_d2.getValue()[i]; else if(f_d2.getValue().size()) d2=f_d2.getValue()[0];
            if(i<f_d3.getValue().size()) d3=f_d3.getValue()[i]; else if(f_d3.getValue().size()) d3=f_d3.getValue()[0];
            this->material[i].init( mu1,mu2,mu3,alpha1,alpha2,alpha3,d1,d2,d3,f_PSDStabilization.getValue() );
        }
        Inherit::reinit();
    }

    void handleEvent(sofa::core::objectmodel::Event *event)
    {
        if ( dynamic_cast<simulation::AnimateEndEvent*>(event))
        {
            if( f_mu1.isDirty() || f_mu2.isDirty() || f_mu3.isDirty() ||
                 f_alpha1.isDirty() || f_alpha2.isDirty() || f_alpha3.isDirty() ||
                    f_d1.isDirty() || f_d2.isDirty() || f_d3.isDirty() ||
                    f_PSDStabilization.isDirty() ) reinit();
        }
    }


protected:
    OgdenForceField(core::behavior::MechanicalState<_DataTypes> *mm = NULL)
        : Inherit(mm)
        , f_mu1(initData(&f_mu1,vector<Real>((int)1,(Real)1000),"mu1",""))
        , f_mu2(initData(&f_mu2,vector<Real>((int)1,(Real)1000),"mu2",""))
        , f_mu3(initData(&f_mu3,vector<Real>((int)1,(Real)1000),"mu3",""))
        , f_alpha1(initData(&f_alpha1,vector<Real>((int)1,(Real)1000),"alpha1",""))
        , f_alpha2(initData(&f_alpha2,vector<Real>((int)1,(Real)1000),"alpha2",""))
        , f_alpha3(initData(&f_alpha3,vector<Real>((int)1,(Real)1000),"alpha3",""))
        , f_d1(initData(&f_d1,vector<Real>((int)1,(Real)1000),"d1",""))
        , f_d2(initData(&f_d2,vector<Real>((int)1,(Real)1000),"d2",""))
        , f_d3(initData(&f_d3,vector<Real>((int)1,(Real)1000),"d3",""))
        , f_PSDStabilization(initData(&f_PSDStabilization,false,"PSDStabilization","project stiffness matrix to its nearest symmetric, positive semi-definite matrix"))
    {
        this->f_listening.setValue(true);
    }

    virtual ~OgdenForceField()     {    }

};


}
}
}

#endif
