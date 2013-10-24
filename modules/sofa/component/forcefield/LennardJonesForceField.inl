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
#ifndef SOFA_COMPONENT_FORCEFIELD_LENNARDJONESFORCEFIELD_INL
#define SOFA_COMPONENT_FORCEFIELD_LENNARDJONESFORCEFIELD_INL

#include <sofa/component/forcefield/LennardJonesForceField.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/simulation/common/Simulation.h>
#include <sofa/helper/system/config.h>
#include <sofa/helper/gl/template.h>
#include <math.h>
#include <iostream>




namespace sofa
{

namespace component
{

namespace forcefield
{


template<class DataTypes>
void LennardJonesForceField<DataTypes>::init()
{
    this->Inherit::init();

    assert( this->mstate );

    a = (p0.getValue() * (Real)pow(d0.getValue(),alpha.getValue())) / (1-alpha.getValue()/beta.getValue());
    b = (p0.getValue() * (Real)pow(d0.getValue(),beta.getValue())) / (beta.getValue()/alpha.getValue()-1);
    sout << "Lennard-Jones initialized: alpha="<<alpha.getValue()<<" beta="<<beta.getValue()<<" d0="<<d0.getValue()<<" p0="<<p0.getValue()<<" a="<<a<<" b="<<b<<sendl;
    // Validity check: compute force and potential at d0
    Real f0 = a*alpha.getValue()*(Real)pow(d0.getValue(),-alpha.getValue()-1)-b*beta.getValue()*(Real)pow(d0.getValue(),-beta.getValue()-1);
    if (fabs(f0)>0.001)
        serr << "Lennard-Jones initialization failed: f0="<<f0<<sendl;
    Real cp0 = (a*(Real)pow(d0.getValue(),-alpha.getValue())-b*(Real)pow(d0.getValue(),-beta.getValue()));
    if (fabs(cp0/p0.getValue()-1)>0.001)
        serr << "Lennard-Jones initialization failed: cp0="<<cp0<<sendl;
    // Debug
    for (Real d = 0; d<dmax.getValue(); d+= dmax.getValue()/60)
    {
        Real f = a*alpha.getValue()*(Real)pow(d,-alpha.getValue()-1)-b*beta.getValue()*(Real)pow(d,-beta.getValue()-1);
        sout << "f("<<d<<")="<<f<<sendl;
    }
}

template<class DataTypes>
void LennardJonesForceField<DataTypes>::addForce(const core::MechanicalParams* /* mparams */ /* PARAMS FIRST */, DataVecDeriv& d_f, const DataVecCoord& d_x, const DataVecDeriv& d_v)
{
    VecDeriv& f1 = *d_f.beginEdit();
    const VecCoord& p1 = d_x.getValue();
    const VecDeriv& v1 = d_v.getValue();

    Real dmax2 = dmax.getValue()*dmax.getValue();
    this->dforces.clear();
    f1.resize(p1.size());
    for (unsigned int ib=1; ib<p1.size(); ib++)
    {
        const Coord pb = p1[ib];
        for (unsigned int ia=0; ia<ib; ia++)
        {
            const Coord pa = p1[ia];
            const Deriv u = pb-pa;
            const Real d2 = u.norm2();
            if (d2 >= dmax2) continue;
            const Real d = (Real)sqrt(d2);
            const Real fa = a*alpha.getValue()*(Real)pow(d,-alpha.getValue()-1);
            const Real fb = b*beta.getValue()*(Real)pow(d,-beta.getValue()-1);
            Real forceIntensity = fa - fb;
            //sout << ia<<"-"<<ib<<" d="<<d<<" f="<<forceIntensity<<sendl;
            DForce df;
            df.a = ia;
            df.b = ib;
            if (forceIntensity > fmax.getValue())
            {
                forceIntensity = fmax.getValue();
                df.df = 0;
            }
            else
            {
                df.df = ((-alpha.getValue()-1)*fa - (-beta.getValue()-1)*fb)/(d*d2);
            }
            this->dforces.push_back(df);
            Deriv force = u*(forceIntensity/d);

            // Add damping
            force += (v1[ib]-v1[ia])*damping.getValue();

            f1[ia]+=force;
            f1[ib]-=force;
        }
    }
    d_f.endEdit();
}

template<class DataTypes>
void LennardJonesForceField<DataTypes>::addDForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& d_df, const DataVecDeriv& d_dx)
{
    VecDeriv& df1 = *d_df.beginEdit();
    const VecDeriv& dx1 = d_dx.getValue();
    Real kFactor = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());

    const VecCoord& p1 = *this->mstate->getX();
    df1.resize(dx1.size());
    for (unsigned int i=0; i<this->dforces.size(); i++)
    {
        const DForce& df = this->dforces[i];
        const unsigned int ia = df.a;
        const unsigned int ib = df.b;
        const Deriv u = p1[ib]-p1[ia];
        const Deriv du = dx1[ib]-dx1[ia];
        const Deriv dforce = (u * (df.df * (du*u))) * kFactor;
        df1[ia] += dforce;
        df1[ib] -= dforce;
    }
    d_df.endEdit();
}

template<class DataTypes>
void LennardJonesForceField<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    if (!vparams->displayFlags().getShowForceFields()) return;
    const VecCoord& p1 = *this->mstate->getX();

    std::vector< defaulttype::Vector3 > points[2];

    const Real d02 = this->d0.getValue()*this->d0.getValue();
    for (unsigned int i=0; i<this->dforces.size(); i++)
    {
        const DForce& df = this->dforces[i];
        if ((p1[df.b]-p1[df.a]).norm2() < d02)
        {
            points[0].push_back(p1[df.a]);
            points[0].push_back(p1[df.b]);
        }
        else
        {
            points[1].push_back(p1[df.a]);
            points[1].push_back(p1[df.b]);
        }
    }
    vparams->drawTool()->drawLines(points[0], 1, defaulttype::Vec<4,float>(1,1,1,1));
    vparams->drawTool()->drawLines(points[1], 1, defaulttype::Vec<4,float>(0,0,1,1));

}

} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FORCEFIELD_LENNARDJONESFORCEFIELD_INL
