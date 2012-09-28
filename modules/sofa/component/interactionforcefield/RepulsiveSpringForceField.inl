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
#ifndef SOFA_COMPONENT_INTERACTIONFORCEFIELD_REPULSIVESPRINGFORCEFIELD_INL
#define SOFA_COMPONENT_INTERACTIONFORCEFIELD_REPULSIVESPRINGFORCEFIELD_INL

#include <sofa/component/interactionforcefield/RepulsiveSpringForceField.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/component/interactionforcefield/StiffSpringForceField.inl>

namespace sofa
{

namespace component
{

namespace interactionforcefield
{



template<class DataTypes>
void RepulsiveSpringForceField<DataTypes>::addForce(const MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecDeriv& data_f1, DataVecDeriv& data_f2, const DataVecCoord& data_x1, const DataVecCoord& data_x2, const DataVecDeriv& data_v1, const DataVecDeriv& data_v2 )
{

    VecDeriv&       f1 = *data_f1.beginEdit();
    const VecCoord& x1 =  data_x1.getValue();
    const VecDeriv& v1 =  data_v1.getValue();
    VecDeriv&       f2 = *data_f2.beginEdit();
    const VecCoord& x2 =  data_x2.getValue();
    const VecDeriv& v2 =  data_v2.getValue();

    const helper::vector<Spring>& springs= this->springs.getValue();
    this->dfdx.resize(springs.size());
    f1.resize(x1.size());
    f2.resize(x2.size());
    for (unsigned int i=0; i<springs.size(); i++)
    {
#if 1
        int a = springs[i].m1;
        int b = springs[i].m2;
        Coord u = x2[b]-x1[a];
        Real d = u.norm();
        if (d < springs[i].initpos)
        {
            Real inverseLength = 1.0f/d;
            u *= inverseLength;
            Real elongation = (Real)(d - springs[i].initpos);
            Deriv relativeVelocity = v2[b]-v1[a];
            Real elongationVelocity = dot(u,relativeVelocity);
            Real forceIntensity = (Real)(springs[i].ks*elongation+springs[i].kd*elongationVelocity);
            Deriv force = u*forceIntensity;
            f1[a]+=force;
            f2[b]-=force;
            if (this->maskInUse)
            {
                this->mstate1->forceMask.insertEntry(a);
                this->mstate2->forceMask.insertEntry(b);
            }

            Mat& m = this->dfdx[i];
            Real tgt = forceIntensity * inverseLength;
            for( int j=0; j<N; ++j )
            {
                for( int k=0; k<N; ++k )
                {
                    m[j][k] = ((Real)springs[i].ks-tgt) * u[j] * u[k];
                }
                m[j][j] += tgt;
            }
        }
        else
#endif
        {
            Mat& m = this->dfdx[i];
            for( int j=0; j<N; ++j )
                for( int k=0; k<N; ++k )
                    m[j][k] = 0.0;
        }
    }

    data_f1.endEdit();
    data_f2.endEdit();
}

template <class DataTypes>
double RepulsiveSpringForceField<DataTypes>::getPotentialEnergy(const core::MechanicalParams* /* PARAMS FIRST */, const DataVecCoord&, const DataVecCoord& ) const
{
    serr<<"RepulsiveSpringForceField::getPotentialEnergy-not-implemented !!!"<<sendl;
    return 0;
}


} // namespace interactionforcefield

} // namespace component

} // namespace sofa

#endif
