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
#ifndef FLEXIBLE_PlasticStrainJacobianBlock_H
#define FLEXIBLE_PlasticStrainJacobianBlock_H

#include "../BaseJacobian.h"

#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/Mat.h>
#include "../types/StrainTypes.h"
#include "../helper.h"

namespace sofa
{

namespace defaulttype
{


/** Template class used to implement one jacobian block for PlasticStrainMapping */
template<class TStrain>
class PlasticStrainJacobianBlock : public BaseJacobianBlock<TStrain,TStrain>
{

public:

    typedef BaseJacobianBlock<TStrain,TStrain> Inherit;
    typedef typename Inherit::InCoord InCoord;
    typedef typename Inherit::InDeriv InDeriv;
    typedef typename Inherit::OutCoord OutCoord;
    typedef typename Inherit::OutDeriv OutDeriv;
    typedef typename Inherit::MatBlock MatBlock;
    typedef typename Inherit::KBlock KBlock;
    typedef typename Inherit::Real Real;

    typedef typename TStrain::StrainMat StrainMat;  ///< Matrix representing a strain
    typedef typename TStrain::StrainVec StrainVec;  ///< Vec representing a strain (Voigt notation)
    enum { strain_size = TStrain::strain_size };
    enum { order = TStrain::order };
    enum { spatial_dimensions = TStrain::spatial_dimensions };

    static const bool constantJ = true;


    /**
    Mapping:   ADDITION -> \f$ E_elastic = E_total - E_plastic \f ; MULTIPLICATION -> \f$ E_elastic = E_total * E_plastic^-1 \f$
    Jacobian:    \f$  dE = Id \f$
    */


    InCoord _plasticStrain;

    void reset()
    {
        _plasticStrain.clear();
    }


    void addapply( OutCoord& /*result*/, const InCoord& /*data*/ ) {}

    void addapply_multiplication( OutCoord& result, const InCoord& data, Real max, Real squaredYield, Real creep )
    {
        // eventually remove a part of the strain to simulate plasticity

        // could be optimized by storing the computation of the previous time step
        StrainMat plasticStrainMat = StrainVoigtToMat( _plasticStrain.getStrain() ) + StrainMat::Identity();
        StrainMat plasticStrainMatInverse; plasticStrainMatInverse.invert( plasticStrainMat );

        // elasticStrain = totalStrain * plasticStrain^-1
        StrainMat elasticStrainMat = ( StrainVoigtToMat( data.getStrain() ) + StrainMat::Identity() ) * plasticStrainMatInverse;
        StrainVec elasticStrainVec = StrainMatToVoigt( elasticStrainMat - StrainMat::Identity() );

        // if( ||elasticStrain||  > c_yield ) plasticStrain += dt * c_creep * dt * elasticStrain
        if( elasticStrainVec.norm2() > squaredYield )
            _plasticStrain.getStrain() += creep * elasticStrainVec;

        // if( ||plasticStrain|| > c_max ) plasticStrain *= c_max / ||plasticStrain||
        Real plasticStrainNorm2 = _plasticStrain.getStrain().norm2();
        if( plasticStrainNorm2 > max*max )
            _plasticStrain.getStrain() *= max / helper::rsqrt( plasticStrainNorm2 );

        plasticStrainMat = StrainVoigtToMat( _plasticStrain.getStrain() ) + StrainMat::Identity();

        // remaining elasticStrain = totalStrain * plasticStrain^-1
        plasticStrainMatInverse.invert( plasticStrainMat );
        elasticStrainMat = ( StrainVoigtToMat( data.getStrain() ) + StrainMat::Identity() ) * plasticStrainMatInverse;
        elasticStrainVec = StrainMatToVoigt( elasticStrainMat - StrainMat::Identity() );

        result.getStrain() += elasticStrainVec;
    }

    void addapply_addition( OutCoord& result, const InCoord& data, Real max, Real squaredYield, Real creep )
    {
        // eventually remove a part of the strain to simulate plasticity

        // elasticStrain = totalStrain - plasticStrain
        InCoord elasticStrain = data - _plasticStrain;

        // if( ||elasticStrain||  > c_yield ) plasticStrain += dt * c_creep * dt * elasticStrain
        if( elasticStrain.getStrain().norm2() > squaredYield )
            _plasticStrain += elasticStrain * creep;

        // if( ||plasticStrain|| > c_max ) plasticStrain *= c_max / ||plasticStrain||
        Real plasticStrainNorm2 = _plasticStrain.getStrain().norm2();
        if( plasticStrainNorm2 > max*max )
            _plasticStrain.getStrain() *= max / helper::rsqrt( plasticStrainNorm2 );

        // remaining elasticStrain = totatStrain - plasticStrain
        elasticStrain = data - _plasticStrain;

        result += elasticStrain;
    }

    void addmult( OutDeriv& result,const InDeriv& data )
    {
        result += data;
    }

    void addMultTranspose( InDeriv& result, const OutDeriv& data )
    {
        result += data;
    }

    MatBlock getJ()
    {
        return MatBlock::Identity();
    }

    // Not Yet implemented..
    KBlock getK(const OutDeriv& /*childForce*/)
    {
        return KBlock();
    }

    void addDForce( InDeriv& /*df*/, const InDeriv& /*dx*/, const OutDeriv& /*childForce*/, const double& /*kfactor */)
    {
    }

}; // class PlasticStrainJacobianBlock

} // namespace defaulttype
} // namespace sofa



#endif // FLEXIBLE_PlasticStrainJacobianBlock_H
