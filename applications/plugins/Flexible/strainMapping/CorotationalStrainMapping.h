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
#ifndef SOFA_COMPONENT_MAPPING_CorotationalStrainMAPPING_H
#define SOFA_COMPONENT_MAPPING_CorotationalStrainMAPPING_H

#include "../initFlexible.h"
#include "../strainMapping/BaseStrainMapping.h"
#include "../strainMapping/CorotationalStrainJacobianBlock.inl"

#include <sofa/helper/OptionsGroup.h>

namespace sofa
{
namespace component
{
namespace mapping
{

using helper::vector;

/** Deformation Gradient to Corotational Lagrangian Strain mapping
*/

template <class TIn, class TOut>
class CorotationalStrainMapping : public BaseStrainMappingT<defaulttype::CorotationalStrainJacobianBlock<TIn,TOut> >
{
public:
    typedef defaulttype::CorotationalStrainJacobianBlock<TIn,TOut> BlockType;
    typedef BaseStrainMappingT<BlockType > Inherit;

    SOFA_CLASS(SOFA_TEMPLATE2(CorotationalStrainMapping,TIn,TOut), SOFA_TEMPLATE(BaseStrainMappingT,BlockType ));

    /** @name  Corotational methods */
    //@{
    enum DecompositionMethod { POLAR=0, QR, SMALL, SVD, NB_DecompositionMethod };
    Data<helper::OptionsGroup> f_method;
    //@}


    Data<bool> f_geometricStiffness; ///< should geometricStiffness be considered?

    virtual void reinit()
    {
        Inherit::reinit();

        switch( f_method.getValue().getSelectedId() )
        {
        case SMALL:
        {
            for( size_t i=0 ; i<this->jacobian.size() ; i++ )
            {
                this->jacobian[i].init_small();
            }
            break;
        }
        case QR:
        {
            for( size_t i=0 ; i<this->jacobian.size() ; i++ )
            {
                this->jacobian[i].init_qr( f_geometricStiffness.getValue() );
            }
            break;
        }
        case POLAR:
        {
            for( size_t i=0 ; i<this->jacobian.size() ; i++ )
            {
                this->jacobian[i].init_polar( f_geometricStiffness.getValue() );
            }
            break;
        }
        case SVD:
        {
            for( size_t i=0 ; i<this->jacobian.size() ; i++ )
            {
                this->jacobian[i].init_svd( f_geometricStiffness.getValue() );
            }
            break;
        }
        }
    }


protected:
    CorotationalStrainMapping (core::State<TIn>* from = NULL, core::State<TOut>* to= NULL)
        : Inherit ( from, to )
        , f_method( initData( &f_method, "method", "Decomposition method" ) )
        , f_geometricStiffness( initData( &f_geometricStiffness, false, "geometricStiffness", "Should geometricStiffness be considered?" ) )
    {
        helper::OptionsGroup Options;
        Options.setNbItems( NB_DecompositionMethod );
        Options.setItemName( SMALL, "small" );
        Options.setItemName( QR,    "qr"    );
        Options.setItemName( POLAR, "polar" );
        Options.setItemName( SVD,   "svd"   );
        Options.setSelectedItem( SVD );
        f_method.setValue( Options );
    }

    virtual ~CorotationalStrainMapping() { }

    virtual void apply( const core::MechanicalParams * /*mparams*/ , Data<typename Inherit::OutVecCoord>& dOut, const Data<typename Inherit::InVecCoord>& dIn )
    {
        if(this->f_printLog.getValue()) std::cout<<this->getName()<<":apply"<<std::endl;

        helper::ReadAccessor<Data<typename Inherit::InVecCoord> > inpos (*this->fromModel->read(core::ConstVecCoordId::position()));
        helper::ReadAccessor<Data<typename Inherit::OutVecCoord> > outpos (*this->toModel->read(core::ConstVecCoordId::position()));
        if(inpos.size()!=outpos.size()) this->resizeOut();

        typename Inherit::OutVecCoord& out = *dOut.beginEdit();
        const typename Inherit::InVecCoord&  in  =  dIn.getValue();

        switch( f_method.getValue().getSelectedId() )
        {
        case SMALL:
        {
#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
            for( int i=0 ; i < static_cast<int>(this->jacobian.size()) ; i++ )
            {
                out[i] = typename Inherit::OutCoord();
                this->jacobian[i].addapply_small( out[i], in[i] );
            }
            break;
        }
        case QR:
        {
#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
            for( int i=0 ; i < static_cast<int>(this->jacobian.size()) ; i++ )
            {
                out[i] = typename Inherit::OutCoord();
                this->jacobian[i].addapply_qr( out[i], in[i] );
            }
            break;
        }
        case POLAR:
        {
#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
            for( int i=0 ; i < static_cast<int>(this->jacobian.size()) ; i++ )
            {
                out[i] = typename Inherit::OutCoord();
                this->jacobian[i].addapply_polar( out[i], in[i] );
            }
            break;
        }
        case SVD:
        {
#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
            for( int i=0 ; i < static_cast<int>(this->jacobian.size()) ; i++ )
            {
                out[i] = typename Inherit::OutCoord();
                this->jacobian[i].addapply_svd( out[i], in[i] );
            }
            break;
        }
        }

        dOut.endEdit();

        /*if(!BlockType::constant)*/ if(this->assemble.getValue()) this->updateJ();
    }

    virtual void applyDJT(const core::MechanicalParams* mparams, core::MultiVecDerivId parentDfId, core::ConstMultiVecDerivId )
    {
        if( !f_geometricStiffness.getValue() ) return;
        if(BlockType::constant) return;

        Data<typename Inherit::InVecDeriv>& parentForceData = *parentDfId[this->fromModel.get(mparams)].write();
        const Data<typename Inherit::InVecDeriv>& parentDisplacementData = *mparams->readDx(this->fromModel);
        const Data<typename Inherit::OutVecDeriv>& childForceData = *mparams->readF(this->toModel);

        helper::WriteAccessor<Data<typename Inherit::InVecDeriv> > parentForce (parentForceData);
        helper::ReadAccessor<Data<typename Inherit::InVecDeriv> > parentDisplacement (parentDisplacementData);
        helper::ReadAccessor<Data<typename Inherit::OutVecDeriv> > childForce (childForceData);

        if(this->assemble.getValue())
        {
            this->updateK(childForce.ref());
            this->K.addMult(parentForceData,parentDisplacementData,mparams->kFactor());
        }
        else
        {
            switch( f_method.getValue().getSelectedId() )
            {
            case SMALL:
            {
                break;
            }
            case QR:
            {
#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
                for( int i=0 ; i < static_cast<int>(this->jacobian.size()) ; i++ )
                {
                    this->jacobian[i].addDForce_qr( parentForce[i], parentDisplacement[i], childForce[i], mparams->kFactor() );
                }
                break;
            }
            case POLAR:
            {
#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
                for( int i=0 ; i < static_cast<int>(this->jacobian.size()) ; i++ )
                {
                    this->jacobian[i].addDForce_polar( parentForce[i], parentDisplacement[i], childForce[i], mparams->kFactor() );
                }
                break;
            }
            case SVD:
            {
#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
                for( int i=0 ; i < static_cast<int>(this->jacobian.size()) ; i++ )
                {
                    this->jacobian[i].addDForce_svd( parentForce[i], parentDisplacement[i], childForce[i], mparams->kFactor() );
                }
                break;
            }
            }
        }
    }

};


} // namespace mapping
} // namespace component
} // namespace sofa

#endif
