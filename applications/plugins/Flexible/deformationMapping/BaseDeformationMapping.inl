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
#ifndef SOFA_COMPONENT_MAPPING_BaseDeformationMAPPING_INL
#define SOFA_COMPONENT_MAPPING_BaseDeformationMAPPING_INL

#include "../deformationMapping/BaseDeformationMapping.h"
#include <sofa/component/visualmodel/VisualModelImpl.h>
#include "../quadrature/BaseGaussPointSampler.h"
#include <sofa/helper/gl/Color.h>

#include <omp.h>
#include <limits>

#include <Eigen/Core>
#include <Eigen/Dense>

namespace sofa
{

/** determinant for 1x1 matrix  to complement 2x2 and 3x3 implementations (used for visualization of  det F ) **/
namespace defaulttype { template<class real> inline real determinant(const Mat<1,1,real>& m) { return m(0,0);} }

namespace component
{
namespace mapping
{


template <class JacobianBlockType>
BaseDeformationMappingT<JacobianBlockType>::BaseDeformationMappingT (core::State<In>* from , core::State<Out>* to)
    : Inherit ( from, to )
    , _shapeFunction(NULL)
    , f_index ( initData ( &f_index,"indices","parent indices for each child" ) )
    , f_w ( initData ( &f_w,"weights","influence weights of the Dofs" ) )
    , f_dw ( initData ( &f_dw,"weightGradients","weight gradients" ) )
    , f_ddw ( initData ( &f_ddw,"weightHessians","weight Hessians" ) )
    , f_F0 ( initData ( &f_F0,"M","Transformations from material to 3d space (linear for now..)" ) )
    , assembleJ ( initData ( &assembleJ,false, "assembleJ","Assemble the Jacobian matrix or use optimized Jacobian/vector multiplications" ) )
    , assembleK ( initData ( &assembleK,false, "assembleK","Assemble the geometric stiffness matrix or use optimized Jacobian/vector multiplications" ) )
    , f_pos0 ( initData ( &f_pos0,"restPosition","initial spatial positions of children" ) )
    , missingInformationDirty(true)
    , KdTreeDirty(true)
    , maskFrom(NULL)
    , maskTo(NULL)
    , triangles(0)
    , extTriangles(0)
    , extvertPosIdx(0)
    , showDeformationGradientScale(initData(&showDeformationGradientScale, (float)0.0, "showDeformationGradientScale", "Scale for deformation gradient display"))
    , showDeformationGradientStyle ( initData ( &showDeformationGradientStyle,"showDeformationGradientStyle","Visualization style for deformation gradients" ) )
    , showColorOnTopology ( initData ( &showColorOnTopology,"showColorOnTopology","Color mapping method" ) )
    , showColorScale(initData(&showColorScale, (float)1.0, "showColorScale", "Color mapping scale"))
{
    helper::OptionsGroup methodOptions(3,"0 - None"
            ,"1 - trace(F^T.F)-3"
            ,"2 - sqrt(det(F^T.F))-1");
    methodOptions.setSelectedItem(0);
    showColorOnTopology.setValue(methodOptions);

    helper::OptionsGroup styleOptions(6,"0 - All axis"
            ,"1 - First axis"
            ,"2 - Second axis"
            ,"3 - Third axis"
            ,"4 - deformation"
            ,"5 - 1st piola stress" );
    styleOptions.setSelectedItem(0);
    showDeformationGradientStyle.setValue(styleOptions);
}



template <class JacobianBlockType>
void BaseDeformationMappingT<JacobianBlockType>::resizeOut()
{
    if(this->f_printLog.getValue()) std::cout<<this->getName()<<"::resizeOut()"<<std::endl;

    helper::ReadAccessor<Data<OutVecCoord> > out (*this->toModel->read(core::ConstVecCoordId::position()));

    helper::WriteAccessor<Data<VecCoord> > pos0 (this->f_pos0);
    this->missingInformationDirty=true; this->KdTreeDirty=true; // need to update mapped spatial positions if needed for visualization

    unsigned int size;

    engine::BaseGaussPointSampler* sampler;
    this->getContext()->get(sampler,core::objectmodel::BaseContext::Local);
    bool restPositionSet=false;
    helper::ReadAccessor<Data< OutVecCoord > >  rest(*this->toModel->read(core::ConstVecCoordId::restPosition()));

    if(sampler) // retrieve initial positions from gauss point sampler (deformation gradient types)
    {
        size = sampler->getNbSamples();
        if(rest.size()==size && size!=1) restPositionSet=true;
        this->toModel->resize(size);
        pos0.resize(size);  for(unsigned int i=0; i<size; i++) pos0[i]=sampler->getSample(i);
        if(this->f_printLog.getValue())  std::cout<<this->getName()<<" : "<< size <<" gauss points imported"<<std::endl;
    }
    else  // retrieve initial positions from children dofs (vec types)
    {
        size = out.size();
        pos0.resize(size);  for(unsigned int i=0; i<size; i++ )  Out::get(pos0[i][0],pos0[i][1],pos0[i][2],out[i]);
    }

    // init shape function
    if( !_shapeFunction ) this->getContext()->get(_shapeFunction,core::objectmodel::BaseContext::SearchUp);
    if ( !_shapeFunction ) serr << "ShapeFunction<"<<ShapeFunctionType::Name()<<"> component not found" << sendl;
    else
    {
        vector<mCoord> mpos0;
        mpos0.resize(pos0.size());
        for(unsigned int i=0; i<pos0.size(); ++i)  StdVectorTypes<mCoord,mCoord>::set( mpos0[i], pos0[i][0] , pos0[i][1] , pos0[i][2]);
        // interpolate weights at sample positions
        _shapeFunction->computeShapeFunction(mpos0,*this->f_F0.beginEdit(),*this->f_index.beginEdit(),*this->f_w.beginEdit(),*this->f_dw.beginEdit(),*this->f_ddw.beginEdit());
        this->f_index.endEdit();     this->f_F0.endEdit();    this->f_w.endEdit();        this->f_dw.endEdit();        this->f_ddw.endEdit();

        // use custom rest positions (to set material directions or set residual deformations)
        if(restPositionSet)
        {
            helper::WriteAccessor<Data< VMaterialToSpatial > >  F0(this->f_F0);
            for(unsigned int i=0; i<rest.size(); ++i) F0[i]=OutDataTypesInfo<Out>::getF(rest[i]);
            if(this->f_printLog.getValue())  std::cout<<this->getName()<<" : "<<rest.size()<<" rest positions imported "<<std::endl;
        }
    }

    // init jacobians
    initJacobianBlocks();

    // clear forces
    if(this->toModel->write(core::VecDerivId::force())) { helper::WriteAccessor<Data< OutVecDeriv > >  f(*this->toModel->write(core::VecDerivId::force())); for(unsigned int i=0;i<f.size();i++) f[i].clear(); }
    // clear velocities
    if(this->toModel->write(core::VecDerivId::velocity())) { helper::WriteAccessor<Data< OutVecDeriv > >  vel(*this->toModel->write(core::VecDerivId::velocity())); for(unsigned int i=0;i<vel.size();i++) vel[i].clear(); }

    reinit();
}



template <class JacobianBlockType>
void BaseDeformationMappingT<JacobianBlockType>::resizeOut(const vector<Coord>& position0, vector<vector<unsigned int> > index,vector<vector<Real> > w, vector<vector<Vec<spatial_dimensions,Real> > > dw, vector<vector<Mat<spatial_dimensions,spatial_dimensions,Real> > > ddw, vector<Mat<spatial_dimensions,spatial_dimensions,Real> > F0)
{
    if(this->f_printLog.getValue()) std::cout<<this->getName()<<"::resizeOut()"<<std::endl;

    helper::ReadAccessor<Data<OutVecCoord> > out (*this->toModel->read(core::ConstVecCoordId::position()));

    helper::WriteAccessor<Data<VecCoord> > pos0 (this->f_pos0);
    this->missingInformationDirty=true; this->KdTreeDirty=true; // need to update mapped spatial positions if needed for visualization

    unsigned int size = position0.size();

    // paste input values
    this->toModel->resize(size);
    pos0.resize(size);  for(unsigned int i=0; i<size; i++ )        pos0[i]=position0[i];

    helper::WriteAccessor<Data<vector<VRef> > > wa_index (this->f_index);   wa_index.resize(size);  for(unsigned int i=0; i<size; i++ )    wa_index[i].assign(index[i].begin(), index[i].end());
    helper::WriteAccessor<Data<vector<VReal> > > wa_w (this->f_w);          wa_w.resize(size);  for(unsigned int i=0; i<size; i++ )    wa_w[i].assign(w[i].begin(), w[i].end());
    helper::WriteAccessor<Data<vector<VGradient> > > wa_dw (this->f_dw);    wa_dw.resize(size);  for(unsigned int i=0; i<size; i++ )    wa_dw[i].assign(dw[i].begin(), dw[i].end());
    helper::WriteAccessor<Data<vector<VHessian> > > wa_ddw (this->f_ddw);   wa_ddw.resize(size);  for(unsigned int i=0; i<size; i++ )    wa_ddw[i].assign(ddw[i].begin(), ddw[i].end());
    helper::WriteAccessor<Data<VMaterialToSpatial> > wa_F0 (this->f_F0);    wa_F0.resize(size);  for(unsigned int i=0; i<size; i++ )    for(unsigned int j=0; j<spatial_dimensions; j++ ) for(unsigned int k=0; k<material_dimensions; k++ )   wa_F0[i][j][k]=F0[i][j][k];

    if(this->f_printLog.getValue())  std::cout<<this->getName()<<" : "<< size <<" custom gauss points imported"<<std::endl;

    // init jacobians
    initJacobianBlocks();

    // clear forces
    if(this->toModel->write(core::VecDerivId::force())) { helper::WriteAccessor<Data< OutVecDeriv > >  f(*this->toModel->write(core::VecDerivId::force())); for(unsigned int i=0;i<f.size();i++) f[i].clear(); }
    // clear velocities
    if(this->toModel->write(core::VecDerivId::velocity())) { helper::WriteAccessor<Data< OutVecDeriv > >  vel(*this->toModel->write(core::VecDerivId::velocity())); for(unsigned int i=0;i<vel.size();i++) vel[i].clear(); }

    reinit();
}


template <class JacobianBlockType>
void BaseDeformationMappingT<JacobianBlockType>::init()
{
    if (core::behavior::BaseMechanicalState* stateFrom = dynamic_cast<core::behavior::BaseMechanicalState*>(this->fromModel.get()))
        maskFrom = &stateFrom->forceMask;
    if (core::behavior::BaseMechanicalState* stateTo = dynamic_cast<core::behavior::BaseMechanicalState*>(this->toModel.get()))
        maskTo = &stateTo->forceMask;

    component::visualmodel::VisualModelImpl *visual;
    this->getContext()->get( visual, core::objectmodel::BaseContext::Local);
    if(visual) {this->extTriangles = &visual->getTriangles(); this->extvertPosIdx = &visual->m_vertPosIdx.getValue(); this->triangles=0; }
    else
    {
        core::topology::BaseMeshTopology *topo;
        this->getContext()->get( topo, core::objectmodel::BaseContext::Local);
        if(topo) {this->triangles = &topo->getTriangles();  this->extTriangles=0; }
    }


    baseMatrices.resize( 1 ); // just a wrapping for getJs()
    baseMatrices[0] = &eigenJacobian;

    resizeOut();

    Inherit::init();
}

template <class JacobianBlockType>
void BaseDeformationMappingT<JacobianBlockType>::reinit()
{
    if(this->assembleJ.getValue()) updateJ();

    apply(NULL, *this->toModel->write(core::VecCoordId::position()), *this->fromModel->read(core::ConstVecCoordId::position()));
    if(this->toModel->write(core::VecDerivId::velocity())) applyJ(NULL, *this->toModel->write(core::VecDerivId::velocity()), *this->fromModel->read(core::ConstVecDerivId::velocity()));

    Inherit::reinit();
}



template <class JacobianBlockType>
void BaseDeformationMappingT<JacobianBlockType>::updateJ()
{
    helper::ReadAccessor<Data<InVecCoord> > in (*this->fromModel->read(core::ConstVecCoordId::position()));
    //helper::ReadAccessor<Data<OutVecCoord> > out (*this->toModel->read(core::ConstVecCoordId::position()));
    eigenJacobian.resizeBlocks(jacobian.size(),in.size());
    for(unsigned int i=0; i<jacobian.size(); i++)
    {
        //            vector<MatBlock> blocks;
        //            vector<unsigned> columns;
        eigenJacobian.beginBlockRow(i);
        for(unsigned int j=0; j<jacobian[i].size(); j++)
        {
            //                columns.push_back( this->f_index.getValue()[i][j] );
            //                blocks.push_back( jacobian[i][j].getJ() );
            eigenJacobian.createBlock( this->f_index.getValue()[i][j], jacobian[i][j].getJ());
        }
        //            eigenJacobian.appendBlockRow( i, columns, blocks );
        eigenJacobian.endBlockRow();
    }
    //        eigenJacobian.endEdit();
    eigenJacobian.compress();
}

template <class JacobianBlockType>
void BaseDeformationMappingT<JacobianBlockType>::updateK(const OutVecDeriv& childForce)
{
    helper::ReadAccessor<Data<InVecCoord> > in (*this->fromModel->read(core::ConstVecCoordId::position()));
    K.resizeBlocks(in.size(),in.size());
    vector<KBlock> diagonalBlocks; diagonalBlocks.resize(in.size());

    for(unsigned int i=0; i<jacobian.size(); i++)
    {
        for(unsigned int j=0; j<jacobian[i].size(); j++)
        {
            unsigned int index=this->f_index.getValue()[i][j];
            diagonalBlocks[index] += jacobian[i][j].getK(childForce[i]);
        }
    }
    for(unsigned int i=0; i<in.size(); i++)
    {
        //            vector<KBlock> blocks;
        //            vector<unsigned> columns;
        //            columns.push_back( i );
        //            blocks.push_back( diagonalBlocks[i] );
        //            K.appendBlockRow( i, columns, blocks );
        K.beginBlockRow(i);
        K.createBlock(i,diagonalBlocks[i]);
        K.endBlockRow();
    }
    //        K.endEdit();
    K.compress();
}



template <class JacobianBlockType>
void BaseDeformationMappingT<JacobianBlockType>::apply(const core::MechanicalParams * /*mparams*/ , Data<OutVecCoord>& dOut, const Data<InVecCoord>& dIn)
{
    if(this->f_printLog.getValue()) std::cout<<this->getName()<<":apply"<<std::endl;

    helper::ReadAccessor<Data<OutVecCoord> > outpos (*this->toModel->read(core::ConstVecCoordId::position()));
//    if(_sampler) if(_sampler->getNbSamples()!=outpos.size()) resizeOut();

    OutVecCoord&  out = *dOut.beginEdit();
    const InVecCoord&  in = dIn.getValue();

#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
    for(unsigned int i=0; i<jacobian.size(); i++)
    {
        out[i]=OutCoord();
        for(unsigned int j=0; j<jacobian[i].size(); j++)
        {
            unsigned int index=this->f_index.getValue()[i][j];
            jacobian[i][j].addapply(out[i],in[index]);
        }
    }
    dOut.endEdit();

    if(!BlockType::constantJ) if(this->assembleJ.getValue()) updateJ();

    this->missingInformationDirty=true; this->KdTreeDirty=true; // need to update spatial positions of defo grads if needed for visualization
}



template <class JacobianBlockType>
void BaseDeformationMappingT<JacobianBlockType>::applyJ(const core::MechanicalParams * /*mparams*/ , Data<OutVecDeriv>& dOut, const Data<InVecDeriv>& dIn)
{
    if(this->assembleJ.getValue())  eigenJacobian.mult(dOut,dIn);
    else
    {
        OutVecDeriv&  out = *dOut.beginEdit();
        const InVecDeriv&  in = dIn.getValue();

        if ((!this->maskTo)||(this->maskTo&& !(this->maskTo->isInUse())) )
        {
#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
            for(unsigned int i=0; i<jacobian.size(); i++)
            {
                out[i]=OutDeriv();
                for(unsigned int j=0; j<jacobian[i].size(); j++)
                {
                    unsigned int index=this->f_index.getValue()[i][j];
                    jacobian[i][j].addmult(out[i],in[index]);
                }
            }
        }
        else
        {
            typedef helper::ParticleMask ParticleMask;
            const ParticleMask::InternalStorage &indices=this->maskTo->getEntries();
            for (ParticleMask::InternalStorage::const_iterator  it=indices.begin(); it!=indices.end(); it++ )
            {
                unsigned int i= ( unsigned int ) ( *it );
                out[i]=OutDeriv();
                for(unsigned int j=0; j<jacobian[i].size(); j++)
                {
                    unsigned int index=this->f_index.getValue()[i][j];
                    jacobian[i][j].addmult(out[i],in[index]);
                }
            }
        }

        dOut.endEdit();
    }
}

template <class JacobianBlockType>
void BaseDeformationMappingT<JacobianBlockType>::applyJT(const core::MechanicalParams * /*mparams*/ , Data<InVecDeriv>& dIn, const Data<OutVecDeriv>& dOut)
{
    if(this->assembleJ.getValue())  eigenJacobian.addMultTranspose(dIn,dOut);
    else
    {
        InVecDeriv&  in = *dIn.beginEdit();
        const OutVecDeriv&  out = dOut.getValue();

        if((!this->maskTo)||(this->maskTo&& !(this->maskTo->isInUse())) )
        {
            //#pragma omp parallel for
            for(unsigned int i=0; i<jacobian.size(); i++)
            {
                for(unsigned int j=0; j<jacobian[i].size(); j++)
                {
                    unsigned int index=this->f_index.getValue()[i][j];
                    jacobian[i][j].addMultTranspose(in[index],out[i]);
                }
            }
        }
        else
        {
            typedef helper::ParticleMask ParticleMask;
            const ParticleMask::InternalStorage &indices=this->maskTo->getEntries();
            for (ParticleMask::InternalStorage::const_iterator  it=indices.begin(); it!=indices.end(); it++ )
            {
                const int i= ( int ) ( *it );
                for(unsigned int j=0; j<jacobian[i].size(); j++)
                {
                    unsigned int index=this->f_index.getValue()[i][j];
                    jacobian[i][j].addMultTranspose(in[index],out[i]);
                }
            }
        }

        dIn.endEdit();
    }
}

template <class JacobianBlockType>
void BaseDeformationMappingT<JacobianBlockType>::applyDJT(const core::MechanicalParams* mparams, core::MultiVecDerivId parentDfId, core::ConstMultiVecDerivId )
{
    if(BlockType::constantJ) return;

    Data<InVecDeriv>& parentForceData = *parentDfId[this->fromModel.get(mparams)].write();
    const Data<InVecDeriv>& parentDisplacementData = *mparams->readDx(this->fromModel);
    const Data<OutVecDeriv>& childForceData = *mparams->readF(this->toModel);

    helper::WriteAccessor<Data<InVecDeriv> > parentForce (parentForceData);
    helper::ReadAccessor<Data<InVecDeriv> > parentDisplacement (parentDisplacementData);
    helper::ReadAccessor<Data<OutVecDeriv> > childForce (childForceData);

    if(this->assembleK.getValue())
    {
        updateK(childForce.ref());
        K.addMult(parentForceData,parentDisplacementData,mparams->kFactor());
    }
    else
    {
        if((!this->maskTo)||(this->maskTo&& !(this->maskTo->isInUse())) )
        {
            //#pragma omp parallel for
            for(unsigned int i=0; i<jacobian.size(); i++)
            {
                for(unsigned int j=0; j<jacobian[i].size(); j++)
                {
                    unsigned int index=this->f_index.getValue()[i][j];
                    jacobian[i][j].addDForce( parentForce[index], parentDisplacement[index], childForce[i], mparams->kFactor() );
                }
            }
        }
        else
        {
            typedef helper::ParticleMask ParticleMask;
            const ParticleMask::InternalStorage &indices=this->maskTo->getEntries();
            for (ParticleMask::InternalStorage::const_iterator  it=indices.begin(); it!=indices.end(); it++ )
            {
                const int i= ( int ) ( *it );
                for(unsigned int j=0; j<jacobian[i].size(); j++)
                {
                    unsigned int index=this->f_index.getValue()[i][j];
                    jacobian[i][j].addDForce( parentForce[index], parentDisplacement[index], childForce[i], mparams->kFactor() );
                }
            }
        }
    }
}



template <class JacobianBlockType>
void BaseDeformationMappingT<JacobianBlockType>::applyJT( const core::ConstraintParams * /*cparams*/, Data<InMatrixDeriv>& _out, const Data<OutMatrixDeriv>& _in )
{
    // TODO handle mask

    InMatrixDeriv& out = *_out.beginEdit();
    const OutMatrixDeriv& in = _in.getValue();

    typename OutMatrixDeriv::RowConstIterator rowItEnd = in.end();

    for (typename OutMatrixDeriv::RowConstIterator rowIt = in.begin(); rowIt != rowItEnd; ++rowIt)
    {
        typename OutMatrixDeriv::ColConstIterator colItEnd = rowIt.end();
        typename OutMatrixDeriv::ColConstIterator colIt = rowIt.begin();

        if (colIt != colItEnd)
        {
            typename InMatrixDeriv::RowIterator o = out.writeLine(rowIt.index());

            for ( ; colIt != colItEnd; ++colIt)
            {
                unsigned int indexIn = colIt.index();

                for(unsigned int j=0; j<jacobian[indexIn].size(); j++)
                {
                    unsigned int indexOut = this->f_index.getValue()[indexIn][j];

                    InDeriv tmp;
                    jacobian[indexIn][j].addMultTranspose( tmp, colIt.val() );

                    o.addCol( indexOut, tmp );
                }
            }
        }
    }

    _out.endEdit();
}




/** abstract implementation of BasePointMapper functions
    they call mapPosition/mapDefoGradient specialized functions (templated on specific jacobianBlockType)
**/

template <class JacobianBlockType>
void BaseDeformationMappingT<JacobianBlockType>::ForwardMapping(Coord& p,const Coord& p0)
{
    if ( !_shapeFunction ) return;

    // interpolate weights at sample positions
    mCoord mp0;        StdVectorTypes<mCoord,mCoord>::set( mp0, p0[0] , p0[1] , p0[2]);
    MaterialToSpatial M;  VRef ref; VReal w;
    _shapeFunction->computeShapeFunction(mp0,M,ref,w);

    // map using specific instanciation
    this->mapPosition(p,p0,ref,w);
}



/** inversion of rectangular deformation gradients (used in backward mapping) **/
template <int L,typename Real>
inline static void invert(defaulttype::Mat<L,L,Real> &Minv, const defaulttype::Mat<L,L,Real> &M)
{
//    Eigen::Map<const Eigen::Matrix<Real,L,L,Eigen::RowMajor> >  eM(&M[0][0]);
//    Eigen::Map<Eigen::Matrix<Real,L,L,Eigen::RowMajor> >  eMinv(&Minv[0][0]);
//    eMinv=eM.inverse();
    Minv.invert(M);
}

template <int L,typename Real>
inline static void invert(defaulttype::Mat<1,L,Real> &Minv, const defaulttype::Mat<L,1,Real> &M)
{
    Real n2inv=0; for(unsigned int i=0; i<L; i++) n2inv+=M[i][0]*M[i][0];
    n2inv=1./n2inv;
    for(unsigned int i=0; i<L; i++)  Minv[0][i]=M[i][0]*n2inv;
}

template <typename Real>
inline static void invert(defaulttype::Mat<2,3,Real> &Minv, const defaulttype::Mat<3,2,Real> &M)
{
    defaulttype::Vec<3,Real> u=M.transposed()[0],v=M.transposed()[1],w=cross(u,v);
    w.normalize();
    defaulttype::Mat<3,3,Real> Mc; for(unsigned int i=0; i<3; i++) {Mc[i][0]=M[i][0]; Mc[i][1]=M[i][1]; Mc[i][2]=w[i];}
    defaulttype::Mat<3,3,Real> Mcinv; invert(Mcinv,Mc);
    for(unsigned int i=0; i<2; i++) for(unsigned int j=0; j<3; j++) Minv[i][j]=Mcinv[i][0];
}



template <class JacobianBlockType>
void BaseDeformationMappingT<JacobianBlockType>::BackwardMapping(Coord& p0,const Coord& p,const Real Thresh, const unsigned int NbMaxIt)
{
    if ( !_shapeFunction ) return;

    // iterate: p0(n+1) = F0.F^-1 (p-p(n)) + p0(n)
    unsigned int count=0;
    mCoord mp0;
    MaterialToSpatial F0;  VRef ref; VReal w; VGradient dw;
    Coord pnew;
    MaterialToSpatial F;
    Mat<material_dimensions,spatial_dimensions,Real> Finv;
    Mat<spatial_dimensions,spatial_dimensions,Real> F0Finv;

    while(count<NbMaxIt)
    {
        StdVectorTypes<mCoord,mCoord>::set( mp0, p0[0] , p0[1] , p0[2]);
        _shapeFunction->computeShapeFunction(mp0,F0,ref,w,&dw);
        if(!w[0]) { p0=Coord(); return; } // outside object

        this->mapPosition(pnew,p0,ref,w);
        if((p-pnew).norm2()<Thresh) return; // has converged
        this->mapDeformationGradient(F,p0,F0,ref,w,dw);

        invert(Finv,F);
        F0Finv=F0*Finv;
        p0+=F0Finv*(p-pnew);
        count++;
    }
}


template <class JacobianBlockType>
unsigned int BaseDeformationMappingT<JacobianBlockType>::getClosestMappedPoint(const Coord& p, Coord& x0,Coord& x, bool useKdTree)
{
    helper::ReadAccessor<Data<OutVecCoord> > out (*this->toModel->read(core::ConstVecCoordId::position()));
    if(this->missingInformationDirty)
    {
        if(!OutDataTypesInfo<Out>::positionMapped) mapPositions();
        this->missingInformationDirty=false;
    }

    unsigned int index=0;
    if(useKdTree)
    {
        if(this->KdTreeDirty)
        {
            if(OutDataTypesInfo<Out>::positionMapped) { f_pos.resize(out.size()); for(unsigned int i=0; i<out.size(); i++ )  Out::get(f_pos[i][0],f_pos[i][1],f_pos[i][2],out[i]);  } // copy to f_pos
            this->f_KdTree.build(f_pos);
            this->KdTreeDirty=false;
        }
        index=this->f_KdTree.getClosest(p);
        x=f_pos[index];
    }
    else
    {
        Real dmin=std::numeric_limits<Real>::max();
        for(unsigned int i=0; i<out.size(); i++)
        {
            Coord P; if(OutDataTypesInfo<Out>::positionMapped) Out::get(P[0],P[1],P[2],out[i]); else P=f_pos[i];
            Real d=(p-P).norm2();
            if(d<dmin) {dmin=d; index=i; x=P;}
        }
    }
    x0=f_pos0.getValue()[index];
    return index;
}

template<int matdim,typename Real>
void drawEllipsoid(const Mat<3,matdim,Real> & F, const Vec<3,Real> &p, const float& scale)
{
    glPushMatrix();

    GLdouble transformMatrix[16];
    for(int i=0; i<3; i++) for(int j=0; j<matdim; j++) transformMatrix[4*j+i] = (double)F(i,j)*scale;

    if(matdim==1)
    {
        for(int i=0; i<3; i++) for(int j=1; j<3; j++) transformMatrix[4*j+i] = 0;
    }
    else if(matdim==2)
    {
        Vec<3,Real> w=cross(F.transposed()[0],F.transposed()[1]); w.normalize();
        for(int i=0; i<3; i++)  transformMatrix[8+i]=(double)w[i]*scale*0.01; // arbitrarily small thickness
    }

    for(int i=0; i<3; i++)  transformMatrix[i+12]=p[i];
    for(int i=0; i<3; i++)  transformMatrix[4*i+3]=0; transformMatrix[15] = 1;
    glMultMatrixd(transformMatrix);

    GLUquadricObj* ellipsoid = gluNewQuadric();
    gluSphere(ellipsoid, 1.0, 10, 10);
    gluDeleteQuadric(ellipsoid);

    glPopMatrix();
}



template <class JacobianBlockType>
void BaseDeformationMappingT<JacobianBlockType>::draw(const core::visual::VisualParams* vparams)
{
    if (!vparams->displayFlags().getShowMechanicalMappings() && !showDeformationGradientScale.getValue() && showColorOnTopology.getValue().getSelectedId()==0) return;

    helper::ReadAccessor<Data<InVecCoord> > in (*this->fromModel->read(core::ConstVecCoordId::position()));
    helper::ReadAccessor<Data<OutVecCoord> > out (*this->toModel->read(core::ConstVecCoordId::position()));
    helper::ReadAccessor<Data<OutVecDeriv> > outf (*this->toModel->read(core::ConstVecDerivId::force()));
    helper::ReadAccessor<Data<vector<VRef> > > ref (this->f_index);
    helper::ReadAccessor<Data<vector<VReal> > > w (this->f_w);

    if(this->missingInformationDirty)
    {
        if(!OutDataTypesInfo<Out>::positionMapped) mapPositions();
        if(!OutDataTypesInfo<Out>::FMapped) if(showDeformationGradientScale.getValue() || showColorOnTopology.getValue().getSelectedId()!=0) mapDeformationGradients();
        this->missingInformationDirty=false;
    }

    if (vparams->displayFlags().getShowMechanicalMappings())
    {
        vector< defaulttype::Vec3d > edge;     edge.resize(2);
        Vec<4,float> col;

        for(unsigned i=0; i<out.size(); i++ )
        {
            if(OutDataTypesInfo<Out>::positionMapped) Out::get(edge[1][0],edge[1][1],edge[1][2],out[i]);
            else edge[1]=f_pos[i];
            for(unsigned j=0; j<ref[i].size(); j++ )
                if(w[i][j])
                {
                    In::get(edge[0][0],edge[0][1],edge[0][2],in[ref[i][j]]);
                    sofa::helper::gl::Color::getHSVA(&col[0],240.*w[i][j],1.,.8,1.);
                    vparams->drawTool()->drawLines ( edge, 1, col );
                }
        }
    }
    if (showDeformationGradientScale.getValue())
    {
        //                glPushAttrib ( GL_LIGHTING_BIT );
        //                glDisable ( GL_LIGHTING );
        float scale=showDeformationGradientScale.getValue();
        Vec<4,float> col( 0.5, 0.5, 0.0, 1.0 );
        MaterialToSpatial F;
        Coord p;
        for(unsigned i=0; i<out.size(); i++ )
        {
            if(OutDataTypesInfo<Out>::FMapped) F=OutDataTypesInfo<Out>::getF(out[i]); else F=f_F[i];
            if(OutDataTypesInfo<Out>::positionMapped) Out::get(p[0],p[1],p[2],out[i]); else p=f_pos[i];

            if(showDeformationGradientStyle.getValue().getSelectedId()==0)
            for(int j=0; j<material_dimensions; j++)
            {
                Coord u=F.transposed()(j)*0.5*scale;
                vparams->drawTool()->drawCylinder(p-u,p+u,0.05*scale,col,3);
            }
            else if(showDeformationGradientStyle.getValue().getSelectedId()==1)
                {
                    Coord u=F.transposed()(0)*0.5*scale;
                    vparams->drawTool()->drawCylinder(p-u,p+u,0.05*scale,col,3);
                }
            else if(showDeformationGradientStyle.getValue().getSelectedId()==2)
                {
                    Coord u=F.transposed()(1)*0.5*scale;
                    vparams->drawTool()->drawCylinder(p-u,p+u,0.05*scale,col,3);
                }
            else if(showDeformationGradientStyle.getValue().getSelectedId()==3)
                {
                    Coord u=F.transposed()(2)*0.5*scale;
                    vparams->drawTool()->drawCylinder(p-u,p+u,0.05*scale,col,3);
                }
            else if(showDeformationGradientStyle.getValue().getSelectedId()==4) // strain
                {
                    vparams->drawTool()->setMaterial(col);
                    drawEllipsoid(F,p,0.5*scale);
                }
            else if(showDeformationGradientStyle.getValue().getSelectedId()==5) // stress
                if(OutDataTypesInfo<Out>::FMapped)
                {
                    F=OutDataTypesInfo<Out>::getF(outf[i]);
                    vparams->drawTool()->setMaterial(col);
                    drawEllipsoid(F,p,0.5*scale);
                }

        }
        //                glPopAttrib();
    }

    if(showColorOnTopology.getValue().getSelectedId() && (this->extTriangles || this->triangles))
    {
        std::vector< Real > val(out.size());
        MaterialToSpatial F;
        for(unsigned i=0; i<out.size(); i++ )
        {
            if(OutDataTypesInfo<Out>::FMapped) F=OutDataTypesInfo<Out>::getF(out[i]); else F=f_F[i];

            if(showColorOnTopology.getValue().getSelectedId()==1) val[i]=(defaulttype::trace(F.transposed()*F)-3.);
            else  val[i]=sqrt(defaulttype::determinant(F.transposed()*F))-1.;

            //if (val[i]<0) val[i]=2*val[i]/(val[i]+1.);
            val[i]*=240 * this->showColorScale.getValue();
            val[i]+=120;
            if (val[i]<0) val[i]=0;
            if (val[i]>240) val[i]=240;
        }

        int nb =0;
        if(triangles) nb+=triangles->size();
        if(extTriangles) nb+=extTriangles->size();

        std::vector< defaulttype::Vector3 > points(3*nb),normals;
        std::vector< Vec<4,float> > colors(3*nb);
        unsigned int count=0;

        if(triangles)
            for ( unsigned int i = 0; i < triangles->size(); i++)
                for ( unsigned int j = 0; j < 3; j++)
                {
                    const unsigned int index = (*triangles)[i][j];
                    if(OutDataTypesInfo<Out>::positionMapped) Out::get(points[count][0],points[count][1],points[count][2],out[index]); else points[count]=f_pos[index];
                    sofa::helper::gl::Color::getHSVA(&colors[count][0],val[index],1.,.8,1.);
                    count++;
                }
        if(extTriangles)
            for ( unsigned int i = 0; i < extTriangles->size(); i++)
                for ( unsigned int j = 0; j < 3; j++)
                {
                    unsigned int index = (*extTriangles)[i][j];
                    if(this->extvertPosIdx) index=(*extvertPosIdx)[index];
                    if(OutDataTypesInfo<Out>::positionMapped) Out::get(points[count][0],points[count][1],points[count][2],out[index]); else points[count]=f_pos[index];
                    sofa::helper::gl::Color::getHSVA(&colors[count][0],val[index],1.,.8,1.);
                    count++;
                }
        glPushAttrib( GL_LIGHTING_BIT );
        glDisable( GL_LIGHTING);
        vparams->drawTool()->drawTriangles(points, normals, colors);
        glPopAttrib();
    }
}




} // namespace mapping
} // namespace component
} // namespace sofa

#endif
