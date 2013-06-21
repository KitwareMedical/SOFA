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
#ifndef SOFA_COMPONENT_MAPPING_DEFORMABLEONRIGIDFRAME_H
#define SOFA_COMPONENT_MAPPING_DEFORMABLEONRIGIDFRAME_H

#include <sofa/core/Multi2Mapping.h>

#include <sofa/component/component.h>

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/topology/Topology.h>
#include <vector>



namespace sofa
{

namespace component
{

namespace mapping
{

/// This class can be overridden if needed for additionnal storage within template specializations.
template<class InDataTypes, class OutDataTypes>
class DeformableOnRigidFrameMappingInternalData
{
public:
};

template <class TIn, class TInRoot, class TOut>
class DeformableOnRigidFrameMapping : public core::Multi2Mapping<TIn, TInRoot, TOut>
{
 public:
    SOFA_CLASS(SOFA_TEMPLATE3(DeformableOnRigidFrameMapping, TIn, TInRoot, TOut), SOFA_TEMPLATE3(core::Multi2Mapping, TIn, TInRoot, TOut) );

    typedef core::Multi2Mapping<TIn, TInRoot, TOut> Inherit;

    typedef TIn In;
    typedef TInRoot InRoot;
    typedef TOut Out;

    typedef typename Out::VecCoord OutVecCoord;
    typedef typename Out::VecDeriv OutVecDeriv;
    typedef typename Out::MatrixDeriv OutMatrixDeriv;
    typedef typename Out::Coord OutCoord;
    typedef typename Out::Deriv OutDeriv;
    typedef typename OutCoord::value_type OutReal;
    typedef Data<OutVecCoord> OutDataVecCoord;
    typedef Data<OutVecDeriv> OutDataVecDeriv;
    typedef Data<OutMatrixDeriv> OutDataMatrixDeriv;

    typedef typename In::VecCoord InVecCoord;
    typedef typename In::VecDeriv InVecDeriv;
    typedef typename In::MatrixDeriv InMatrixDeriv;
    typedef typename In::Coord InCoord;
    typedef typename In::Deriv InDeriv;
    typedef typename In::Real InReal;
    typedef Data<InVecCoord> InDataVecCoord;
    typedef Data<InVecDeriv> InDataVecDeriv;
    typedef Data<InMatrixDeriv> InDataMatrixDeriv;

    typedef typename InRoot::VecCoord InRootVecCoord;
    typedef typename InRoot::VecDeriv InRootVecDeriv;
    typedef typename InRoot::MatrixDeriv InRootMatrixDeriv;
    typedef typename InRoot::Coord InRootCoord;
    typedef typename InRoot::Deriv InRootDeriv;
    typedef typename InRoot::Real InRootReal;
    typedef Data<InRootVecCoord> InRootDataVecCoord;
    typedef Data<InRootVecDeriv> InRootDataVecDeriv;
    typedef Data<InRootMatrixDeriv> InRootDataMatrixDeriv;

    typedef typename OutCoord::value_type Real;
    typedef OutCoord Coord;
    typedef OutDeriv Deriv;
    enum { N=Out::spatial_dimensions };
    typedef defaulttype::Mat<N,N,Real> Mat;
    typedef defaulttype::Vec<N,Real> Vector ;

    OutVecCoord rotatedPoints;
    DeformableOnRigidFrameMappingInternalData<In, Out> data;
    Data<unsigned int> index;
    Data< bool > indexFromEnd;
    Data<sofa::helper::vector<unsigned int> >  repartition;
    Data< bool > globalToLocalCoords;

    Data< Real > m_rootAngularForceScaleFactor;
    Data< Real > m_rootLinearForceScaleFactor;

    helper::ParticleMask* maskFrom;
    helper::ParticleMask* maskTo;

    int addPoint ( const OutCoord& c );
    int addPoint ( const OutCoord& c, int indexFrom );

    void init();

	void handleTopologyChange(core::topology::Topology* t);

    /// Return true if the destination model has the same topology as the source model.
    ///
    /// This is the case for mapping keeping a one-to-one correspondance between
    /// input and output DOFs (mostly identity or data-conversion mappings).
    virtual bool sameTopology() const { return true; }

    //Apply
    void apply( OutVecCoord& out, const InVecCoord& in, const InRootVecCoord* inroot  );
    void apply(
        const core::MechanicalParams* /* mparams */ /* PARAMS FIRST */, const helper::vector<OutDataVecCoord*>& dataVecOutPos,
        const helper::vector<const InDataVecCoord*>& dataVecInPos ,
        const helper::vector<const InRootDataVecCoord*>& dataVecInRootPos)
    {
        if(dataVecOutPos.empty() || dataVecInPos.empty())
            return;

        const InRootVecCoord* inroot = NULL;

        //We need only one input In model and input Root model (if present)
        OutVecCoord& out = *dataVecOutPos[0]->beginEdit();
        const InVecCoord& in = dataVecInPos[0]->getValue();

        if (!dataVecInRootPos.empty())
            inroot = &dataVecInRootPos[0]->getValue();

        apply(out, in, inroot);

        dataVecOutPos[0]->endEdit();
    }

    //ApplyJ
    void applyJ( OutVecDeriv& out, const InVecDeriv& in, const InRootVecDeriv* inroot );
    void applyJ(
        const core::MechanicalParams* /* mparams */ /* PARAMS FIRST */, const helper::vector< OutDataVecDeriv*>& dataVecOutVel,
        const helper::vector<const InDataVecDeriv*>& dataVecInVel,
        const helper::vector<const InRootDataVecDeriv*>& dataVecInRootVel)
    {
        if(dataVecOutVel.empty() || dataVecInVel.empty())
            return;

        const InRootVecDeriv* inroot = NULL;

        //We need only one input In model and input Root model (if present)
        OutVecDeriv& out = *dataVecOutVel[0]->beginEdit();
        const InVecDeriv& in = dataVecInVel[0]->getValue();

        if (!dataVecInRootVel.empty())
            inroot = &dataVecInRootVel[0]->getValue();

        applyJ(out,in, inroot);

        dataVecOutVel[0]->endEdit();
    }

    //ApplyJT Force
    void applyJT( InVecDeriv& out, const OutVecDeriv& in, InRootVecDeriv* outroot );
    void applyJT(
        const core::MechanicalParams* /* mparams */ /* PARAMS FIRST */, const helper::vector< InDataVecDeriv*>& dataVecOutForce,
        const helper::vector< InRootDataVecDeriv*>& dataVecOutRootForce,
        const helper::vector<const OutDataVecDeriv*>& dataVecInForce)
    {
        if(dataVecOutForce.empty() || dataVecInForce.empty())
            return;

        InRootVecDeriv* outroot = NULL;

        //We need only one input In model and input Root model (if present)
        InVecDeriv& out = *dataVecOutForce[0]->beginEdit();
        const OutVecDeriv& in = dataVecInForce[0]->getValue();

        if (!dataVecOutRootForce.empty())
            outroot = dataVecOutRootForce[0]->beginEdit();

        applyJT(out,in, outroot);

        dataVecOutForce[0]->endEdit();
        if (outroot != NULL)
            dataVecOutRootForce[0]->endEdit();
    }

    virtual void applyDJT(const core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, core::MultiVecDerivId /*inForce*/, core::ConstMultiVecDerivId /*outForce*/)
    {
        //serr<<"Warning ! DeformableOnRigidFrameMapping::applyDJT not implemented"<<sendl;
    }


    //ApplyJT Constraint
    void applyJT( InMatrixDeriv& out, const OutMatrixDeriv& in, InRootMatrixDeriv* outroot );
    void applyJT(
        const core::ConstraintParams* /* cparams */ /* PARAMS FIRST */, const helper::vector< InDataMatrixDeriv*>& dataMatOutConst ,
        const helper::vector< InRootDataMatrixDeriv*>&  dataMatOutRootConst ,
        const helper::vector<const OutDataMatrixDeriv*>& dataMatInConst)
    {
        if(dataMatOutConst.empty() || dataMatInConst.empty())
            return;

        InRootMatrixDeriv* outroot = NULL;

        //We need only one input In model and input Root model (if present)
        InMatrixDeriv& out = *dataMatOutConst[0]->beginEdit();
        const OutMatrixDeriv& in = dataMatInConst[0]->getValue();

        if (!dataMatOutRootConst.empty())
            outroot = dataMatOutRootConst[0]->beginEdit();

        applyJT(out,in, outroot);

        dataMatOutConst[0]->endEdit();
        if (outroot != NULL)
            dataMatOutRootConst[0]->endEdit();
    }

    /**
      * @brief
      MAP the mass: this function recompute the rigid mass (gravity center position and inertia) of the object
          based on its deformed shape
      */
    void recomputeRigidMass() {}

    //@}

    void draw(const core::visual::VisualParams* vparams);

    void clear ( int reserve=0 );

    void setRepartition ( unsigned int value );
    void setRepartition ( sofa::helper::vector<unsigned int> values );

protected:
    DeformableOnRigidFrameMapping();

    virtual ~DeformableOnRigidFrameMapping()
    {}

    core::State<In>* m_fromModel;
    core::State<Out>* m_toModel;
    core::State<InRoot>* m_fromRootModel;

    InRootCoord rootX;
};

#ifndef SOFA_FLOAT
using sofa::defaulttype::Vec3dTypes;
#endif

#ifndef SOFA_DOUBLE
using sofa::defaulttype::Vec3fTypes;
using sofa::defaulttype::ExtVec3fTypes;
#endif

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_MAPPING_DEFORMABLEONRIGIDFRAME_CPP)  //// ATTENTION PB COMPIL WIN3Z
#ifndef SOFA_FLOAT
extern template class SOFA_MISC_MAPPING_API DeformableOnRigidFrameMapping< Vec3dTypes, Rigid3dTypes, Vec3dTypes >;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_MISC_MAPPING_API DeformableOnRigidFrameMapping< Vec3fTypes, Rigid3fTypes, Vec3fTypes >;
#endif
#endif

} // namespace mapping

} // namespace component

} // namespace sofa

#endif
