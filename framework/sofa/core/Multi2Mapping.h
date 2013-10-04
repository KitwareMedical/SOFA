/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *
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
*                              SOFA :: Framework                              *
*                                                                             *
* Authors: The SOFA Team (see Authors.txt)                                    *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_CORE_MULTI2MAPPING_H
#define SOFA_CORE_MULTI2MAPPING_H

#include <sofa/core/BaseMapping.h>
#include <sofa/core/core.h>
#include <sofa/core/VecId.h>


namespace sofa
{

namespace core
{

/**
 *  \brief Specialized interface to describe many to many mapping.
 *   The inputs can be of two different types, while all the outputs must be of the same type.
 */

template <class TIn1, class TIn2, class TOut>
class Multi2Mapping : public BaseMapping
{
public:
    SOFA_ABSTRACT_CLASS(SOFA_TEMPLATE3(Multi2Mapping,TIn1, TIn2,TOut), BaseMapping);

    /// Input Model Type
    typedef TIn1 In1;
    typedef TIn2 In2;
    /// Output Model Type
    typedef TOut Out;

    typedef typename In1::VecCoord In1VecCoord;
    typedef typename In1::VecDeriv In1VecDeriv;
    typedef typename In1::MatrixDeriv In1MatrixDeriv;
    typedef Data<In1VecCoord> In1DataVecCoord;
    typedef Data<In1VecDeriv> In1DataVecDeriv;
    typedef Data<In1MatrixDeriv> In1DataMatrixDeriv;
    typedef typename In2::VecCoord In2VecCoord;
    typedef typename In2::VecDeriv In2VecDeriv;
    typedef typename In2::MatrixDeriv In2MatrixDeriv;
    typedef Data<In2VecCoord> In2DataVecCoord;
    typedef Data<In2VecDeriv> In2DataVecDeriv;
    typedef Data<In2MatrixDeriv> In2DataMatrixDeriv;

    typedef typename Out::VecCoord OutVecCoord;
    typedef typename Out::VecDeriv OutVecDeriv;
    typedef typename Out::MatrixDeriv OutMatrixDeriv;
    typedef Data<OutVecCoord> OutDataVecCoord;
    typedef Data<OutVecDeriv> OutDataVecDeriv;
    typedef Data<OutMatrixDeriv> OutDataMatrixDeriv;

    typedef MultiLink<Multi2Mapping<In1,In2,Out>, State< In1 >, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> LinkFromModels1;
    typedef typename LinkFromModels1::Container VecFromModels1;
    typedef MultiLink<Multi2Mapping<In1,In2,Out>, State< In2 >, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> LinkFromModels2;
    typedef typename LinkFromModels2::Container VecFromModels2;
    typedef MultiLink<Multi2Mapping<In1,In2,Out>, State< Out >, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> LinkToModels;
    typedef typename LinkToModels::Container VecToModels;

protected:
    /// Input Models container. New inputs are added through addInputModel(In* ).
    LinkFromModels1 fromModels1;
    LinkFromModels2 fromModels2;
    LinkToModels toModels;

public:

    Data<bool> f_applyRestPosition; ///< @todo document this
protected:

    /// Constructor
    Multi2Mapping();
    /// Destructor
    virtual ~Multi2Mapping() {};
public:

    virtual void addInputModel1(State<In1>*, const std::string& path = "");
    virtual void addInputModel2(State<In2>*, const std::string& path = "");
    virtual void addOutputModel(State<Out>*, const std::string& path = "");

    /// Return the reference to fromModels (In1).
    const VecFromModels1& getFromModels1();
    /// Return the reference to fromModels (In2).
    const VecFromModels2& getFromModels2();
    /// Return reference to toModels.
    const VecToModels& getToModels();

    /// Return a container of input models statically casted as BaseObject*
    helper::vector<BaseState*> getFrom();
    /// Return container of output model statically casted as BaseObject*.
    helper::vector<BaseState*> getTo();

    /// Get the source (upper) model.
    virtual helper::vector<behavior::BaseMechanicalState*> getMechFrom();

    /// Get the destination (lower, mapped) model.
    virtual helper::vector<behavior::BaseMechanicalState*> getMechTo();

    /// Apply ///
    /// Apply the mapping to position vectors.
    ///
    /// If the Mapping can be represented as a matrix J, this method computes
    /// $ out = J in $
    virtual void apply (const MechanicalParams* mparams /* PARAMS FIRST  = MechanicalParams::defaultInstance()*/, MultiVecCoordId outPos, ConstMultiVecCoordId inPos )
    {
        helper::vector<OutDataVecCoord*> vecOutPos;
        getVecOutCoord(outPos, vecOutPos);
        helper::vector<const In1DataVecCoord*> vecIn1Pos;
        getConstVecIn1Coord(inPos, vecIn1Pos);
        helper::vector<const In2DataVecCoord*> vecIn2Pos;
        getConstVecIn2Coord(inPos, vecIn2Pos);

        this->apply(mparams /* PARAMS FIRST */, vecOutPos, vecIn1Pos, vecIn2Pos);
    }
    /// This method must be reimplemented by all mappings.
    /// InPos and OutPos by default contains VecIds of type V_COORD.
    /// The size of InPos vector is the same as the number of fromModels.
    /// The size of OutPos vector is the same as the number of OutModels.
    virtual void apply(
        const MechanicalParams* mparams /* PARAMS FIRST */, const helper::vector<OutDataVecCoord*>& dataVecOutPos,
        const helper::vector<const In1DataVecCoord*>& dataVecIn1Pos ,
        const helper::vector<const In2DataVecCoord*>& dataVecIn2Pos)
#ifdef SOFA_DEPRECATE_OLD_API
        = 0;
#else
    {
        //Not optimized at all...
        helper::vector<OutVecCoord*> vecOutPos;
        for(unsigned int i=0; i<dataVecOutPos.size(); i++)
            vecOutPos.push_back(dataVecOutPos[i]->beginEdit(mparams));

        helper::vector<const In1VecCoord*> vecIn1Pos;
        for(unsigned int i=0; i<dataVecIn1Pos.size(); i++)
            vecIn1Pos.push_back(&dataVecIn1Pos[i]->getValue(mparams));
        helper::vector<const In2VecCoord*> vecIn2Pos;
        for(unsigned int i=0; i<dataVecIn2Pos.size(); i++)
            vecIn2Pos.push_back(&dataVecIn2Pos[i]->getValue(mparams));

        this->apply(vecOutPos, vecIn1Pos, vecIn2Pos);

        //Really Not optimized at all...
        for(unsigned int i=0; i<dataVecOutPos.size(); i++)
            dataVecOutPos[i]->endEdit(mparams);
    }
    /// Compat Method
    /// @deprecated
    virtual void apply(const helper::vector<OutVecCoord*>& /* outPos */,
            const helper::vector<const In1VecCoord*>& /* inPos1 */,
            const helper::vector<const In2VecCoord*>& /* inPos2 */) { };
#endif //SOFA_DEPRECATE_OLD_API

    /// ApplyJ ///
    /// This method computes
    /// $ out = J in $
    /// where J is the tangent operator (the linear approximation) of the mapping
    virtual void applyJ (const MechanicalParams* mparams /* PARAMS FIRST  = MechanicalParams::defaultInstance()*/, MultiVecDerivId outVel, ConstMultiVecDerivId inVel )
    {
        helper::vector<OutDataVecDeriv*> vecOutVel;
        getVecOutDeriv(outVel, vecOutVel);
        helper::vector<const In1DataVecDeriv*> vecIn1Vel;
        getConstVecIn1Deriv(inVel, vecIn1Vel);
        helper::vector<const In2DataVecDeriv*> vecIn2Vel;
        getConstVecIn2Deriv(inVel, vecIn2Vel);
        this->applyJ(mparams /* PARAMS FIRST */, vecOutVel, vecIn1Vel, vecIn2Vel);
    }
    /// This method must be reimplemented by all mappings.
    /// InDeriv and OutDeriv by default contains VecIds of type V_DERIV.
    /// The size of InDeriv vector is the same as the number of fromModels.
    /// The size of OutDeriv vector is the same as the number of OutModels.
    virtual void applyJ(
        const MechanicalParams* mparams /* PARAMS FIRST */, const helper::vector< OutDataVecDeriv*>& dataVecOutVel,
        const helper::vector<const In1DataVecDeriv*>& dataVecIn1Vel,
        const helper::vector<const In2DataVecDeriv*>& dataVecIn2Vel)
#ifdef SOFA_DEPRECATE_OLD_API
        = 0;
#else
    {
        //Not optimized at all...
        helper::vector<OutVecDeriv*> vecOutVel;
        for(unsigned int i=0; i<dataVecOutVel.size(); i++)
            vecOutVel.push_back(dataVecOutVel[i]->beginEdit(mparams));

        helper::vector<const In1VecDeriv*> vecIn1Vel;
        for(unsigned int i=0; i<dataVecIn1Vel.size(); i++)
            vecIn1Vel.push_back(&dataVecIn1Vel[i]->getValue(mparams));
        helper::vector<const In2VecDeriv*> vecIn2Vel;
        for(unsigned int i=0; i<dataVecIn2Vel.size(); i++)
            vecIn2Vel.push_back(&dataVecIn2Vel[i]->getValue(mparams));
        this->applyJ(vecOutVel, vecIn1Vel, vecIn2Vel);

        //Really Not optimized at all...
        for(unsigned int i=0; i<dataVecOutVel.size(); i++)
            dataVecOutVel[i]->endEdit(mparams);
    }
    /// Compat Method
    /// @deprecated
    virtual void applyJ(const helper::vector< OutVecDeriv*>& /* outDeriv */,
            const helper::vector<const In1VecDeriv*>& /* inDeriv1 */,
            const helper::vector<const In2VecDeriv*>& /* inDeriv2 */) { };
#endif //SOFA_DEPRECATE_OLD_API

    /// ApplyJT (Force)///
    /// Apply the mapping to Force vectors.
    virtual void applyJT (const MechanicalParams* mparams /* PARAMS FIRST  = MechanicalParams::defaultInstance()*/, MultiVecDerivId inForce, ConstMultiVecDerivId outForce )
    {
        helper::vector<In1DataVecDeriv*> vecOut1Force;
        getVecIn1Deriv(inForce, vecOut1Force);
        helper::vector<In2DataVecDeriv*> vecOut2Force;
        getVecIn2Deriv(inForce, vecOut2Force);

        helper::vector<const OutDataVecDeriv*> vecInForce;
        getConstVecOutDeriv(outForce, vecInForce);
        this->applyJT(mparams /* PARAMS FIRST */, vecOut1Force, vecOut2Force, vecInForce);
    }
    /// This method must be reimplemented by all mappings.
    /// InDeriv and OutDeriv by default contains VecIds of type V_DERIV.
    /// The size of InDeriv vector is the same as the number of fromModels.
    /// The size of OutDeriv vector is the same as the number of OutModels.
    virtual void applyJT(
        const MechanicalParams* mparams /* PARAMS FIRST */, const helper::vector< In1DataVecDeriv*>& dataVecOut1Force,
        const helper::vector< In2DataVecDeriv*>& dataVecOut2Force,
        const helper::vector<const OutDataVecDeriv*>& dataVecInForce)
#ifdef SOFA_DEPRECATE_OLD_API
        = 0;
#else
    {
        //Not optimized at all...
        helper::vector<In1VecDeriv*> vecOut1Force;
        for(unsigned int i=0; i<dataVecOut1Force.size(); i++)
            vecOut1Force.push_back(dataVecOut1Force[i]->beginEdit(mparams));
        helper::vector<In2VecDeriv*> vecOut2Force;
        for(unsigned int i=0; i<dataVecOut2Force.size(); i++)
            vecOut2Force.push_back(dataVecOut2Force[i]->beginEdit(mparams));

        helper::vector<const OutVecDeriv*> vecInForce;
        for(unsigned int i=0; i<dataVecInForce.size(); i++)
            vecInForce.push_back(&dataVecInForce[i]->getValue(mparams));

        this->applyJT(vecOut1Force, vecOut2Force, vecInForce);

        //Really Not optimized at all...
        for(unsigned int i=0; i<dataVecOut1Force.size(); i++)
        {
            dataVecOut1Force[i]->endEdit(mparams);
        }
        for(unsigned int i=0; i<dataVecOut2Force.size(); i++)
        {
            dataVecOut2Force[i]->endEdit(mparams);
        }
    }
    /// Compat Method
    /// @deprecated
    virtual void applyJT(const helper::vector< In1VecDeriv*>& /* out1Deriv */,
            const helper::vector< In2VecDeriv*>& /* out2Deriv */,
            const helper::vector<const OutVecDeriv*>& /* inDeriv */) { };
#endif //SOFA_DEPRECATE_OLD_API

    /// ApplyJT (Constraint)///
    virtual void applyJT(const ConstraintParams* cparams /* PARAMS FIRST  = ConstraintParams::defaultInstance()*/, MultiMatrixDerivId inConst, ConstMultiMatrixDerivId outConst )
    {
        helper::vector<In1DataMatrixDeriv*> matOut1Const;
        getMatIn1Deriv(inConst, matOut1Const);
        helper::vector<In2DataMatrixDeriv*> matOut2Const;
        getMatIn2Deriv(inConst, matOut2Const);

        helper::vector<const OutDataMatrixDeriv*> matInConst;
        getConstMatOutDeriv(outConst, matInConst);
        this->applyJT(cparams /* PARAMS FIRST */, matOut1Const, matOut2Const, matInConst);
    }
    /// This method must be reimplemented by all mappings if they need to support constraints.
    virtual void applyJT(
        const ConstraintParams* cparams /* PARAMS FIRST */, const helper::vector< In1DataMatrixDeriv*>& dataMatOut1Const ,
        const helper::vector< In2DataMatrixDeriv*>&  dataMatOut2Const ,
        const helper::vector<const OutDataMatrixDeriv*>& dataMatInConst)
#ifdef SOFA_DEPRECATE_OLD_API
    {
        serr << "This mapping does not support constraint because Multi2Mapping::applyJT(const ConstraintParams*, const helper::vector< In1DataMatrixDeriv*>&, const helper::vector< In2DataMatrixDeriv*>&, const helper::vector<const OutDataMatrixDeriv*>&) is not overloaded." << sendl;
    }
#else
    {
        //Not optimized at all...
        helper::vector<In1MatrixDeriv*> matOut1Const;
        for(unsigned int i=0; i<dataMatOut1Const.size(); i++)
            matOut1Const.push_back(dataMatOut1Const[i]->beginEdit(cparams));
        helper::vector<In2MatrixDeriv*> matOut2Const;
        for(unsigned int i=0; i<dataMatOut2Const.size(); i++)
            matOut2Const.push_back(dataMatOut2Const[i]->beginEdit(cparams));

        helper::vector<const OutMatrixDeriv*> matInConst;
        for(unsigned int i=0; i<dataMatInConst.size(); i++)
            matInConst.push_back(&dataMatInConst[i]->getValue(cparams));

        this->applyJT(matOut1Const, matOut2Const, matInConst);

        //Really Not optimized at all...
        for(unsigned int i=0; i<dataMatOut1Const.size(); i++)
        {
            dataMatOut1Const[i]->endEdit(cparams);
        }
        for(unsigned int i=0; i<dataMatOut2Const.size(); i++)
        {
            dataMatOut2Const[i]->endEdit(cparams);
        }
    }
    /// Compat Method
    /// @deprecated
    virtual void applyJT( const helper::vector< In1MatrixDeriv*>& /*outConstraint1*/ ,
            const helper::vector< In2MatrixDeriv*>& /*outConstraint2*/ ,
            const helper::vector<const OutMatrixDeriv*>& /*inConstraint*/ )
    {
        serr << "This mapping does not support constraints since Multi2Mapping::applyJT(const helper::vector< In1MatrixDeriv*>&, const helper::vector< In2MatrixDeriv*>&,helper::vector<const OutMatrixDeriv*>&) is not overloaded" << sendl;
    }
#endif //SOFA_DEPRECATE_OLD_API

    /// computeAccFromMapping
    virtual void computeAccFromMapping(const MechanicalParams* mparams /* PARAMS FIRST  = MechanicalParams::defaultInstance()*/, MultiVecDerivId outAcc, ConstMultiVecDerivId inVel, ConstMultiVecDerivId inAcc )
    {
        helper::vector<OutDataVecDeriv*> vecOutAcc;
        getVecOutDeriv(outAcc, vecOutAcc);

        helper::vector<const In1DataVecDeriv*> vecIn1Vel;
        getConstVecIn1Deriv(inVel, vecIn1Vel);
        helper::vector<const In1DataVecDeriv*> vecIn1Acc;
        getConstVecIn1Deriv(inAcc, vecIn1Acc);

        helper::vector<const In2DataVecDeriv*> vecIn2Vel;
        getConstVecIn2Deriv(inVel, vecIn2Vel);
        helper::vector<const In2DataVecDeriv*> vecIn2Acc;
        getConstVecIn2Deriv(inAcc, vecIn2Acc);

        this->computeAccFromMapping(mparams /* PARAMS FIRST */, vecOutAcc, vecIn1Vel, vecIn2Vel,vecIn1Acc, vecIn2Acc);
    }
    /// This method must be reimplemented by all mappings if they need to support composite accelerations
    virtual void computeAccFromMapping(
        const MechanicalParams* mparams /* PARAMS FIRST */, const helper::vector< OutDataVecDeriv*>& dataVecOutAcc,
        const helper::vector<const In1DataVecDeriv*>& dataVecIn1Vel,
        const helper::vector<const In2DataVecDeriv*>& dataVecIn2Vel,
        const helper::vector<const In1DataVecDeriv*>& dataVecIn1Acc,
        const helper::vector<const In2DataVecDeriv*>& dataVecIn2Acc)
#ifdef SOFA_DEPRECATE_OLD_API
    {
    }
#else
    {
        //Not optimized at all...
        helper::vector<OutVecDeriv*> vecOutAcc;
        for(unsigned int i=0; i<dataVecOutAcc.size(); i++)
            vecOutAcc.push_back(dataVecOutAcc[i]->beginEdit(mparams));

        helper::vector<const In1VecDeriv*> vecIn1Vel;
        for(unsigned int i=0; i<dataVecIn1Vel.size(); i++)
            vecIn1Vel.push_back(&dataVecIn1Vel[i]->getValue(mparams));
        helper::vector<const In1VecDeriv*> vecIn1Acc;
        for(unsigned int i=0; i<dataVecIn1Acc.size(); i++)
            vecIn1Acc.push_back(&dataVecIn1Acc[i]->getValue(mparams));

        helper::vector<const In2VecDeriv*> vecIn2Vel;
        for(unsigned int i=0; i<dataVecIn2Vel.size(); i++)
            vecIn2Vel.push_back(&dataVecIn2Vel[i]->getValue(mparams));
        helper::vector<const In2VecDeriv*> vecIn2Acc;
        for(unsigned int i=0; i<dataVecIn2Acc.size(); i++)
            vecIn2Acc.push_back(&dataVecIn2Acc[i]->getValue(mparams));

        this->computeAccFromMapping(vecOutAcc, vecIn1Vel, vecIn2Vel, vecIn1Acc, vecIn2Acc);

        //Really Not optimized at all...
        for(unsigned int i=0; i<dataVecOutAcc.size(); i++)
            dataVecOutAcc[i]->endEdit(mparams);
    }
    /// Compat Method
    /// @deprecated
    virtual void computeAccFromMapping( const helper::vector< OutVecDeriv*>& /*outDx*/,
            const helper::vector<const In1VecDeriv*>& /*inV1 */,
            const helper::vector<const In2VecDeriv*>& /*inV2 */,
            const helper::vector<const In1VecDeriv*>& /*inDx1 */,
            const helper::vector<const In2VecDeriv*>& /*inDx2 */ )
    {
    }
#endif //SOFA_DEPRECATE_OLD_API

    virtual void init();

    ///<TO REMOVE>
    /// Apply the mapping to position and velocity vectors.
    ///
    /// This method call the internal apply(helper::vector<VecId>& inPos, helper::vector<VecId>& outPos)
    /// and applyJ(helper::vector<VecId>& inDeriv, helper::vector<VecId>& outDeriv) methods.
    //virtual void updateMapping();

    /// Disable the mapping to get the original coordinates of the mapped model.
    ///
    /// It is for instance used in RigidMapping to get the local coordinates of the object.
    virtual void disable();

    virtual std::string getTemplateName() const
    {
        return templateName(this);
    }

    static std::string templateName(const Multi2Mapping<TIn1,TIn2, TOut>* = NULL);

    /// Pre-construction check method called by ObjectFactory.
    ///
    /// This implementation read the object1 and object2 attributes and check
    /// if they are compatible with the input and output models types of this
    /// mapping.
    template<class T>
    static bool canCreate(T*& obj, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg)
    {
        std::string input1 = arg->getAttribute("input1","");
        std::string input2 = arg->getAttribute("input2","");
        std::string output = arg->getAttribute("output","");
        if (!input1.empty() && !LinkFromModels1::CheckPaths(input1, context))
            return false;
        if (!input2.empty() && !LinkFromModels2::CheckPaths(input2, context))
            return false;
        if (output.empty() || !LinkToModels::CheckPaths(output, context))
            return false;

        return BaseMapping::canCreate(obj, context, arg);
    }

    /// Construction method called by ObjectFactory.
    ///
    /// This implementation read the input and output attributes to
    /// find the input and output models of this mapping.
    template<class T>
    static typename T::SPtr create(T*, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg)
    {
        typename T::SPtr obj = sofa::core::objectmodel::New<T>();

        if (context)
            context->addObject(obj);

        if (arg)
            obj->parse(arg);

        return obj;
    }


protected:
    void getVecIn1Coord     (const MultiVecCoordId id,         helper::vector<      In1DataVecCoord*> &v) const
    {   for (unsigned int i=0; i<fromModels1.size(); ++i) v.push_back(id[fromModels1[i]].write()); }
    void getConstVecIn1Coord(const ConstMultiVecCoordId id,    helper::vector<const In1DataVecCoord*> &v) const
    {   for (unsigned int i=0; i<fromModels1.size(); ++i) v.push_back(id[fromModels1[i]].read());  }
    void getVecIn1Deriv     (const MultiVecDerivId id,         helper::vector<      In1DataVecDeriv*> &v) const
    {   for (unsigned int i=0; i<fromModels1.size(); ++i) v.push_back(id[fromModels1[i]].write()); }
    void getConstVecIn1Deriv(const ConstMultiVecDerivId id,    helper::vector<const In1DataVecDeriv*> &v) const
    {   for (unsigned int i=0; i<fromModels1.size(); ++i) v.push_back(id[fromModels1[i]].read());  }
    void getMatIn1Deriv     (const MultiMatrixDerivId id,      helper::vector<      In1DataMatrixDeriv*> &v) const
    {   for (unsigned int i=0; i<fromModels1.size(); ++i) v.push_back(id[fromModels1[i]].write()); }
    void getConstMatIn1Deriv(const ConstMultiMatrixDerivId id, helper::vector<const In1DataMatrixDeriv*> &v) const
    {   for (unsigned int i=0; i<fromModels1.size(); ++i) v.push_back(id[fromModels1[i]].read());  }

    void getVecIn2Coord     (const MultiVecCoordId id,         helper::vector<      In2DataVecCoord*> &v) const
    {   for (unsigned int i=0; i<fromModels2.size(); ++i) v.push_back(id[fromModels2[i]].write()); }
    void getConstVecIn2Coord(const ConstMultiVecCoordId id,    helper::vector<const In2DataVecCoord*> &v) const
    {   for (unsigned int i=0; i<fromModels2.size(); ++i) v.push_back(id[fromModels2[i]].read());  }
    void getVecIn2Deriv     (const MultiVecDerivId id,         helper::vector<      In2DataVecDeriv*> &v) const
    {   for (unsigned int i=0; i<fromModels2.size(); ++i) v.push_back(id[fromModels2[i]].write()); }
    void getConstVecIn2Deriv(const ConstMultiVecDerivId id,    helper::vector<const In2DataVecDeriv*> &v) const
    {   for (unsigned int i=0; i<fromModels2.size(); ++i) v.push_back(id[fromModels2[i]].read());  }
    void getMatIn2Deriv     (const MultiMatrixDerivId id,      helper::vector<      In2DataMatrixDeriv*> &v) const
    {   for (unsigned int i=0; i<fromModels2.size(); ++i) v.push_back(id[fromModels2[i]].write()); }
    void getConstMatIn2Deriv(const ConstMultiMatrixDerivId id, helper::vector<const In2DataMatrixDeriv*> &v) const
    {   for (unsigned int i=0; i<fromModels2.size(); ++i) v.push_back(id[fromModels2[i]].read());  }

    void getVecOutCoord     (const MultiVecCoordId id,         helper::vector<      OutDataVecCoord*> &v) const
    {   for (unsigned int i=0; i<toModels.size(); ++i)  v.push_back(id[toModels[i]].write());      }
    void getConstVecOutCoord(const ConstMultiVecCoordId id,    helper::vector<const OutDataVecCoord*> &v) const
    {   for (unsigned int i=0; i<toModels.size(); ++i)  v.push_back(id[toModels[i]].read());       }
    void getVecOutDeriv     (const MultiVecDerivId id,         helper::vector<      OutDataVecDeriv*> &v) const
    {   for (unsigned int i=0; i<toModels.size(); ++i)  v.push_back(id[toModels[i]].write());      }
    void getConstVecOutDeriv(const ConstMultiVecDerivId id,    helper::vector<const OutDataVecDeriv*> &v) const
    {   for (unsigned int i=0; i<toModels.size(); ++i)  v.push_back(id[toModels[i]].read());       }
    void getMatOutDeriv     (const MultiMatrixDerivId id,      helper::vector<      OutDataMatrixDeriv*> &v) const
    {   for (unsigned int i=0; i<toModels.size(); ++i)  v.push_back(id[toModels[i]].write());      }
    void getConstMatOutDeriv(const ConstMultiMatrixDerivId id, helper::vector<const OutDataMatrixDeriv*> &v) const
    {   for (unsigned int i=0; i<toModels.size(); ++i)  v.push_back(id[toModels[i]].read());       }
};

} // namespace core

} // namespace sofa

#endif

