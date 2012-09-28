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
#ifndef SOFA_CORE_MULTIMAPPING_INL
#define SOFA_CORE_MULTIMAPPING_INL

#include <sofa/core/MultiMapping.h>

namespace sofa
{

namespace core
{

template< class In, class Out>
MultiMapping<In,Out>::MultiMapping()
    : fromModels(initLink("input", "Input Object(s)"))
    , toModels(initLink("output", "Output Object(s)"))
    , f_applyRestPosition( initData( &f_applyRestPosition, false, "applyRestPosition", "set to true to apply this mapping to restPosition at init"))
{

}

template < class In, class Out >
void MultiMapping<In,Out>::addInputModel(BaseState* fromModel, const std::string& path)
{
    State<In>* from = dynamic_cast<State<In>*>(fromModel);
    assert(from && "MultiMapping needs a State of the appropriate type to add as input model");
    this->fromModels.add(from, path);
}

template< class In, class Out >
void MultiMapping<In,Out>::addOutputModel(BaseState* toModel, const std::string& path)
{
    State<Out>* to = dynamic_cast<State<Out>*>(toModel);
    assert(to);
    this->toModels.add(to, path);
    if (isMechanical())
    {
        if(to != NULL && !testMechanicalState(to))
            setNonMechanical();
    }
}

template< class In, class Out>
const typename MultiMapping<In,Out>::VecFromModels& MultiMapping<In,Out>::getFromModels()
{
    return this->fromModels.getValue();
}

template< class In, class Out>
const typename MultiMapping<In,Out>::VecToModels& MultiMapping<In,Out>::getToModels()
{
    return this->toModels.getValue();
}

template< class In, class Out >
helper::vector<BaseState*> MultiMapping<In,Out>::getFrom()
{
    const VecFromModels& models = getFromModels();
    unsigned int size = models.size();
    helper::vector<BaseState*> baseModels(size);
    for (unsigned int i=0; i<size; ++i) baseModels[i] = models[i].ptr.get();
    return baseModels;
}

template< class In, class Out >
helper::vector<BaseState* > MultiMapping<In,Out>::getTo()
{
    const VecToModels& models = getToModels();
    unsigned int size = models.size();
    helper::vector<BaseState*> baseModels(size);
    for (unsigned int i=0; i<size; ++i) baseModels[i] = models[i].ptr.get();
    return baseModels;
}

template <class In, class Out>
helper::vector<behavior::BaseMechanicalState*> MultiMapping<In,Out>::getMechFrom()
{
    helper::vector<behavior::BaseMechanicalState*> mechFromVec;
    for (unsigned int i=0 ; i<this->fromModels.size() ; i++)
    {
        behavior::BaseMechanicalState* meshFrom = dynamic_cast<behavior::BaseMechanicalState*> (this->fromModels.get(i));
        if(meshFrom)
            mechFromVec.push_back(meshFrom);
    }
    return mechFromVec;
}

template <class In, class Out>
helper::vector<behavior::BaseMechanicalState*> MultiMapping<In,Out>::getMechTo()
{
    helper::vector<behavior::BaseMechanicalState*> mechToVec;
    for (unsigned int i=0 ; i<this->toModels.size() ; i++)
    {
        behavior::BaseMechanicalState* meshTo = dynamic_cast<behavior::BaseMechanicalState*> (this->toModels.get(i));
        if(meshTo)
            mechToVec.push_back(meshTo);
    }
    return mechToVec;
}

template <class In, class Out>
void MultiMapping<In,Out>::init()
{
    apply(MechanicalParams::defaultInstance()  /* PARAMS FIRST */, VecCoordId::position(), ConstVecCoordId::position());
    applyJ(MechanicalParams::defaultInstance()  /* PARAMS FIRST */, VecDerivId::velocity(), ConstVecDerivId::velocity());
    if (f_applyRestPosition.getValue())
        apply(MechanicalParams::defaultInstance() /* PARAMS FIRST */, VecCoordId::restPosition(), ConstVecCoordId::restPosition());
}

#ifdef SOFA_SMP
template<class T>
struct ParallelMultiMappingApply
{
    void operator()(const MechanicalParams* mparams  /* PARAMS FIRST */, void *m, Shared_rw<defaulttype::SharedVector<typename T::Out::VecCoord*> > out, Shared_r<defaulttype::SharedVector<const typename T::In::VecCoord*> > in)
    {
        ((T *)m)->apply(mparams  /* PARAMS FIRST */, out.access(), in.read());
    }
};

template<class T>
struct ParallelMultiMappingApplyJ
{
    void operator()(void *m, Shared_rw<defaulttype::SharedVector<typename T::Out::VecDeriv*> > out, Shared_r<defaulttype::SharedVector<const typename T::In::VecDeriv*> > in)
    {
        ((T *)m)->applyJ(out.access(), in.read());
    }
};

template<class T>
struct accessOutPos
{
    void operator()(void *m, Shared_rw<typename T::Out::VecCoord> out)
    {
        out.access();
    }
};

template<class T>
struct ParallelMultiMappingApply3
{
    void operator()(void *m, Shared_rw<typename T::Out::VecCoord> out, Shared_r<typename T::In::VecCoord> in1, Shared_r<typename T::In::VecCoord> in2)
    {
        out.access();
        in1.read();
        in2.read();
        ((T *)m)->apply(((T *)m)->VecOutPos,((T *)m)->VecInPos);
    }
};

template<class T>
struct ParallelMultiMappingApplyJ3
{
    void operator()(void *m, Shared_rw<typename T::Out::VecDeriv> out, Shared_r<typename T::In::VecDeriv> in1,Shared_r<typename T::In::VecDeriv> in2)
    {
        out.access();
        in1.read();
        in2.read();
        ((T *)m)->applyJ(((T *)m)->VecOutVel,((T *)m)->VecInVel);
    }
};
#endif /* SOFA_SMP */

template <class In, class Out>
void MultiMapping<In,Out>::apply(const MechanicalParams* mparams /* PARAMS FIRST */, MultiVecCoordId outPos, ConstMultiVecCoordId inPos)
{
    helper::vector<OutDataVecCoord*> vecOutPos;
    getVecOutCoord(outPos, vecOutPos);
    helper::vector<const InDataVecCoord*> vecInPos;
    getConstVecInCoord(inPos, vecInPos);

#ifdef SOFA_SMP
//		if (mparams->execMode() == ExecParams::EXEC_KAAPI)
//			Task<ParallelMultiMappingApply< MultiMapping<In,Out> > >(mparams /* PARAMS FIRST */, this,
//					**defaulttype::getShared(*out), **defaulttype::getShared(*in));
//		else
#endif /* SOFA_SMP */
    this->apply(mparams /* PARAMS FIRST */, vecOutPos, vecInPos);
}// MultiMapping::apply

template <class In, class Out>
void MultiMapping<In,Out>::applyJ(const MechanicalParams* mparams /* PARAMS FIRST */, MultiVecDerivId outVel, ConstMultiVecDerivId inVel)
{
    helper::vector<OutDataVecDeriv*> vecOutVel;
    getVecOutDeriv(outVel, vecOutVel);
    helper::vector<const InDataVecDeriv*> vecInVel;
    getConstVecInDeriv(inVel, vecInVel);

#ifdef SOFA_SMP
//		if (mparams->execMode() == ExecParams::EXEC_KAAPI)
//			Task<ParallelMultiMappingApplyJ< MultiMapping<In,Out> > >(mparams /* PARAMS FIRST */, this,
//					**defaulttype::getShared(*out), **defaulttype::getShared(*in));
//		else
#endif /* SOFA_SMP */
    this->applyJ(mparams /* PARAMS FIRST */, vecOutVel, vecInVel);
}// MultiMapping::applyJ

template <class In, class Out>
void MultiMapping<In,Out>::applyJT(const MechanicalParams* mparams /* PARAMS FIRST */, MultiVecDerivId inForce, ConstMultiVecDerivId outForce)
{
    helper::vector<InDataVecDeriv*> vecOutForce;
    getVecInDeriv(inForce, vecOutForce);
    helper::vector<const OutDataVecDeriv*> vecInForce;
    getConstVecOutDeriv(outForce, vecInForce);

    this->applyJT(mparams /* PARAMS FIRST */, vecOutForce, vecInForce);
}// MultiMapping::applyJT

template <class In, class Out>
std::string MultiMapping<In,Out>::templateName(const MultiMapping<In, Out>* /*mapping*/)
{
    //return std::string("MultiMapping<") + In::Name() + std::string(",") + Out::Name() + std::string(">");
    return In::Name() + std::string(",") + Out::Name();
}

template <class In, class Out>
void MultiMapping<In,Out>::disable()
{
}


} // namespace core

} // namespace sofa

#endif
