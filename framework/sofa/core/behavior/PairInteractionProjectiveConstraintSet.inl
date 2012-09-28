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
#ifndef SOFA_CORE_BEHAVIOR_PAIRINTERACTIONPROJECTIVECONSTRAINTSET_INL
#define SOFA_CORE_BEHAVIOR_PAIRINTERACTIONPROJECTIVECONSTRAINTSET_INL

#include <sofa/core/behavior/PairInteractionProjectiveConstraintSet.h>


namespace sofa
{

namespace core
{

namespace behavior
{

template<class DataTypes>
PairInteractionProjectiveConstraintSet<DataTypes>::PairInteractionProjectiveConstraintSet(MechanicalState<DataTypes> *mm1, MechanicalState<DataTypes> *mm2)
    : endTime( initData(&endTime,(double)-1,"endTime","The constraint stops acting after the given value.\nUse a negative value for infinite constraints") )
    , mstate1(initLink("object1", "First object to constrain"), mm1)
    , mstate2(initLink("object2", "Second object to constrain"), mm2)
{
    if (!mm1)
        mstate1.setPath("@./"); // default to state of the current node
    if (!mm2)
        mstate2.setPath("@./"); // default to state of the current node
}

template<class DataTypes>
PairInteractionProjectiveConstraintSet<DataTypes>::~PairInteractionProjectiveConstraintSet()
{
}

template<class DataTypes>
void PairInteractionProjectiveConstraintSet<DataTypes>::init()
{
    BaseInteractionProjectiveConstraintSet::init();
    if (mstate1 == NULL || mstate2 == NULL)
    {
        mstate1 = mstate2 = dynamic_cast< MechanicalState<DataTypes>* >(getContext()->getMechanicalState());
    }

    this->mask1 = &mstate1->forceMask;
    this->mask2 = &mstate2->forceMask;
}

template<class DataTypes>
bool PairInteractionProjectiveConstraintSet<DataTypes>::isActive() const
{
    if( endTime.getValue()<0 ) return true;
    return endTime.getValue()>getContext()->getTime();
}

template<class DataTypes>
void PairInteractionProjectiveConstraintSet<DataTypes>::projectJacobianMatrix(const MechanicalParams* /*mparams*/ /* PARAMS FIRST */, MultiMatrixDerivId /*cId*/)
{
    serr << "NOT IMPLEMENTED YET" << sendl;
}

#ifdef SOFA_SMP
template<class DataTypes>
struct PairConstraintProjectResponseTask
{
    void operator()(const MechanicalParams* mparams /* PARAMS FIRST */, PairInteractionProjectiveConstraintSet<DataTypes>  *c, Shared_rw< objectmodel::Data< typename DataTypes::VecDeriv> > dx1,Shared_rw< objectmodel::Data< typename DataTypes::VecDeriv> > dx2)
    {
        c->projectResponse(mparams /* PARAMS FIRST */, dx1.access(), dx2.access());
    }
};

template<class DataTypes>
struct PairConstraintProjectVelocityTask
{
    void operator()(const MechanicalParams* mparams /* PARAMS FIRST */, PairInteractionProjectiveConstraintSet<DataTypes>  *c, Shared_rw< objectmodel::Data< typename DataTypes::VecDeriv> > v1, Shared_rw< objectmodel::Data< typename DataTypes::VecDeriv> > v2)
    {
        c->projectVelocity(mparams /* PARAMS FIRST */, v1.access(), v2.access());
    }
};

template<class DataTypes>
struct PairConstraintProjectPositionTask
{
    void operator()(const MechanicalParams* mparams /* PARAMS FIRST */, PairInteractionProjectiveConstraintSet<DataTypes>  *c, Shared_rw< objectmodel::Data< typename DataTypes::VecCoord> > x1, Shared_rw< objectmodel::Data< typename DataTypes::VecCoord> > x2)
    {
        c->projectPosition(mparams /* PARAMS FIRST */, x1.access(), x2.access());
    }
};
#endif /* SOFA_SMP */

template<class DataTypes>
void PairInteractionProjectiveConstraintSet<DataTypes>::projectResponse(const MechanicalParams* mparams /* PARAMS FIRST */, MultiVecDerivId dxId)
{
    if( !isActive() ) return;
    if (mstate1 && mstate2)
    {
        this->mask1 = &mstate1->forceMask;
        this->mask2 = &mstate2->forceMask;
#ifdef SOFA_SMP
        if (mparams->execMode() == ExecParams::EXEC_KAAPI)
            Task<PairConstraintProjectResponseTask<DataTypes> >(mparams /* PARAMS FIRST */, this, **defaulttype::getShared(*dxId[mstate1.get(mparams)].write()), **defaulttype::getShared(*dxId[mstate2.get(mparams)].write()));
        else
#endif /* SOFA_SMP */
            projectResponse(mparams /* PARAMS FIRST */, *dxId[mstate1.get(mparams)].write(), *dxId[mstate2.get(mparams)].write());
    }
}

template<class DataTypes>
void PairInteractionProjectiveConstraintSet<DataTypes>::projectVelocity(const MechanicalParams* mparams /* PARAMS FIRST */, MultiVecDerivId vId)
{
    if( !isActive() ) return;
    if (mstate1 && mstate2)
    {
        this->mask1 = &mstate1->forceMask;
        this->mask2 = &mstate2->forceMask;
#ifdef SOFA_SMP
        if (mparams->execMode() == ExecParams::EXEC_KAAPI)
            Task<PairConstraintProjectVelocityTask<DataTypes> >(mparams /* PARAMS FIRST */, this, **defaulttype::getShared(*vId[mstate1.get(mparams)].write()), **defaulttype::getShared(*vId[mstate2.get(mparams)].write()));
        else
#endif /* SOFA_SMP */
            projectVelocity(mparams /* PARAMS FIRST */, *vId[mstate1.get(mparams)].write(), *vId[mstate2.get(mparams)].write());
    }
}

template<class DataTypes>
void PairInteractionProjectiveConstraintSet<DataTypes>::projectPosition(const MechanicalParams* mparams /* PARAMS FIRST */, MultiVecCoordId xId)
{
    if( !isActive() ) return;
    if (mstate1 && mstate2)
    {
        this->mask1 = &mstate1->forceMask;
        this->mask2 = &mstate2->forceMask;
#ifdef SOFA_SMP
        if (mparams->execMode() == ExecParams::EXEC_KAAPI)
            Task<PairConstraintProjectPositionTask<DataTypes> >(mparams /* PARAMS FIRST */, this, **defaulttype::getShared(*xId[mstate1.get(mparams)].write()), **defaulttype::getShared(*xId[mstate2.get(mparams)].write()));
        else
#endif /* SOFA_SMP */
            projectPosition(mparams /* PARAMS FIRST */, *xId[mstate1.get(mparams)].write(), *xId[mstate2.get(mparams)].write());
    }
}

} // namespace behavior

} // namespace core

} // namespace sofa

#endif
