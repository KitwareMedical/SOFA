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
#ifndef SOFA_CORE_BEHAVIOR_MIXEDINTERACTIONCONSTRAINT_INL
#define SOFA_CORE_BEHAVIOR_MIXEDINTERACTIONCONSTRAINT_INL

#include "MixedInteractionConstraint.h"
#include <sofa/core/behavior/Constraint.inl>

namespace sofa
{

namespace core
{

namespace behavior
{

template<class DataTypes1, class DataTypes2>
MixedInteractionConstraint<DataTypes1, DataTypes2>::MixedInteractionConstraint(MechanicalState<DataTypes1> *mm1, MechanicalState<DataTypes2> *mm2)
    : endTime( initData(&endTime,(double)-1,"endTime","The constraint stops acting after the given value.\nUse a negative value for infinite constraints") )
    , mstate1(initLink("object1", "First object to constrain"), mm1)
    , mstate2(initLink("object2", "Second object to constrain"), mm2)
{
}

template<class DataTypes1, class DataTypes2>
MixedInteractionConstraint<DataTypes1, DataTypes2>::~MixedInteractionConstraint()
{
}

template<class DataTypes1, class DataTypes2>
void MixedInteractionConstraint<DataTypes1, DataTypes2>::init()
{
    BaseInteractionConstraint::init();
    this->mask1 = &mstate1->forceMask;
    this->mask2 = &mstate2->forceMask;
}

template<class DataTypes1, class DataTypes2>
bool MixedInteractionConstraint<DataTypes1, DataTypes2>::isActive() const
{
    if( endTime.getValue()<0 ) return true;
    return endTime.getValue()>getContext()->getTime();
}

template<class DataTypes1, class DataTypes2>
void MixedInteractionConstraint<DataTypes1, DataTypes2>::getConstraintViolation(const ConstraintParams* cParams /* PARAMS FIRST */, defaulttype::BaseVector *v)
{
    if (cParams)
    {
        getConstraintViolation(cParams /* PARAMS FIRST */, v, *cParams->readX(mstate1), *cParams->readX(mstate2), *cParams->readV(mstate1), *cParams->readV(mstate2));
    }
}

template<class DataTypes1, class DataTypes2>
void MixedInteractionConstraint<DataTypes1, DataTypes2>::buildConstraintMatrix(const ConstraintParams* cParams /* PARAMS FIRST */, MultiMatrixDerivId cId, unsigned int &cIndex)
{
    if (cParams)
    {
        buildConstraintMatrix(cParams /* PARAMS FIRST */, *cId[mstate1.get(cParams)].write(), *cId[mstate2.get(cParams)].write(), cIndex, *cParams->readX(mstate1), *cParams->readX(mstate2));
    }
}

} // namespace behavior

} // namespace core

} // namespace sofa

#endif
