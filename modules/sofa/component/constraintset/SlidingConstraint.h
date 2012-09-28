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
#ifndef SOFA_COMPONENT_CONSTRAINTSET_SLIDINGCONSTRAINT_H
#define SOFA_COMPONENT_CONSTRAINTSET_SLIDINGCONSTRAINT_H

#include <sofa/component/component.h>

#include <sofa/core/behavior/PairInteractionConstraint.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <iostream>

namespace sofa
{

namespace component
{

namespace constraintset
{

template<class DataTypes>
class SlidingConstraint : public core::behavior::PairInteractionConstraint<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(SlidingConstraint,DataTypes), SOFA_TEMPLATE(core::behavior::PairInteractionConstraint,DataTypes));

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::MatrixDeriv MatrixDeriv;
    typedef typename DataTypes::MatrixDeriv::RowIterator MatrixDerivRowIterator;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename Coord::value_type Real;
    typedef typename core::behavior::MechanicalState<DataTypes> MechanicalState;
    typedef typename core::behavior::PairInteractionConstraint<DataTypes> Inherit;

    typedef core::objectmodel::Data<VecCoord>		DataVecCoord;
    typedef core::objectmodel::Data<VecDeriv>		DataVecDeriv;
    typedef core::objectmodel::Data<MatrixDeriv>    DataMatrixDeriv;

protected:

    bool yetIntegrated;

    unsigned int cid;

    Data<int> m1;
    Data<int> m2a;
    Data<int> m2b;

    Real m_dist;	// constraint violation
    Real thirdConstraint; // 0 if A<proj<B, -1 if proj<A, 1 if B<proj



    SlidingConstraint(MechanicalState* object1, MechanicalState* object2)
        : Inherit(object1, object2)
        , yetIntegrated(false)
        , m1(initData(&m1, 0, "sliding_point","index of the spliding point on the first model"))
        , m2a(initData(&m2a, 0, "axis_1","index of one end of the sliding axis"))
        , m2b(initData(&m2b, 0, "axis_2","index of the other end of the sliding axis"))
    {
    }

    SlidingConstraint(MechanicalState* object)
        : Inherit(object, object)
        , yetIntegrated(false)
        , m1(initData(&m1, 0, "sliding_point","index of the spliding point on the first model"))
        , m2a(initData(&m2a, 0, "axis_1","index of one end of the sliding axis"))
        , m2b(initData(&m2b, 0, "axis_2","index of the other end of the sliding axis"))
    {
    }

    SlidingConstraint()
        : yetIntegrated(false)
        , m1(initData(&m1, 0, "sliding_point","index of the spliding point on the first model"))
        , m2a(initData(&m2a, 0, "axis_1","index of one end of the sliding axis"))
        , m2b(initData(&m2b, 0, "axis_2","index of the other end of the sliding axis"))
    {
    }

    virtual ~SlidingConstraint()
    {
    }
public:
    virtual void init();

    virtual void buildConstraintMatrix(const core::ConstraintParams* cParams /* PARAMS FIRST =core::ConstraintParams::defaultInstance()*/, DataMatrixDeriv &c1, DataMatrixDeriv &c2, unsigned int &cIndex
            , const DataVecCoord &x1, const DataVecCoord &x2);

    void getConstraintViolation(const core::ConstraintParams* cParams /* PARAMS FIRST =core::ConstraintParams::defaultInstance()*/, defaulttype::BaseVector *v, const DataVecCoord &x1, const DataVecCoord &x2
            , const DataVecDeriv &v1, const DataVecDeriv &v2);

    virtual void getConstraintResolution(std::vector<core::behavior::ConstraintResolution*>& resTab, unsigned int& offset);

    void draw(const core::visual::VisualParams* vparams);
};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_CONSTRAINTSET_SLIDINGCONSTRAINT_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_CONSTRAINT_API SlidingConstraint< defaulttype::Vec3dTypes >;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_CONSTRAINT_API SlidingConstraint< defaulttype::Vec3fTypes >;
#endif
#endif

} // namespace constraintset

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_CONSTRAINTSET_SLIDINGCONSTRAINT_H
