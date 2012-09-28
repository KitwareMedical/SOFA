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
#ifndef SOFA_COMPONENT_CONSTRAINTSET_STOPPERCONSTRAINT_H
#define SOFA_COMPONENT_CONSTRAINTSET_STOPPERCONSTRAINT_H

#include <sofa/core/behavior/Constraint.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/behavior/OdeSolver.h>

namespace sofa
{

namespace component
{

namespace constraintset
{

class StopperConstraintResolution1Dof : public core::behavior::ConstraintResolution
{
protected:
    double _invW, _w, _min, _max ;

public:

    StopperConstraintResolution1Dof(const double &min, const double &max) { nbLines=1; _min=min; _max=max; }

    virtual void init(int line, double** w, double *force)
    {
        _w = w[line][line];
        _invW = 1.0/_w;
        force[line  ] = 0.0;
    }

    virtual void resolution(int line, double** /*w*/, double* d, double* force)
    {
        double dfree = d[line] - _w * force[line];

        if (dfree > _max)
            force[line] = (_max - dfree) * _invW;
        else if (dfree < _min)
            force[line] = (_min - dfree) * _invW;
        else
            force[line] = 0;
    }
};

template< class DataTypes >
class StopperConstraint : public core::behavior::Constraint<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(StopperConstraint,DataTypes), SOFA_TEMPLATE(core::behavior::Constraint,DataTypes));

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::MatrixDeriv MatrixDeriv;
    typedef typename Coord::value_type Real;
    typedef typename core::behavior::MechanicalState<DataTypes> MechanicalState;
    typedef typename core::behavior::Constraint<DataTypes> Inherit;

    typedef typename DataTypes::MatrixDeriv::RowIterator MatrixDerivRowIterator;
    typedef core::objectmodel::Data<VecCoord>		DataVecCoord;
    typedef core::objectmodel::Data<VecDeriv>		DataVecDeriv;
    typedef core::objectmodel::Data<MatrixDeriv>    DataMatrixDeriv;

protected:

    unsigned int cid;

    Data<int> index;
    Data<double> min, max;



    StopperConstraint(MechanicalState* object)
        : Inherit(object)
        , index(initData(&index, 0, "index", "index of the stop constraint"))
        , min(initData(&min, -100.0, "min", "minimum value accepted"))
        , max(initData(&max, 100.0, "max", "maximum value accepted"))
    {
    }


    StopperConstraint()
        : index(initData(&index, 0, "index", "index of the stop constraint"))
        , min(initData(&min, -100.0, "min", "minimum value accepted"))
        , max(initData(&max, 100.0, "max", "maximum value accepted"))
    {
    }

    virtual ~StopperConstraint() {}
public:
    virtual void init();
    virtual void buildConstraintMatrix(const core::ConstraintParams* cParams /* PARAMS FIRST =core::ConstraintParams::defaultInstance()*/, DataMatrixDeriv &c_d, unsigned int &cIndex, const DataVecCoord &x);
    virtual void getConstraintViolation(const core::ConstraintParams* cParams /* PARAMS FIRST =core::ConstraintParams::defaultInstance()*/, defaulttype::BaseVector *resV, const DataVecCoord &x, const DataVecDeriv &v);

    virtual void getConstraintResolution(std::vector<core::behavior::ConstraintResolution*>& resTab, unsigned int& offset);
};
} // namespace constraintset

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_CONSTRAINTSET_BILATERALINTERACTIONCONSTRAINT_H
