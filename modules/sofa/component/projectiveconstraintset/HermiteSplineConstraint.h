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
#ifndef SOFA_COMPONENT_PROJECTIVECONSTRAINTSET_HERMITESPLINECONSTRAINT_H
#define SOFA_COMPONENT_PROJECTIVECONSTRAINTSET_HERMITESPLINECONSTRAINT_H

#include <sofa/core/behavior/ProjectiveConstraintSet.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/helper/vector.h>
#include <sofa/component/topology/TopologySubsetData.h>


namespace sofa
{

namespace component
{

namespace projectiveconstraintset
{

using core::objectmodel::Data;
using namespace sofa::core::objectmodel;
using namespace sofa::defaulttype;

/**
	Impose a trajectory to given Dofs following a Hermite cubic spline constraint.
	Control parameters are :
	  - begin and end points
	  - derivates at this points
	  - acceleration curve on the trajectory
	*/
template <class DataTypes>
class HermiteSplineConstraint : public core::behavior::ProjectiveConstraintSet<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(HermiteSplineConstraint,DataTypes),SOFA_TEMPLATE(sofa::core::behavior::ProjectiveConstraintSet, DataTypes));

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::MatrixDeriv MatrixDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename Coord::value_type Real;
    typedef typename MatrixDeriv::RowIterator MatrixDerivRowIterator;
    typedef typename MatrixDeriv::RowType MatrixDerivRowType;
    typedef Data<VecCoord> DataVecCoord;
    typedef Data<VecDeriv> DataVecDeriv;
    typedef Data<MatrixDeriv> DataMatrixDeriv;
    typedef helper::vector<unsigned int> SetIndexArray;
    typedef sofa::component::topology::PointSubsetData< SetIndexArray > SetIndex;
    typedef typename defaulttype::Vec<3, Real> Vec3R;
    typedef typename defaulttype::Vec<2, Real> Vec2R;
    typedef typename helper::Quater<Real> QuatR;

public:
    ///indices of the DOFs constraints
    SetIndex m_indices;

    /// the time steps defining the duration of the constraint
    Data<Real> m_tBegin;
    Data<Real> m_tEnd;

    /// control parameters :
    /// first control point
    Data<Vec3R> m_x0;
    /// first derivated control point
    Data<Vec3R> m_dx0;
    /// second control point
    Data<Vec3R> m_x1;
    /// second derivated control point
    Data<Vec3R> m_dx1;
    /// acceleration parameters : the accaleration curve is itself a hermite spline, with first point at (0,0) and second at (1,1)
    /// and derivated on this points are :
    Data<Vec2R> m_sx0;
    Data<Vec2R> m_sx1;



protected:
    HermiteSplineConstraint();

    HermiteSplineConstraint(core::behavior::MechanicalState<DataTypes>* mstate);

    ~HermiteSplineConstraint();
public:
    void clearConstraints();
    void addConstraint(unsigned index );

    void setBeginTime(const Real &t) {m_tBegin.setValue(t);}
    void setEndTime(const Real &t) {m_tEnd.setValue(t);}

    Real getBeginTime() {return m_tBegin.getValue();}
    Real getEndTime() {return m_tEnd.getValue();}

    void computeHermiteCoefs( const Real u, Real &H00, Real &H10, Real &H01, Real &H11);
    void computeDerivateHermiteCoefs( const Real u, Real &dH00, Real &dH10, Real &dH01, Real &dH11);

    /// -- Constraint interface
    void init();
    void reinit();


    void projectResponse(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& resData);
    void projectVelocity(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& vData);
    void projectPosition(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecCoord& xData);
    void projectJacobianMatrix(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataMatrixDeriv& cData);

    void draw(const core::visual::VisualParams* vparams);

protected:
    template <class DataDeriv>
    void projectResponseT(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataDeriv& dx);

    /// Pointer to the current topology
    sofa::core::topology::BaseMeshTopology* topology;

};

} // namespace projectiveconstraintset

} // namespace component

} // namespace sofa

#endif
