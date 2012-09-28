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
#ifndef SOFA_COMPONENT_CONSTRAINTSET_LMCONSTRAINTSOLVER_H
#define SOFA_COMPONENT_CONSTRAINTSET_LMCONSTRAINTSOLVER_H

#include <sofa/core/behavior/ConstraintSolver.h>
#include <sofa/core/behavior/BaseLMConstraint.h>
#include <sofa/core/behavior/BaseConstraintCorrection.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/simulation/common/MechanicalVisitor.h>
#include <sofa/component/component.h>


#include <sofa/component/linearsolver/EigenMatrixManipulator.h>


namespace sofa
{

namespace component
{

namespace constraintset
{

using core::behavior::BaseLMConstraint;
using core::behavior::ConstraintGroup;

class SOFA_CONSTRAINT_API LMConstraintSolver : public sofa::core::behavior::ConstraintSolver
{
protected:
    typedef sofa::core::VecId VecId;
    typedef sofa::core::MultiVecId MultiVecId;
    typedef sofa::core::ConstraintParams::ConstOrder ConstOrder;

    typedef Eigen::Matrix<SReal, Eigen::Dynamic, Eigen::Dynamic> MatrixEigen;
    typedef linearsolver::VectorEigen          VectorEigen;
    typedef linearsolver::SparseMatrixEigen    SparseMatrixEigen;
    typedef linearsolver::SparseVectorEigen    SparseVectorEigen;

    typedef helper::set< sofa::core::behavior::BaseMechanicalState* > SetDof;
    typedef std::map< const sofa::core::behavior::BaseMechanicalState *, SparseMatrixEigen > DofToMatrix;
    typedef std::map< const sofa::core::behavior::BaseMechanicalState *, helper::set<unsigned int> > DofToMask;
    typedef std::map< const sofa::core::behavior::BaseMechanicalState *, core::behavior::BaseConstraintCorrection* > DofToConstraintCorrection;

public:
    SOFA_CLASS(LMConstraintSolver, sofa::core::behavior::ConstraintSolver);
protected:
    LMConstraintSolver();
public:
    virtual void init();
    virtual void reinit() {graphKineticEnergy.setDisplayed(traceKineticEnergy.getValue());};


    virtual bool prepareStates(const core::ConstraintParams *, MultiVecId res1, MultiVecId res2=MultiVecId::null());
    virtual bool buildSystem(const core::ConstraintParams *, MultiVecId res1, MultiVecId res2=MultiVecId::null());
    virtual bool solveSystem(const core::ConstraintParams *, MultiVecId res1, MultiVecId res2=MultiVecId::null());
    virtual bool applyCorrection(const core::ConstraintParams *, MultiVecId res1, MultiVecId res2=MultiVecId::null());

    virtual void handleEvent( core::objectmodel::Event *e);



    Data<bool> constraintAcc;
    Data<bool> constraintVel;
    Data<bool> constraintPos;
    Data<unsigned int> numIterations;
    Data<double> maxError;
    mutable Data<std::map < std::string, sofa::helper::vector<double> > > graphGSError;
    Data< bool > traceKineticEnergy;
    mutable Data<std::map < std::string, sofa::helper::vector<double> > > graphKineticEnergy;



    template <class T>
    std::string printDimension(const T& m)
    {
        std::ostringstream out;
        out << "(" << m.rows() << "," << m.cols() << ")";
        return out.str();
    }

    void convertSparseToDense(const SparseMatrixEigen& sparseM, MatrixEigen& out) const;
protected:
    /// Explore the graph, looking for LMConstraints: each LMConstraint can tell if they need State Propagation in order to compute the right hand term of the system
    virtual bool needPriorStatePropagation(core::ConstraintParams::ConstOrder order) const;

    /// Construct the Right hand term of the system
    virtual void buildRightHandTerm      ( const helper::vector< core::behavior::BaseLMConstraint* > &LMConstraints,
            VectorEigen &c, MultiVecId, ConstOrder Order ) const;
    /// Construct the Inverse of the mass matrix for a set of Dofs
    virtual void buildInverseMassMatrices( const SetDof &setDofs,
            DofToMatrix& invMassMatrices);
    /// Construct the L matrices: write the constraint equations, and use dofUsed to remember the particles used in order to speed up the constraint correction
    virtual void buildLMatrices          ( ConstOrder Order, const helper::vector< core::behavior::BaseLMConstraint* > &LMConstraints,
            DofToMatrix &LMatrices, DofToMask &dofUsed) const;
    /// Construct the Left Matrix A=sum_i(L_i.M^{-1}_i.L_i^T), store the matrices M^{-1}_i.L_i^T in order to compute later the constraint correction to apply
    virtual void buildLeftMatrix         ( const DofToMatrix& invMassMatrices,
            DofToMatrix& LMatrices, SparseMatrixEigen &LeftMatrix, DofToMatrix &invMass_Ltrans) const;
    /// Solve the System using a projective Gauss-Seidel algorithm: compute the Lagrange Multipliers Lambda
    virtual bool solveConstraintSystemUsingGaussSeidel(MultiVecId id, ConstOrder Order,
            const helper::vector< core::behavior::BaseLMConstraint* > &LMConstraints,
            const MatrixEigen &W,
            const VectorEigen &c,
            VectorEigen &Lambda);

    /// Compute Kinetic Energy
    virtual void computeKineticEnergy(MultiVecId id);

    /** Apply the correction to the state corresponding
     * @param id nature of the constraint, and correction to apply
     * @param dof MechanicalState to correct
     * @param invM_Ltrans matrix M^-1.L^T to apply the correction from the independant dofs through the mapping
     * @param c correction vector
     * @param propageVelocityChange need to propagate the correction done to the velocity for the position
     **/
    virtual void constraintStateCorrection(VecId id, core::ConstraintParams::ConstOrder order,
            bool isPositionChangesUpdateVelocity,
            const SparseMatrixEigen  &invM_Ltrans,
            const VectorEigen  &Lambda,
            const sofa::helper::set< unsigned int > &dofUsed,
            sofa::core::behavior::BaseMechanicalState* dof) const;


    ///
    virtual void buildLMatrix          ( const sofa::core::behavior::BaseMechanicalState *dof,
            const std::list<unsigned int> &idxEquations,unsigned int constraintOffset,
            SparseMatrixEigen& L, sofa::helper::set< unsigned int > &dofUsed ) const;
    virtual void buildInverseMassMatrix( const sofa::core::behavior::BaseMechanicalState* mstate,
            const core::behavior::BaseConstraintCorrection* constraintCorrection,
            SparseMatrixEigen& matrix) const;
    virtual void buildInverseMassMatrix( const sofa::core::behavior::BaseMechanicalState* mstate,
            const core::behavior::BaseMass* mass,
            SparseMatrixEigen& matrix) const;

    core::ConstraintParams::ConstOrder orderState;
    unsigned int numConstraint;

    //Variables used to do the computation
    sofa::simulation::MechanicalWriteLMConstraint LMConstraintVisitor;
    SetDof setDofs;
    DofToMask dofUsed;
    DofToMatrix invMass_Ltrans;
    DofToMatrix LMatrices;

    MatrixEigen W;
    VectorEigen c;
    VectorEigen Lambda;

    //Persitent datas
    DofToMatrix invMassMatrix;
    DofToConstraintCorrection constraintCorrections;

};

} // namespace constraintset

} // namespace component

} // namespace sofa

#endif
