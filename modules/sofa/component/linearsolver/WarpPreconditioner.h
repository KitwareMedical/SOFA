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
#ifndef SOFA_COMPONENT_LINEARSOLVER_WARPPRECONDITIONER_H
#define SOFA_COMPONENT_LINEARSOLVER_WARPPRECONDITIONER_H

#include <sofa/simulation/common/MechanicalVisitor.h>
#include <sofa/core/behavior/LinearSolver.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/component/forcefield/TetrahedronFEMForceField.h>
#include <sofa/defaulttype/Mat.h>
#include <sofa/component/linearsolver/FullVector.h>
#include <math.h>
#include <sofa/core/behavior/RotationMatrix.h>
#include <sofa/core/behavior/BaseRotationFinder.h>
#include <sofa/core/behavior/RotationMatrix.h>
#include <sofa/component/linearsolver/MatrixLinearSolver.h>

#include <map>

namespace sofa
{

namespace component
{

namespace linearsolver
{

/// Linear system solver wrapping another (precomputed) linear solver by a per-node rotation matrix
template<class TMatrix, class TVector, class ThreadManager = NoThreadManager>
class WarpPreconditioner : public sofa::component::linearsolver::MatrixLinearSolver<TMatrix,TVector,ThreadManager>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE3(WarpPreconditioner,TMatrix,TVector,ThreadManager),SOFA_TEMPLATE3(sofa::component::linearsolver::MatrixLinearSolver,TMatrix,TVector,ThreadManager));
    typedef TMatrix Matrix;
    typedef TVector Vector;
    typedef typename TMatrix::Real Real;
    typedef sofa::component::linearsolver::MatrixLinearSolver<TMatrix,TVector,ThreadManager> Inherit;
    typedef sofa::defaulttype::MatNoInit<3, 3, Real> Transformation;
    typedef TMatrix TRotationMatrix;
    typedef typename Inherit::JMatrixType JMatrixType;

    Data <std::string> solverName;
    Data<unsigned> f_useRotationFinder;

protected:
    WarpPreconditioner();

public:

    ~WarpPreconditioner();

    void bwdInit();

    void setSystemMBKMatrix(const core::MechanicalParams* mparams);

    virtual void invert(Matrix& M);

    virtual void solve(Matrix& M, Vector& solution, Vector& rh);

    virtual bool addJMInvJt(defaulttype::BaseMatrix* result, defaulttype::BaseMatrix* J, double fact);

    virtual bool addMInvJt(defaulttype::BaseMatrix* result, defaulttype::BaseMatrix* J, double fact);

    unsigned getSystemDimention(const sofa::core::MechanicalParams* mparams);

    void computeResidual(const core::ExecParams* params, defaulttype::BaseVector* /*f*/);

private :

    core::behavior::LinearSolver* realSolver;

    int updateSystemSize,currentSystemSize;

    int indexwork;
    bool first;

    TRotationMatrix * rotationWork[2];
    std::vector<sofa::core::behavior::BaseRotationFinder *> rotationFinders;

    JMatrixType j_local;
};


} // namespace linearsolver

} // namespace component

} // namespace sofa

#endif
