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
#ifndef SOFA_COMPONENT_LINEARSOLVER_SparseLDLSolver_H
#define SOFA_COMPONENT_LINEARSOLVER_SparseLDLSolver_H

#include <sofa/core/behavior/LinearSolver.h>
#include <sofa/simulation/common/MechanicalVisitor.h>
#include <sofa/component/linearsolver/FullMatrix.h>
#include <sofa/component/linearsolver/SparseMatrix.h>
#include <sofa/component/linearsolver/CompressedRowSparseMatrix.h>
#include <sofa/helper/map.h>
#include <math.h>

#include <sofa/component/linearsolver/SparseLDLSolverImpl.h>
#include <sofa/defaulttype/BaseMatrix.h>

namespace sofa
{

namespace component
{

namespace linearsolver
{

/// Direct linear solver based on Sparse LDL^T factorization, implemented with the CSPARSE library
template<class TMatrix, class TVector, class TThreadManager = NoThreadManager>
class SparseLDLSolver : public sofa::component::linearsolver::SparseLDLSolverImpl<TMatrix,TVector, TThreadManager>
{
public :
    SOFA_CLASS(SOFA_TEMPLATE3(SparseLDLSolver,TMatrix,TVector,TThreadManager),SOFA_TEMPLATE3(sofa::component::linearsolver::SparseLDLSolverImpl,TMatrix,TVector,TThreadManager));
    typedef sofa::component::linearsolver::SparseLDLSolverImpl<TMatrix,TVector,TThreadManager> Inherit;

public:
    typedef TMatrix Matrix;
    typedef TVector Vector;
    typedef TThreadManager ThreadManager;
    typedef typename Matrix::Real Real;

    SparseLDLSolver();
    void solve (Matrix& M, Vector& x, Vector& b);
    void invert(Matrix& M);

    MatrixInvertData * createInvertData()
    {
        return new SparseLDLSolverInvertData();
    }

protected :
    helper::vector<Real> B;

    class SparseLDLSolverInvertData : public MatrixInvertData
    {
    public :
        sofa::component::linearsolver::CompressedRowSparseMatrix<Real> Mfiltered;
        int n,nnz;

        helper::vector<Real> values,D;
        helper::vector<int> rowind,colptr;
        helper::vector<int> perm, invperm; //premutation inverse
    };

};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_BUILD_SPARSE_SOLVER)
extern template class SOFA_SPARSE_SOLVER_API SparseLDLSolver< CompressedRowSparseMatrix< double>,FullVector<double> >;
extern template class SOFA_SPARSE_SOLVER_API SparseLDLSolver< CompressedRowSparseMatrix< defaulttype::Mat<3,3,double> >,FullVector<double> >;
extern template class SOFA_SPARSE_SOLVER_API SparseLDLSolver< CompressedRowSparseMatrix< float>,FullVector<float> >;
extern template class SOFA_SPARSE_SOLVER_API SparseLDLSolver< CompressedRowSparseMatrix< defaulttype::Mat<3,3,float> >,FullVector<float> >;
#endif


} // namespace linearsolver

} // namespace component

} // namespace sofa

#endif
