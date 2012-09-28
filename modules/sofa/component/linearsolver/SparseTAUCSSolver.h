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
#ifndef SOFA_COMPONENT_LINEARSOLVER_SparseTAUCSSolver_H
#define SOFA_COMPONENT_LINEARSOLVER_SparseTAUCSSolver_H

#include <sofa/core/behavior/LinearSolver.h>
#include <sofa/component/linearsolver/MatrixLinearSolver.h>
#include <sofa/simulation/common/MechanicalVisitor.h>
#include <sofa/component/linearsolver/SparseMatrix.h>
#include <sofa/component/linearsolver/FullMatrix.h>
#include <sofa/component/linearsolver/CompressedRowSparseMatrix.h>
#include <sofa/helper/map.h>
#include <math.h>

// include all headers included in taucs.h to fix errors on macx
#ifndef WIN32
#include <complex.h>
#endif

#include <assert.h>
#include <float.h>
#include <stdlib.h>

#include <taucs_lib.h>

namespace sofa
{

namespace component
{

namespace linearsolver
{



/// Direct linear solvers implemented with the TAUCS library
template<class TMatrix, class TVector>
class SparseTAUCSSolver : public sofa::component::linearsolver::MatrixLinearSolver<TMatrix,TVector>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(SparseTAUCSSolver,TMatrix,TVector),SOFA_TEMPLATE2(sofa::component::linearsolver::MatrixLinearSolver,TMatrix,TVector));

    typedef TMatrix Matrix;
    typedef TVector Vector;
    typedef typename Matrix::Real Real;
    typedef sofa::component::linearsolver::MatrixLinearSolver<TMatrix,TVector> Inherit;

    Data< helper::vector<std::string> > f_options;
    Data<bool> f_symmetric;
    Data<bool> f_verbose;
#ifdef SOFA_HAVE_CILK
    Data<unsigned> f_nproc;
#endif

    SparseTAUCSSolver();
    void solve (Matrix& M, Vector& x, Vector& b);
    void invert(Matrix& M);

    MatrixInvertData * createInvertData()
    {
        return new SparseTAUCSSolverInvertData();
    }

public:
    class SparseTAUCSSolverInvertData : public MatrixInvertData
    {
    public :
        CompressedRowSparseMatrix<Real> Mfiltered;
        void* factorization;
        taucs_ccs_matrix matrix_taucs;

        SparseTAUCSSolverInvertData()
        {
            factorization = NULL;
        }

        ~SparseTAUCSSolverInvertData()
        {
            if (factorization) taucs_linsolve(NULL, &factorization, 0, NULL, NULL, NULL, NULL);
        }
    };


};

} // namespace linearsolver

} // namespace component

} // namespace sofa

#endif
