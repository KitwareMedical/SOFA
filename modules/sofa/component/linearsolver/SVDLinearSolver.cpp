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
// Author: François Faure, INRIA-UJF, (C) 2011
//
// Copyright: See COPYING file that comes with this distribution
#include <sofa/component/linearsolver/SVDLinearSolver.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/component/linearsolver/NewMatMatrix.h>
#include <sofa/component/linearsolver/FullMatrix.h>
#include <sofa/component/linearsolver/SparseMatrix.h>
#include <sofa/component/linearsolver/CompressedRowSparseMatrix.h>
#include <sofa/simulation/common/MechanicalVisitor.h>
#include <sofa/helper/system/thread/CTime.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/core/ObjectFactory.h>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>

namespace sofa
{
namespace component
{
namespace linearsolver
{
using core::VecId;
using namespace sofa::defaulttype;
using namespace sofa::core::behavior;
using namespace sofa::simulation;
#ifdef DISPLAY_TIME
using sofa::helper::system::thread::CTime;
#endif

template<class TMatrix, class TVector>
SVDLinearSolver<TMatrix,TVector>::SVDLinearSolver()
    : f_verbose( initData(&f_verbose,false,"verbose","Dump system state at each iteration") )
    , f_minSingularValue( initData(&f_minSingularValue,(Real)1.0e-6,"minSingularValue","Thershold under which a singular value is set to 0, for the stabilization of ill-conditioned system.") )
    , f_conditionNumber( initData(&f_conditionNumber,(Real)0.0,"conditionNumber","Condition number of the matrix: ratio between the largest and smallest singular values. Computed in method solve.") )
{
#ifdef DISPLAY_TIME
    timeStamp = 1.0 / (double)CTime::getRefTicksPerSec();
#endif
}


/// Solve Mx=b
template<class TMatrix, class TVector>
void SVDLinearSolver<TMatrix,TVector>::solve(Matrix& M, Vector& x, Vector& b)
{
#ifdef SOFA_DUMP_VISITOR_INFO
    simulation::Visitor::printComment("SVD");
#endif
#ifdef DISPLAY_TIME
    CTime * timer;
    double time1 = (double) timer->getTime();
#endif
    const bool printLog = this->f_printLog.getValue();
    const bool verbose  = f_verbose.getValue();

    /// Convert the matrix and the right-hand vector to Eigen objects
    Eigen::MatrixXd m(M.rowSize(),M.colSize());
    Eigen::VectorXd rhs(M.rowSize());
    for(unsigned i=0; i<M.rowSize(); i++ )
    {
        for( unsigned j=0; j<M.colSize(); j++ )
            m(i,j) = M[i][j];
        rhs(i) = b[i];
    }
    if(verbose)
    {
        serr << "SVDLinearSolver<TMatrix,TVector>::solve, Here is the matrix m:" << sendl << m << sendl;
    }

    /// Compute the SVD decomposition and the condition number
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(m, Eigen::ComputeThinU | Eigen::ComputeThinV);
    f_conditionNumber.setValue( (Real)(svd.singularValues()(0) / svd.singularValues()(M.rowSize()-1)) );
    if(printLog)
    {
        serr << "SVDLinearSolver<TMatrix,TVector>::solve, the singular values are:" << sendl << svd.singularValues() << sendl;
    }
    if(verbose)
    {
        serr << "Its left singular vectors are the columns of the thin U matrix:" << sendl << svd.matrixU() << sendl;
        serr << "Its right singular vectors are the columns of the thin V matrix:" << sendl << svd.matrixV() << sendl;
    }

    /// Solve the equation system and copy the solution to the SOFA vector
//    Eigen::VectorXd solution = svd.solve(rhs);
//    for(unsigned i=0; i<M.rowSize(); i++ ){
//        x[i] = solution(i);
//    }
    Eigen::VectorXd Ut_b = svd.matrixU().transpose() *  rhs;
    Eigen::VectorXd S_Ut_b(M.colSize());
    for( unsigned i=0; i<M.colSize(); i++ )   /// product with the diagonal matrix, using the threshold for near-null values
    {
        if( svd.singularValues()[i] > f_minSingularValue.getValue() )
            S_Ut_b[i] = Ut_b[i]/svd.singularValues()[i];
        else
            S_Ut_b[i] = (Real)0.0 ;
    }
    Eigen::VectorXd solution = svd.matrixV() * S_Ut_b;
    for(unsigned i=0; i<M.rowSize(); i++ )
    {
        x[i] = (Real) solution(i);
    }

    if( printLog )
    {
#ifdef DISPLAY_TIME
        time1 = (double)(((double) timer->getTime() - time1) * timeStamp / (nb_iter-1));
        std::cerr<<"SVDLinearSolver::solve, SVD = "<<time1<<std::endl;
#endif
        serr << "SVDLinearSolver<TMatrix,TVector>::solve, rhs vector = " << sendl << rhs.transpose() << sendl;
        serr << "SVDLinearSolver<TMatrix,TVector>::solve, solution = " << sendl << x << sendl;
        serr << "SVDLinearSolver<TMatrix,TVector>::solve, verification, mx - b = " << sendl << (m * solution - rhs ).transpose() << sendl;
    }
}


SOFA_DECL_CLASS(SVDLinearSolver)

int SVDLinearSolverClass = core::RegisterObject("Linear system solver using the conjugate gradient iterative algorithm")
        .add< SVDLinearSolver< FullMatrix<double>, FullVector<double> > >()
        .add< SVDLinearSolver< FullMatrix<float>, FullVector<float> > >()
        .addAlias("SVDLinear")
        .addAlias("SVD")
        ;

} // namespace linearsolver

} // namespace component

} // namespace sofa

