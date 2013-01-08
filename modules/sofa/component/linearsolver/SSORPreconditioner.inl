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
// Author: François Faure, INRIA-UJF, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
#ifndef SOFA_COMPONENT_LINEARSOLVER_SSORPRECONDITIONER_INL
#define SOFA_COMPONENT_LINEARSOLVER_SSORPRECONDITIONER_INL
#include <sofa/component/linearsolver/SSORPreconditioner.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/component/linearsolver/NewMatMatrix.h>
#include <sofa/component/linearsolver/FullMatrix.h>
#include <sofa/component/linearsolver/SparseMatrix.h>
#include <sofa/component/linearsolver/CompressedRowSparseMatrix.h>
#include <iostream>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/behavior/LinearSolver.h>
#include <math.h>
#include <sofa/helper/system/thread/CTime.h>

namespace sofa
{

namespace component
{

namespace linearsolver
{

using namespace sofa::defaulttype;
using namespace sofa::core::behavior;
using namespace sofa::simulation;
using namespace sofa::core::objectmodel;

template<class TMatrix, class TVector, class TThreadManager>
SSORPreconditioner<TMatrix,TVector,TThreadManager>::SSORPreconditioner()
    : f_verbose( initData(&f_verbose,false,"verbose","Dump system state at each iteration") )
    , f_omega( initData(&f_omega,1.0, "omega","Omega coefficient") )
{
}

// solve (D+U) * D^-1 * ( D + U)
template<class TMatrix, class TVector, class TThreadManager>
void SSORPreconditioner<TMatrix,TVector,TThreadManager>::solve (Matrix& M, Vector& z, Vector& r)
{
    SSORPreconditionerInvertData * data = (SSORPreconditionerInvertData *) getMatrixInvertData(&M);

    const int n = M.rowSize();
    const Real w = (Real)f_omega.getValue();
    //Solve (D/w+u) * u3 = r;
    for (int j=n-1; j>=0; j--)
    {
        double temp = 0.0;
        for (int i=j+1; i<n; i++)
        {
            temp += z[i] * M.element(i,j);
        }
        z[j] = (r[j] - temp) * w * data->inv_diag[j];
    }

    //Solve (I + w D^-1 * L) * z = u3
    for (int j=0; j<n; j++)
    {
        double temp = 0.0;
        for (int i=0; i<j; i++)
        {
            temp += z[i] * M.element(i,j);
        }
        z[j] = z[j] - temp * w * data->inv_diag[j];
        // we can reuse z because all values that we read are updated
    }

    if (w != (Real)1.0)
        for (unsigned j=0; j<M.rowSize(); j++)
            z[j] *= 2-w;

}

template<>
void SSORPreconditioner<SparseMatrix<double>, FullVector<double> >::solve (Matrix& M, Vector& z, Vector& r)
{
    SSORPreconditionerInvertData * data = (SSORPreconditionerInvertData *) getMatrixInvertData(&M);

    const int n = M.rowSize();
    const Real w = (Real)f_omega.getValue();

    //Solve (D/w+U) * t = r;
    for (int j=n-1; j>=0; j--)
    {
        double temp = 0.0;
        for (Matrix::LElementConstIterator it = ++M[j].find(j), end = M[j].end(); it != end; ++it)
        {
            int i = it->first;
            double e = it->second;
            temp += z[i] * e;
        }
        z[j] = (r[j] - temp) * w * data->inv_diag[j];
    }

    //Solve (I + w * D^-1 * L) * z = t
    for (int j=0; j<n; j++)
    {
        double temp = 0.0;
        for (Matrix::LElementConstIterator it = M[j].begin(), end = M[j].find(j); it != end; ++it)
        {
            int i = it->first;
            double e = it->second;
            temp += z[i] * e;
        }
        z[j] -= temp * w * data->inv_diag[j];
        // we can reuse z because all values that we read are updated
    }

    if (w != (Real)1.0)
        for (unsigned j=0; j<M.rowSize(); j++)
            z[j] *= 2-w;
}

template<>
void SSORPreconditioner<CompressedRowSparseMatrix<double>, FullVector<double> >::solve (Matrix& M, Vector& z, Vector& r)
{
    SSORPreconditionerInvertData * data = (SSORPreconditionerInvertData *) getMatrixInvertData(&M);

    const int n = M.rowSize();
    const Real w = (Real)f_omega.getValue();

    //const Matrix::VecIndex& rowIndex = M.getRowIndex();
    const Matrix::VecIndex& colsIndex = M.getColsIndex();
    const Matrix::VecBloc& colsValue = M.getColsValue();
    //Solve (D/w+U) * t = r;
    for (int j=n-1; j>=0; j--)
    {
        double temp = 0.0;
        Matrix::Range rowRange = M.getRowRange(j);
        int xi = rowRange.begin();
        while (xi < rowRange.end() && colsIndex[xi] <= j) ++xi;
        for (; xi < rowRange.end(); ++xi)
        {
            int i = colsIndex[xi];
            double e = colsValue[xi];
            temp += z[i] * e;
        }
        z[j] = (r[j] - temp) * w * data->inv_diag[j];
    }

    //Solve (I + w D^-1 * L) * z = t
    for (int j=0; j<n; j++)
    {
        double temp = 0.0;
        Matrix::Range rowRange = M.getRowRange(j);
        int xi = rowRange.begin();
        for (; xi < rowRange.end() && colsIndex[xi] < j; ++xi)
        {
            int i = colsIndex[xi];
            double e = colsValue[xi];
            temp += z[i] * e;
        }
        z[j] -= temp * w * data->inv_diag[j];
        // we can reuse z because all values that we read are updated
    }

    if (w != (Real)1.0)
        for (unsigned j=0; j<M.rowSize(); j++)
            z[j] *= 2-w;
}

#define B 3
#define Real double
#define typename
//template<int B, class Real>
template<>
void SSORPreconditioner< CompressedRowSparseMatrix< defaulttype::Mat<B,B,Real> >, FullVector<Real> >::solve(Matrix& M, Vector& z, Vector& r)
{
    SSORPreconditionerInvertData * data = (SSORPreconditionerInvertData *) getMatrixInvertData(&M);

    //const int n = M.rowSize();
    const int nb = M.rowBSize();
    const Real w = (Real)f_omega.getValue();

    //const Matrix::VecIndex& rowIndex = M.getRowIndex();
    const typename Matrix::VecIndex& colsIndex = M.getColsIndex();
    const typename Matrix::VecBloc& colsValue = M.getColsValue();
    //Solve (D+U) * t = r;
    for (int jb=nb-1; jb>=0; jb--)
    {
        int j0 = jb*B;
        defaulttype::Vec<B,Real> temp;
        typename Matrix::Range rowRange = M.getRowRange(jb);
        int xi = rowRange.begin();
        while (xi < rowRange.end() && colsIndex[xi] < jb) ++xi;
        // bloc on the diagonal
        const typename Matrix::Bloc& bdiag = colsValue[xi];
        // upper triangle matrix
        for (++xi; xi < rowRange.end(); ++xi)
        {
            int i0 = colsIndex[xi]*B;
            const typename Matrix::Bloc& b = colsValue[xi];
            for (int j1=0; j1<B; ++j1)
            {
                //int j = j0+j1;
                for (int i1=0; i1<B; ++i1)
                {
                    int i = i0+i1;
                    temp[j1] += z[i] * b[j1][i1];
                }
            }
        }
        // then the diagonal
        {
            const typename Matrix::Bloc& b = bdiag;
            for (int j1=B-1; j1>=0; j1--)
            {
                int j = j0+j1;
                for (int i1=j1+1; i1<B; ++i1)
                {
                    int i = j0+i1;
                    temp[j1]+= z[i] * b[j1][i1];
                }
                z[j] = (r[j] - temp[j1]) * w * data->inv_diag[j];
            }
        }
    }

    //Solve (I + D^-1 * L) * z = t
    for (int jb=0; jb<nb; jb++)
    {
        int j0 = jb*B;
        defaulttype::Vec<B,Real> temp;
        typename Matrix::Range rowRange = M.getRowRange(jb);
        int xi = rowRange.begin();
        // lower triangle matrix
        for (; xi < rowRange.end() && colsIndex[xi] < jb; ++xi)
        {
            int i0 = colsIndex[xi]*B;
            const typename Matrix::Bloc& b = colsValue[xi];
            for (int j1=0; j1<B; ++j1)
            {
                //int j = j0+j1;
                for (int i1=0; i1<B; ++i1)
                {
                    int i = i0+i1;
                    temp[j1] += z[i] * b[j1][i1];
                }
            }
        }
        // then the diagonal
        {
            const typename Matrix::Bloc& b = colsValue[xi];
            for (int j1=0; j1<B; ++j1)
            {
                int j = j0+j1;
                for (int i1=0; i1<j1; ++i1)
                {
                    int i = j0+i1;
                    temp[j1] += z[i] * b[j1][i1];
                }
                // we can reuse z because all values that we read are updated
                z[j] -= temp[j1] * w * data->inv_diag[j];
            }
        }
    }
}

#undef B
#undef Real
#undef typename


template<class TMatrix, class TVector, class TThreadManager>
void SSORPreconditioner<TMatrix,TVector,TThreadManager>::invert(Matrix& M)
{
    SSORPreconditionerInvertData * data = (SSORPreconditionerInvertData *) this->getMatrixInvertData(&M);

    int n = M.rowSize();
    data->inv_diag.resize(n);
    for (int j=0; j<n; j++) data->inv_diag[j] = 1.0 / M.element(j,j);
}

} // namespace linearsolver

} // namespace component

} // namespace sofa

#endif
