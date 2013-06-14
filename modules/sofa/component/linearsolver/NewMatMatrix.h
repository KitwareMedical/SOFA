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
#ifndef SOFA_COMPONENT_LINEARSOLVER_NEWMATMATRIX_H
#define SOFA_COMPONENT_LINEARSOLVER_NEWMATMATRIX_H

#include <sofa/defaulttype/BaseMatrix.h>
#include <sofa/component/linearsolver/MatrixLinearSolver.h>
#include <sofa/component/linearsolver/SparseMatrix.h>
#include "NewMatVector.h"

namespace sofa
{

namespace component
{

namespace linearsolver
{

//#define NEWMAT_CHECK
//#define NEWMAT_VERBOSE

template<class Mat>
class TNewMatMatrix : public Mat, public defaulttype::BaseMatrix
{
public:
    typedef Mat M;
    //typedef NEWMAT::Matrix SubMatrixType;
    typedef TNewMatMatrix<NEWMAT::Matrix> SubMatrixType;
    typedef TNewMatMatrix<NEWMAT::Matrix> InvMatrixType;
    // return the dimension of submatrices when requesting a given size
    static int getSubMatrixDim(int n) { return n; }
    typedef NEWMAT::LinearEquationSolver LUSolver;
    explicit TNewMatMatrix(int defaultBandWidth = 11)
        : bandWidth(defaultBandWidth)
    {
    }

    void resize(int nbRow, int nbCol)
    {
#ifdef NEWMAT_VERBOSE
        std::cout << /* this->Name()  <<  */": resize("<<nbRow<<","<<nbCol<<")"<<std::endl;
#endif
        M::ReSize(nbRow, nbCol);
        (*this) = 0.0;
    }

    unsigned int rowSize(void) const
    {
        return M::Nrows();
    }

    unsigned int colSize(void) const
    {
        return M::Ncols();
    }

    SReal element(int i, int j) const
    {
#ifdef NEWMAT_CHECK
        if ((unsigned)i >= (unsigned)rowSize() || (unsigned)j >= (unsigned)colSize())
        {
            std::cerr << "ERROR: invalid read access to element ("<<i<<","<<j<<") in "<</* this->Name() <<*/" of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
            return 0.0;
        }
#endif
        return M::element(i,j);
    }

    void set(int i, int j, SReal v)
    {
#ifdef NEWMAT_VERBOSE
        std::cout << /* this->Name()  <<  */"("<<rowSize()<<","<<colSize()<<"): element("<<i<<","<<j<<") = "<<v<<std::endl;
#endif
#ifdef NEWMAT_CHECK
        if ((unsigned)i >= (unsigned)rowSize() || (unsigned)j >= (unsigned)colSize())
        {
            std::cerr << "ERROR: invalid write access to element ("<<i<<","<<j<<") in "<</* this->Name() <<*/" of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
            return;
        }
#endif
        M::element(i,j) = v;
    }

    void add(int i, int j, SReal v)
    {
#ifdef NEWMAT_VERBOSE
        std::cout << /* this->Name()  <<  */"("<<rowSize()<<","<<colSize()<<"): element("<<i<<","<<j<<") += "<<v<<std::endl;
#endif
#ifdef NEWMAT_CHECK
        if ((unsigned)i >= (unsigned)rowSize() || (unsigned)j >= (unsigned)colSize())
        {
            std::cerr << "ERROR: invalid write access to element ("<<i<<","<<j<<") in "<</* this->Name() <<*/" of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
            return;
        }
#endif
        M::element(i,j) += v;
    }

    void clear(int i, int j)
    {
#ifdef NEWMAT_VERBOSE
        std::cout << /* this->Name()  <<  */"("<<rowSize()<<","<<colSize()<<"): element("<<i<<","<<j<<") = 0"<<std::endl;
#endif
#ifdef NEWMAT_CHECK
        if ((unsigned)i >= (unsigned)rowSize() || (unsigned)j >= (unsigned)colSize())
        {
            std::cerr << "ERROR: invalid write access to element ("<<i<<","<<j<<") in "<</* this->Name() <<*/" of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
            return;
        }
#endif
        M::element(i,j) = 0.0;
    }

    void clearRow(int i)
    {
#ifdef NEWMAT_VERBOSE
        std::cout << /* this->Name()  <<  */"("<<rowSize()<<","<<colSize()<<"): row("<<i<<") = 0"<<std::endl;
#endif
#ifdef NEWMAT_CHECK
        if ((unsigned)i >= (unsigned)rowSize())
        {
            std::cerr << "ERROR: invalid write access to row "<<i<<" in "<</* this->Name() <<*/" of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
            return;
        }
#endif
        M::Row(1+i) = 0.0;
    }

    void clearCol(int j)
    {
#ifdef NEWMAT_VERBOSE
        std::cout << /* this->Name()  <<  */"("<<rowSize()<<","<<colSize()<<"): col("<<j<<") = 0"<<std::endl;
#endif
#ifdef NEWMAT_CHECK
        if ((unsigned)j >= (unsigned)colSize())
        {
            std::cerr << "ERROR: invalid write access to column "<<j<<" in "<</* this->Name() <<*/" of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
            return;
        }
#endif
        M::Column(1+j) = 0.0;
    }

    void clearRowCol(int i)
    {
#ifdef NEWMAT_VERBOSE
        std::cout << /* this->Name()  <<  */"("<<rowSize()<<","<<colSize()<<"): row("<<i<<") = 0 and col("<<i<<") = 0"<<std::endl;
#endif
#ifdef NEWMAT_CHECK
        if ((unsigned)i >= (unsigned)rowSize() || (unsigned)i >= (unsigned)colSize())
        {
            std::cerr << "ERROR: invalid write access to row and column "<<i<<" in "<</* this->Name() <<*/" of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
            return;
        }
#endif
        M::Row(1+i) = 0.0;
        M::Column(1+i) = 0.0;
    }

    NEWMAT::GetSubMatrix sub(int i, int j, int nrow, int ncol)
    {
        return M::SubMatrix(i+1,i+nrow,j+1,j+ncol);
    }

    template<class T>
    void getSubMatrix(int i, int j, int nrow, int ncol, T& m)
    {
        m = M::SubMatrix(i+1,i+nrow,j+1,j+ncol);
    }

    template<class T>
    void setSubMatrix(int i, int j, int nrow, int ncol, const T& m)
    {
        M::SubMatrix(i+1,i+nrow,j+1,j+ncol) = m;
    }

    NEWMAT::GetSubMatrix asub(int bi, int bj, int nrow, int ncol)
    {
        return M::SubMatrix(bi*nrow+1,bi*nrow+nrow,bj*ncol+1,bj*ncol+ncol);
    }

    template<class T>
    void getAlignedSubMatrix(int bi, int bj, int nrow, int ncol, T& m)
    {
        m = M::SubMatrix(bi*nrow+1,bi*nrow+nrow,bj*ncol+1,bj*ncol+ncol);
    }

    template<class T>
    void setAlignedSubMatrix(int bi, int bj, int nrow, int ncol, const T& m)
    {
        M::SubMatrix(bi*nrow+1,bi*nrow+nrow,bj*ncol+1,bj*ncol+ncol) = m;
    }

    void solve(NewMatVector *rv, NewMatVector *ov)
    {
#ifdef NEWMAT_VERBOSE
        std::cout << /* this->Name()  <<  */"("<<rowSize()<<","<<colSize()<<"): solve("<<*ov<<") = "<<std::endl;
#endif
        *rv = this->i() * *ov;
#ifdef NEWMAT_VERBOSE
        std::cout << /* this->Name()  <<  */"("<<rowSize()<<","<<colSize()<<"): solve("<<*ov<<") = "<<*rv<<std::endl;
#endif
    }

    virtual void solve(defaulttype::BaseVector *op, defaulttype::BaseVector *res)
    {
        NewMatVector *rv = dynamic_cast<NewMatVector *>(res);
        NewMatVector *ov = dynamic_cast<NewMatVector *>(op);

        assert((ov!=NULL) && (rv!=NULL));
        solve(rv,ov);
    }

    LUSolver* makeLUSolver()
    {
        return new LUSolver(*this);
    }

    void solve(NewMatVector *rv, NewMatVector *ov, LUSolver* solver)
    {
#ifdef NEWMAT_VERBOSE
        std::cout << /* this->Name()  <<  */"("<<rowSize()<<","<<colSize()<<"): solve("<<*ov<<") = "<<std::endl;
#endif
        *rv = solver->i() * *ov;
#ifdef NEWMAT_VERBOSE
        std::cout << /* this->Name()  <<  */"("<<rowSize()<<","<<colSize()<<"): solve("<<*ov<<") = "<<*rv<<std::endl;
#endif
    }

    template<class T>
    void operator=(const T& m) { M::operator=(m); }

    void clear() { (*this) = 0.0; }

    friend std::ostream& operator << (std::ostream& out, const TNewMatMatrix& v )
    {
        int nx = v.Ncols();
        int ny = v.Nrows();
        out << "[";
        for (int y=0; y<ny; ++y)
        {
            out << "\n[";
            for (int x=0; x<nx; ++x)
            {
                out << " " << v.element(y,x);
            }
            out << " ]";
        }
        out << " ]";
        return out;
    }

    static const char* Name();

    int bandWidth;
};

typedef TNewMatMatrix<NEWMAT::Matrix> NewMatMatrix;
typedef TNewMatMatrix<NEWMAT::SymmetricMatrix> NewMatSymmetricMatrix;
typedef TNewMatMatrix<NEWMAT::BandMatrix> NewMatBandMatrix;
typedef TNewMatMatrix<NEWMAT::SymmetricBandMatrix> NewMatSymmetricBandMatrix;

template<>
inline const char* TNewMatMatrix<NEWMAT::Matrix>::Name() { return "NewMat"; }

template<>
inline const char* TNewMatMatrix<NEWMAT::SymmetricMatrix>::Name() { return "NewMatSymmetric"; }

template<>
inline const char* TNewMatMatrix<NEWMAT::BandMatrix>::Name() { return "NewMatBand"; }

template<>
inline const char* TNewMatMatrix<NEWMAT::SymmetricBandMatrix>::Name() { return "NewMatSymmetricBand"; }

template<>
inline void TNewMatMatrix<NEWMAT::SymmetricMatrix>::resize(int nbRow, int nbCol)
{
    if (nbCol != nbRow)
        std::cerr << "ERROR: NEWMAT::SymmetricMatrix must be square, size "<<nbRow<<"x"<<nbCol<<" not supported."<<std::endl;
    M::ReSize(nbRow);
}

template<>
inline void TNewMatMatrix<NEWMAT::BandMatrix>::resize(int nbRow, int nbCol)
{
    if (nbCol != nbRow)
        std::cerr << "ERROR: NEWMAT::BandMatrix must be square, size "<<nbRow<<"x"<<nbCol<<" not supported."<<std::endl;
    M::ReSize(nbRow, bandWidth, bandWidth);
}

template<>
inline void TNewMatMatrix<NEWMAT::SymmetricBandMatrix>::resize(int nbRow, int nbCol)
{
    if (nbCol != nbRow)
        std::cerr << "ERROR: NEWMAT::SymmetricBandMatrix must be square, size "<<nbRow<<"x"<<nbCol<<" not supported."<<std::endl;
    M::ReSize(nbRow, bandWidth);
}

template<>
inline void TNewMatMatrix<NEWMAT::SymmetricMatrix>::set(int i, int j, SReal v)
{
#ifdef NEWMAT_VERBOSE
    std::cout << /* this->Name()  <<  */"("<<rowSize()<<","<<colSize()<<"): element("<<i<<","<<j<<") = "<<v<<std::endl;
#endif
#ifdef NEWMAT_CHECK
    if ((unsigned)i >= (unsigned)rowSize() || (unsigned)j >= (unsigned)colSize())
    {
        std::cerr << "ERROR: invalid write access to element ("<<i<<","<<j<<") in "<</* this->Name() <<*/" of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
        return;
    }
#endif
    if (j <= i)
        M::element(i,j) = v;
}

template<>
inline void TNewMatMatrix<NEWMAT::SymmetricMatrix>::add(int i, int j, SReal v)
{
#ifdef NEWMAT_VERBOSE
    std::cout << /* this->Name()  <<  */"("<<rowSize()<<","<<colSize()<<"): element("<<i<<","<<j<<") += "<<v<<std::endl;
#endif
#ifdef NEWMAT_CHECK
    if ((unsigned)i >= (unsigned)rowSize() || (unsigned)j >= (unsigned)colSize())
    {
        std::cerr << "ERROR: invalid write access to element ("<<i<<","<<j<<") in "<</* this->Name() <<*/" of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
        return;
    }
#endif
    if (j <= i)
        M::element(i,j) += v;
}

template<>
inline SReal TNewMatMatrix<NEWMAT::BandMatrix>::element(int i, int j) const
{
#ifdef NEWMAT_CHECK
    if ((unsigned)i >= (unsigned)rowSize() || (unsigned)j >= (unsigned)colSize())
    {
        std::cerr << "ERROR: invalid read access to element ("<<i<<","<<j<<") in "<</* this->Name() <<*/" of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
        return 0.0;
    }
#endif
    if (j < i-bandWidth || j > i+bandWidth)
        return 0.0;
    else
        return M::element(i,j);
}


template<>
inline void TNewMatMatrix<NEWMAT::BandMatrix>::set(int i, int j, SReal v)
{
#ifdef NEWMAT_VERBOSE
    std::cout << /* this->Name()  <<  */"("<<rowSize()<<","<<colSize()<<"): element("<<i<<","<<j<<") = "<<v<<std::endl;
#endif
#ifdef NEWMAT_CHECK
    if ((unsigned)i >= (unsigned)rowSize() || (unsigned)j >= (unsigned)colSize())
    {
        std::cerr << "ERROR: invalid write access to element ("<<i<<","<<j<<") in "<</* this->Name() <<*/" of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
        return;
    }
#endif
    if (j >= i-bandWidth && j <= i+bandWidth)
    {
        M::element(i,j) = v;
    }
}

template<>
inline void TNewMatMatrix<NEWMAT::BandMatrix>::add(int i, int j, SReal v)
{
#ifdef NEWMAT_VERBOSE
    std::cout << /* this->Name()  <<  */"("<<rowSize()<<","<<colSize()<<"): element("<<i<<","<<j<<") += "<<v<<std::endl;
#endif
#ifdef NEWMAT_CHECK
    if ((unsigned)i >= (unsigned)rowSize() || (unsigned)j >= (unsigned)colSize())
    {
        std::cerr << "ERROR: invalid write access to element ("<<i<<","<<j<<") in "<</* this->Name() <<*/" of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
        return;
    }
#endif
    if (j >= i-bandWidth && j <= i+bandWidth)
    {
        M::element(i,j) += v;
    }
}

template<>
inline SReal TNewMatMatrix<NEWMAT::SymmetricBandMatrix>::element(int i, int j) const
{
#ifdef NEWMAT_CHECK
    if ((unsigned)i >= (unsigned)rowSize() || (unsigned)j >= (unsigned)colSize())
    {
        std::cerr << "ERROR: invalid read access to element ("<<i<<","<<j<<") in "<</* this->Name() <<*/" of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
        return 0.0;
    }
#endif
    if (j < i-bandWidth || j > i+bandWidth)
        return 0.0;
    else
        return M::element(i,j);
}

template<>
inline void TNewMatMatrix<NEWMAT::SymmetricBandMatrix>::set(int i, int j, SReal v)
{
#ifdef NEWMAT_VERBOSE
    std::cout << /* this->Name()  <<  */"("<<rowSize()<<","<<colSize()<<"): element("<<i<<","<<j<<") = "<<v<<std::endl;
#endif
#ifdef NEWMAT_CHECK
    if ((unsigned)i >= (unsigned)rowSize() || (unsigned)j >= (unsigned)colSize())
    {
        std::cerr << "ERROR: invalid write access to element ("<<i<<","<<j<<") in "<</* this->Name() <<*/" of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
        return;
    }
#endif
    if (j < i-bandWidth || j > i+bandWidth)
    {
        std::cerr << "ERROR: trying to set "<<v<<" to element ("<<i<<","<<j<<") in NEWMAT::SymmetricBandMatrix of bandwidth "<<bandWidth<<std::endl;
        return;
    }
    if (j <= i)
        M::element(i,j) = v;
}

template<>
inline void TNewMatMatrix<NEWMAT::SymmetricBandMatrix>::add(int i, int j, SReal v)
{
#ifdef NEWMAT_VERBOSE
    std::cout << /* this->Name()  <<  */"("<<rowSize()<<","<<colSize()<<"): element("<<i<<","<<j<<") += "<<v<<std::endl;
#endif
#ifdef NEWMAT_CHECK
    if ((unsigned)i >= (unsigned)rowSize() || (unsigned)j >= (unsigned)colSize())
    {
        std::cerr << "ERROR: invalid write access to element ("<<i<<","<<j<<") in "<</* this->Name() <<*/" of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
        return;
    }
#endif
    if (j < i-bandWidth || j > i+bandWidth)
    {
        std::cerr << "ERROR: trying to set "<<v<<" to element ("<<i<<","<<j<<") in NEWMAT::SymmetricBandMatrix of bandwidth "<<bandWidth<<std::endl;
        return;
    }
    if (j <= i)
        M::element(i,j) += v;
}


//template<>
//class MatrixLinearSolverInternalData< component::linearsolver::NewMatVector >
//{
//public:
//	typedef NewMatVector::Real Real;
//    typedef SparseMatrix<Real> JMatrixType;
//    typedef defaulttype::BaseMatrix ResMatrixType;

//    Data<int> bandWidth;
//    MatrixLinearSolverInternalData(core::objectmodel::BaseObject* o)
//    : bandWidth( o->initData(&bandWidth, 11, "bandWidth", "width of the band on each side of the diagonal (i.e. total values per lines is 2*bandWidth+1)"))
//    {}
//};

//template<>
//inline component::linearsolver::TNewMatMatrix<NEWMAT::SymmetricBandMatrix>* MatrixLinearSolver< component::linearsolver::TNewMatMatrix<NEWMAT::SymmetricBandMatrix> , component::linearsolver::NewMatVector >::createMatrix()
//{
//    return new component::linearsolver::TNewMatMatrix<NEWMAT::SymmetricBandMatrix>(this->internalData.bandWidth.getValue());
//}


} // namespace linearsolver

} // namespace component

} // namespace sofa

#endif
