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
#ifndef SOFA_COMPONENT_LINEARSOLVER_BTDLINEARSOLVER_H
#define SOFA_COMPONENT_LINEARSOLVER_BTDLINEARSOLVER_H

#include <sofa/core/behavior/LinearSolver.h>
#include <sofa/component/linearsolver/MatrixLinearSolver.h>
#include <sofa/component/linearsolver/SparseMatrix.h>
#include <sofa/component/linearsolver/FullMatrix.h>
#include <math.h>
#include <sofa/defaulttype/Mat.h>

namespace sofa
{

namespace component
{

namespace linearsolver
{

/// Simple bloc full matrix container (used for InvMatrixType)
template<int N, typename T>
class BlocFullMatrix : public defaulttype::BaseMatrix
{
public:

    enum { BSIZE = N };
    typedef T Real;
    typedef int Index;

    class TransposedBloc
    {
    public:
        const defaulttype::Mat<BSIZE,BSIZE,Real>& m;
        TransposedBloc(const defaulttype::Mat<BSIZE,BSIZE,Real>& m) : m(m) {}
        defaulttype::Vec<BSIZE,Real> operator*(const defaulttype::Vec<BSIZE,Real>& v)
        {
            return m.multTranspose(v);
        }
        defaulttype::Mat<BSIZE,BSIZE,Real> operator-() const
        {
            return -m.transposed();
        }
    };

    class Bloc : public defaulttype::Mat<BSIZE,BSIZE,Real>
    {
    public:
        int Nrows() const { return BSIZE; }
        int Ncols() const { return BSIZE; }
        void resize(int, int)
        {
            clear();
        }
        const T& element(int i, int j) const { return (*this)[i][j]; }
        void set(int i, int j, const T& v) { (*this)[i][j] = v; }
        void add(int i, int j, const T& v) { (*this)[i][j] += v; }
        void operator=(const defaulttype::Mat<BSIZE,BSIZE,Real>& v)
        {
            defaulttype::Mat<BSIZE,BSIZE,Real>::operator=(v);
        }
        defaulttype::Mat<BSIZE,BSIZE,Real> operator-() const
        {
            return defaulttype::Mat<BSIZE,BSIZE,Real>::operator-();
        }
        defaulttype::Mat<BSIZE,BSIZE,Real> operator-(const defaulttype::Mat<BSIZE,BSIZE,Real>& m) const
        {
            return defaulttype::Mat<BSIZE,BSIZE,Real>::operator-(m);
        }
        defaulttype::Vec<BSIZE,Real> operator*(const defaulttype::Vec<BSIZE,Real>& v)
        {
            return defaulttype::Mat<BSIZE,BSIZE,Real>::operator*(v);
        }
        defaulttype::Mat<BSIZE,BSIZE,Real> operator*(const defaulttype::Mat<BSIZE,BSIZE,Real>& m)
        {
            return defaulttype::Mat<BSIZE,BSIZE,Real>::operator*(m);
        }
        defaulttype::Mat<BSIZE,BSIZE,Real> operator*(const Bloc& m)
        {
            return defaulttype::Mat<BSIZE,BSIZE,Real>::operator*(m);
        }
        defaulttype::Mat<BSIZE,BSIZE,Real> operator*(const TransposedBloc& mt)
        {
            return defaulttype::Mat<BSIZE,BSIZE,Real>::operator*(mt.m.transposed());
        }
        TransposedBloc t() const
        {
            return TransposedBloc(*this);
        }
        Bloc i() const
        {
            Bloc r;
            r.invert(*this);
            return r;
        }
    };
    typedef Bloc SubMatrixType;
    typedef FullMatrix<T> InvMatrixType;
    // return the dimension of submatrices when requesting a given size
    static int getSubMatrixDim(int) { return BSIZE; }

protected:
    Bloc* data;
    Index nTRow,nTCol;
    Index nBRow,nBCol;
    Index allocsize;

public:

    BlocFullMatrix()
        : data(NULL), nTRow(0), nTCol(0), nBRow(0), nBCol(0), allocsize(0)
    {
    }

    BlocFullMatrix(int nbRow, int nbCol)
        : data(new T[nbRow*nbCol]), nTRow(nbRow), nTCol(nbCol), nBRow(nbRow/BSIZE), nBCol(nbCol/BSIZE), allocsize((nbCol/BSIZE)*(nbRow/BSIZE))
    {
    }

    ~BlocFullMatrix()
    {
        if (allocsize>0)
            delete data;
    }

    Bloc* ptr() { return data; }
    const Bloc* ptr() const { return data; }

    const Bloc& bloc(int bi, int bj) const
    {
        return data[bi*nBCol + bj];
    }
    Bloc& bloc(int bi, int bj)
    {
        return data[bi*nBCol + bj];
    }

    void resize(int nbRow, int nbCol)
    {
        if (nbCol != nTCol || nbRow != nTRow)
        {
            if (allocsize < 0)
            {
                if ((nbCol/BSIZE)*(nbRow/BSIZE) > -allocsize)
                {
                    std::cerr << "ERROR: cannot resize preallocated matrix to size ("<<nbRow<<","<<nbCol<<")"<<std::endl;
                    return;
                }
            }
            else
            {
                if ((nbCol/BSIZE)*(nbRow/BSIZE) > allocsize)
                {
                    if (allocsize > 0)
                        delete[] data;
                    allocsize = (nbCol/BSIZE)*(nbRow/BSIZE);
                    data = new Bloc[allocsize];
                }
            }
            nTCol = nbCol;
            nTRow = nbRow;
            nBCol = nbCol/BSIZE;
            nBRow = nbRow/BSIZE;
        }
        clear();
    }

    unsigned int rowSize(void) const
    {
        return nTRow;
    }

    unsigned int colSize(void) const
    {
        return nTCol;
    }

    SReal element(int i, int j) const
    {
        int bi = i / BSIZE; i = i % BSIZE;
        int bj = j / BSIZE; j = j % BSIZE;
        return bloc(bi,bj)[i][j];
    }

    const Bloc& asub(int bi, int bj, int, int) const
    {
        return bloc(bi,bj);
    }

    const Bloc& sub(int i, int j, int, int) const
    {
        return asub(i/BSIZE,j/BSIZE);
    }

    Bloc& asub(int bi, int bj, int, int)
    {
        return bloc(bi,bj);
    }

    Bloc& sub(int i, int j, int, int)
    {
        return asub(i/BSIZE,j/BSIZE);
    }

    template<class B>
    void getSubMatrix(int i, int j, int nrow, int ncol, B& m)
    {
        m = sub(i,j, nrow, ncol);
    }

    template<class B>
    void getAlignedSubMatrix(int bi, int bj, int nrow, int ncol, B& m)
    {
        m = asub(bi, bj, nrow, ncol);
    }

    template<class B>
    void setSubMatrix(int i, int j, int nrow, int ncol, const B& m)
    {
        sub(i,j, nrow, ncol) = m;
    }

    template<class B>
    void setAlignedSubMatrix(int bi, int bj, int nrow, int ncol, const B& m)
    {
        asub(bi, bj, nrow, ncol) = m;
    }

    void set(int i, int j, double v)
    {
        int bi = i / BSIZE; i = i % BSIZE;
        int bj = j / BSIZE; j = j % BSIZE;
        bloc(bi,bj)[i][j] = (Real)v;
    }

    void add(int i, int j, double v)
    {
        int bi = i / BSIZE; i = i % BSIZE;
        int bj = j / BSIZE; j = j % BSIZE;
        bloc(bi,bj)[i][j] += (Real)v;
    }

    void clear(int i, int j)
    {
        int bi = i / BSIZE; i = i % BSIZE;
        int bj = j / BSIZE; j = j % BSIZE;
        bloc(bi,bj)[i][j] = (Real)0;
    }

    void clearRow(int i)
    {
        int bi = i / BSIZE; i = i % BSIZE;
        for (int bj = 0; bj < nBCol; ++bj)
            for (int j=0; j<BSIZE; ++j)
                bloc(bi,bj)[i][j] = (Real)0;
    }

    void clearCol(int j)
    {
        int bj = j / BSIZE; j = j % BSIZE;
        for (int bi = 0; bi < nBRow; ++bi)
            for (int i=0; i<BSIZE; ++i)
                bloc(bi,bj)[i][j] = (Real)0;
    }

    void clearRowCol(int i)
    {
        clearRow(i);
        clearCol(i);
    }

    void clear()
    {
        for (Index i=0; i<3*nBRow; ++i)
            data[i].clear();
    }

    template<class Real2>
    FullVector<Real2> operator*(const FullVector<Real2>& v) const
    {
        FullVector<Real2> res(rowSize());
        for (int bi=0; bi<nBRow; ++bi)
        {
            int bj = 0;
            for (int i=0; i<BSIZE; ++i)
            {
                Real r = 0;
                for (int j=0; j<BSIZE; ++j)
                {
                    r += bloc(bi,bj)[i][j] * v[(bi + bj - 1)*BSIZE + j];
                }
                res[bi*BSIZE + i] = r;
            }
            for (++bj; bj<nBCol; ++bj)
            {
                for (int i=0; i<BSIZE; ++i)
                {
                    Real r = 0;
                    for (int j=0; j<BSIZE; ++j)
                    {
                        r += bloc(bi,bj)[i][j] * v[(bi + bj - 1)*BSIZE + j];
                    }
                    res[bi*BSIZE + i] += r;
                }
            }
        }
        return res;
    }

    friend std::ostream& operator << (std::ostream& out, const BlocFullMatrix<N,T>& v)
    {
        int nx = v.colSize();
        int ny = v.rowSize();
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
};

template<int N, typename T>
class BlockVector : public FullVector<T>
{
public:
    typedef FullVector<T> Inherit;
    typedef T Real;
    typedef typename Inherit::Index Index;

    typedef typename Inherit::value_type value_type;
    typedef typename Inherit::size_type size_type;
    typedef typename Inherit::iterator iterator;
    typedef typename Inherit::const_iterator const_iterator;

    class Bloc : public defaulttype::Vec<N,T>
    {
    public:
        int Nrows() const { return N; }
        void resize(int) { this->clear(); }
        void operator=(const defaulttype::Vec<N,T>& v)
        {
            defaulttype::Vec<N,T>::operator=(v);
        }
        void operator=(int v)
        {
            defaulttype::Vec<N,T>::fill((float)v);
        }
        void operator=(float v)
        {
            defaulttype::Vec<N,T>::fill(v);
        }
        void operator=(double v)
        {
            defaulttype::Vec<N,T>::fill(v);
        }
    };

    typedef Bloc SubVectorType;

public:

    BlockVector()
    {
    }

    explicit BlockVector(Index n)
        : Inherit(n)
    {
    }

    virtual ~BlockVector()
    {
    }

    const Bloc& sub(int i, int) const
    {
        return (const Bloc&)*(this->ptr()+i);
    }

    Bloc& sub(int i, int)
    {
        return (Bloc&)*(this->ptr()+i);
    }

    const Bloc& asub(int bi, int) const
    {
        return (const Bloc&)*(this->ptr()+bi*N);
    }

    Bloc& asub(int bi, int)
    {
        return (Bloc&)*(this->ptr()+bi*N);
    }
};

/// Simple BTD matrix container
template<int N, typename T>
class BTDMatrix : public defaulttype::BaseMatrix
{
public:
    enum { BSIZE = N };
    typedef T Real;
    typedef int Index;


    class TransposedBloc
    {
    public:
        const defaulttype::Mat<BSIZE,BSIZE,Real>& m;
        TransposedBloc(const defaulttype::Mat<BSIZE,BSIZE,Real>& m) : m(m) {}
        defaulttype::Vec<BSIZE,Real> operator*(const defaulttype::Vec<BSIZE,Real>& v)
        {
            return m.multTranspose(v);
        }
        defaulttype::Mat<BSIZE,BSIZE,Real> operator-() const
        {
            defaulttype::Mat<BSIZE,BSIZE,Real> r;
            for (int i=0; i<BSIZE; i++)
                for (int j=0; j<BSIZE; j++)
                    r[i][j]=-m[j][i];
            return r;
        }
    };

    class Bloc : public defaulttype::Mat<BSIZE,BSIZE,Real>
    {
    public:
        int Nrows() const { return BSIZE; }
        int Ncols() const { return BSIZE; }
        void resize(int, int)
        {
            clear();
        }
        const T& element(int i, int j) const { return (*this)[i][j]; }
        void set(int i, int j, const T& v) { (*this)[i][j] = v; }
        void add(int i, int j, const T& v) { (*this)[i][j] += v; }
        void operator=(const defaulttype::Mat<BSIZE,BSIZE,Real>& v)
        {
            defaulttype::Mat<BSIZE,BSIZE,Real>::operator=(v);
        }
        defaulttype::Mat<BSIZE,BSIZE,Real> operator-() const
        {
            defaulttype::Mat<BSIZE,BSIZE,Real> r;
            for (int i=0; i<BSIZE; i++)
                for (int j=0; j<BSIZE; j++)
                    r[i][j]=-(*this)[i][j];
            return r;
        }
        defaulttype::Mat<BSIZE,BSIZE,Real> operator-(const defaulttype::Mat<BSIZE,BSIZE,Real>& m) const
        {
            return defaulttype::Mat<BSIZE,BSIZE,Real>::operator-(m);
        }
        defaulttype::Vec<BSIZE,Real> operator*(const defaulttype::Vec<BSIZE,Real>& v)
        {
            return defaulttype::Mat<BSIZE,BSIZE,Real>::operator*(v);
        }
        defaulttype::Mat<BSIZE,BSIZE,Real> operator*(const defaulttype::Mat<BSIZE,BSIZE,Real>& m)
        {
            return defaulttype::Mat<BSIZE,BSIZE,Real>::operator*(m);
        }
        defaulttype::Mat<BSIZE,BSIZE,Real> operator*(const Bloc& m)
        {
            return defaulttype::Mat<BSIZE,BSIZE,Real>::operator*(m);
        }
        defaulttype::Mat<BSIZE,BSIZE,Real> operator*(const TransposedBloc& mt)
        {
            return defaulttype::Mat<BSIZE,BSIZE,Real>::operator*(mt.m.transposed());
        }
        TransposedBloc t() const
        {
            return TransposedBloc(*this);
        }
        Bloc i() const
        {
            Bloc r;
            r.invert(*this);
            return r;
        }
    };

    typedef Bloc SubMatrixType;
    typedef sofa::defaulttype::Mat<N,N,Real> BlocType;
    typedef BlocFullMatrix<N,T> InvMatrixType;
    // return the dimension of submatrices when requesting a given size
    static int getSubMatrixDim(int) { return BSIZE; }

protected:
    Bloc* data;
    Index nTRow,nTCol;
    Index nBRow,nBCol;
    Index allocsize;

public:

    BTDMatrix()
        : data(NULL), nTRow(0), nTCol(0), nBRow(0), nBCol(0), allocsize(0)
    {
    }

    BTDMatrix(int nbRow, int nbCol)
        : data(new T[3*(nbRow/BSIZE)]), nTRow(nbRow), nTCol(nbCol), nBRow(nbRow/BSIZE), nBCol(nbCol/BSIZE), allocsize(3*(nbRow/BSIZE))
    {
    }

    ~BTDMatrix()
    {
        if (allocsize>0)
            delete data;
    }

    Bloc* ptr() { return data; }
    const Bloc* ptr() const { return data; }

    //Real* operator[](Index i)
    //{
    //    return data+i*pitch;
    //}
    const Bloc& bloc(int bi, int bj) const
    {
        return data[3*bi + (bj - bi + 1)];
    }
    Bloc& bloc(int bi, int bj)
    {
        return data[3*bi + (bj - bi + 1)];
    }

    void resize(int nbRow, int nbCol)
    {
        if (nbCol != nTCol || nbRow != nTRow)
        {
            if (allocsize < 0)
            {
                if ((nbRow/BSIZE)*3 > -allocsize)
                {
                    std::cerr << "ERROR: cannot resize preallocated matrix to size ("<<nbRow<<","<<nbCol<<")"<<std::endl;
                    return;
                }
            }
            else
            {
                if ((nbRow/BSIZE)*3 > allocsize)
                {
                    if (allocsize > 0)
                        delete[] data;
                    allocsize = (nbRow/BSIZE)*3;
                    data = new Bloc[allocsize];
                }
            }
            nTCol = nbCol;
            nTRow = nbRow;
            nBCol = nbCol/BSIZE;
            nBRow = nbRow/BSIZE;
        }
        clear();
    }

    unsigned int rowSize(void) const
    {
        return nTRow;
    }

    unsigned int colSize(void) const
    {
        return nTCol;
    }

    SReal element(int i, int j) const
    {
        int bi = i / BSIZE; i = i % BSIZE;
        int bj = j / BSIZE; j = j % BSIZE;
        int bindex = bj - bi + 1;
        if ((unsigned)bindex >= 3) return (SReal)0;
        return data[bi*3+bindex][i][j];
    }

    const Bloc& asub(int bi, int bj, int, int) const
    {
        static Bloc b;
        int bindex = bj - bi + 1;
        if ((unsigned)bindex >= 3) return b;
        return data[bi*3+bindex];
    }

    const Bloc& sub(int i, int j, int, int) const
    {
        return asub(i/BSIZE,j/BSIZE);
    }

    Bloc& asub(int bi, int bj, int, int)
    {
        static Bloc b;
        int bindex = bj - bi + 1;
        if ((unsigned)bindex >= 3) return b;
        return data[bi*3+bindex];
    }

    Bloc& sub(int i, int j, int, int)
    {
        return asub(i/BSIZE,j/BSIZE);
    }

    template<class B>
    void getSubMatrix(int i, int j, int nrow, int ncol, B& m)
    {
        m = sub(i,j, nrow, ncol);
    }

    template<class B>
    void getAlignedSubMatrix(int bi, int bj, int nrow, int ncol, B& m)
    {
        m = asub(bi, bj, nrow, ncol);
    }

    template<class B>
    void setSubMatrix(int i, int j, int nrow, int ncol, const B& m)
    {
        sub(i,j, nrow, ncol) = m;
    }

    template<class B>
    void setAlignedSubMatrix(int bi, int bj, int nrow, int ncol, const B& m)
    {
        asub(bi, bj, nrow, ncol) = m;
    }

    void set(int i, int j, double v)
    {
        int bi = i / BSIZE; i = i % BSIZE;
        int bj = j / BSIZE; j = j % BSIZE;
        int bindex = bj - bi + 1;
        if ((unsigned)bindex >= 3) return;
        data[bi*3+bindex][i][j] = (Real)v;
    }

    void add(int i, int j, double v)
    {
        int bi = i / BSIZE; i = i % BSIZE;
        int bj = j / BSIZE; j = j % BSIZE;
        int bindex = bj - bi + 1;
        if ((unsigned)bindex >= 3) return;
        data[bi*3+bindex][i][j] += (Real)v;
    }

    void clear(int i, int j)
    {
        int bi = i / BSIZE; i = i % BSIZE;
        int bj = j / BSIZE; j = j % BSIZE;
        int bindex = bj - bi + 1;
        if ((unsigned)bindex >= 3) return;
        data[bi*3+bindex][i][j] = (Real)0;
    }

    void clearRow(int i)
    {
        int bi = i / BSIZE; i = i % BSIZE;
        for (int bj = 0; bj < 3; ++bj)
            for (int j=0; j<BSIZE; ++j)
                data[bi*3+bj][i][j] = (Real)0;
    }

    void clearCol(int j)
    {
        int bj = j / BSIZE; j = j % BSIZE;
        if (bj > 0)
            for (int i=0; i<BSIZE; ++i)
                data[(bj-1)*3+2][i][j] = (Real)0;
        for (int i=0; i<BSIZE; ++i)
            data[bj*3+1][i][j] = (Real)0;
        if (bj < nBRow-1)
            for (int i=0; i<BSIZE; ++i)
                data[(bj+1)*3+0][i][j] = (Real)0;
    }

    void clearRowCol(int i)
    {
        clearRow(i);
        clearCol(i);
    }

    void clear()
    {
        for (Index i=0; i<3*nBRow; ++i)
            data[i].clear();
    }

    template<class Real2>
    FullVector<Real2> operator*(const FullVector<Real2>& v) const
    {
        FullVector<Real2> res(rowSize());
        for (int bi=0; bi<nBRow; ++bi)
        {
            int b0 = (bi > 0) ? 0 : 1;
            int b1 = ((bi < nBRow - 1) ? 3 : 2);
            for (int i=0; i<BSIZE; ++i)
            {
                Real r = 0;
                for (int bj = b0; bj < b1; ++bj)
                {
                    for (int j=0; j<BSIZE; ++j)
                    {
                        r += data[bi*3+bj][i][j] * v[(bi + bj - 1)*BSIZE + j];
                    }
                }
                res[bi*BSIZE + i] = r;
            }
        }
        return res;
    }

    friend std::ostream& operator << (std::ostream& out, const BTDMatrix<N,T>& v)
    {
        int nx = v.colSize();
        int ny = v.rowSize();
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
};


/// Linear system solver using Thomas Algorithm for Block Tridiagonal matrices
///
/// References:
/// Conte, S.D., and deBoor, C. (1972). Elementary Numerical Analysis. McGraw-Hill, New York
/// http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
/// http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)
/// http://www4.ncsu.edu/eos/users/w/white/www/white/ma580/chap2.5.PDF
template<class Matrix, class Vector>
class BTDLinearSolver : public sofa::component::linearsolver::MatrixLinearSolver<Matrix,Vector>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(BTDLinearSolver, Matrix, Vector), SOFA_TEMPLATE2(sofa::component::linearsolver::MatrixLinearSolver, Matrix, Vector));

    Data<bool> f_verbose;
    Data<bool> problem;
    Data<bool> subpartSolve;

    Data<bool> verification;
    Data<bool> test_perf;

    typedef typename Vector::SubVectorType SubVector;
    typedef typename Matrix::SubMatrixType SubMatrix;
    typedef typename Vector::Real Real;
    typedef typename Matrix::BlocType BlocType;
    typedef std::list<int> ListIndex;
    typedef std::pair<int,int> IndexPair;
    typedef std::map<IndexPair, SubMatrix> MysparseM;
    typedef typename std::map<IndexPair, SubMatrix>::iterator MysparseMit;

    //helper::vector<SubMatrix> alpha;
    helper::vector<SubMatrix> alpha_inv;
    helper::vector<SubMatrix> lambda;
    helper::vector<SubMatrix> B;
    typename Matrix::InvMatrixType Minv;  //inverse matrix


    //////////////////////////// for subpartSolve
    MysparseM H; // force transfer
    MysparseMit H_it;
    Vector bwdContributionOnLH;  //
    Vector fwdContributionOnRH;

    Vector _rh_buf;		 //				// buf the right hand term
    //Vector _df_buf;		 //
    SubVector _acc_rh_bloc;		// accumulation of rh through the browsing of the structure
    SubVector _acc_lh_bloc;		// accumulation of lh through the browsing of the strucutre
    int	current_bloc, first_block;
    std::vector<SubVector> Vec_dRH;			// buf the dRH on block that are not current_bloc...
    ////////////////////////////

    helper::vector<int> nBlockComputedMinv;
    Vector Y;

    Data<int> f_blockSize;
protected:
    BTDLinearSolver()
        : f_verbose( initData(&f_verbose,false,"verbose","Dump system state at each iteration") )
        , problem(initData(&problem, false,"showProblem", "display debug informations about subpartSolve computation") )
        , subpartSolve(initData(&subpartSolve, false,"subpartSolve", "Allows for the computation of a subpart of the system") )
        , verification(initData(&verification, false,"verification", "verification of the subpartSolve"))
        , test_perf(initData(&test_perf, false,"test_perf", "verification of performance"))
        , f_blockSize( initData(&f_blockSize,6,"blockSize","dimension of the blocks in the matrix") )
    {
        int bsize = Matrix::getSubMatrixDim(0);
        if (bsize > 0)
        {
            // the template uses fixed bloc size
            f_blockSize.setValue(bsize);
            f_blockSize.setReadOnly(true);
        }
    }
public:
    void my_identity(SubMatrix& Id, const int size_id);

    void invert(SubMatrix& Inv, const BlocType& m);

    void invert(Matrix& M);

    void computeMinvBlock(int i, int j);

    double getMinvElement(int i, int j);

    /// Solve Mx=b
    void solve (Matrix& /*M*/, Vector& x, Vector& b);



    /// Multiply the inverse of the system matrix by the transpose of the given matrix, and multiply the result with the given matrix J
    ///
    /// @param result the variable where the result will be added
    /// @param J the matrix J to use
    /// @return false if the solver does not support this operation, of it the system matrix is not invertible
    bool addJMInvJt(defaulttype::BaseMatrix* result, defaulttype::BaseMatrix* J, double fact);


    /// Init the partial solve
    void init_partial_solve();

    /// partial solve :
    /// b is accumulated
    /// db is a sparse vector that is added to b
    /// partial_x is a sparse vector (with sparse map given) that provide the result of M x = b+db
    /// Solve Mx=b
    //void partial_solve_old(ListIndex&  Iout, ListIndex&  Iin , bool NewIn);
    void partial_solve(ListIndex&  Iout, ListIndex&  Iin , bool NewIn);



    void init_partial_inverse(const int &nb, const int &bsize);



    template<class RMatrix, class JMatrix>
    bool addJMInvJt(RMatrix& result, JMatrix& J, double fact);



private:


    int _indMaxNonNullForce; // point with non null force which index is the greatest and for which globalAccumulate was not proceed

    int _indMaxFwdLHComputed;  // indice of node from which bwdLH is accurate

    /// private functions for partial solve
    /// step1=> accumulate RH locally for the InBloc (only if a new force is detected on RH)
    void bwdAccumulateRHinBloc(int indMaxBloc);   // indMaxBloc should be equal to _indMaxNonNullForce


    /// step2=> accumulate LH globally to step down the value of current_bloc to 0
    void bwdAccumulateLHGlobal( );


    /// step3=> accumulate RH globally to step up the value of current_bloc to the smallest value needed in OutBloc
    void fwdAccumulateRHGlobal(int indMinBloc);


    /// step4=> compute solution for the indices in the bloc
    /// (and accumulate the potential local dRH (set in Vec_dRH) [set in step1] that have not been yet taken into account by the global bwd and fwd
    void fwdComputeLHinBloc(int indMaxBloc);














};

} // namespace linearsolver

} // namespace component

} // namespace sofa

#endif
