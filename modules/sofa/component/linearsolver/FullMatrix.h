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
#ifndef SOFA_COMPONENT_LINEARSOLVER_FULLMATRIX_H
#define SOFA_COMPONENT_LINEARSOLVER_FULLMATRIX_H

#include <sofa/defaulttype/BaseMatrix.h>
#include "FullVector.h"

#include <map>

namespace sofa
{

namespace component
{

namespace linearsolver
{

//#define FULLMATRIX_CHECK
//#define FULLMATRIX_VERBOSE


/// Simple full matrix container
template<typename T>
class FullMatrix : public defaulttype::BaseMatrix
{
public:
    typedef T Real;
    typedef int Index;
    typedef FullVector<Real> Line;
    class LineConstIterator
    {
    public:
        Index first;
        Line second;
        LineConstIterator(Real* p, Index i, Index size, Index pitch) : first(i), second(p+i*pitch,size,pitch) { }
        const Line& operator*() const { return second; }
        const Line* operator->() const { return &second; }
        void operator++() { ++first; second.setptr(second.ptr() + second.capacity()); }
        void operator++(int) { ++first; second.setptr(second.ptr() + second.capacity()); }
        void operator--() { --first; second.setptr(second.ptr() - second.capacity()); }
        void operator--(int) { --first; second.setptr(second.ptr() - second.capacity()); }
        bool operator==(const LineConstIterator& i) const { return i.second.ptr() == second.ptr(); }
        bool operator!=(const LineConstIterator& i) const { return i.second.ptr() != second.ptr(); }
    };
    class LineIterator : public LineConstIterator
    {
    public:
        LineIterator(Real* p, Index i, Index size, Index pitch) : LineConstIterator(p,i,size,pitch) {}
        Line& operator*() { return this->second; }
        Line* operator->() { return &this->second; }
    };
    typedef typename Line::iterator LElementIterator;
    typedef typename Line::const_iterator LElementConstIterator;

protected:
    Real* data;
    Index nRow,nCol;
    Index pitch;
    Index allocsize;

public:

    FullMatrix()
        : data(NULL), nRow(0), nCol(0), pitch(0), allocsize(0)
    {
    }

    FullMatrix(int nbRow, int nbCol)
        : data(new T[nbRow*nbCol]), nRow(nbRow), nCol(nbCol), pitch(nbCol), allocsize(nbRow*nbCol)
    {
    }

    FullMatrix(Real* p, int nbRow, int nbCol)
        : data(p), nRow(nbRow), nCol(nbCol), pitch(nbCol), allocsize(-nbRow*nbCol)
    {
    }

    FullMatrix(Real* p, int nbRow, int nbCol, int pitch)
        : data(p), nRow(nbRow), nCol(nbCol), pitch(pitch), allocsize(-nbRow*pitch)
    {
    }

    ~FullMatrix()
    {
        if (allocsize>0)
            delete[] data;
    }

    Real* ptr() { return data; }
    const Real* ptr() const { return data; }

    LineIterator begin() { return LineIterator(data, 0, nCol, pitch); }
    LineIterator end()   { return LineIterator(data, nRow, nCol, pitch);   }
    LineConstIterator begin() const { return LineConstIterator(data, 0, nCol, pitch); }
    LineConstIterator end()   const { return LineConstIterator(data, nRow, nCol, pitch);   }

    //Line operator[](Index i)
    //{
    //    return Line(data+i*pitch, nCol, pitch);
    //}
    Real* operator[](Index i)
    {
        return data+i*pitch;
    }


    //const Line operator[](Index i) const
    //{
    //    return Line(data+i*pitch, nCol, pitch);
    //}
    const Real* operator[](Index i) const
    {
        return data+i*pitch;
    }

    void resize(int nbRow, int nbCol)
    {
#ifdef FULLMATRIX_VERBOSE
        if (nbRow != rowSize() || nbCol != colSize())
            std::cout << /*this->Name() << */": resize("<<nbRow<<","<<nbCol<<")"<<std::endl;
#endif
        if (nbCol != nCol || nbRow != nRow)
        {
            if (allocsize < 0)
            {
                if (nbRow*nbCol > -allocsize)
                {
                    std::cerr << "ERROR: cannot resize preallocated matrix to size ("<<nbRow<<","<<nbCol<<")"<<std::endl;
                    return;
                }
            }
            else
            {
                if (nbRow*nbCol > allocsize)
                {
                    if (allocsize > 0)
                        delete[] data;
                    allocsize = nbRow*nbCol;
                    data = new Real[allocsize];
                }
            }
            pitch = nbCol;
            nCol = nbCol;
            nRow = nbRow;
        }
        clear();
    }


    unsigned int rowSize(void) const
    {
        return nRow;
    }

    unsigned int colSize(void) const
    {
        return nCol;
    }

    SReal element(int i, int j) const
    {
#ifdef FULLMATRIX_CHECK
        if ((unsigned)i >= (unsigned)rowSize() || (unsigned)j >= (unsigned)colSize())
        {
            std::cerr << "ERROR: invalid read access to element ("<<i<<","<<j<<") in "<</*this->Name()<<*/" of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
            return 0.0;
        }
#endif
        return data[i*pitch+j];
    }

    void set(int i, int j, double v)
    {
#ifdef FULLMATRIX_VERBOSE
        std::cout << /*this->Name() <<*/ "("<<rowSize()<<","<<colSize()<<"): element("<<i<<","<<j<<") = "<<v<<std::endl;
#endif
#ifdef FULLMATRIX_CHECK
        if ((unsigned)i >= (unsigned)rowSize() || (unsigned)j >= (unsigned)colSize())
        {
            std::cerr << "ERROR: invalid write access to element ("<<i<<","<<j<<") in "<</*this->Name()<<*/" of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
            return;
        }
#endif
        data[i*pitch+j] = (Real)v;
    }

    void add(int i, int j, double v)
    {
#ifdef FULLMATRIX_VERBOSE
        std::cout << /*this->Name() << */"("<<rowSize()<<","<<colSize()<<"): element("<<i<<","<<j<<") += "<<v<<std::endl;
#endif
#ifdef FULLMATRIX_CHECK
        if ((unsigned)i >= (unsigned)rowSize() || (unsigned)j >= (unsigned)colSize())
        {
            std::cerr << "ERROR: invalid write access to element ("<<i<<","<<j<<") in "/*<<this->Name()*/<<" of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
            return;
        }
#endif
        data[i*pitch+j] += (Real)v;
    }

    void clear(int i, int j)
    {
#ifdef FULLMATRIX_VERBOSE
        std::cout << /*this->Name() <<*/ "("<<rowSize()<<","<<colSize()<<"): element("<<i<<","<<j<<") = 0"<<std::endl;
#endif
#ifdef FULLMATRIX_CHECK
        if ((unsigned)i >= (unsigned)rowSize() || (unsigned)j >= (unsigned)colSize())
        {
            std::cerr << "ERROR: invalid write access to element ("<<i<<","<<j<<") in "<</*this->Name()<<*/" of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
            return;
        }
#endif
        data[i*pitch+j] = (Real)0;
    }

    void clearRow(int i)
    {
#ifdef FULLMATRIX_VERBOSE
        std::cout << /*this->Name() <<*/ "("<<rowSize()<<","<<colSize()<<"): row("<<i<<") = 0"<<std::endl;
#endif
#ifdef FULLMATRIX_CHECK
        if ((unsigned)i >= (unsigned)rowSize())
        {
            std::cerr << "ERROR: invalid write access to row "<<i<<" in "<</*this->Name()<<*/" of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
            return;
        }
#endif
        for (Index j=0; j<nCol; ++j)
            data[i*pitch+j] = (Real)0;
    }

    void clearCol(int j)
    {
#ifdef FULLMATRIX_VERBOSE
        std::cout <</* this->Name() << */"("<<rowSize()<<","<<colSize()<<"): col("<<j<<") = 0"<<std::endl;
#endif
#ifdef FULLMATRIX_CHECK
        if ((unsigned)j >= (unsigned)colSize())
        {
            std::cerr << "ERROR: invalid write access to column "<<j<<" in "<</*this->Name()<<*/" of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
            return;
        }
#endif
        for (Index i=0; i<nRow; ++i)
            data[i*pitch+j] = (Real)0;
    }

    void clearRowCol(int i)
    {
#ifdef FULLMATRIX_VERBOSE
        std::cout << /*this->Name() << */"("<<rowSize()<<","<<colSize()<<"): row("<<i<<") = 0 and col("<<i<<") = 0"<<std::endl;
#endif
#ifdef FULLMATRIX_CHECK
        if ((unsigned)i >= (unsigned)rowSize() || (unsigned)i >= (unsigned)colSize())
        {
            std::cerr << "ERROR: invalid write access to row and column "<<i<<" in "<</*this->Name()<<*/" of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
            return;
        }
#endif
        clearRow(i);
        clearCol(i);
    }

    void clear()
    {
        //if (pitch == nCol)
        //    std::fill(data, data+nRow*pitch, (Real)0);
        //else
        for (Index i=0; i<nRow; ++i)
            for (Index j=0; j<nCol; ++j)
                data[i*pitch+j] = (Real)0;
    }

    /// matrix-vector product
    /// @returns this * v
    template<class Real2>
    FullVector<Real2> operator*( const FullVector<Real2>& v ) const
    {
        FullVector<Real2> res( rowSize() );
        mul( res, v );
        return res;
    }

    /// matrix-vector product
    /// res = this * v
    template<class Real2>
    void mul( FullVector<Real2>& res,const FullVector<Real2>& b ) const
    {
        for( Index i=0 ; i<nRow ; ++i )
        {
            Real r = 0;
            for( Index j=0 ; j<nCol ; ++j )
                r += data[i*pitch+j] * b[j];
            res[i] = r;
        }
    }

    /// transposed matrix-vector product
    /// res = this^T * v
    template<class Real2>
    void mulT( FullVector<Real2>& res, const FullVector<Real2>& b ) const
    {
        for( Index i=0 ; i<nCol ; ++i )
        {
            Real r = 0;
            for( Index j=0 ; j<nRow ; ++j )
                r += data[j*pitch+i] * b[j];
            res[i] = r;
        }
    }




    /// matrix multiplication
    /// @returns this * m
    template<class Real2>
    FullMatrix<Real2> operator*( const FullMatrix<Real2>& m ) const
    {
        FullMatrix<Real2> res( rowSize(), colSize() );
        mul( res, m );
        return res;
    }

    /// matrix multiplication
    /// res = this * m
    template<class Real2>
    void mul( FullMatrix<Real2>& res, const FullMatrix<Real2>& m ) const
    {
        assert( m.rowSize() == (unsigned)nCol );

        res.resize( nRow, m.colSize() );
        for( Index i=0 ; i<nRow ; ++i )
        {
            for( unsigned j=0 ; j<m.colSize() ; ++j )
            {
                res.set( i, j, element(i,0)*m.element(0,j) );
                for( Index k=1 ; k<nCol; ++k )
                    res.add( i, j, element(i,k)*m.element(k,j) );
            }
        }
    }

    /// transposed matrix multiplication
    /// res = this^T * m
    template<class Real2>
    void mulT( FullMatrix<Real2>& res, const FullMatrix<Real2>& m ) const
    {
        assert( m.rowSize() == (unsigned)nRow );

        res.resize( nCol, m.colSize() );
        for( Index i=0 ; i<nCol ; ++i )
        {
            for( unsigned j=0 ; j<m.colSize() ; ++j )
            {
                res.set( i, j, element(0,i)*m.element(0,j) );
                for( Index k=1 ; k<nRow ; ++k )
                    res.add( i, j, element(k,i)*m.element(k,j) );
            }
        }
    }




    friend std::ostream& operator << (std::ostream& out, const FullMatrix<T>& v )
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

template<> inline const char* FullMatrix<double>::Name() { return "FullMatrix"; }
template<> inline const char* FullMatrix<float>::Name() { return "FullMatrixf"; }


/// Simple full matrix container, with an additionnal pointer per line, to be able do get a T** pointer and use [i][j] directly
template<typename T>
class LPtrFullMatrix : public FullMatrix<T>
{
protected:
    T** ldata;
    int lallocsize;
public:
    LPtrFullMatrix()
        : ldata(NULL), lallocsize(0)
    {
    }

    virtual ~LPtrFullMatrix()
    {
        if (lallocsize > 0)
            delete[] ldata;
    }

    void resize(int nbRow, int nbCol)
    {
        if (nbRow == this->nRow && nbCol == this->nCol)
            this->clear();
        else
        {
            this->FullMatrix<T>::resize(nbRow, nbCol);
            if (nbRow > lallocsize)
            {
                if (lallocsize > 0)
                    delete[] ldata;
                ldata = new T*[nbRow];
                lallocsize = nbRow;
            }
            for (int i=0; i<nbRow; ++i)
                ldata[i] = this->data + i*this->pitch;
        }
    }

    T** lptr() { return ldata; }
    const T** lptr() const { return ldata; }
    //operator T**() { return ldata; }
    //operator const T**() const { return ldata; }
};



} // namespace linearsolver

} // namespace component

} // namespace sofa

#endif
