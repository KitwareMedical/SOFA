/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *
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
*                              SOFA :: Framework                              *
*                                                                             *
* Authors: The SOFA Team (see Authors.txt)                                    *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#include <sofa/helper/system/config.h>
#include <sofa/defaulttype/BaseMatrix.h>
#include <sofa/defaulttype/Mat.h>
#include <sofa/defaulttype/Vec.h>

namespace sofa
{

namespace defaulttype
{

BaseMatrix::BaseMatrix() {}

BaseMatrix::~BaseMatrix()
{}

void BaseMatrix::compress() {
//    std::cerr << "BaseMatrix::compress() : Do Nothing" << std::endl;
}


static inline void opVresize(defaulttype::BaseVector& vec, int n) { vec.resize(n); }
template<class Real2> static inline void opVresize(Real2* vec, int n) { for (const Real2* end=vec+n; vec != end; ++vec) *vec = (Real2)0; }
static inline SReal opVget(const defaulttype::BaseVector& vec, int i) { return (SReal)vec.element(i); }
template<class Real2> static inline SReal opVget(const Real2* vec, int i) { return (SReal)vec[i]; }
static inline void opVset(defaulttype::BaseVector& vec, int i, SReal v) { vec.set(i, v); }
template<class Real2> static inline void opVset(Real2* vec, int i, SReal v) { vec[i] = (Real2)v; }
static inline void opVadd(defaulttype::BaseVector& vec, int i, double v) { vec.add(i, (SReal)v); }
static inline void opVadd(defaulttype::BaseVector& vec, int i, float v) { vec.add(i, (SReal)v); }
template<class Real2> static inline void opVadd(Real2* vec, int i, double v) { vec[i] += (Real2)v; }
template<class Real2> static inline void opVadd(Real2* vec, int i, float v) { vec[i] += (Real2)v; }

template <class Real, int NL, int NC, bool add, bool transpose, class M, class V1, class V2>
struct BaseMatrixLinearOpMV_BlockDiagonal
{
    typedef typename M::RowBlockConstIterator RowBlockConstIterator;
    typedef typename M::ColBlockConstIterator ColBlockConstIterator;
    typedef typename M::BlockConstAccessor BlockConstAccessor;
    typedef Mat<NL,NC,Real> BlockData;
    void operator()(const M* mat, V1& result, const V2& v)
    {
        const int rowSize = mat->rowSize();
        const int colSize = mat->colSize();
        BlockData buffer;
        Vec<NC,Real> vtmpj;
        Vec<NL,Real> vtmpi;
        if (!add)
            opVresize(result, (transpose ? colSize : rowSize));
        for (std::pair<RowBlockConstIterator, RowBlockConstIterator> rowRange = mat->bRowsRange();
                rowRange.first != rowRange.second;
                ++rowRange.first)
        {
            std::pair<ColBlockConstIterator,ColBlockConstIterator> colRange = rowRange.first.range();
            if (colRange.first != colRange.second) // diagonal bloc exists
            {
                BlockConstAccessor block = colRange.first.bloc();
                const BlockData& bdata = *(const BlockData*)block.elements(buffer.ptr());
                const int i = block.getRow() * NL;
                const int j = block.getCol() * NC;
                if (!transpose)
                {
                    VecNoInit<NC,Real> vj;
                    for (int bj = 0; bj < NC; ++bj)
                        vj[bj] = (Real)opVget(v, j+bj);
                    Vec<NL,Real> resi;
                    for (int bi = 0; bi < NL; ++bi)
                        for (int bj = 0; bj < NC; ++bj)
                            resi[bi] += bdata[bi][bj] * vj[bj];
                    for (int bi = 0; bi < NL; ++bi)
                        opVadd(result, i+bi, resi[bi]);
                }
                else
                {
                    VecNoInit<NL,Real> vi;
                    for (int bi = 0; bi < NL; ++bi)
                        vi[bi] = (Real)opVget(v, i+bi);
                    Vec<NC,Real> resj;
                    for (int bi = 0; bi < NL; ++bi)
                        for (int bj = 0; bj < NC; ++bj)
                            resj[bj] += bdata[bi][bj] * vi[bi];
                    for (int bj = 0; bj < NC; ++bj)
                        opVadd(result, j+bj, resj[bj]);
                }
            }
        }
    }
};

// specialication for 1x1 blocs
template <class Real, bool add, bool transpose, class M, class V1, class V2>
struct BaseMatrixLinearOpMV_BlockDiagonal<Real, 1, 1, add, transpose, M, V1, V2>
{
    void operator()(const M* mat, V1& result, const V2& v)
    {
        const int rowSize = mat->rowSize();
        const int colSize = mat->colSize();
        if (!add)
            opVresize(result, (transpose ? colSize : rowSize));
        const int size = (rowSize < colSize) ? rowSize : colSize;
        for (int i=0; i<size; ++i)
        {
            opVadd(result, i, opVget(v, i));
        }
    }
};


template <class Real, int NL, int NC, bool add, bool transpose, class M, class V1, class V2>
struct BaseMatrixLinearOpMV_BlockSparse
{
    typedef typename M::RowBlockConstIterator RowBlockConstIterator;
    typedef typename M::ColBlockConstIterator ColBlockConstIterator;
    typedef typename M::BlockConstAccessor BlockConstAccessor;
    typedef Mat<NL,NC,Real> BlockData;
    void operator()(const M* mat, V1& result, const V2& v)
    {
//         std::cout << "BaseMatrixLinearOpMV_BlockSparse: " << mat->bRowSize() << "x" << mat->bColSize() << " " << NL << "x" << NC << " blocks, " << (add ? "add" : "write") << " to result vector, use " << (transpose ? "transposed " : "") << "matrix." << std::endl;
        const int rowSize = mat->rowSize();
        const int colSize = mat->colSize();
        BlockData buffer;
        Vec<NC,Real> vtmpj;
        Vec<NL,Real> vtmpi;
        if (!add)
            opVresize(result, (transpose ? colSize : rowSize));
        for (std::pair<RowBlockConstIterator, RowBlockConstIterator> rowRange = mat->bRowsRange();
                rowRange.first != rowRange.second;
                ++rowRange.first)
        {
            const int i = rowRange.first.row() * NL;
            if (!transpose)
            {
                for (int bi = 0; bi < NL; ++bi)
                    vtmpi[bi] = (Real)0;
            }
            else
            {
                for (int bi = 0; bi < NL; ++bi)
                    vtmpi[bi] = (Real)opVget(v, i+bi);
            }
            for (std::pair<ColBlockConstIterator,ColBlockConstIterator> colRange = rowRange.first.range();
                    colRange.first != colRange.second;
                    ++colRange.first)
            {
                BlockConstAccessor block = colRange.first.bloc();
                const BlockData& bdata = *(const BlockData*)block.elements(buffer.ptr());
                const int j = block.getCol() * NC;
                if (!transpose)
                {
                    for (int bj = 0; bj < NC; ++bj)
                        vtmpj[bj] = (Real)opVget(v, j+bj);
                    for (int bi = 0; bi < NL; ++bi)
                        for (int bj = 0; bj < NC; ++bj)
                            vtmpi[bi] += bdata[bi][bj] * vtmpj[bj];
                }
                else
                {
                    for (int bj = 0; bj < NC; ++bj)
                        vtmpj[bj] = (Real)0;
                    for (int bi = 0; bi < NL; ++bi)
                        for (int bj = 0; bj < NC; ++bj)
                            vtmpj[bj] += bdata[bi][bj] * vtmpi[bi];
                    for (int bj = 0; bj < NC; ++bj)
                        opVadd(result, j+bj, vtmpj[bj]);
                }
            }
            if (!transpose)
            {
                for (int bi = 0; bi < NL; ++bi)
                    opVadd(result, i+bi, vtmpi[bi]);
            }
            else
            {
            }
        }
    }
};

template<bool add, bool transpose>
class BaseMatrixLinearOpMV
{
public:
    template <class M, class V1, class V2>
    static inline void opFull(const M* mat, V1& result, const V2& v)
    {
        const int rowSize = mat->rowSize();
        const int colSize = mat->colSize();
        if (!add)
            opVresize(result, (transpose ? colSize : rowSize));
        if (!transpose)
        {
            for (int i=0; i<rowSize; ++i)
            {
                double r = 0;
                for (int j=0; j<colSize; ++j)
                {
                    r += mat->element(i,j) * opVget(v, j);
                }
                opVadd(result, i, r);
            }
        }
        else
        {
            for (int i=0; i<rowSize; ++i)
            {
                const double val = opVget(v, i);
                for (int j=0; j<colSize; ++j)
                {
                    opVadd(result, j, mat->element(i,j) * val);
                }
            }
        }
    }

    template <class M, class V1, class V2>
    static inline void opIdentity(const M* mat, V1& result, const V2& v)
    {
        const int rowSize = mat->rowSize();
        const int colSize = mat->colSize();
        if (!add)
            opVresize(result, (transpose ? colSize : rowSize));
        const int size = (rowSize < colSize) ? rowSize : colSize;
        for (int i=0; i<size; ++i)
        {
            opVadd(result, i, opVget(v, i));
        }
    }


    template <int NL, int NC, class M, class V1, class V2>
    static inline void opDiagonal(const M* mat, V1& result, const V2& v)
    {
        const int rowSize = mat->rowSize();
        const int colSize = mat->colSize();
        if (!add)
            opVresize(result, (transpose ? colSize : rowSize));
        const int size = (rowSize < colSize) ? rowSize : colSize;
        for (int i=0; i<size; ++i)
        {
            opVadd(result, i, mat->element(i,i) * opVget(v, i));
        }
    }

    template <class Real, class M, class V1, class V2>
    static inline void opDynamicRealDefault(const M* mat, V1& result, const V2& v, int NL, int NC, BaseMatrix::MatrixCategory /*category*/)
    {
        std::cout << "PERFORMANCE WARNING: multiplication by matric with block size " << NL << "x" << NC << " not optimized." << std::endl;
        opFull(mat, result, v);
    }

    template <class Real, int NL, int NC, class M, class V1, class V2>
    static inline void opDynamicRealNLNC(const M* mat, V1& result, const V2& v, BaseMatrix::MatrixCategory category)
    {
        switch(category)
        {
        case BaseMatrix::MATRIX_DIAGONAL:
        {
            BaseMatrixLinearOpMV_BlockDiagonal<Real, NL, NC, add, transpose, M, V1, V2> op;
            op(mat, result, v);
            break;
        }
        default: // default to sparse
        {
            BaseMatrixLinearOpMV_BlockSparse<Real, NL, NC, add, transpose, M, V1, V2> op;
            op(mat, result, v);
            break;
        }
        }
    }

    template <class Real, int NL, class M, class V1, class V2>
    static inline void opDynamicRealNL(const M* mat, V1& result, const V2& v, int NC, BaseMatrix::MatrixCategory category)
    {
        switch(NC)
        {
        case 1: opDynamicRealNLNC<Real, NL, 1, M, V1, V2>(mat, result, v, category); break;
        case 2: opDynamicRealNLNC<Real, NL, 2, M, V1, V2>(mat, result, v, category); break;
        case 3: opDynamicRealNLNC<Real, NL, 3, M, V1, V2>(mat, result, v, category); break;
        case 4: opDynamicRealNLNC<Real, NL, 4, M, V1, V2>(mat, result, v, category); break;
            //case 5: opDynamicRealNLNC<Real, NL, 5, M, V1, V2>(mat, result, v, category); break;
        case 6: opDynamicRealNLNC<Real, NL, 6, M, V1, V2>(mat, result, v, category); break;
        default: opDynamicRealDefault<Real, M, V1, V2>(mat, result, v, NL, NC, category); break;
        }
    }

    template <class Real, class M, class V1, class V2>
    static inline void opDynamicReal(const M* mat, V1& result, const V2& v, int NL, int NC, BaseMatrix::MatrixCategory category)
    {
        switch(NL)
        {
        case 1: opDynamicRealNL<Real, 1, M, V1, V2>(mat, result, v, NC, category); break;
        case 2: opDynamicRealNL<Real, 2, M, V1, V2>(mat, result, v, NC, category); break;
        case 3: opDynamicRealNL<Real, 3, M, V1, V2>(mat, result, v, NC, category); break;
        case 4: opDynamicRealNL<Real, 4, M, V1, V2>(mat, result, v, NC, category); break;
            //case 5: opDynamicRealNL<Real, 5, M, V1, V2>(mat, result, v, NC, category); break;
        case 6: opDynamicRealNL<Real, 6, M, V1, V2>(mat, result, v, NC, category); break;
        default: opDynamicRealDefault<Real, M, V1, V2>(mat, result, v, NL, NC, category); break;
        }
    }

    template <class M, class V1, class V2>
    static inline void opDynamic(const M* mat, V1& result, const V2& v)
    {
        const int NL = mat->getBlockRows();
        const int NC = mat->getBlockCols();
        const BaseMatrix::MatrixCategory category = mat->getCategory();
        const BaseMatrix::ElementType elementType = mat->getElementType();
        const unsigned int elementSize = mat->getElementSize();
//         std::cout << "BaseMatrixLinearOpMV: " << mat->bRowSize() << "x" << mat->bColSize() << " " << NL << "x" << NC << " blocks, ";
//         switch (category)
//         {
// 	  case BaseMatrix::MATRIX_IDENTITY: std::cout << "identity"; break;
// 	  case BaseMatrix::MATRIX_DIAGONAL: std::cout << "diagonal"; break;
// 	  case BaseMatrix::MATRIX_BAND: std::cout << "band"; break;
// 	  case BaseMatrix::MATRIX_SPARSE: std::cout << "sparse"; break;
// 	  case BaseMatrix::MATRIX_FULL: std::cout << "full"; break;
// 	  case BaseMatrix::MATRIX_UNKNOWN: std::cout << "unknown"; break;
//         }
//         if (elementType == BaseMatrix::ELEMENT_INT) std::cout << " int" << elementSize*8;
//         else if (elementType == BaseMatrix::ELEMENT_FLOAT) std::cout << " float" << elementSize*8;
//         std::cout << " matrix, ";
//         std::cout << (add ? "add" : "write") << " to result vector, use " << (transpose ? "transposed " : "") << "matrix." << std::endl;
        if (category == BaseMatrix::MATRIX_IDENTITY)
        {
            opIdentity(mat, result, v);
        }
        else if (category == BaseMatrix::MATRIX_FULL)
        {
            opFull(mat, result, v);
        }
        //else if (category == BaseMatrix::MATRIX_DIAGONAL && NL == 1 && NC == 1)
        //{
        //    opDiagonal(mat, result, v);
        //}
        //else if (elementType == BaseMatrix::ELEMENT_INT && elementSize == sizeof(int))
        //{
        //    opDynamicReal<int, M, V1, V2>(mat, result, v, NL, NC, category);
        //}
        else if (elementType == BaseMatrix::ELEMENT_FLOAT && elementSize == sizeof(float))
        {
            opDynamicReal<float, M, V1, V2>(mat, result, v, NL, NC, category);
        }
        else // default to double
        {
            opDynamicReal<double, M, V1, V2>(mat, result, v, NL, NC, category);
        }
    }
};


class BaseMatrixLinearOpMulV : public BaseMatrixLinearOpMV<false, false> {};
class BaseMatrixLinearOpPMulV : public BaseMatrixLinearOpMV<true, false> {};
class BaseMatrixLinearOpMulTV : public BaseMatrixLinearOpMV<false, true> {};
class BaseMatrixLinearOpPMulTV : public BaseMatrixLinearOpMV<true, true> {};

/// Multiply the matrix by vector v and put the result in vector result
void BaseMatrix::opMulV(defaulttype::BaseVector* result, const defaulttype::BaseVector* v) const
{
    //std::cout << "BaseMatrix::opMulV(BaseVector)" << std::endl;
    BaseMatrixLinearOpMulV::opDynamic(this, *result, *v);
}

/// Multiply the matrix by float vector v and put the result in vector result
void BaseMatrix::opMulV(float* result, const float* v) const
{
    //std::cout << "BaseMatrix::opMulV(float)" << std::endl;
    BaseMatrixLinearOpMulV::opDynamic(this, result, v);
}

/// Multiply the matrix by double vector v and put the result in vector result
void BaseMatrix::opMulV(double* result, const double* v) const
{
    //std::cout << "BaseMatrix::opMulV(double)" << std::endl;
    BaseMatrixLinearOpMulV::opDynamic(this, result, v);
}

/// Multiply the matrix by vector v and add the result in vector result
void BaseMatrix::opPMulV(defaulttype::BaseVector* result, const defaulttype::BaseVector* v) const
{
    BaseMatrixLinearOpPMulV::opDynamic(this, *result, *v);
}

/// Multiply the matrix by float vector v and add the result in vector result
void BaseMatrix::opPMulV(float* result, const float* v) const
{
    BaseMatrixLinearOpPMulV::opDynamic(this, result, v);
}

/// Multiply the matrix by double vector v and add the result in vector result
void BaseMatrix::opPMulV(double* result, const double* v) const
{
    BaseMatrixLinearOpPMulV::opDynamic(this, result, v);
}


/// Multiply the transposed matrix by vector v and put the result in vector result
void BaseMatrix::opMulTV(defaulttype::BaseVector* result, const defaulttype::BaseVector* v) const
{
    BaseMatrixLinearOpMulTV::opDynamic(this, *result, *v);
}

/// Multiply the transposed matrix by float vector v and put the result in vector result
void BaseMatrix::opMulTV(float* result, const float* v) const
{
    BaseMatrixLinearOpMulTV::opDynamic(this, result, v);
}

/// Multiply the transposed matrix by double vector v and put the result in vector result
void BaseMatrix::opMulTV(double* result, const double* v) const
{
    BaseMatrixLinearOpMulTV::opDynamic(this, result, v);
}

/// Multiply the transposed matrix by vector v and add the result in vector result
void BaseMatrix::opPMulTV(defaulttype::BaseVector* result, const defaulttype::BaseVector* v) const
{
    BaseMatrixLinearOpPMulTV::opDynamic(this, *result, *v);
}

/// Multiply the transposed matrix by float vector v and add the result in vector result
void BaseMatrix::opPMulTV(float* result, const float* v) const
{
    BaseMatrixLinearOpPMulTV::opDynamic(this, result, v);
}

/// Multiply the transposed matrix by double vector v and add the result in vector result
void BaseMatrix::opPMulTV(double* result, const double* v) const
{
    BaseMatrixLinearOpPMulTV::opDynamic(this, result, v);
}

/// Multiply the transposed matrix by matrix m and store the result in matrix result
void BaseMatrix::opMulTM(BaseMatrix * /*result*/,BaseMatrix * /*m*/) const
{
    std::cerr <<"WARNING : BASEMATRIX::opMulTM not yet implemented"<<std::endl;
}








template <class Real, int NL, int NC, bool transpose, class M1, class M2 >
struct BaseMatrixLinearOpAM_BlockSparse
{
    typedef typename M1::RowBlockConstIterator RowBlockConstIterator;
    typedef typename M1::ColBlockConstIterator ColBlockConstIterator;
    typedef typename M1::BlockConstAccessor BlockConstAccessor;
    typedef Mat<NL,NC,Real> BlockData;

    void operator()(const M1* m1, M2* m2, double & fact)
    {
        //const int rowSize = m1->rowSize();
        //const int colSize = m1->colSize();
        BlockData buffer;

        for (std::pair<RowBlockConstIterator, RowBlockConstIterator> rowRange = m1->bRowsRange();
                rowRange.first != rowRange.second;
                ++rowRange.first)
        {
            const int i = rowRange.first.row() * NL;

            for (std::pair<ColBlockConstIterator,ColBlockConstIterator> colRange = rowRange.first.range();
                    colRange.first != colRange.second;
                    ++colRange.first)
            {

                BlockConstAccessor block = colRange.first.bloc();
                const BlockData& bdata = *(const BlockData*)block.elements(buffer.ptr());
                const int j = block.getCol() * NC;

                if (!transpose)
                {
                    for (int bi = 0; bi < NL; ++bi)
                        for (int bj = 0; bj < NC; ++bj)
                            m2->add(i+bi,j+bj,bdata[bi][bj] * fact);
                    //m2->add(i,j,bdata * fact);

                }
                else
                {
                    for (int bi = 0; bi < NL; ++bi)
                        for (int bj = 0; bj < NC; ++bj)
                            m2->add(j+bj,i+bi,bdata[bi][bj] * fact);
                    //m2->add(j,i,bdata * fact);

                }
            }
        }
    }
};


template <class Real, int NL, int NC, bool transpose, class M1 , class M2 >
struct BaseMatrixLinearOpAMS_BlockSparse
{
    typedef typename M1::RowBlockConstIterator RowBlockConstIterator;
    typedef typename M1::ColBlockConstIterator ColBlockConstIterator;
    typedef typename M1::BlockConstAccessor BlockConstAccessor;
    typedef Mat<NL,NC,Real> BlockData;

    void operator()(const M1* m1, M1* m2, double & fact)
    {
        //const int rowSize = m1->rowSize();
        //const int colSize = m1->colSize();
        BlockData buffer;

        for (std::pair<RowBlockConstIterator, RowBlockConstIterator> rowRange = m1->bRowsRange();
                rowRange.first != rowRange.second;
                ++rowRange.first)
        {
            const int i = rowRange.first.row() * NL;

            for (std::pair<ColBlockConstIterator,ColBlockConstIterator> colRange = rowRange.first.range();
                    colRange.first != colRange.second;
                    ++colRange.first)
            {

                BlockConstAccessor block = colRange.first.bloc();
                const BlockData& bdata = *(const BlockData*)block.elements(buffer.ptr());
                const int j = block.getCol() * NC;

//                 if (!transpose) {
// 		    m2->add(i,j,bdata);
//                 } else {
// 		    m2->add(j,i,bdata);
//                 }
                if (!transpose)
                {
                    for (int bi = 0; bi < NL; ++bi)
                        for (int bj = 0; bj < NC; ++bj)
                            m2->add(i+bi,j+bj,bdata[bi][bj] * fact);
                    //m2->add(i,j,bdata * fact);

                }
                else
                {
                    for (int bi = 0; bi < NL; ++bi)
                        for (int bj = 0; bj < NC; ++bj)
                            m2->add(j+bj,i+bi,bdata[bi][bj] * fact);
                    //m2->add(j,i,bdata * fact);
                }
            }
        }
    }
};

template <class Real, bool transpose, class M1, class M2 >
struct BaseMatrixLinearOpAM1_BlockSparse
{
    typedef typename M1::RowBlockConstIterator RowBlockConstIterator;
    typedef typename M1::ColBlockConstIterator ColBlockConstIterator;
    typedef typename M1::BlockConstAccessor BlockConstAccessor;
    typedef Real BlockData;

    void operator()(const M1* m1, M2* m2, double & fact)
    {
        BlockData buffer;

        for (std::pair<RowBlockConstIterator, RowBlockConstIterator> rowRange = m1->bRowsRange();
                rowRange.first != rowRange.second;
                ++rowRange.first)
        {
            const int i = rowRange.first.row();

            for (std::pair<ColBlockConstIterator,ColBlockConstIterator> colRange = rowRange.first.range();
                    colRange.first != colRange.second;
                    ++colRange.first)
            {

                BlockConstAccessor block = colRange.first.bloc();
                const BlockData& bdata = *(const BlockData*)block.elements(&buffer);
                const int j = block.getCol();

                if (!transpose)
                {
                    m2->add(i,j,bdata * fact);
                }
                else
                {
                    m2->add(j,i,bdata * fact);
                }
            }
        }
    }
};

template< bool transpose>
class BaseMatrixLinearOpAM
{
public:
    template <class M1, class M2 >
    static inline void opFull(const M1 * m1, M2 * m2, double fact)
    {
        const int rowSize = m1->rowSize();
        const int colSize = m2->colSize();
        if (!transpose)
        {
            for (int j=0; j<rowSize; ++j)
            {
                for (int i=0; i<colSize; ++i)
                {
                    m2->add(j,i,m1->element(i,j)*fact);
                }
            }
        }
        else
        {
            for (int j=0; j<rowSize; ++j)
            {
                for (int i=0; i<colSize; ++i)
                {
                    m2->add(j,i,m1->element(j,i)*fact);
                }
            }
        }
    }

    template <class M1, class M2 >
    static inline void opIdentity(const M1 * m1, M2 * m2, double fact)
    {
        const int colSize = m1->colSize();
        for (int j=0; j<colSize; ++j)
        {
            m2->add(j,j,fact);
        }
    }


    template <int NL, int NC, class M1, class M2 >
    static inline void opDiagonal(const M1 * m1, M2 * m2, double fact)
    {
        const int colSize = m1->colSize();
        for (int j=0; j<colSize; ++j)
        {
            m2->add(j,j,m1->element(j,j) * fact);
        }
    }

    template <class Real, class M1, class M2 >
    static inline void opDynamicRealDefault(const M1 * m1, M2 * m2, double fact, int NL, int NC, BaseMatrix::MatrixCategory /*category*/)
    {
        std::cout << "PERFORMANCE WARNING: multiplication by matric with block size " << NL << "x" << NC << " not optimized." << std::endl;
        opFull(m1, m2, fact);
    }

    template <class Real, int NL, int NC, class M1, class M2>
    static inline void opDynamicRealNLNC(const M1 * m1, M2 * m2, double fact, BaseMatrix::MatrixCategory /* category */)
    {
//        switch(category)
//        {
// 	  case BaseMatrix::MATRIX_DIAGONAL:
// 	  {
// 	      BaseMatrixLinearOpAM_BlockDiagonal<Real, NL, NC, add, transpose, M, V1, V2> op;
// 	      op(m1, m2, fact);
// 	      break;
// 	  }
//	  default: // default to sparse
        {
            const int NL1 = m2->getBlockRows();
            const int NC2 = m2->getBlockCols();
            if ((NL1==NL) && (NC==NC2))
            {
                BaseMatrixLinearOpAMS_BlockSparse<Real, NL, NC, transpose, M1 , M2 > op;
                op(m1, m2, fact);
            }
            else
            {
                BaseMatrixLinearOpAM_BlockSparse<Real, NL, NC, transpose, M1 , M2 > op;
                op(m1, m2, fact);
            }
//	      break;
//	  }
        }
    }

    template <class Real, class M1, class M2>
    static inline void opDynamicReal1(const M1 * m1, M2 * m2, double fact, BaseMatrix::MatrixCategory /* category */)
    {
//        switch(category)
//        {
// 	  case BaseMatrix::MATRIX_DIAGONAL:
// 	  {
// 	      BaseMatrixLinearOpAM_BlockDiagonal<Real, NL, NC, add, transpose, M, V1, V2> op;
// 	      op(m1, m2, fact);
// 	      break;
// 	  }
//	  default: // default to sparse
        {
            BaseMatrixLinearOpAM1_BlockSparse<Real , transpose, M1, M2 > op;
            op(m1, m2, fact);
//	      break;
//	  }
        }
    }

    template <class Real, int NL, class M1, class M2 >
    static inline void opDynamicRealNL(const M1 * m1, M2 * m2, double fact, int NC, BaseMatrix::MatrixCategory category)
    {
        switch(NC)
        {
        case 1: if (NL==1) opDynamicReal1<Real, M1, M2 >(m1, m2, fact, category);
            else opDynamicRealNLNC<Real, NL, 1, M1, M2 >(m1, m2, fact, category);
            break;
        case 2: opDynamicRealNLNC<Real, NL, 2, M1, M2 >(m1, m2, fact, category); break;
        case 3: opDynamicRealNLNC<Real, NL, 3, M1, M2 >(m1, m2, fact, category); break;
        case 4: opDynamicRealNLNC<Real, NL, 4, M1, M2 >(m1, m2, fact, category); break;
        case 5: opDynamicRealNLNC<Real, NL, 5, M1, M2 >(m1, m2, fact, category); break;
        case 6: opDynamicRealNLNC<Real, NL, 6, M1, M2 >(m1, m2, fact, category); break;
        default: opDynamicRealDefault<Real, M1, M2 >(m1, m2, fact, NL, NC, category); break;
        }
    }

    template <class Real, class M1, class M2 >
    static inline void opDynamicReal(const M1 * m1, M2 * m2, double fact, int NL, int NC, BaseMatrix::MatrixCategory category)
    {
        switch(NL)
        {
            //case 1: opDynamicRealNL<Real, 1, M1, M2 >(m1, m2, fact, NC, category); break;
        case 2: opDynamicRealNL<Real, 2, M1, M2 >(m1, m2, fact, NC, category); break;
        case 3: opDynamicRealNL<Real, 3, M1, M2 >(m1, m2, fact, NC, category); break;
        case 4: opDynamicRealNL<Real, 4, M1, M2 >(m1, m2, fact, NC, category); break;
        case 5: opDynamicRealNL<Real, 5, M1, M2 >(m1, m2, fact, NC, category); break;
        case 6: opDynamicRealNL<Real, 6, M1, M2 >(m1, m2, fact, NC, category); break;
        default: opDynamicRealDefault<Real, M1, M2 >(m1, m2, fact, NL, NC, category); break;
        }
    }

    template <class M1, class M2 >
    static inline void opDynamic(const M1 * m1, M2 * m2, double fact)
    {
        const int NL = m1->getBlockRows();
        const int NC = m1->getBlockCols();
        const BaseMatrix::MatrixCategory category = m1->getCategory();
        const BaseMatrix::ElementType elementType = m1->getElementType();
        const unsigned int elementSize = m1->getElementSize();

        if (category == BaseMatrix::MATRIX_IDENTITY)
        {
            opIdentity(m1, m2, fact);
        }
        else if (category == BaseMatrix::MATRIX_FULL)
        {
            opFull(m1, m2, fact);
        }
        //else if (category == BaseMatrix::MATRIX_DIAGONAL && NL == 1 && NC == 1)
        //{
        //    opDiagonal(mat, result, v);
        //}
        //else if (elementType == BaseMatrix::ELEMENT_INT && elementSize == sizeof(int))
        //{
        //    opDynamicReal<int, M, V1, V2>(mat, result, v, NL, NC, category);
        //}
        else if (elementType == BaseMatrix::ELEMENT_FLOAT && elementSize == sizeof(float))
        {
            opDynamicReal<float, M1, M2 >(m1, m2, fact, NL, NC, category);
        }
        else // default to double
        {
            opDynamicReal<double, M1, M2 >(m1, m2, fact, NL, NC, category);
        }
    }
};

class BaseMatrixLinearOpAddM : public BaseMatrixLinearOpAM< false> {};
class BaseMatrixLinearOpAddMT : public BaseMatrixLinearOpAM< true> {};

/// Multiply the matrix by vector v and put the result in vector result
void BaseMatrix::opAddM(defaulttype::BaseMatrix* result,double fact) const
{
    BaseMatrixLinearOpAddM::opDynamic(this, result, fact);
}

/// Multiply the matrix by vector v and put the result in vector result
void BaseMatrix::opAddMT(defaulttype::BaseMatrix* result,double fact) const
{
    BaseMatrixLinearOpAddMT::opDynamic(this, result, fact);
}


} // nampespace defaulttype

} // nampespace sofa
