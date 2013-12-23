/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU General Public License as published by the Free  *
* Software Foundation; either version 2 of the License, or (at your option)   *
* any later version.                                                          *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    *
* more details.                                                               *
*                                                                             *
* You should have received a copy of the GNU General Public License along     *
* with this program; if not, write to the Free Software Foundation, Inc., 51  *
* Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.                   *
*******************************************************************************
*                            SOFA :: Applications                             *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
/** Sparse matrix test suite.
  The same suite is instanciated using different parameters: entry types (float/double) and BlockMN size in CompressedRowSparse.
  */


#include <gtest/gtest.h>
#include "Sofa_test.h"

#include <sofa/component/linearsolver/EigenSparseMatrix.h>
#include <sofa/component/linearsolver/SparseMatrix.h>
#include <sofa/component/linearsolver/CompressedRowSparseMatrix.h>
#include <sofa/component/linearsolver/FullMatrix.h>
#include <sofa/component/linearsolver/FullVector.h>
#include <sofa/defaulttype/Mat.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/helper/RandomGenerator.h>

#include <ctime>

#define BENCHMARK_MATRIX_MULTIPLICATION 1

#if BENCHMARK_MATRIX_MULTIPLICATION
#include <sys/time.h>
double get_time() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (double) tv.tv_sec + tv.tv_usec*1e-6;
}
#endif


namespace sofa {
using std::cout;
using std::cerr;
using std::endl;



/** Sparse matrix test suite.

  Creates matrices of different types, sets their entries and checks that all the matrices are the same.

  Perform matrix-vector products and compare the results.
  */
template <typename _Real, unsigned NumRows, unsigned NumCols, unsigned BlockRows, unsigned BlockCols>
struct TestSparseMatrices : public Sofa_test<_Real>
{
    // Scalar type and dimensions of the matrices to test
    typedef _Real Real;
    static const unsigned NROWS=NumRows;   // matrix size
    static const unsigned NCOLS=NumCols;
    static const unsigned BROWS=BlockRows; // block size used for matrices with block-wise storage
    static const unsigned BCOLS=BlockCols;


    // Dense implementation
    typedef sofa::component::linearsolver::FullMatrix<Real> FullMatrix;
    typedef sofa::component::linearsolver::FullVector<Real> FullVector;

    // Simple sparse matrix implemented using map< map< > >
    typedef sofa::component::linearsolver::SparseMatrix<Real> MapMatrix;

    // Blockwise Compressed Sparse Row format
    typedef sofa::defaulttype::Mat<BROWS,BCOLS,Real> BlockMN;
    typedef sofa::component::linearsolver::CompressedRowSparseMatrix<BlockMN> CRSMatrixMN;
    typedef sofa::defaulttype::Mat<BCOLS,BROWS,Real> BlockNM;
    typedef sofa::component::linearsolver::CompressedRowSparseMatrix<BlockNM> CRSMatrixNM;
    typedef sofa::defaulttype::Mat<BROWS,BROWS,Real> BlockMM;
    typedef sofa::component::linearsolver::CompressedRowSparseMatrix<BlockMM> CRSMatrixMM;
    typedef sofa::defaulttype::Mat<BCOLS,BCOLS,Real> BlockNN;
    typedef sofa::component::linearsolver::CompressedRowSparseMatrix<BlockNN> CRSMatrixNN;

    // Implementation based on Eigen
    typedef sofa::defaulttype::StdVectorTypes< sofa::defaulttype::Vec<BCOLS,Real>, sofa::defaulttype::Vec<BCOLS,Real> > InTypes;
    typedef sofa::defaulttype::StdVectorTypes< sofa::defaulttype::Vec<BROWS,Real>, sofa::defaulttype::Vec<BROWS,Real> > OutTypes;
    typedef sofa::component::linearsolver::EigenSparseMatrix<InTypes,OutTypes> EigenBlockMatrix;
    typedef sofa::component::linearsolver::EigenBaseSparseMatrix<Real> EigenBaseMatrix;


    // Regular Mat implementation
    typedef sofa::defaulttype::Mat<NROWS,NCOLS,Real> MatMN;
    typedef sofa::defaulttype::Mat<NCOLS,NROWS,Real> MatNM;
    typedef sofa::defaulttype::Mat<NROWS,NROWS,Real> MatMM;
    typedef sofa::defaulttype::Mat<NCOLS,NCOLS,Real> MatNN;
    typedef sofa::defaulttype::Vec<NROWS,Real> VecM;
    typedef sofa::defaulttype::Vec<NCOLS,Real> VecN;


    // The matrices used in the tests
    CRSMatrixMN crs1,crs2;
    FullMatrix fullMat;
    MapMatrix mapMat;
    EigenBlockMatrix eiBlock1,eiBlock2,eiBlock3;
    EigenBaseMatrix eiBase;
    // matrices for multiplication test
    CRSMatrixNM crsMultiplier;
    CRSMatrixMM crsMultiplication;
    CRSMatrixNN crsTransposeMultiplication;
    FullMatrix fullMultiplier;
    FullMatrix fullMultiplication;
    FullMatrix fullTransposeMultiplication;
    MatMN mat;
    MatNM matMultiplier;
    MatMM matMultiplication;
    MatNN matTransposeMultiplication;
    EigenBaseMatrix eiBaseMultiplier;
    EigenBaseMatrix eiBaseMultiplication;


    // The vectors used in the tests
    FullVector fullVec_ncols;
    FullVector fullVec_nrows_reference,fullVec_nrows_result;
    VecM vecM;
    VecN vecN;


    /// generating a random Mat
    /// if sparse=0 a lot a null entries will be created
    template<int nbrows,int nbcols>
    static void generateRandomMat( defaulttype::Mat<nbrows,nbcols,Real>& mat, bool sparse=false )
    {
        sofa::helper::RandomGenerator randomGenerator;
        randomGenerator.initSeed( (long)time(0) );

        for( int j=0; j<mat.nbCols; j++)
        {
            for( int i=0; i<mat.nbLines; i++)
            {
                Real random = randomGenerator.random<Real>( (Real) -1, (Real) 1 );
                if( sparse && random > -0.5 && random < 0.5 ) random = 0;
                mat(i,j)=random;
            }
        }
    }

    /// filling a BaseMatrix from a given Mat
    template<int nbrows,int nbcols>
    static void copyFromMat( defaulttype::BaseMatrix& baseMat, const defaulttype::Mat<nbrows,nbcols,Real>& mat )
    {
        baseMat.resize(mat.nbLines,mat.nbCols);

        for( int j=0; j<mat.nbCols; j++)
        {
            for( int i=0; i<mat.nbLines; i++)
            {
                if( !baseMat.isSparse() || mat(i,j)!=0 ) baseMat.add( i, j, mat(i,j) );
            }
        }

        baseMat.compress();
    }

    /// Create the context for the matrix tests.
    TestSparseMatrices()
    {
        //std::cout<<"Matrix_test "<<NumRows<<" "<<NumCols<<" "<<BlockRows<<" "<<BlockCols<<std::endl;

        // resize and fill the matrices
        generateRandomMat( mat, true );
        copyFromMat( crs1, mat );
        copyFromMat( crs2, mat );
        copyFromMat( fullMat, mat );
        copyFromMat( mapMat, mat );
        copyFromMat( eiBlock1, mat );
        copyFromMat( eiBlock2, mat );
//        copyFromMat( eiBlock3, mat );
        copyFromMat( eiBase, mat );
        eiBlock3.copyFrom(crs1);




        // resize and fill the vectors
        fullVec_ncols.resize(NCOLS);
        fullVec_nrows_reference.resize(NROWS);
        fullVec_nrows_result.resize(NROWS);
        for( unsigned i=0; i<NCOLS; i++)
        {
            fullVec_ncols[i] = i;
            vecN[i] = i;
        }
        fullMat.mul(fullVec_nrows_reference,fullVec_ncols); //    cerr<<"MatrixTest: vref = " << vref << endl;

        vecM = mat * vecN;




        // matrix multiplication

        generateRandomMat( matMultiplier, true );
        copyFromMat( crsMultiplier, matMultiplier );
        copyFromMat( fullMultiplier, matMultiplier );
        copyFromMat( eiBaseMultiplier, matMultiplier );


        matMultiplication = mat * matMultiplier;
        crs1.mul( crsMultiplication, crsMultiplier );
        fullMat.mul( fullMultiplication, fullMultiplier );
        eiBase.mul_MT( eiBaseMultiplication, eiBaseMultiplier );

        matTransposeMultiplication = mat.multTranspose( mat );
        crs1.mulTranspose( crsTransposeMultiplication, crs1 );
        fullMat.mulT( fullTransposeMultiplication, fullMat );

    }


    /** Check that EigenMatrix update works as well as direct init. Return true if the test succeeds.*/
    bool checkEigenMatrixUpdate()
    {
        // fill two matrices with the same values, one directly, one it two passes, then compare their values
        EigenBlockMatrix a,b;
        a.resize(NROWS,NCOLS);
        b.resize(NROWS,NCOLS);
        for( unsigned j=0; j<NCOLS; j++)
        {
            for( unsigned i=0; i<NROWS; i++)
            {
                double valij = i*NCOLS+j;
                a.add(i,j,valij);
                if( i==j )
                    b.add(i,j,valij);
            }
        }
        a.compress();
        b.compress();

        //    cerr<<"MatrixTest<Real,RN,CN>::checkEigenMatrixUpdate, a = " << a << endl;
        //    cerr<<"MatrixTest<Real,RN,CN>::checkEigenMatrixUpdate, b incomplete = " << b << endl;

        // second pass for b. Some values are set with the right value, some with the wrong value, some are not set
        for( unsigned j=0; j<NCOLS; j++)
        {
            for( unsigned i=0; i<NROWS; i++)
            {
                double valij = i*NCOLS+j;
                if( i!=j )
                    b.add(i,j,valij);
            }
        }
        b.compress();
        //    cerr<<"MatrixTest<Real,RN,CN>::checkEigenMatrixUpdate, b complete = " << b << endl;
        return Sofa_test<_Real>::matricesAreEqual(a,b);
    }

    /** Check the filling of EigenMatrix per rows of blocks. Return true if the test succeeds.*/
    bool checkEigenMatrixBlockRowFilling()
    {
        EigenBlockMatrix mb;
        FullMatrix ma;
        unsigned br=3, bc=3;
        ma.resize(br*BROWS,bc*BCOLS);
        mb.resizeBlocks(br,bc);
        for( unsigned i=0; i<br; i++ )
        {
            mb.beginBlockRow(i);
            if( i%2==0 ) // leave some rows empty
            {
                for( int j=bc-1; j>=0; j--) // set the blocs in reverse order, for fun.
                {
                    // create a block and give it some value
                    BlockMN b;
                    for( unsigned k=0; k<BROWS && k<BCOLS; k++ ){
                        b[k][k] = i+j;
                        ma.set(i*BROWS+k, j*BCOLS+k, i+j);
                    }

                    // insert the block in the matrix
                    mb.createBlock(j,b);
                }
            }
            mb.endBlockRow();
        }
        mb.compress();
        //    serr()<<"MatrixTest<Real,RN,CN>::checkEigenMatrixBlockRowFilling, ma = " << ma << endl;
        //    serr()<<"MatrixTest<Real,RN,CN>::checkEigenMatrixBlockRowFilling, mb = " << mb << endl;
        return Sofa_test<_Real>::matricesAreEqual(ma,mb);
    }

    bool checkEigenMatrixBlockFromCompressedRowSparseMatrix()
    {
//        if( !matricesAreEqual(crs1,eiBlock3)){
//            cout<<"heckEigenMatrixBlockFromCompressedRowSparseMatrix, crs1 = " << crs1 << endl;
//            cout<<"heckEigenMatrixBlockFromCompressedRowSparseMatrix, eiBlock3 = " << eiBlock3 << endl;
//        }
        return Sofa_test<_Real>::matricesAreEqual(crs1,eiBlock3);
    }
};

#ifndef SOFA_FLOAT
///////////////////
// double precision
///////////////////
// trivial blocs
typedef TestSparseMatrices<double,4,8,4,8> Ts4848;
#define TestMatrix Ts4848
#include "Matrix_test.inl"
#undef TestMatrix

// semi-trivial blocs
typedef TestSparseMatrices<double,4,8,4,2> Ts4842;
#define TestMatrix Ts4842
#include "Matrix_test.inl"
#undef TestMatrix

typedef TestSparseMatrices<double,4,8,1,8> Ts4818;
#define TestMatrix Ts4818
#include "Matrix_test.inl"
#undef TestMatrix

// well-fitted blocs
typedef TestSparseMatrices<double,4,8,2,2> Ts4822;
#define TestMatrix Ts4822
#include "Matrix_test.inl"
#undef TestMatrix

/// not fitted blocs
//typedef TestSparseMatrices<double,4,8,2,3> Ts4823;
//#define TestMatrix Ts4823
//#include "Matrix_test.inl"
//#undef TestMatrix


#if BENCHMARK_MATRIX_MULTIPLICATION
///// multiplication timing
typedef TestSparseMatrices<double,360,900,3,3> TsMultiplicationTimings;
TEST_F(TsMultiplicationTimings, benchmark )
{
    matMultiplication.clear();
    double start = get_time();
    matMultiplication = mat * matMultiplier;
    double stop = get_time();
    std::cerr<<"Mat:\t\t"<<stop-start<<" (ms)"<<std::endl;

    fullMultiplication.clear();
    start = get_time();
    fullMat.mul( fullMultiplication, fullMultiplier );
    stop = get_time();
    std::cerr<<"Full:\t\t"<<stop-start<<" (ms)"<<std::endl;

    crsMultiplication.clear();
    start = get_time();
    crs1.mul( crsMultiplication, crsMultiplier );
    stop = get_time();
    std::cerr<<"CRS:\t\t"<<stop-start<<" (ms)"<<std::endl;

    eiBaseMultiplication.clear();
    start = get_time();
    eiBase.mul( eiBaseMultiplication, eiBaseMultiplier );
    stop = get_time();
    std::cerr<<"Eigen Base ST:\t\t"<<stop-start<<" (ms)"<<std::endl;

#ifdef USING_OMP_PRAGMAS
    eiBaseMultiplication.clear();
    start = get_time();
    eiBase.mul_MT( eiBaseMultiplication, eiBaseMultiplier );
    stop = get_time();
    std::cerr<<"Eigen Base MT:\t\t"<<stop-start<<" (ms)"<<std::endl;
#endif

    ASSERT_TRUE( true );
}

#endif

#endif

#ifndef SOFA_DOUBLE
///////////////////
// simple precision
///////////////////
// trivial blocs
typedef TestSparseMatrices<float,4,8,4,8> Ts4848f;
#define TestMatrix Ts4848f
#include "Matrix_test.inl"
#undef TestMatrix

// semi-trivial blocs
typedef TestSparseMatrices<float,4,8,4,2> Ts4842f;
#define TestMatrix Ts4842f
#include "Matrix_test.inl"
#undef TestMatrix

typedef TestSparseMatrices<float,4,8,1,8> Ts4818f;
#define TestMatrix Ts4818f
#include "Matrix_test.inl"
#undef TestMatrix

/// well-fitted blocs
typedef TestSparseMatrices<float,4,8,2,2> Ts4822f;
#define TestMatrix Ts4822f
#include "Matrix_test.inl"
#undef TestMatrix

#endif

}// namespace sofa

