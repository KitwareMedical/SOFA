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
#ifndef SOFA_COMPONENT_LINEARSOLVER_EigenSparseMatrix_H
#define SOFA_COMPONENT_LINEARSOLVER_EigenSparseMatrix_H

#include "EigenBaseSparseMatrix.h"
#include <sofa/defaulttype/Mat.h>
#include <sofa/component/linearsolver/CompressedRowSparseMatrix.h>
#include <sofa/helper/SortedPermutation.h>
#include <sofa/helper/vector.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
using std::cerr;
using std::endl;

namespace sofa
{

namespace component
{

namespace linearsolver
{
using helper::vector;


/** Variant of EigenBaseSparseMatrix, capable of block-view access.
  The blocks correspond to matrix blocks of the size of the DataTypes Deriv.

  There are two ways of filling the matrix:
  - Random block access is provided by method wBlock. Use compress() after the last insertion.
  - Block rows can be efficiently appended using methods beginBlockRow, createBlock, endBlockRow. Use compress() after the last insertion. The rows must be created in increasing index order.

  The two ways of filling the matrix can not be used at the same time.
  */
template<class InDataTypes, class OutDataTypes>
class EigenSparseMatrix : public EigenBaseSparseMatrix<typename OutDataTypes::Real>
{
public:
    typedef EigenBaseSparseMatrix<typename OutDataTypes::Real> Inherit;
    typedef typename OutDataTypes::Real OutReal;

    typedef typename OutDataTypes::Real Real; // what's inside the eigen matrix

    typedef Eigen::SparseMatrix<OutReal,Eigen::RowMajor> CompressedMatrix;
    typedef Eigen::Matrix<OutReal,Eigen::Dynamic,1>  VectorEigenOut;

    typedef typename InDataTypes::Deriv InDeriv;
    typedef typename InDataTypes::VecDeriv InVecDeriv;
    typedef typename InDataTypes::Real InReal;
    typedef typename OutDataTypes::Deriv OutDeriv;
    typedef typename OutDataTypes::VecDeriv OutVecDeriv;
    enum { Nin=InDataTypes::deriv_total_size, Nout=OutDataTypes::deriv_total_size };
    typedef defaulttype::Mat<Nout,Nin, OutReal> Block;  ///< block relating an OutDeriv to an InDeriv. This is used for input only, not for internal storage.

protected:
    typedef std::map<int,Block> BlockRowMap;        ///< Map which represents one block-view row of the matrix. The index represents the block-view column index of an entry.
    typedef std::map<int,BlockRowMap> BlockMatMap;  ///< Map which represents a block-view matrix. The index represents the block-view index of a block-view row.
    BlockMatMap incomingBlocks;                     ///< To store block-view data before it is compressed in optimized format.
    typedef Eigen::Matrix<InReal,Eigen::Dynamic,1>  VectorEigenIn;


  	// some helpers to ease mapping from/to eigen types, lvalue and rvalue
	  template<class VecDeriv>
	  struct map_traits {
		  typedef typename VecDeriv::value_type value_type;
		  typedef typename defaulttype::DataTypeInfo<value_type>::ValueType real_type;
		  
		  static const unsigned size = defaulttype::DataTypeInfo<value_type>::Size;
		  
		  typedef Eigen::Matrix< real_type, Eigen::Dynamic, 1 > matrix_type;

		  typedef Eigen::Map< matrix_type > map_type;
		  typedef Eigen::Map< const matrix_type > const_map_type;
		  
          static map_type map(real_type* data, unsigned k) {
              return map_type( data, k * size );
		  }
		  
		  static const_map_type const_map(const real_type* data, unsigned k) {
			  return const_map_type( data, k * size );
		  }
		  

	  };
	  
      template<class VecDeriv>
      static typename map_traits<VecDeriv>::const_map_type
      map(const helper::ReadAccessor< Data<VecDeriv> >& data) {
          return map_traits<VecDeriv>::const_map(&data[0][0], data.size());
      }
	
      template<class VecDeriv>
      static typename map_traits<VecDeriv>::map_type
      map(helper::WriteAccessor< Data<VecDeriv> >& data) {
          return map_traits<VecDeriv>::map(&data[0][0], data.size());
      }

	  template<class VecDeriv>
	  static typename map_traits<VecDeriv>::const_map_type
	  map(const VecDeriv& data) {
		  return map_traits<VecDeriv>::const_map(&data[0][0], data.size());
	  }

	  template<class VecDeriv>
	  static typename map_traits<VecDeriv>::map_type
	  map( VecDeriv& data) {
		   return map_traits<VecDeriv>::map(&data[0][0], data.size());
	  }
	

	template<class LHS, class RHS>
	static bool alias(const LHS& lhs, const RHS& rhs) {
		return (void*)(&lhs[0][0]) == (void*)(&rhs[0][0]);
	}

public:

    EigenSparseMatrix(int nbRow=0, int nbCol=0):Inherit(nbRow,nbCol) {}

    /// Resize the matrix without preserving the data (the matrix is set to zero), with the size given in number of blocks
    void resizeBlocks(int nbBlockRows, int nbBlockCols)
    {
//        this->compressedMatrix.resize(nbBlockRows * Nout, nbBlockCols * Nin);
        this->resize(nbBlockRows * Nout, nbBlockCols * Nin);
    }


//    /// Finalize the matrix after a series of insertions. Add the values from the temporary list to the compressed matrix, and clears the list.
//    virtual void compress()
//    {
//        Inherit::compress();

//        if( incomingBlocks.empty() ) return;
//        compress_incomingBlocks();
//        //        cerr<<"compress, before incoming blocks " << this->eigenMatrix << endl;
//        //        cerr<<"compress, incoming blocks " << this->compressedIncoming << endl;
//        this->compressedMatrix += this->compressedIncoming;
//        //        cerr<<"compress, final value " << this->eigenMatrix << endl;
//        this->compressedMatrix.finalize();
//    }

//    /** Return write access to an incoming block.
//    Note that this does not give access to the compressed matrix.
//    The block belongs to a temporary list which will be added to the compressed matrix using method compress().
//    */
//    Block& wBlock( int i, int j )
//    {
//        return incomingBlocks[i][j];
//    }


    /** Prepare the insertion of a new row of blocks in the matrix.
       Then create blocks using createBlock( unsigned column,  const Block& b ).
        Then finally use endBlockRow() to validate the row insertion.
        @sa createBlock( unsigned column,  const Block& b )
        @sa endBlockRow()
        */
    void beginBlockRow(unsigned row)
    {
        bRow = row;
        bColumns.clear();
        blocks.clear();
    }

    /** Create a block in the current row, which must be previously
        initialized using beginBlockRow(unsigned row).  The blocks
        need not be created in column order. The blocks are not
        actually created in the matrix until method endBlockRow() is
        called.  If the block already exists, they will be added.
        */
    void createBlock( unsigned column,  const Block& b )
    {
        blocks.push_back(b);
        bColumns.push_back(column);
    }

    /** Finalize the creation of the current block row.
      @sa beginBlockRow(unsigned row)
      @sa createBlock( unsigned column,  const Block& b )
      */
    void endBlockRow()
    {
        vector<unsigned> p = helper::sortedPermutation(bColumns); // indices in ascending column order

        for( unsigned r=0; r<Nout; r++ )   // process one scalar row after another
        {
//           this->beginRow(r+ bRow*Nout);
            for(unsigned i=0; i<p.size(); i++ )  // process the blocks in ascending order
            {
                const Block& b = blocks[p[i]];
                for( unsigned c=0; c<Nin; c++ )
                {
                    if( b[r][c]!=0.0 ){
//                        this->insertBack( r + bRow*Nout, c + bColumns[p[i]] * Nin, b[r][c]);
                        this->add( r + bRow*Nout, c + bColumns[p[i]] * Nin, b[r][c]);
                    }
                }
            }
        }
    }


    // max: added template real type to work around the mess between
    // mapping ::Real member types (sometimes it's In::Real, sometimes it's
    // Out::Real, see e.g. BarycentricMapping/RigidMapping)

    /** Set from a CompressedRowSparseMatrix. @pre crs must be compressed
      */
    template<class AnyReal>
    void copyFrom( const CompressedRowSparseMatrix< defaulttype::Mat<Nout,Nin, AnyReal> >& crs )
    {
        this->resize( crs.rowSize(), crs.colSize() );
//        cerr<<"copyFrom, size " << crs.rowSize() << ", " << crs.colSize()<< ", block rows: " << crs.rowIndex.size() << endl;
//        cerr<<"copyFrom, crs = " << crs << endl;

//        int rowStarted = 0;
        for (unsigned int xi = 0; xi < crs.rowIndex.size(); ++xi)  // for each non-null block row
        {
            int blRow = crs.rowIndex[xi];      // block row

//            while( rowStarted<blRow*Nout )   // make sure all the rows are started, even the empty ones
//            {
//                this->compressedMatrix.startVec(rowStarted);
//                rowStarted++;
//            }

            typename CompressedRowSparseMatrix<Block>::Range rowRange(crs.rowBegin[xi], crs.rowBegin[xi+1]);

            for( unsigned r=0; r<Nout; r++ )   // process one scalar row after another
            {
                if(r+ blRow*Nout >= this->rowSize() ) break;
//                cerr<<"copyFrom,  startVec " << rowStarted << endl;
//                this->compressedMatrix.startVec(rowStarted++);


                for (int xj = rowRange.begin(); xj < rowRange.end(); ++xj)  // for each non-null block
                {
                    int blCol = crs.colsIndex[xj];     // block column
                    const Block& b = crs.colsValue[xj]; // block value
                    for( unsigned c=0; c<Nin; c++ ) if( c+ blCol*Nin < this->colSize() )
                        {
                        this->add(r + blRow*Nout, c + blCol*Nin, b[r][c]);
//                        this->compressedMatrix.insertBack(r + blRow*Nout, c + blCol*Nin) = b[r][c];
//                        cerr<<"copyFrom,  insert at " << r + blRow*Nout << ", " << c + blCol*Nin << endl;
                        }

                }
            }
        }
        this->compress();

    }



protected:

	// max: factored out the two exact same implementations for
	// VecDeriv/Data<VecDeriv>
	
	template<class OutType, class InType>
	void mult_impl(OutType& result, const InType& data) const {
		 
		// use optimized product if possible
		if(canCast(data)) {

            if( alias(result, data) ) {
				this->map(result) = (this->compressedMatrix * 
				                     this->map(data).template cast<Real>()).template cast<OutReal>();
            } else {
                this->map(result).noalias() = (this->compressedMatrix *
                                               this->map(data).template cast<Real>()).template cast<OutReal>();
            }
			
			return;
		}

		// convert the data to Eigen type
		VectorEigenOut aux1(this->colSize(),1), aux2(this->rowSize(),1);
		for(unsigned i = 0, n = data.size(); i < n; ++i) {
			for(unsigned j = 0; j < Nin; ++j) {
				aux1[Nin * i + j] = data[i][j];
			}
		}
		
		// compute the product
		aux2.noalias() = this->compressedMatrix * aux1;
		// convert the result back to the Sofa type
		for(unsigned i = 0, n = result.size(); i < n; ++i) {
			for(unsigned j = 0; j < Nout; ++j) {
				result[i][j] = aux2[Nout * i + j];
			}
		}
	}


	template<class OutType, class InType>
	void addMult_impl( OutType& result, const InType& data, Real fact) const {
		
		// use optimized product if possible
		if( canCast(data) ) {

			// TODO multiply only the smallest dimension by fact 

            if( alias(result, data) ) {
                map(result) += (this->compressedMatrix * (map(data).template cast<Real>() * fact)).template cast<OutReal>();
            } else {
                typename map_traits<OutType>::map_type r = map(result);
                r.noalias() = r + (this->compressedMatrix * (map(data).template cast<Real>() * fact)).template cast<OutReal>();
            }
			
			return;
		}

		// convert the data to Eigen type
		VectorEigenOut aux1(this->colSize()),aux2(this->rowSize());
		for(unsigned i = 0, n = data.size(); i < n; ++i) {
			for(unsigned j = 0; j < Nin; ++j) {
				aux1[Nin * i + j] = data[i][j];
			}
		}
        
		// compute the product
		aux2.noalias() = this->compressedMatrix * aux1;
        
		// convert the result back to the Sofa type
		for(unsigned i = 0, n = result.size(); i < n; ++i) {
			for(unsigned j = 0; j < Nout; ++j) {
				result[i][j] += aux2[Nout * i + j] * fact;
			}
		}
	}

	template<class InType, class OutType>
    void addMultTranspose_impl( InType& result, const OutType& data, Real fact) const {

		// use optimized product if possible
		if(canCast(result)) {

            // TODO multiply only the smallest dimension by fact

            if( alias(result, data) ) {
                map(result) += (this->compressedMatrix.transpose() * (map(data).template cast<Real>() * fact)).template cast<InReal>();
            } else {
                typename map_traits<InType>::map_type r = map(result);
                r.noalias() = r + (this->compressedMatrix.transpose() * (map(data).template cast<Real>() * fact)).template cast<InReal>();
            }
			
			return;
		}

		// convert the data to Eigen type
		VectorEigenOut aux1(this->rowSize()), aux2(this->colSize());

		for(unsigned i = 0, n = data.size(); i < n; ++i) {
			for(unsigned j = 0; j < Nout; ++j) {
				aux1[Nout * i + j] = data[i][j];
			}
		}
		
		// compute the product
		aux2.noalias() = this->compressedMatrix.transpose() * aux1;
		
		// convert the result back to the Sofa type
		for(unsigned i = 0, n = result.size(); i < n; ++i) {
			for(unsigned j = 0; j < Nin; ++j) {
                result[i][j] += aux2[Nin * i + j] * fact;
			}
		}
	}

	
public:

    /// compute result = A * data
    void mult( OutVecDeriv& result, const InVecDeriv& data ) const {
	    mult_impl(result, data);
    }


    /// compute result = A * data
    void mult( Data<OutVecDeriv>& _result, const Data<InVecDeriv>& _data ) const {
        helper::WriteAccessor<Data<OutVecDeriv> > result (_result);
        helper::ReadAccessor<Data<InVecDeriv> > data (_data);

        mult_impl(result, data);
    }

    /// compute result += A * data
    void addMult( OutVecDeriv& result, const InVecDeriv& data ) const {
	    addMult_impl(result, data, 1.0);
    }
      
    /// compute result += A * data * fact
    void addMult( OutVecDeriv& result, const InVecDeriv& data, const OutReal fact ) const {
        addMult_impl(result, data, fact);
    }

    /// compute result += A * data
    void addMult( Data<OutVecDeriv>& result, const Data<InVecDeriv>& data) const {
        helper::WriteAccessor<Data<OutVecDeriv> > res (result);
        helper::ReadAccessor<Data<InVecDeriv> > dat (data);

        addMult_impl(res, dat, 1.0);
    }
	
   
    /// compute result += A * data * fact
    void addMult( Data<OutVecDeriv>& result, const Data<InVecDeriv>& data, const OutReal fact) const {
        helper::WriteAccessor<Data<OutVecDeriv> > res (result);
        helper::ReadAccessor<Data<InVecDeriv> > dat (data);

        addMult_impl(res, dat, fact);
    }
    
	

    /// compute result += A^T * data
    void addMultTranspose( InVecDeriv& result, const OutVecDeriv& data ) const {
        addMultTranspose_impl(result, data, 1.0);
    }

    /// compute result += A^T * data * fact
    void addMultTranspose( InVecDeriv& result, const OutVecDeriv& data, const OutReal fact ) const {
        addMultTranspose_impl(result, data, fact);
    }

    /// compute result += A^T * data
    void addMultTranspose( Data<InVecDeriv>& result, const Data<OutVecDeriv>& data ) const {
        helper::WriteAccessor<Data<InVecDeriv> > res (result);
        helper::ReadAccessor<Data<OutVecDeriv> > dat (data);

        addMultTranspose_impl(res, dat, 1.0);
    }

    /// compute result += A^T * data * fact
    void addMultTranspose( Data<InVecDeriv>& result, const Data<OutVecDeriv>& data, const OutReal fact ) const {
        helper::WriteAccessor<Data<InVecDeriv> > res (result);
        helper::ReadAccessor<Data<OutVecDeriv> > dat (data);

        addMultTranspose_impl(res, dat, fact);
    }

    static const char* Name();

private:
    //@{
    /** Auxiliary variables for methods beginBlockRow(unsigned row),
     * createBlock( unsigned column, const Block& b ) and
     * endBlockRow() */

    unsigned bRow;
    vector<unsigned> bColumns;
    vector<Block> blocks;
    //@}

	
	template<class T>
	bool canCast( const T& v ) const { 
		// is it contiguous ?
		typedef typename T::value_type value_type;
		typedef typename value_type::value_type scalar_type;
		return  (v.size() - 1) * sizeof(value_type) == (&v[v.size() - 1][0] - &v[0][0]) * sizeof(scalar_type); 
	}
	
	


};

#ifndef SOFA_FLOAT
template<> inline const char* EigenSparseMatrix<defaulttype::Vec3dTypes, defaulttype::Vec1dTypes >::Name() { return "EigenSparseMatrix3d1d"; }
#endif

#ifndef SOFA_DOUBLE
template<> inline const char* EigenSparseMatrix<defaulttype::Vec3fTypes, defaulttype::Vec1fTypes >::Name() { return "EigenSparseMatrix3f1f"; }
#endif

// max: much cleaner like this :)


} // namespace linearsolver

} // namespace component

} // namespace sofa

#endif
