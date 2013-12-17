#ifndef PROJECTIONMAPPING_H
#define PROJECTIONMAPPING_H

#include "AssembledMapping.h"
#include "initCompliant.h"

#include "../utils/map.h"
#include "../utils/pair.h"

namespace sofa {

namespace component {

namespace mapping {

/*
  
  A simple projection-mapping: projects input orthogonally (canonical
  dot-product) onto given (input dof index, deriv) pairs. This is a kind of
  general-purpose mapping. Input: any vector dof. Output: 1d dofs.

  @author: maxime.tournier@inria.fr

*/


template <class TIn, class TOut >
class ProjectionMapping : public AssembledMapping<TIn, TOut> {
public:
	SOFA_CLASS(SOFA_TEMPLATE2(ProjectionMapping,TIn,TOut), SOFA_TEMPLATE2(AssembledMapping,TIn,TOut));
	
	typedef typename TIn::Real in_real;
	typedef typename TOut::Real out_real;

	typedef defaulttype::SerializablePair<unsigned, typename TIn::Coord> set_type;
    Data< vector< set_type > > set;
	
	typedef AssembledMapping<TIn, TOut> base;
	typedef ProjectionMapping self;

	ProjectionMapping()
		: set(initData(&set, "set", 
						"(index, coord) vector describing projection space"))  {
		assert( base::Nout == 1 );
		assert( TIn::deriv_total_size == TIn::coord_total_size && "not vector input dofs !" );
	}
	
protected:
	
	virtual void assemble( const typename self::in_pos_type& in) {
		
		const vector<set_type>& s = set.getValue();
		
		// resize/clear jacobian
		typename self::jacobian_type::CompressedMatrix& J = this->jacobian.compressedMatrix;
		J.resize( base::Nout * s.size(), 
				  base::Nin * in.size() );
		
		J.setZero();

		for( unsigned i = 0, n = s.size(); i < n; ++i) {
			unsigned row = i;

			J.startVec( row );
			
			for( unsigned j = 0, m = self::Nin; j < m; ++j) {
				unsigned col = self::Nin * s[i].pair.first + j;

				J.insertBack(row, col) = s[i].pair.second[j];
			}
	
		}
		
		J.finalize();
		
	}
	
	virtual void apply(typename self::out_pos_type& out, 
					   const typename self::in_pos_type& in ) {
		const vector<set_type>& s = set.getValue();
		assert(s.size() == out.size());	

		for( unsigned i = 0, n = s.size(); i < n; ++i) {
			map(out[i])(0) = map(in[s[i].pair.first]).dot( map(s[i].pair.second ) );
		}
		
	}
	
};

}
}
}



#endif
