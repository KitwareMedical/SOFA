#include "MinresSolver.h"

#include <sofa/core/ObjectFactory.h>

#include "utils/scoped.h"
#include "utils/minres.h"

#include "utils/kkt.h"
#include "utils/schur.h"

namespace sofa {
namespace component {
namespace linearsolver {
 
SOFA_DECL_CLASS(MinresSolver);
int MinresSolverClass = core::RegisterObject("Sparse Minres linear solver").add< MinresSolver >();


MinresSolver::MinresSolver() 
	: // use_schur(initData(&use_schur, false, "use_schur", "use Schur complement when solving. *warning* this will invert the response matrix at each time step unless fast_schur is true")),
	  // fast_schur(initData(&fast_schur, false, "fast_schur", "only invert response matrix once")),
	parallel(initData(&parallel, false, "parallel", "use openmp to parallelize matrix-vector products when use_schur is false"))
{
	
}

			

void MinresSolver::solve_schur(AssembledSystem::vec& x,
                               const AssembledSystem& sys,
                               const AssembledSystem::vec& b) const {
	// unconstrained velocity
	vec tmp(sys.m);
	response->solve(tmp, b.head(sys.m));
	x.head( sys.m ) = tmp;
	
	if( sys.n ) {
		
		::schur<response_type> A(sys, *response);
		
		vec rhs = b.tail(sys.n) - sys.J * x.head(sys.m);
		
		vec lambda = x.tail(sys.n);

		typedef ::minres<SReal> solver_type;		
		
		solver_type::params p = params(rhs);
		solver_type::solve(lambda, A, rhs, p);
		
		// constraint velocity correction
		response->solve(tmp, sys.J.transpose() * lambda );

		x.head( sys.m ) += tmp;
		x.tail( sys.n ) = lambda;
		
		report("minres (schur)", p );
	}


}


void MinresSolver::solve_kkt(AssembledSystem::vec& x,
                             const AssembledSystem& system,
                             const AssembledSystem::vec& b) const {
	params_type p = params(b);
			
	vec rhs = b;
	if( system.n ) rhs.tail(system.n) = -rhs.tail(system.n);
	
	kkt A(system, parallel.getValue() );
	
	typedef minres<real> solver_type;
	solver_type::solve(x, A, rhs, p);
	
	report("minres (kkt)", p );
}



			
}
}
}


