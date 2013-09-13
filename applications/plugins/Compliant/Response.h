#ifndef COMPLIANT_RESPONSE_H
#define COMPLIANT_RESPONSE_H

#include "initCompliant.h"

#include "AssembledSystem.h"
#include <sofa/core/objectmodel/BaseObject.h>

namespace sofa {
namespace component {
namespace linearsolver {

struct Response : core::objectmodel::BaseObject {

	SOFA_CLASS(Response, core::objectmodel::BaseObject);

	typedef AssembledSystem system_type;
	
	typedef system_type::real real;
	typedef system_type::vec vec;

	typedef system_type::mat mat;
	typedef system_type::cmat cmat;
	
	virtual void factor(const mat& ) = 0;
	
	virtual void solve(cmat&, const cmat& ) const = 0;
	virtual void solve(vec&, const vec& ) const = 0;
	
};


}
}
}

#endif
