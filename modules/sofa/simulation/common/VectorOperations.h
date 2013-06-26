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
#ifndef SOFA_SIMULATION_COMMON_VECTOROPERATIONS_H
#define SOFA_SIMULATION_COMMON_VECTOROPERATIONS_H

#include <sofa/core/behavior/BaseVectorOperations.h>
#include <sofa/simulation/common/common.h>
#include <sofa/simulation/common/VisitorExecuteFunc.h>

namespace sofa
{

namespace core
{
class ExecParams;
}

namespace simulation
{

namespace common
{

class SOFA_SIMULATION_COMMON_API VectorOperations : public sofa::core::behavior::BaseVectorOperations
{
public:

    VectorOperations(const sofa::core::ExecParams* params, sofa::core::objectmodel::BaseContext* ctx);

    /// Allocate a temporary vector
    void v_alloc(sofa::core::MultiVecCoordId& v);
    void v_alloc(sofa::core::MultiVecDerivId& v);
    /// Free a previously allocated temporary vector
    void v_free(sofa::core::MultiVecCoordId& id);
    void v_free(sofa::core::MultiVecDerivId& id, bool interactionForceField=false);

    void v_realloc(sofa::core::MultiVecCoordId& id, bool interactionForceField=false);
    void v_realloc(sofa::core::MultiVecDerivId& id, bool interactionForceField=false);

    void v_clear(core::MultiVecId v); ///< v=0
    void v_eq(core::MultiVecId v, core::MultiVecId a); ///< v=a
    void v_eq(core::MultiVecId v, core::MultiVecId a, double f); ///< v=f*a
    void v_peq(core::MultiVecId v, core::MultiVecId a, double f=1.0); ///< v+=f*a
#ifdef SOFA_SMP
    void v_peq(core::MultiVecId v, core::MultiVecId a, Shared<double> &fSh, double f=1.0) ; ///< v+=f*a
    void v_meq(core::MultiVecId v, core::MultiVecId a, Shared<double> &fSh) ; ///< v+=f*a
#endif
    void v_teq(core::MultiVecId v, double f) ; ///< v*=f
    void v_op(core::MultiVecId v, core::ConstMultiVecId a, core::ConstMultiVecId  b, double f=1.0) ; ///< v=a+b*f
#ifdef SOFA_SMP
    void v_op(core::MultiVecId v, core::MultiVecId a, core::MultiVecId b, Shared<double> &f) ; ///< v=a+b*f
#endif
    void v_multiop(const core::behavior::BaseMechanicalState::VMultiOp& o);
    void v_dot(core::ConstMultiVecId a, core::ConstMultiVecId  b); ///< a dot b ( get result using finish )
#ifdef SOFA_SMP
    void v_dot(Shared<double> &result,core::MultiVecId a, core::MultiVecId b) ; ///< a dot b
#endif
    void v_threshold(core::MultiVecId a, double threshold); ///< nullify the values below the given threshold

    double finish();
    void print( core::MultiVecId v, std::ostream& out );
protected:
    VisitorExecuteFunc executeVisitor;
    /// Result of latest v_dot operation
    double result;
};

}
}
}


#endif //SOFA_SIMULATION_COMMON_VECTOROPERATIONS_H
