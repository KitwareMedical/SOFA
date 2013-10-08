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
#ifndef SOFA_CORE_COLLISION_GENERICCONTACTCORRECTION_H
#define SOFA_CORE_COLLISION_GENERICCONTACTCORRECTION_H

#include <sofa/core/behavior/ConstraintCorrection.h>

#include <sofa/core/behavior/OdeSolver.h>
#include <sofa/core/behavior/LinearSolver.h>

#include <sofa/defaulttype/Mat.h>
#include <sofa/defaulttype/Vec.h>

#include <sofa/component/linearsolver/SparseMatrix.h>
#include <sofa/component/linearsolver/FullMatrix.h>

namespace sofa {

namespace component {

namespace constraintset {

using namespace sofa::core;
using namespace sofa::core::behavior;
using namespace sofa::defaulttype;

class GenericConstraintCorrection : public BaseConstraintCorrection {
public:
    SOFA_CLASS(GenericConstraintCorrection, BaseConstraintCorrection);

protected:
    GenericConstraintCorrection();
    virtual ~GenericConstraintCorrection();

public:
    virtual void bwdInit();

    virtual void addComplianceInConstraintSpace(const ConstraintParams *cparams, defaulttype::BaseMatrix* W);

    virtual void getComplianceMatrix(defaulttype::BaseMatrix* ) const;

    virtual void computeAndApplyMotionCorrection(const ConstraintParams *cparams, MultiVecCoordId x, MultiVecDerivId v, MultiVecDerivId f, const BaseVector * lambda);

    virtual void computeAndApplyPositionCorrection(const ConstraintParams *cparams, MultiVecCoordId x, MultiVecDerivId f, const BaseVector *lambda);

    virtual void computeAndApplyVelocityCorrection(const ConstraintParams *cparams, MultiVecDerivId v, MultiVecDerivId f, const BaseVector *lambda);

    virtual void applyPredictiveConstraintForce(const core::ConstraintParams * /*cparams*/, core::MultiVecDerivId /*f*/, const defaulttype::BaseVector *lambda);

    virtual void rebuildSystem(double massFactor, double forceFactor);

    virtual void applyContactForce(const defaulttype::BaseVector *f);

    virtual void resetContactForce();

    virtual void computeResidual(const core::ExecParams* /*params*/ /* PARAMS FIRST */, BaseVector *lambda);

    Data< helper::vector< std::string > >  solverName;

    /// Pre-construction check method called by ObjectFactory.
    template<class T>
    static bool canCreate(T*& obj, objectmodel::BaseContext* context, objectmodel::BaseObjectDescription* arg) {
        return BaseConstraintCorrection::canCreate(obj, context, arg);
    }

protected:

    behavior::OdeSolver* odesolver;
    std::vector<sofa::core::behavior::LinearSolver*> linearsolvers;
};

} // namespace collision

} // namespace component

} // namespace sofa

#endif
