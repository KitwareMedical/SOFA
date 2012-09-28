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
#ifndef SOFA_COMPONENT_CONSTRAINTSET_CONSTRAINTSOLVERIMPL_H
#define SOFA_COMPONENT_CONSTRAINTSET_CONSTRAINTSOLVERIMPL_H

#include <sofa/core/behavior/ConstraintSolver.h>

#include <sofa/simulation/common/MechanicalVisitor.h>

#include <sofa/component/linearsolver/FullMatrix.h>

namespace sofa
{

namespace component
{

namespace constraintset
{

using namespace sofa::defaulttype;
using namespace sofa::component::linearsolver;

class SOFA_CONSTRAINT_API ConstraintProblem
{
public:
    LPtrFullMatrix<double> W;
    FullVector<double> dFree, f;

    ConstraintProblem() : tolerance(0.00001), maxIterations(1000), dimension(0) {}
    virtual ~ConstraintProblem() {};

    double tolerance;
    int maxIterations;

    virtual void clear(int nbConstraints);
    int getDimension()	{ return dimension; }
    double** getW()		{ return W.lptr(); }
    double* getDfree()	{ return dFree.ptr(); }
    double* getF()		{ return f.ptr(); }

    virtual void solveTimed(double tolerance, int maxIt, double timeout) = 0;

protected:
    int dimension;
};


class SOFA_CONSTRAINT_API ConstraintSolverImpl : public sofa::core::behavior::ConstraintSolver
{
public:
    SOFA_CLASS(ConstraintSolverImpl, sofa::core::behavior::ConstraintSolver);

    virtual ConstraintProblem* getConstraintProblem() = 0;

    /// Do not use the following LCPs until the next call to this function.
    /// This is used to prevent concurent access to the LCP when using a LCPForceFeedback through an haptic thread.
    virtual void lockConstraintProblem(ConstraintProblem* p1, ConstraintProblem* p2=NULL) = 0;
};



/// Gets the vector of constraint violation values
class MechanicalGetConstraintViolationVisitor : public simulation::BaseMechanicalVisitor
{
public:

    MechanicalGetConstraintViolationVisitor(const core::ConstraintParams* params /* PARAMS FIRST */, BaseVector *v)
        : simulation::BaseMechanicalVisitor(params)
        , cparams(params)
        , m_v(v)
    {
#ifdef SOFA_DUMP_VISITOR_INFO
        setReadWriteVectors();
#endif
    }

    virtual Result fwdConstraintSet(simulation::Node* node, core::behavior::BaseConstraintSet* cSet)
    {
        if (core::behavior::BaseConstraintSet *c=dynamic_cast<core::behavior::BaseConstraintSet*>(cSet))
        {
            ctime_t t0 = begin(node, c);
            c->getConstraintViolation(cparams /* PARAMS FIRST */, m_v);
            end(node, c, t0);
        }
        return RESULT_CONTINUE;
    }

    /// This visitor must go through all mechanical mappings, even if isMechanical flag is disabled
    virtual bool stopAtMechanicalMapping(simulation::Node* /*node*/, core::BaseMapping* /*map*/)
    {
        return false; // !map->isMechanical();
    }

    /// Return a class name for this visitor
    /// Only used for debugging / profiling purposes
    virtual const char* getClassName() const { return "MechanicalGetConstraintViolationVisitor";}

#ifdef SOFA_DUMP_VISITOR_INFO
    void setReadWriteVectors()
    {
    }
#endif

private:
    /// Constraint parameters
    const sofa::core::ConstraintParams *cparams;

    /// Vector for constraint values
    BaseVector* m_v;
};

} // namespace constraintset

} // namespace component

} // namespace sofa

#endif
