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
#include <sofa/component/linearsolver/GraphScatteredTypes.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/simulation/common/MechanicalVisitor.h>
#include <sofa/simulation/common/MechanicalMatrixVisitor.h>
#include <sofa/simulation/common/MechanicalVPrintVisitor.h>
#include <sofa/simulation/common/VelocityThresholdVisitor.h>
#include <sofa/core/behavior/LinearSolver.h>

#include <stdlib.h>
#include <math.h>

namespace sofa
{

namespace component
{

namespace linearsolver
{

using sofa::core::behavior::LinearSolver;
using sofa::core::objectmodel::BaseContext;

void GraphScatteredMatrix::apply(GraphScatteredVector& res, GraphScatteredVector& x)
{
    // matrix-vector product
#if 1
    // new more powerful visitors
    parent->propagateDxAndResetDf(x,res);
    parent->addMBKdx(res,parent->mparams.mFactor(),parent->mparams.bFactor(),parent->mparams.kFactor(), false); // df = (m M + b B + k K) dx

#else
    parent->propagateDx(x);          // dx = p
    parent->computeDf(res);            // q = K p

    if (parent->mparams.kFactor() != 1.0)
        res *= parent->mparams.kFactor(); // q = k K p

    // apply global Rayleigh damping
    if (parent->mparams.mFactor() == 1.0)
    {
        parent->addMdx(res); // no need to propagate p as dx again
    }
    else if (parent->mparams.mFactor() != 0.0)
    {
        parent->addMdx(res,core::MultiVecDerivId(),parent->mparams.mFactor()); // no need to propagate p as dx again
    }
    // q = (m M + k K) p

    /// @TODO: non-rayleigh damping (i.e. the B factor)
#endif

    // filter the product to take the constraints into account
    //
    parent->projectResponse(res);     // q is projected to the constrained space
}

#ifdef SOFA_SMP
void GraphScatteredMatrix::apply(ParallelGraphScatteredVector& res, ParallelGraphScatteredVector& x)
{
    // matrix-vector product
#if 0
    // may not have their SMP version
    // new more powerful visitors
    parent->propagateDxAndResetDf(x,res);
    parent->addMBKdx(res,mFact,bFact,kFact, false); // df = (m M + b B + k K) dx

#else
    parent->propagateDx(x);          // dx = p
    parent->computeDf(res);            // q = K p

    if (parent->mparams.kFactor() != 1.0)
        res *= parent->mparams.kFactor(); // q = k K p

    // apply global Rayleigh damping
    // TODO : check
//    if (parent->mparams.kFactor() == 1.0)
//    {
//        parent->addMdx(res); // no need to propagate p as dx again
//    }
    /* else */ if (parent->mparams.mFactor() != 0.0)
    {
        parent->addMdx(res,core::MultiVecDerivId::null(),parent->mparams.mFactor()); // no need to propagate p as dx again
    }
    // q = (m M + k K) p

    /// @TODO: non-rayleigh damping (i.e. the B factor)
#endif

    // filter the product to take the constraints into account
    //
    parent->projectResponse(res);     // q is projected to the constrained space
}
#endif

} // namespace linearsolver

} // namespace component

} // namespace sofa
