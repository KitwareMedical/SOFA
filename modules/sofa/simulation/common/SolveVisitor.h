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
#ifndef SOFA_SIMULATION_SOLVEACTION_H
#define SOFA_SIMULATION_SOLVEACTION_H

#if !defined(__GNUC__) || (__GNUC__ > 3 || (_GNUC__ == 3 && __GNUC_MINOR__ > 3))
#pragma once
#endif

#include <sofa/simulation/common/Visitor.h>
#include <sofa/core/behavior/OdeSolver.h>

namespace sofa
{

namespace simulation
{

/** Used by the animation loop: send the solve signal to the others solvers

 */
class SOFA_SIMULATION_COMMON_API SolveVisitor : public Visitor
{
public:
    SolveVisitor(const sofa::core::ExecParams* params /* PARAMS FIRST */, double _dt) : Visitor(params), dt(_dt), freeMotion(false) {}
    SolveVisitor(const sofa::core::ExecParams* params /* PARAMS FIRST */, double _dt, bool free) : Visitor(params), dt(_dt), freeMotion(free) {}
    virtual void processSolver(simulation::Node* node, core::behavior::OdeSolver* b);
    virtual Result processNodeTopDown(simulation::Node* node);

    /// Specify whether this action can be parallelized.
    virtual bool isThreadSafe() const { return true; }

    /// Return a category name for this action.
    /// Only used for debugging / profiling purposes
    virtual const char* getCategoryName() const { return "behavior update position"; }
    virtual const char* getClassName() const { return "SolveVisitor"; }

    void setDt(double _dt) {dt = _dt;}
    double getDt() {return dt;}
protected:
    double dt;
    bool freeMotion;
};

} // namespace simulation

} // namespace sofa

#endif
