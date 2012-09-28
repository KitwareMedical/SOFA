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
#ifndef SOFA_SIMULATION_ANIMATEACTION_H
#define SOFA_SIMULATION_ANIMATEACTION_H

#include <sofa/simulation/common/common.h>
#include <sofa/simulation/common/Visitor.h>
#include <sofa/core/VecId.h>
#include <sofa/core/MultiVecId.h>
#include <sofa/core/ExecParams.h>
#include <sofa/core/MechanicalParams.h>
#include <sofa/core/BehaviorModel.h>
#include <sofa/core/behavior/BaseInteractionForceField.h>
#include <sofa/core/behavior/OdeSolver.h>
#include <sofa/core/behavior/BaseAnimationLoop.h>
#include <sofa/core/collision/Pipeline.h>

using namespace sofa::core;


namespace sofa
{

namespace simulation
{

class SOFA_SIMULATION_COMMON_API AnimateVisitor : public Visitor
{

protected :
    double dt;
#ifdef SOFA_HAVE_EIGEN2
    bool firstNodeVisited;
#endif
public:
    AnimateVisitor(const core::ExecParams* params = ExecParams::defaultInstance());
    AnimateVisitor(const core::ExecParams* params /* PARAMS FIRST  = ExecParams::defaultInstance()*/, double dt);

    void setDt(double v) { dt = v; }
    double getDt() const { return dt; }

    virtual void processCollisionPipeline(simulation::Node* node, core::collision::Pipeline* obj);
    virtual void processBehaviorModel(simulation::Node* node, core::BehaviorModel* obj);
    virtual void fwdInteractionForceField(simulation::Node* node, core::behavior::BaseInteractionForceField* obj);
    virtual void processOdeSolver(simulation::Node* node, core::behavior::OdeSolver* obj);

    virtual Result processNodeTopDown(simulation::Node* node);
    //virtual void processNodeBottomUp(simulation::Node* node);

    /// Specify whether this action can be parallelized.
    virtual bool isThreadSafe() const { return true; }

    /// Return a category name for this action.
    /// Only used for debugging / profiling purposes
    virtual const char* getCategoryName() const { return "animate"; }
    virtual const char* getClassName() const { return "AnimateVisitor"; }
};

} // namespace simulation

} // namespace sofa

#endif
