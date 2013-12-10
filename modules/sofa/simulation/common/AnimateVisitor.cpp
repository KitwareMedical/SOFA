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
#include <sofa/simulation/common/AnimateVisitor.h>
#include <sofa/simulation/common/MechanicalVisitor.h>
#include <sofa/simulation/common/CollisionVisitor.h>

#include <sofa/simulation/common/PropagateEventVisitor.h>
#include <sofa/simulation/common/CollisionBeginEvent.h>
#include <sofa/simulation/common/CollisionEndEvent.h>
#include <sofa/simulation/common/IntegrateBeginEvent.h>
#include <sofa/simulation/common/IntegrateEndEvent.h>
#include <sofa/simulation/common/PropagateEventVisitor.h>


#include <sofa/helper/AdvancedTimer.h>

//#include "MechanicalIntegration.h"

namespace sofa
{

namespace simulation
{


AnimateVisitor::AnimateVisitor(const core::ExecParams* params /* PARAMS FIRST */, double dt)
    : Visitor(params)
    , dt(dt)
#ifdef SOFA_HAVE_EIGEN2
    , firstNodeVisited(false)
#endif
{
}

AnimateVisitor::AnimateVisitor(const core::ExecParams* params)
    : Visitor(params)
    , dt(0)
#ifdef SOFA_HAVE_EIGEN2
    , firstNodeVisited(false)
#endif
{
}

void AnimateVisitor::processBehaviorModel(simulation::Node*, core::BehaviorModel* obj)
{
    sofa::helper::AdvancedTimer::stepBegin("BehaviorModel",obj);

    obj->updatePosition(getDt());
    sofa::helper::AdvancedTimer::stepEnd("BehaviorModel",obj);
}

void AnimateVisitor::fwdInteractionForceField(simulation::Node*, core::behavior::BaseInteractionForceField* obj)
{
    //cerr<<"AnimateVisitor::IFF "<<obj->getName()<<endl;
    sofa::helper::AdvancedTimer::stepBegin("InteractionFF",obj);

    MultiVecDerivId   ffId      = VecDerivId::externalForce();
    MechanicalParams mparams; // = MechanicalParams::defaultInstance();
    mparams.setDt(this->dt);
    obj->addForce(&mparams, ffId);

    sofa::helper::AdvancedTimer::stepEnd("InteractionFF",obj);
}

void AnimateVisitor::processCollisionPipeline(simulation::Node* node, core::collision::Pipeline* obj)
{
     sofa::helper::AdvancedTimer::stepBegin("Collision",obj);

    sofa::helper::AdvancedTimer::stepBegin("begin collision",obj);
    {
        CollisionBeginEvent evBegin;
        PropagateEventVisitor eventPropagation( params /* PARAMS FIRST */, &evBegin);
        eventPropagation.execute(node->getContext());
    }
    sofa::helper::AdvancedTimer::stepEnd("begin collision",obj);

    CollisionVisitor act(this->params);
    node->execute(&act);    

    sofa::helper::AdvancedTimer::stepBegin("end collision",obj);
    {
        CollisionEndEvent evEnd;
        PropagateEventVisitor eventPropagation( params /* PARAMS FIRST */, &evEnd);
        eventPropagation.execute(node->getContext());
    }
    sofa::helper::AdvancedTimer::stepEnd("end collision",obj);

    sofa::helper::AdvancedTimer::stepEnd("Collision",obj);
}

void AnimateVisitor::processOdeSolver(simulation::Node* node, core::behavior::OdeSolver* solver)
{
    sofa::helper::AdvancedTimer::stepBegin("Mechanical",node);
    /*    MechanicalIntegrationVisitor act(getDt());
        node->execute(&act);*/
//  cerr<<"AnimateVisitor::processOdeSolver "<<solver->getName()<<endl;
    solver->solve(params /* PARAMS FIRST */, getDt());
    sofa::helper::AdvancedTimer::stepEnd("Mechanical",node);
}

Visitor::Result AnimateVisitor::processNodeTopDown(simulation::Node* node)
{

    //cerr<<"AnimateVisitor::process Node  "<<node->getName()<<endl;
    if (!node->isActive()) return Visitor::RESULT_PRUNE;
#ifdef SOFA_HAVE_EIGEN2
    if (!firstNodeVisited)
    {
        firstNodeVisited=true;

//        core::behavior::BaseAnimationLoop* presenceAnimationManager;
//        node->get(presenceAnimationManager, core::objectmodel::BaseContext::SearchDown);
//        if (!presenceAnimationManager)
//        {
//          std::cerr << "AnimateVisitor::processNodeTopDown, ERROR: no BaseAnimationLoop found while searching down from node: " << node->getName() << std::endl;

//        }
        sofa::core::MechanicalParams mparams(*this->params);
        mparams.setDt(dt);
        MechanicalResetConstraintVisitor resetConstraint(&mparams);
        node->execute(&resetConstraint);
    }
#endif

    if (dt == 0) setDt(node->getDt());

    if (node->collisionPipeline != NULL)
    {

        //ctime_t t0 = begin(node, node->collisionPipeline);
#ifndef SOFA_SMP
        {
            CollisionBeginEvent evBegin;
            PropagateEventVisitor eventPropagation(this->params /* PARAMS FIRST */, &evBegin);
            eventPropagation.execute(node);
        }
        processCollisionPipeline(node, node->collisionPipeline);
        {
            CollisionEndEvent evEnd;
            PropagateEventVisitor eventPropagation(this->params /* PARAMS FIRST */, &evEnd);
            eventPropagation.execute(node);
        }
#endif
        //end(node, node->collisionPipeline, t0);
    }
    /*	if (node->solver != NULL)
    	{
    		ctime_t t0 = begin(node, node->solver);
    		processOdeSolver(node, node->solver);
    		end(node, node->solver, t0);
    		return RESULT_PRUNE;
            }*/
    if (!node->solver.empty() )
    {
        sofa::helper::AdvancedTimer::StepVar timer("Mechanical",node);
        double nextTime = node->getTime() + dt;


        {
            IntegrateBeginEvent evBegin;
            PropagateEventVisitor eventPropagation( this->params /* PARAMS FIRST */, &evBegin);
            eventPropagation.execute(node);
        }

        MechanicalBeginIntegrationVisitor beginVisitor(this->params /* PARAMS FIRST */, dt);
        node->execute(&beginVisitor);

        sofa::core::MechanicalParams m_mparams(*this->params);
        m_mparams.setDt(dt);

#ifdef SOFA_HAVE_EIGEN2
        {
            unsigned int constraintId=0;
            core::ConstraintParams cparams;
            //MechanicalAccumulateConstraint(&m_mparams /* PARAMS FIRST */, constraintId, VecCoordId::position()).execute(node);
            simulation::MechanicalAccumulateConstraint(&cparams /* PARAMS FIRST */, core::MatrixDerivId::holonomicC(),constraintId).execute(node);
        }
#endif

        for( unsigned i=0; i<node->solver.size(); i++ )
        {
            ctime_t t0 = begin(node, node->solver[i]);
            //cerr<<"AnimateVisitor::processNodeTpDown  solver  "<<node->solver[i]->getName()<<endl;
            node->solver[i]->solve(params /* PARAMS FIRST */, getDt());
            end(node, node->solver[i], t0);
        }

        MechanicalPropagatePositionAndVelocityVisitor(&m_mparams /* PARAMS FIRST */, nextTime,VecCoordId::position(),VecDerivId::velocity(),
#ifdef SOFA_SUPPORT_MAPPED_MASS
                VecDerivId::dx(),
#endif
                true).execute( node );

        MechanicalEndIntegrationVisitor endVisitor(this->params /* PARAMS FIRST */, dt);
        node->execute(&endVisitor);

        {
            IntegrateEndEvent evBegin;
            PropagateEventVisitor eventPropagation(this->params /* PARAMS FIRST */, &evBegin);
            eventPropagation.execute(node);
        }

        return RESULT_PRUNE;
    }
    /*
    if (node->mechanicalModel != NULL)
    {
    	std::cerr << "Graph Error: MechanicalState without solver." << std::endl;
    	return RESULT_PRUNE;
    }
    */
    {
        // process InteractionForceFields
        for_each(this, node, node->interactionForceField, &AnimateVisitor::fwdInteractionForceField);
        return RESULT_CONTINUE;
    }
}

// void AnimateVisitor::processNodeBottomUp(simulation::Node* node)
// {
//     node->setTime( node->getTime() + node->getDt() );
// }


} // namespace simulation

} // namespace sofa

