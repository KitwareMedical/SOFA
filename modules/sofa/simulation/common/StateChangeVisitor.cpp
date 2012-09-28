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

#include <sofa/simulation/common/StateChangeVisitor.h>
#include <sofa/simulation/common/Node.h>
#include <sofa/core/topology/TopologicalMapping.h>

namespace sofa
{

namespace simulation
{

StateChangeVisitor::StateChangeVisitor(const sofa::core::ExecParams* params /* PARAMS FIRST */, core::topology::Topology* source)
    : Visitor(params), root(true), source(source)
{
}

void StateChangeVisitor::processStateChange(core::behavior::BaseMechanicalState* obj)
{
    obj->handleStateChange(source);
}

Visitor::Result StateChangeVisitor::processNodeTopDown(simulation::Node* node)
{
    if (!root && node->mechanicalMapping)
    {
        if (!node->mechanicalMapping->sameTopology())
        {
            // stop all topological computations
            return RESULT_PRUNE;
        }
    }
    if (node->mechanicalState)
    {
        this->processStateChange(node->mechanicalState);
    }

    root = false; // now we process child nodes
    return RESULT_CONTINUE; // continue the propagation of state changes
}

} // namespace simulation

} // namespace sofa

