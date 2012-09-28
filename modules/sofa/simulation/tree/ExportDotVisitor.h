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
#ifndef SOFA_SIMULATION_TREE_EXPORTDOTACTION_H
#define SOFA_SIMULATION_TREE_EXPORTDOTACTION_H

#include <sofa/simulation/tree/GNodeVisitor.h>
#include <string>
#include <iostream>


namespace sofa
{

namespace simulation
{

namespace tree
{

class SOFA_SIMULATION_TREE_API ExportDotVisitor : public GNodeVisitor
{
public:
    std::ostream* out;

    bool showNode;
    bool showObject;
    bool showBehaviorModel;
    bool showCollisionModel;
    bool showVisualModel;
    bool showMapping;
    bool showContext;
    bool showCollisionPipeline;
    bool showSolver;
    bool showMechanicalState;
    bool showForceField;
    bool showInteractionForceField;
    bool showConstraint;
    bool showMass;
    bool showTopology;
    bool showMechanicalMapping;

    bool labelNodeName;
    bool labelNodeClass;
    bool labelObjectName;
    bool labelObjectClass;

    ExportDotVisitor(const sofa::core::ExecParams* params /* PARAMS FIRST */, std::ostream* out);
    ~ExportDotVisitor();

    void processObject(GNode* node, core::objectmodel::BaseObject* obj);

    virtual Result processNodeTopDown(GNode* node);
    virtual void processNodeBottomUp(GNode* node);

    virtual const char* getClassName() const { return "ExportDotVisitor"; }

protected:

    /// None names in output
    std::map<core::objectmodel::Base*, std::string> names;
    /// Next indice available for duplicated names
    std::map<std::string, int> nextIndex;

    /// Test if a node should be displayed
    bool display(GNode* node, const char** color=NULL);

    /// Test if an object should be displayed
    bool display(core::objectmodel::BaseObject* obj, const char** color=NULL);

    /// Find the node or object a given object should be attached to.
    /// This is the parent node if it is displayed, otherwise it is the attached MechanicalState or Solver.
    /// Returns an empty string if not found.
    std::string getParentName(core::objectmodel::BaseObject* obj);

    /// Compute the name of a given node or object
    std::string getName(core::objectmodel::Base* o, std::string prefix);

    /// Compute the name of a given node
    std::string getName(GNode* node);

    /// Compute the name of a given object
    std::string getName(core::objectmodel::BaseObject* obj);

};

} // namespace tree

} // namespace simulation

} // namespace sofa

#endif
