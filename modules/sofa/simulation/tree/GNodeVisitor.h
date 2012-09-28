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
#ifndef SOFA_SIMULATION_TREE_GNODEVISITOR_H
#define SOFA_SIMULATION_TREE_GNODEVISITOR_H

#include <sofa/simulation/common/Visitor.h>
#include <sofa/simulation/tree/GNode.h>

namespace sofa
{

namespace simulation
{

namespace tree
{

/**
Base class for the Visitors which deal with GNodes specifically rather than Node.

	@author The SOFA team </www.sofa-framework.org>
*/
class SOFA_SIMULATION_TREE_API GNodeVisitor : public sofa::simulation::Visitor
{
public:
    GNodeVisitor(const sofa::core::ExecParams* params);

    ~GNodeVisitor();

    /// Callback method called when decending to a new node. Recursion will stop if this method returns RESULT_PRUNE
    virtual Result processNodeTopDown(GNode* /*node*/) { return RESULT_CONTINUE; }

    /// Callback method called after child node have been processed and before going back to the parent node.
    virtual void processNodeBottomUp(GNode* /*node*/) { }

    /// Callback method called when decending to a new node. Recursion will stop if this method returns RESULT_PRUNE
    /// This version is offered a LocalStorage to store temporary data
    virtual Result processNodeTopDown(GNode* node, LocalStorage*) { return processNodeTopDown(node); }

    /// Callback method called after child node have been processed and before going back to the parent node.
    /// This version is offered a LocalStorage to store temporary data
    virtual void processNodeBottomUp(GNode* node, LocalStorage*) { processNodeBottomUp(node); }

    /// Callback method called when decending to a new node. Recursion will stop if this method returns RESULT_PRUNE
    virtual Result processNodeTopDown(simulation::Node* node)
    {
        GNode* g = dynamic_cast<GNode*>(node);
        if (!g)
        {
            std::cerr << "GNodeVisitor: node is not a GNode !\n";
            return RESULT_PRUNE;
        }
        else
        {
            return processNodeTopDown(g);
        }
    }

    /// Callback method called after child node have been processed and before going back to the parent node.
    virtual void processNodeBottomUp(simulation::Node* node)
    {
        GNode* g = dynamic_cast<GNode*>(node);
        if (!g)
        {
            std::cerr << "GNodeVisitor: node is not a GNode !\n";
        }
        else
        {
            processNodeBottomUp(g);
        }
    }

    virtual const char* getClassName() const { return "GNodeVisitor"; }
    /// Helper method to enumerate objects in the given list. The callback gets the pointer to node
    template < class Visit, class Container, class Object >
    void for_each(Visit* visitor, GNode* ctx, const Container& list, void (Visit::*fn)(GNode*, Object*))
    {
        for (typename Container::iterator it=list.begin(); it != list.end(); ++it)
        {
            typename Container::pointed_type* ptr = &*(*it);
            if(testTags(ptr))
            {
                debug_write_state_before(ptr);
                ctime_t t=begin(ctx, ptr);
                (visitor->*fn)(ctx, ptr);
                end(ctx, ptr, t);
                debug_write_state_after(ptr);
            }
        }
    }

    /// Helper method to enumerate objects in the given list. The callback gets the pointer to node
    template < class Visit, class Container, class Object >
    Visitor::Result for_each_r(Visit* visitor, GNode* ctx, const Container& list, Visitor::Result (Visit::*fn)(GNode*, Object*))
    {
        Visitor::Result res = Visitor::RESULT_CONTINUE;
        for (typename Container::iterator it=list.begin(); it != list.end(); ++it)
        {
            typename Container::pointed_type* ptr = &*(*it);
            if(testTags(ptr))
            {
                debug_write_state_before(ptr);
                ctime_t t=begin(ctx, ptr);
                res = (visitor->*fn)(ctx, ptr);
                end(ctx, ptr, t);
                debug_write_state_after(ptr);
            }
        }
        return res;
    }

};

}

}

}

#endif
