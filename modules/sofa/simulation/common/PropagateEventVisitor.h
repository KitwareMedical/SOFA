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
//
// C++ Interface: PropagateEventVisitor
//
// Description:
//
//
// Author: The SOFA team </www.sofa-framework.org>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SOFA_SIMULATION_PROPAGATEEVENTACTION_H
#define SOFA_SIMULATION_PROPAGATEEVENTACTION_H

#include <sofa/simulation/common/Visitor.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/core/ExecParams.h>

namespace sofa
{

namespace simulation
{

/**
Visitor used to propagate an event in the the data structure.
Propagation is done top-down until the event is handled.

	@author The SOFA team </www.sofa-framework.org>
*/
class SOFA_SIMULATION_COMMON_API PropagateEventVisitor : public sofa::simulation::Visitor
{
public:
    PropagateEventVisitor(const core::ExecParams* params /* PARAMS FIRST */, sofa::core::objectmodel::Event* e);

    ~PropagateEventVisitor();

    Visitor::Result processNodeTopDown(simulation::Node* node);
    void processObject(simulation::Node*, core::objectmodel::BaseObject* obj);

    virtual const char* getClassName() const { return "PropagateEventVisitor"; }
    virtual std::string getInfos() const { return std::string(m_event->getClassName());  }
protected:
    sofa::core::objectmodel::Event* m_event;
};


}

}

#endif
