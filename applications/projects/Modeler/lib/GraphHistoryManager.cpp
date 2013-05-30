/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU General Public License as published by the Free  *
* Software Foundation; either version 2 of the License, or (at your option)   *
* any later version.                                                          *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    *
* more details.                                                               *
*                                                                             *
* You should have received a copy of the GNU General Public License along     *
* with this program; if not, write to the Free Software Foundation, Inc., 51  *
* Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.                   *
*******************************************************************************
*                            SOFA :: Applications                             *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/

#include "GraphHistoryManager.h"
#include "GraphModeler.h"

#include <sofa/simulation/common/Simulation.h>
#include <sofa/simulation/common/XMLPrintVisitor.h>
#include <sofa/helper/system/SetDirectory.h>

namespace sofa
{

namespace gui
{

namespace qt
{
GraphHistoryManager::GraphHistoryManager(GraphModeler *g): graph(g)
{
};

GraphHistoryManager::~GraphHistoryManager()
{
    for (unsigned int i=0; i<historyOperation.size(); ++i) undo();
}

void GraphHistoryManager::operationPerformed(Operation& o)
{
    clearHistoryUndo();
    historyOperation.push_back(o);

    emit graphModified(true);
    emit undoEnabled(true);
}


void GraphHistoryManager::undo()
{
    if (historyOperation.empty()) return;
    Operation o=historyOperation.back();
    historyOperation.pop_back();

    undoOperation(o);
    historyUndoOperation.push_back(o);

    emit undoEnabled(historyOperation.size());
    emit redoEnabled(true);

    emit graphModified(historyOperation.size());
}

void GraphHistoryManager::redo()
{
    if (historyUndoOperation.empty()) return;
    Operation o=historyUndoOperation.back();
    historyUndoOperation.pop_back();

    undoOperation(o);
    historyOperation.push_back(o);

    emit redoEnabled(historyUndoOperation.size());
    emit undoEnabled(true);

    emit graphModified(historyOperation.size());
}

void GraphHistoryManager::undoOperation(Operation &o)
{
    std::string message;
    switch(o.ID)
    {
    case Operation::DELETE_OBJECT:
        o.parent->addObject(sofa::core::objectmodel::SPtr_dynamic_cast<BaseObject>(o.sofaComponent));
        graph->moveItem(graph->graphListener->items[o.sofaComponent.get()],graph->graphListener->items[o.above.get()]);
        o.ID = Operation::ADD_OBJECT;
        message=std::string("Undo Delete NODE ") + " (" + o.sofaComponent->getClassName() + ") " + o.sofaComponent->getName();
        emit displayMessage(message);
        break;
    case Operation::DELETE_Node:
        o.parent->addChild(sofa::core::objectmodel::SPtr_dynamic_cast<Node>(o.sofaComponent));
        //If the node is not alone below another parent node
        graph->moveItem(graph->graphListener->items[o.sofaComponent.get()],graph->graphListener->items[o.above.get()]);
        o.ID = Operation::ADD_Node;

        message=std::string("Undo Add NODE ") +" (" + o.sofaComponent->getClassName() + ") " + o.sofaComponent->getName();
        emit displayMessage(message);
        break;
    case Operation::ADD_OBJECT:
        o.parent=graph->getNode(graph->graphListener->items[o.sofaComponent.get()]);
        o.above=graph->getComponentAbove(graph->graphListener->items[o.sofaComponent.get()]);
        o.parent->removeObject(sofa::core::objectmodel::SPtr_dynamic_cast<BaseObject>(o.sofaComponent));
        o.ID = Operation::DELETE_OBJECT;

        message=std::string("Undo Delete OBJECT ") +" (" + o.sofaComponent->getClassName() + ") " + o.sofaComponent->getName();
        emit displayMessage(message);
        break;
    case Operation::ADD_Node:
        o.parent=graph->getNode(graph->graphListener->items[o.sofaComponent.get()]->parent());
        o.above=graph->getComponentAbove(graph->graphListener->items[o.sofaComponent.get()]);
        o.parent->removeChild(sofa::core::objectmodel::SPtr_dynamic_cast<Node>(o.sofaComponent));
        o.ID = Operation::DELETE_Node;

        message=std::string("Undo Add OBJECT ") + " (" + o.sofaComponent->getClassName() + ") " + o.sofaComponent->getName();
        emit displayMessage(message);
        break;
    case Operation::COMPONENT_MODIFICATION:
    {
        std::string previousState=componentState(o.sofaComponent.get());
        const std::string &modifDone=setComponentState(o.sofaComponent.get(), o.info);
        QString name=QString(o.sofaComponent->getClassName().c_str()) + QString(" ") + QString(o.sofaComponent->getName().c_str());
        graph->graphListener->items[o.sofaComponent.get()]->setText(0, name);
        graph->setSelected(graph->graphListener->items[o.sofaComponent.get()],true);
        message=std::string("Undo Modifications on OBJECT ") + " (" + o.sofaComponent->getClassName() + ") " + o.sofaComponent->getName() + " | " + modifDone;
        emit displayMessage(message);
        o.info=previousState;
        break;
    }
    case Operation::NODE_MODIFICATION:
    {
        std::string previousState=componentState(o.sofaComponent.get());
        const std::string &modifDone=setComponentState(o.sofaComponent.get(), o.info);
        QString name=QString(o.sofaComponent->getName().c_str());
        graph->graphListener->items[o.sofaComponent.get()]->setText(0, name);
        graph->setSelected(graph->graphListener->items[o.sofaComponent.get()],true);
        message=std::string("Undo Modifications on NODE ") + o.sofaComponent->getName() + " | " + modifDone;
        emit displayMessage(message);
        o.info=previousState;
        break;
    }
    }
}

std::string GraphHistoryManager::componentState(Base *base) const
{
    std::ostringstream out;
    std::map<std::string, std::string*> datas;
    base->writeDatas(datas);

    std::map<std::string, std::string*>::const_iterator it;
    for (it=datas.begin(); it!=datas.end(); ++it)
    {
        if (BaseData *d=base->findField(it->first))
            out << it->first << "=\"" << d->getValueString() << "\" ";
    }
    return out.str();
}

std::string GraphHistoryManager::setComponentState(Base *base, const std::string &datasStr)
{
    std::string modifications;
    std::istringstream in(datasStr);
    while (!in.eof())
    {
        char k;
        std::string dataName;
        while (!in.eof())
        {
            in >> k;
            if (k == '=') break;
            dataName += k;
        }
        std::string dataValue;
        in >> k;//peeking the "
        while (!in.eof())
        {
            k=in.get();
            if (k == '\"') break;
            dataValue += k;
        }
        if (in.eof()) break;
        BaseData *data=base->findField(dataName);
        if (data)
        {
            const std::string &currentValue = data->getValueString();
            if (currentValue != dataValue)
            {
                data->read(dataValue);
                modifications += data->getName() + ", ";
            }
        }

    }
    if (!modifications.empty()) modifications=modifications.substr(0,modifications.size()-2);
    return modifications;
}

void GraphHistoryManager::beginModification(sofa::core::objectmodel::Base* base)
{
    componentPriorModificationState.insert(std::make_pair(base,componentState(base) ));
}

void GraphHistoryManager::endModification(sofa::core::objectmodel::Base* base)
{
    if (componentPriorModificationState[base] != componentState(base))
    {
        if (dynamic_cast<BaseObject*>(base))
        {
            Operation o(base, Operation::COMPONENT_MODIFICATION);
            o.info=componentPriorModificationState[base];
            operationPerformed(o);
        }
        else
        {
            Operation o(base, Operation::NODE_MODIFICATION);
            o.info=componentPriorModificationState[base];
            operationPerformed(o);
        }
    }
    componentPriorModificationState.erase(base);
}

void GraphHistoryManager::graphClean()
{
    clearHistoryUndo();
    clearHistory();
    emit graphModified(false);
}

void GraphHistoryManager::clearHistory()
{
    for ( int i=historyOperation.size()-1; i>=0; --i)
    {
        if (historyOperation[i].ID == Operation::DELETE_OBJECT)
            historyOperation[i].sofaComponent.reset();
        else if (historyOperation[i].ID == Operation::DELETE_Node)
        {
            Node *n=dynamic_cast<Node*>(historyOperation[i].sofaComponent.get());
            simulation::getSimulation()->unload(n);
            historyOperation[i].sofaComponent.reset();
        }
    }
    historyOperation.clear();
    emit( undoEnabled(false) );
}

void GraphHistoryManager::clearHistoryUndo()
{
    for ( int i=historyUndoOperation.size()-1; i>=0; --i)
    {
        if (historyUndoOperation[i].ID == Operation::DELETE_OBJECT)
            historyUndoOperation[i].sofaComponent.reset();
        else if (historyUndoOperation[i].ID == Operation::DELETE_Node)
        {
            Node *n=dynamic_cast<Node*>(historyUndoOperation[i].sofaComponent.get());
            simulation::getSimulation()->unload(n);
            historyUndoOperation[i].sofaComponent.reset();
        }
    }
    historyUndoOperation.clear();
    emit( redoEnabled(false) );
}
}
}
}
