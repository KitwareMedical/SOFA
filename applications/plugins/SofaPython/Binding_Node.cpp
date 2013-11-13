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
*                               SOFA :: Plugins                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#include <sofa/simulation/common/Node.h>
#include <sofa/simulation/common/Simulation.h>
#include <sofa/core/objectmodel/KeypressedEvent.h>
#include <sofa/core/objectmodel/KeyreleasedEvent.h>
#include "ScriptEnvironment.h"
using namespace sofa::simulation;
#include <sofa/core/ExecParams.h>
using namespace sofa::core;

#include "Binding_Node.h"
#include "Binding_Context.h"
#include "PythonVisitor.h"
#include "PythonScriptEvent.h"

extern "C" PyObject * Node_executeVisitor(PyObject *self, PyObject * args)
{
    Node* node=dynamic_cast<Node*>(((PySPtr<Base>*)self)->object.get());

    PyObject* pyVisitor;
    if (!PyArg_ParseTuple(args, "O",&pyVisitor))
        Py_RETURN_NONE;
    PythonVisitor visitor(ExecParams::defaultInstance(),pyVisitor);
    node->executeVisitor(&visitor);

    Py_RETURN_NONE;
}

extern "C" PyObject * Node_getRoot(PyObject *self, PyObject * /*args*/)
{
    Node* node=dynamic_cast<Node*>(((PySPtr<Base>*)self)->object.get());

    // BaseNode is not binded in SofaPython, so getRoot is binded in Node instead of BaseNode
    return SP_BUILD_PYSPTR(node->getRoot());
}

// step the simulation
extern "C" PyObject * Node_simulationStep(PyObject * self, PyObject * args)
{
    Node* node=dynamic_cast<Node*>(((PySPtr<Base>*)self)->object.get());
    double dt;
    if (!PyArg_ParseTuple(args, "d",&dt))
        Py_RETURN_NONE;

    printf("Node_simulationStep node=%s dt=%f\n",node->getName().c_str(),(float)dt);

    getSimulation()->animate ( node, (SReal)dt );
//    simulation::getSimulation()->updateVisual( root );


    Py_RETURN_NONE;
}

extern "C" PyObject * Node_getChild(PyObject * self, PyObject * args)
{
    // BaseNode is not binded in SofaPython, so getChildNode is binded in Node instead of BaseNode
    Node* node=dynamic_cast<Node*>(((PySPtr<Base>*)self)->object.get());
    char *path;
    if (!PyArg_ParseTuple(args, "s",&path))
        Py_RETURN_NONE;
    if (!node || !path)
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }

    const objectmodel::BaseNode::Children& children = node->getChildren();
    Node *childNode = 0;
    // BaseNode ne pouvant pas être bindé en Python, et les BaseNodes des graphes étant toujours des Nodes,
    // on caste directement en Node.
    for (unsigned int i=0; i<children.size(); ++i)
        if (children[i]->getName() == path)
        {
            childNode = dynamic_cast<Node*>(children[i]);
            break;
        }
    if (!childNode)
    {
        printf("<SofaPython> Error: Node.getChildNode(%s) not found.\n",path);
        Py_RETURN_NONE;
    }
    return SP_BUILD_PYSPTR(childNode);
}

extern "C" PyObject * Node_getChildren(PyObject * self, PyObject * /*args*/)
{
    // BaseNode is not binded in SofaPython, so getChildNode is binded in Node instead of BaseNode
    Node* node=dynamic_cast<Node*>(((PySPtr<Base>*)self)->object.get());

    const objectmodel::BaseNode::Children& children = node->getChildren();

    // BaseNode ne pouvant pas être bindé en Python, et les BaseNodes des graphes étant toujours des Nodes,
    // on caste directement en Node.
    PyObject *list = PyList_New(children.size());

    for (unsigned int i=0; i<children.size(); ++i)
        PyList_SetItem(list,i,SP_BUILD_PYSPTR(children[i]));

    return list;
}

extern "C" PyObject * Node_getParents(PyObject * self, PyObject * /*args*/)
{
    // BaseNode is not binded in SofaPython, so getChildNode is binded in Node instead of BaseNode
    Node* node=dynamic_cast<Node*>(((PySPtr<Base>*)self)->object.get());

    const objectmodel::BaseNode::Children& parents = node->getParents();

    // BaseNode ne pouvant pas être bindé en Python, et les BaseNodes des graphes étant toujours des Nodes,
    // on caste directement en Node.
    PyObject *list = PyList_New(parents.size());

    for (unsigned int i=0; i<parents.size(); ++i)
        PyList_SetItem(list,i,SP_BUILD_PYSPTR(parents[i]));

    return list;
}

extern "C" PyObject * Node_createChild(PyObject *self, PyObject * args)
{
    Node* obj=dynamic_cast<Node*>(((PySPtr<Base>*)self)->object.get());
    char *nodeName;
    if (!PyArg_ParseTuple(args, "s",&nodeName))
        Py_RETURN_NONE;
    Node* child = obj->createChild(nodeName).get();
    ScriptEnvironment::nodeCreatedByScript(child);
    return SP_BUILD_PYSPTR(child);
}

extern "C" PyObject * Node_addObject(PyObject *self, PyObject * args)
{
    Node* node=dynamic_cast<Node*>(((PySPtr<Base>*)self)->object.get());
    PyObject* pyChild;
    if (!PyArg_ParseTuple(args, "O",&pyChild))
        Py_RETURN_NONE;
    BaseObject* object=dynamic_cast<BaseObject*>(((PySPtr<Base>*)pyChild)->object.get());
    if (!object)
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }
    node->addObject(object);

    if (node->isInitialized())
        printf("<SofaPython> WARNING Sofa.Node.addObject called on a node(%s) that is already initialized\n",node->getName().c_str());
    if (!ScriptEnvironment::isNodeCreatedByScript(node))
        printf("<SofaPython> WARNING Sofa.Node.addObject called on a node(%s) that is not created by the script\n",node->getName().c_str());

    //object->init();
    // plus besoin !! node->init(sofa::core::ExecParams::defaultInstance());

    Py_RETURN_NONE;
}

extern "C" PyObject * Node_removeObject(PyObject *self, PyObject * args)
{
    Node* node=dynamic_cast<Node*>(((PySPtr<Base>*)self)->object.get());
    PyObject* pyChild;
    if (!PyArg_ParseTuple(args, "O",&pyChild))
        Py_RETURN_NONE;
    BaseObject* object=dynamic_cast<BaseObject*>(((PySPtr<Base>*)pyChild)->object.get());
    if (!object)
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }
    node->removeObject(object);
    node->init(sofa::core::ExecParams::defaultInstance());

    Py_RETURN_NONE;
}

extern "C" PyObject * Node_addChild(PyObject *self, PyObject * args)
{
    Node* obj=dynamic_cast<Node*>(((PySPtr<Base>*)self)->object.get());
    PyObject* pyChild;
    if (!PyArg_ParseTuple(args, "O",&pyChild))
        Py_RETURN_NONE;
    BaseNode* child=dynamic_cast<BaseNode*>(((PySPtr<Base>*)pyChild)->object.get());
    if (!child)
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }
    obj->addChild(child);
    Py_RETURN_NONE;
}

extern "C" PyObject * Node_removeChild(PyObject *self, PyObject * args)
{
    Node* obj=dynamic_cast<Node*>(((PySPtr<Base>*)self)->object.get());
    PyObject* pyChild;
    if (!PyArg_ParseTuple(args, "O",&pyChild))
        Py_RETURN_NONE;
    BaseNode* child=dynamic_cast<BaseNode*>(((PySPtr<Base>*)pyChild)->object.get());
    if (!child)
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }
    obj->removeChild(child);
    Py_RETURN_NONE;
}

extern "C" PyObject * Node_moveChild(PyObject *self, PyObject * args)
{
    Node* obj=dynamic_cast<Node*>(((PySPtr<Base>*)self)->object.get());
    PyObject* pyChild;
    if (!PyArg_ParseTuple(args, "O",&pyChild))
        Py_RETURN_NONE;
    BaseNode* child=dynamic_cast<BaseNode*>(((PySPtr<Base>*)pyChild)->object.get());
    if (!child)
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }
    obj->moveChild(child);
    Py_RETURN_NONE;
}

extern "C" PyObject * Node_detachFromGraph(PyObject *self, PyObject * /*args*/)
{
    Node* node=dynamic_cast<Node*>(((PySPtr<Base>*)self)->object.get());
    node->detachFromGraph();
    Py_RETURN_NONE;
}

extern "C" PyObject * Node_sendScriptEvent(PyObject *self, PyObject * args)
{
    Node* node=dynamic_cast<Node*>(((PySPtr<Base>*)self)->object.get());
    PyObject* pyUserData;
    char* eventName;
    if (!PyArg_ParseTuple(args, "sO",&eventName,&pyUserData))
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }
    PythonScriptEvent event(node,eventName,pyUserData);
    dynamic_cast<Node*>(node->getRoot())->propagateEvent(sofa::core::ExecParams::defaultInstance(), &event);
    Py_RETURN_NONE;
}

extern "C" PyObject * Node_sendKeypressedEvent(PyObject *self, PyObject * args)
{
    Node* node=dynamic_cast<Node*>(((PySPtr<Base>*)self)->object.get());
    char* eventName;
    if (!PyArg_ParseTuple(args, "s",&eventName))
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }
    sofa::core::objectmodel::KeypressedEvent event(eventName ? eventName[0] : '\0');
    dynamic_cast<Node*>(node->getRoot())->propagateEvent(sofa::core::ExecParams::defaultInstance(), &event);
    Py_RETURN_NONE;
}

extern "C" PyObject * Node_sendKeyreleasedEvent(PyObject *self, PyObject * args)
{
    Node* node=dynamic_cast<Node*>(((PySPtr<Base>*)self)->object.get());
    char* eventName;
    if (!PyArg_ParseTuple(args, "s",&eventName))
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }
    sofa::core::objectmodel::KeyreleasedEvent event(eventName ? eventName[0] : '\0');
    dynamic_cast<Node*>(node->getRoot())->propagateEvent(sofa::core::ExecParams::defaultInstance(), &event);
    Py_RETURN_NONE;
}

SP_CLASS_METHODS_BEGIN(Node)
SP_CLASS_METHOD(Node,executeVisitor)
SP_CLASS_METHOD(Node,getRoot)
SP_CLASS_METHOD(Node,simulationStep)
SP_CLASS_METHOD(Node,getChild)
SP_CLASS_METHOD(Node,getChildren)
SP_CLASS_METHOD(Node,getParents)
SP_CLASS_METHOD(Node,createChild)
SP_CLASS_METHOD(Node,addObject)
SP_CLASS_METHOD(Node,removeObject)
SP_CLASS_METHOD(Node,addChild)
SP_CLASS_METHOD(Node,removeChild)
SP_CLASS_METHOD(Node,moveChild)
SP_CLASS_METHOD(Node,detachFromGraph)
SP_CLASS_METHOD(Node,sendScriptEvent)
SP_CLASS_METHOD(Node,sendKeypressedEvent)
SP_CLASS_METHOD(Node,sendKeyreleasedEvent)
SP_CLASS_METHODS_END

SP_CLASS_TYPE_SPTR(Node,Node,Context)
