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
#include "PythonScriptController.h"
using namespace sofa::component::controller;

#include "Binding_PythonScriptController.h"
#include "Binding_Base.h"

#include <sofa/simulation/common/Node.h>
using namespace sofa::simulation;

// These functions are empty ones;
// they are meant to be overriden by real python controller scripts


//#define LOG_UNIMPLEMENTED_METHODS   // prints a message each time a non-implemented (in the script) method is called


extern "C" PyObject * PythonScriptController_onLoaded(PyObject * /*self*/, PyObject * args)
{
    PyObject *pyNode;
    if (!PyArg_ParseTuple(args, "O",&pyNode))
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }

#ifdef LOG_UNIMPLEMENTED_METHODS
    PythonScriptController* obj=dynamic_cast<PythonScriptController*>(((PySPtr<Base>*)self)->object.get());
    Node* node=dynamic_cast<Node*>(((PySPtr<Base>*)pyNode)->object.get());
    std::cerr << obj->m_classname.getValueString() << ".onLoaded not implemented in " << obj->name.getValueString() << std::endl;
#endif

    Py_RETURN_NONE;
}

extern "C" PyObject * PythonScriptController_createGraph(PyObject * /*self*/, PyObject * args)
{
    PyObject *pyNode;
    if (!PyArg_ParseTuple(args, "O",&pyNode))
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }

#ifdef LOG_UNIMPLEMENTED_METHODS
    PythonScriptController* obj=dynamic_cast<PythonScriptController*>(((PySPtr<Base>*)self)->object.get());
    Node* node=dynamic_cast<Node*>(((PySPtr<Base>*)pyNode)->object.get());
    std::cerr << obj->m_classname.getValueString() << ".createGraph not implemented in " << obj->name.getValueString() << std::endl;
#endif

    Py_RETURN_NONE;
}

extern "C" PyObject * PythonScriptController_initGraph(PyObject * /*self*/, PyObject * args)
{
    PyObject *pyNode;
    if (!PyArg_ParseTuple(args, "O",&pyNode))
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }

#ifdef LOG_UNIMPLEMENTED_METHODS
    PythonScriptController* obj=dynamic_cast<PythonScriptController*>(((PySPtr<Base>*)self)->object.get());
    Node* node=dynamic_cast<Node*>(((PySPtr<Base>*)pyNode)->object.get());
    std::cerr << obj->m_classname.getValueString() << ".initGraph not implemented in " << obj->name.getValueString() << std::endl;
#endif

    Py_RETURN_NONE;
}

extern "C" PyObject * PythonScriptController_bwdInitGraph(PyObject * /*self*/, PyObject * args)
{
    PyObject *pyNode;
    if (!PyArg_ParseTuple(args, "O",&pyNode))
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }

#ifdef LOG_UNIMPLEMENTED_METHODS
    PythonScriptController* obj=dynamic_cast<PythonScriptController*>(((PySPtr<Base>*)self)->object.get());
    Node* node=dynamic_cast<Node*>(((PySPtr<Base>*)pyNode)->object.get());
    std::cerr << obj->m_classname.getValueString() << ".bwdInitGraph not implemented in " << obj->name.getValueString() << std::endl;
#endif

    Py_RETURN_NONE;
}

extern "C" PyObject * PythonScriptController_onBeginAnimationStep(PyObject * /*self*/, PyObject * args)
{
    double dt;
    if (!PyArg_ParseTuple(args, "d",&dt))
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }
    else
    {
      dt = 1;
    }

#ifdef LOG_UNIMPLEMENTED_METHODS
    PythonScriptController* obj=dynamic_cast<PythonScriptController*>(((PySPtr<Base>*)self)->object.get());
    std::cerr << obj->m_classname.getValueString() << ".onBeginAnimationStep not implemented in " << obj->name.getValueString() << std::endl;
#endif

    Py_RETURN_NONE;
}

extern "C" PyObject * PythonScriptController_onEndAnimationStep(PyObject * /*self*/, PyObject * args)
{
    double dt;
    if (!PyArg_ParseTuple(args, "d",&dt))
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }

#ifdef LOG_UNIMPLEMENTED_METHODS
    PythonScriptController* obj=dynamic_cast<PythonScriptController*>(((PySPtr<Base>*)self)->object.get());
    std::cerr << obj->m_classname.getValueString() << ".onEndAnimationStep not implemented in " << obj->name.getValueString() << std::endl;
#endif

    Py_RETURN_NONE;
}

extern "C" PyObject * PythonScriptController_storeResetState(PyObject * /*self*/, PyObject * /*args*/)
{

#ifdef LOG_UNIMPLEMENTED_METHODS
    PythonScriptController* obj=dynamic_cast<PythonScriptController*>(((PySPtr<Base>*)self)->object.get());
    std::cerr << obj->m_classname.getValueString() << ".storeresetState not implemented in " << obj->name.getValueString() << std::endl;
#endif

    Py_RETURN_NONE;
}

extern "C" PyObject * PythonScriptController_reset(PyObject * /*self*/, PyObject * /*args*/)
{

#ifdef LOG_UNIMPLEMENTED_METHODS
    PythonScriptController* obj=dynamic_cast<PythonScriptController*>(((PySPtr<Base>*)self)->object.get());
    std::cerr << obj->m_classname.getValueString() << ".reset not implemented in " << obj->name.getValueString() << std::endl;
#endif

    Py_RETURN_NONE;
}

extern "C" PyObject * PythonScriptController_cleanup(PyObject * /*self*/, PyObject * /*args*/)
{

#ifdef LOG_UNIMPLEMENTED_METHODS
    PythonScriptController* obj=dynamic_cast<PythonScriptController*>(((PySPtr<Base>*)self)->object.get());
    std::cerr << obj->m_classname.getValueString() << ".cleanup not implemented in " << obj->name.getValueString() << std::endl;
#endif

    Py_RETURN_NONE;
}

extern "C" PyObject * PythonScriptController_onGUIEvent(PyObject * /*self*/, PyObject * args)
{
    char* controlID;
    char* valueName;
    char* value;
    if (!PyArg_ParseTuple(args, "sss",&controlID,&valueName,&value))
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }

#ifdef LOG_UNIMPLEMENTED_METHODS
    PythonScriptController* obj=dynamic_cast<PythonScriptController*>(((PySPtr<Base>*)self)->object.get());
    std::cerr << obj->m_classname.getValueString() << ".onGUIEvent not implemented in " << obj->name.getValueString() << std::endl;
#endif

    Py_RETURN_NONE;
}

extern "C" PyObject * PythonScriptController_onKeyPressed(PyObject * /*self*/, PyObject * args)
{
    char k;
    if (!PyArg_ParseTuple(args, "c",&k))
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }

#ifdef LOG_UNIMPLEMENTED_METHODS
    PythonScriptController* obj=dynamic_cast<PythonScriptController*>(((PySPtr<Base>*)self)->object.get());
    std::cerr << obj->m_classname.getValueString() << ".onKeyPressed not implemented in " << obj->name.getValueString() << std::endl;
#endif

    Py_RETURN_NONE;
}

extern "C" PyObject * PythonScriptController_onKeyReleased(PyObject * /*self*/, PyObject * args)
{
    char k;
    if (!PyArg_ParseTuple(args, "c",&k))
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }

#ifdef LOG_UNIMPLEMENTED_METHODS
    PythonScriptController* obj=dynamic_cast<PythonScriptController*>(((PySPtr<Base>*)self)->object.get());
    std::cerr << obj->m_classname.getValueString() << ".onKeyReleased not implemented in " << obj->name.getValueString() << std::endl;
#endif

    Py_RETURN_NONE;
}

extern "C" PyObject * PythonScriptController_onMouseButtonLeft(PyObject * /*self*/, PyObject * args)
{
    int x,y;
    bool pressed;
    if (!PyArg_ParseTuple(args, "iib",&x,&y,&pressed))
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }

#ifdef LOG_UNIMPLEMENTED_METHODS
    PythonScriptController* obj=dynamic_cast<PythonScriptController*>(((PySPtr<Base>*)self)->object.get());
    std::cerr << obj->m_classname.getValueString() << ".onMouseButtonLeft not implemented in " << obj->name.getValueString() << std::endl;
#endif

    Py_RETURN_NONE;
}

extern "C" PyObject * PythonScriptController_onMouseButtonMiddle(PyObject * /*self*/, PyObject * args)
{
    int x,y;
    bool pressed;
    if (!PyArg_ParseTuple(args, "iib",&x,&y,&pressed))
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }

#ifdef LOG_UNIMPLEMENTED_METHODS
    PythonScriptController* obj=dynamic_cast<PythonScriptController*>(((PySPtr<Base>*)self)->object.get());
    std::cerr << obj->m_classname.getValueString() << ".onMouseButtonMiddle not implemented in " << obj->name.getValueString() << std::endl;
#endif

    Py_RETURN_NONE;
}

extern "C" PyObject * PythonScriptController_onMouseButtonRight(PyObject * /*self*/, PyObject * args)
{
    int x,y;
    bool pressed;
    if (!PyArg_ParseTuple(args, "iib",&x,&y,&pressed))
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }

#ifdef LOG_UNIMPLEMENTED_METHODS
    PythonScriptController* obj=dynamic_cast<PythonScriptController*>(((PySPtr<Base>*)self)->object.get());
    std::cerr << obj->m_classname.getValueString() << ".onMouseButtonRight not implemented in " << obj->name.getValueString() << std::endl;
#endif

    Py_RETURN_NONE;
}

extern "C" PyObject * PythonScriptController_onMouseWheel(PyObject * /*self*/, PyObject * args)
{
    int x,y,delta;
    if (!PyArg_ParseTuple(args, "iii",&x,&y,&delta))
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }

#ifdef LOG_UNIMPLEMENTED_METHODS
    PythonScriptController* obj=dynamic_cast<PythonScriptController*>(((PySPtr<Base>*)self)->object.get());
    std::cerr << obj->m_classname.getValueString() << ".onMouseWheel not implemented in " << obj->name.getValueString() << std::endl;
#endif

    Py_RETURN_NONE;
}

extern "C" PyObject * PythonScriptController_onScriptEvent(PyObject * /*self*/, PyObject * args)
{
    PyObject *pySenderNode;
    char* eventName;
    PyObject *pyData;
    if (!PyArg_ParseTuple(args, "OsO",&pySenderNode,&eventName,&pyData))
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }
    Node* senderNode=dynamic_cast<Node*>(((PySPtr<Base>*)pySenderNode)->object.get());
    if (!senderNode)
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }
    // TODO check pyData

#ifdef LOG_UNIMPLEMENTED_METHODS
    PythonScriptController* obj=dynamic_cast<PythonScriptController*>(((PySPtr<Base>*)self)->object.get());
    std::cerr << obj->m_classname.getValueString() << ".onScriptEvent not implemented in " << obj->name.getValueString() << std::endl;
#endif

    Py_RETURN_NONE;
}





SP_CLASS_METHODS_BEGIN(PythonScriptController)
SP_CLASS_METHOD(PythonScriptController,onLoaded)
SP_CLASS_METHOD(PythonScriptController,createGraph)
SP_CLASS_METHOD(PythonScriptController,initGraph)
SP_CLASS_METHOD(PythonScriptController,bwdInitGraph)
SP_CLASS_METHOD(PythonScriptController,onKeyPressed)
SP_CLASS_METHOD(PythonScriptController,onKeyReleased)
SP_CLASS_METHOD(PythonScriptController,onMouseButtonLeft)
SP_CLASS_METHOD(PythonScriptController,onMouseButtonRight)
SP_CLASS_METHOD(PythonScriptController,onMouseButtonMiddle)
SP_CLASS_METHOD(PythonScriptController,onMouseWheel)
SP_CLASS_METHOD(PythonScriptController,onBeginAnimationStep)
SP_CLASS_METHOD(PythonScriptController,onEndAnimationStep)
SP_CLASS_METHOD(PythonScriptController,storeResetState)
SP_CLASS_METHOD(PythonScriptController,reset)
SP_CLASS_METHOD(PythonScriptController,cleanup)
SP_CLASS_METHOD(PythonScriptController,onGUIEvent)
SP_CLASS_METHOD(PythonScriptController,onScriptEvent)
SP_CLASS_METHODS_END


SP_CLASS_TYPE_SPTR(PythonScriptController,PythonScriptController,Base)
