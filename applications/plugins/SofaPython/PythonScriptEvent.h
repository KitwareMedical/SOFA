/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *
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
*                              SOFA :: Plugins                                *
*                                                                             *
* Authors: The SOFA Team (see Authors.txt)                                    *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef PYTHONSCRIPTEVENT_H
#define PYTHONSCRIPTEVENT_H

#include "initSofaPython.h"
#include "ScriptEvent.h"

#if defined(_MSC_VER)
// undefine _DEBUG since we want to always link to the release version of
// python and pyconfig.h automatically links debug version if _DEBUG is
// defined.
#ifdef _DEBUG
#define _DEBUG_UNDEFED
#undef _DEBUG
#endif
#endif
//#if defined(__APPLE__) && defined(__MACH__)
//#include <Python/Python.h>
//#else
#include <Python.h>
//#endif
#if defined(_MSC_VER)
// redefine _DEBUG if it was undefed
#ifdef _DEBUG_UNDEFED
#define _DEBUG
#endif
#endif


namespace sofa
{

namespace core
{

namespace objectmodel
{

/**
 * @brief This event notifies about GUI interaction.
 */
class SOFA_SOFAPYTHON_API PythonScriptEvent : public sofa::core::objectmodel::ScriptEvent
{
public:


    /**
     * @brief Constructor.
     */
    PythonScriptEvent(sofa::simulation::Node::SPtr sender, const char* eventName, PyObject* userData);

    /**
     * @brief Destructor.
     */
    virtual ~PythonScriptEvent();

    /**
     * @brief Get the event name
     */
    PyObject* getUserData(void) const {return m_userData;};

    virtual const char* getClassName() const { return "PythonScriptEvent"; }
private:

    PyObject* m_userData;

};

} // namespace objectmodel

} // namespace core

} // namespace sofa

#endif // PYTHONSCRIPTEVENT_H
