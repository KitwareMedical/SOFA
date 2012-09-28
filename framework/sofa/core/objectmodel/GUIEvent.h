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
*                              SOFA :: Framework                              *
*                                                                             *
* Authors: The SOFA Team (see Authors.txt)                                    *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_CORE_OBJECTMODEL_GUIEVENT_H
#define SOFA_CORE_OBJECTMODEL_GUIEVENT_H

#include <sofa/core/objectmodel/Event.h>
#include <string>


namespace sofa
{

namespace core
{

namespace objectmodel
{

/**
 * @brief This event notifies about GUI interaction.
 */
class SOFA_CORE_API GUIEvent : public sofa::core::objectmodel::Event
{
public:


    /**
     * @brief Constructor.
     */
    GUIEvent(const char* controlID, const char* valueName, const char* value);

    /**
     * @brief Destructor.
     */
    virtual ~GUIEvent();

    /**
     * @brief Get the emitter control ID
     */
    const std::string getControlID(void) const {return m_controlID;};

    /**
     * @brief Get the value name
     */
    const std::string getValueName(void) const {return m_valueName;};

    /**
     * @brief Get the value
     */
    const std::string getValue(void) const {return m_value;};


    virtual const char* getClassName() const { return "GUIEvent"; }
private:

    std::string     m_controlID;
    std::string     m_valueName;
    std::string     m_value;

};

} // namespace objectmodel

} // namespace core

} // namespace sofa

#endif // SOFA_CORE_OBJECTMODEL_GUIEVENT_H
