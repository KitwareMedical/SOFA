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
#ifndef SCRIPTCONTROLLER_H
#define SCRIPTCONTROLLER_H

#include <sofa/component/controller/Controller.h>
#include <sofa/core/objectmodel/Context.h>
#include <sofa/core/objectmodel/BaseObjectDescription.h>
#include <sofa/simulation/common/Node.h>
#include "ScriptEvent.h"

namespace sofa
{

namespace component
{

namespace controller
{

class SOFA_SOFAPYTHON_API ScriptController : public Controller
{
public:
    SOFA_CLASS(ScriptController,Controller);

protected:
    ScriptController();


public:
    /// @name control
    ///   Basic control (from BaseObject)
    /// @{

    /// Parse the given description to assign values to this object's fields and potentially other parameters
    virtual void parse ( sofa::core::objectmodel::BaseObjectDescription* arg );

    /// Initialization method called at graph creation and modification, during top-down traversal.
    virtual void init();

    /// Initialization method called at graph creation and modification, during bottom-up traversal.
//    virtual void bwdInit();

    /// Update method called when variables used in precomputation are modified.
//    virtual void reinit();

    /// Save the initial state for later uses in reset()
    virtual void storeResetState();

    /// Reset to initial state
    virtual void reset();

    /// Called just before deleting this object
    /// Any object in the tree bellow this object that are to be removed will be removed only after this call,
    /// so any references this object holds should still be valid.
    virtual void cleanup();

    /// @}



    /// @name Controller notifications
    /// @{

    /**
    * @brief Mouse event callback.
    */
    virtual void onMouseEvent(core::objectmodel::MouseEvent *);

    /**
    * @brief Key Press event callback.
    */
    virtual void onKeyPressedEvent(core::objectmodel::KeypressedEvent *);

    /**
    * @brief Key Release event callback.
    */
    virtual void onKeyReleasedEvent(core::objectmodel::KeyreleasedEvent *);

    /**
    * @brief Begin Animation event callback.
    */
    virtual void onBeginAnimationStep(const double);

    /**
    * @brief End Animation event callback.
    */
    virtual void onEndAnimationStep(const double);

    /// @}

    /**
    * @brief GUI event callback.
    */
    virtual void onGUIEvent(core::objectmodel::GUIEvent *);

    virtual void handleEvent(core::objectmodel::Event *);


protected:

    /// @name Script interface
    ///   Function that need to be implemented for each script language
    /// Typically, all "script_*" functions call the corresponding "*" function of the script, if it exists
    /// @{

    virtual void loadScript() = 0;      // load & bind functions

    virtual void script_onLoaded(sofa::simulation::Node* node) = 0;   // called once, immediately after the script is loaded
    virtual void script_createGraph(sofa::simulation::Node* node) = 0;       // called when the script must create its graph
    virtual void script_initGraph(sofa::simulation::Node* node) = 0;         // called when the script must init its graph, once all the graph has been create

    virtual void script_storeResetState() = 0;
    virtual void script_reset() = 0;

    virtual void script_cleanup() = 0;

    /// keyboard & mouse events
    virtual void script_onKeyPressed(const char c) = 0;
    virtual void script_onKeyReleased(const char c) = 0;

    virtual void script_onMouseButtonLeft(const int posX,const int posY,const bool pressed) = 0;
    virtual void script_onMouseButtonRight(const int posX,const int posY,const bool pressed) = 0;
    virtual void script_onMouseButtonMiddle(const int posX,const int posY,const bool pressed) = 0;
    virtual void script_onMouseWheel(const int posX,const int posY,const int delta) = 0;

    /// called once per frame
    virtual void script_onBeginAnimationStep(const double dt) = 0;
    virtual void script_onEndAnimationStep(const double dt) = 0;

    /// GUI interaction
    virtual void script_onGUIEvent(const char* controlID, const char* valueName, const char* value) = 0;

    /// Script events; user data is implementation-dependant
    virtual void script_onScriptEvent(core::objectmodel::ScriptEvent *) = 0;

    /// @}


};


} // namespace controller

} // namespace component

} // namespace sofa

#endif // SCRIPTCONTROLLER_H
