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
#ifndef SOFA_GUI_BASEGUI_H
#define SOFA_GUI_BASEGUI_H

#include "SofaGUI.h"
#include <sofa/simulation/common/Node.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/component/configurationsetting/ViewerSetting.h>
#include <sofa/component/configurationsetting/MouseButtonSetting.h>

#include <list>


namespace sofa
{

namespace gui
{

class BaseViewer;

class SOFA_SOFAGUI_API BaseGUI
{

public:

    /// @name methods each GUI must implement
    /// @{
    /// Start the GUI loop
    virtual int mainLoop()=0;
    /// Update the GUI
    virtual void redraw()=0;
    /// Close the GUI
    virtual int closeGUI()=0;
    /// Register the scene in our GUI
    virtual void setScene(sofa::simulation::Node::SPtr groot, const char* filename=NULL, bool temporaryFile=false)=0;
    /// Get the rootNode of the sofa scene
    virtual sofa::simulation::Node* currentSimulation() = 0;
    /// @}

    /// Use a component setting to configure our GUI
    virtual void configureGUI(sofa::simulation::Node::SPtr groot);

    /// @name methods to configure the GUI
    /// @{
    virtual void setDumpState(bool) {}
    virtual void setLogTime(bool) {}
    virtual void setExportState(bool) {}
#ifdef SOFA_DUMP_VISITOR_INFO
    virtual void setTraceVisitors(bool) {}
#endif
    virtual void setRecordPath(const std::string & /*path*/) {}
    virtual void setGnuplotPath(const std::string & /*path*/) {}

    virtual void initViewer(BaseViewer* /*viewer*/) {}
    virtual void setViewerConfiguration(sofa::component::configurationsetting::ViewerSetting* /*viewerConf*/) {}
    virtual void setViewerResolution(int /* width */, int /* height */) {}
    virtual void setFullScreen() {}
    virtual void setBackgroundColor(const defaulttype::Vector3& /*color*/) {}
    virtual void setBackgroundImage(const std::string& /*image*/) {}

    virtual BaseViewer* getViewer() {return NULL;}
    virtual void registerViewer(BaseViewer* /*viewer*/) {}

    virtual void setMouseButtonConfiguration(sofa::component::configurationsetting::MouseButtonSetting* /*button*/) {}
    /// @}

    /// @name methods to communicate with the GUI
    /// @{
    virtual void sendMessage(const std::string & /*msgType*/,const std::string & /*msgValue*/) {}
    /// @}

    void exportGnuplot(sofa::simulation::Node* node, std::string gnuplot_directory="");

    static std::string& GetGUIName() { return mGuiName; }

    static const char* GetProgramName() { return mProgramName; }
    static void SetProgramName(const char* argv0) { if(argv0) mProgramName = argv0;}

protected:
    BaseGUI();
    /// The destructor should not be called directly. Use the closeGUI() method instead.
    virtual ~BaseGUI();

    static std::string mGuiName; // would like to make it const but not possible with the current implementation of RealGUI...
    static const char* mProgramName;
};

////// TO declare into BaseViewer
//setScene();
//resetView();
//setBackgroundColour(...)
//setBackgroundImage(...)
//setScene()
//getSceneFileName()

} // namespace gui

} // namespace sofa

#endif
