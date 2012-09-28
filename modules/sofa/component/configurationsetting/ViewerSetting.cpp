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

#include <sofa/component/configurationsetting/ViewerSetting.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace configurationsetting
{

SOFA_DECL_CLASS(ViewerSetting)
int ViewerSettingClass = core::RegisterObject("Configuration for the Viewer of your application")
        .add< ViewerSetting >()
        .addAlias("Viewer")
        ;

ViewerSetting::ViewerSetting():
    resolution(initData(&resolution, defaulttype::Vec<2,int>(800,600), "resolution", "resolution of the Viewer"))
    ,fullscreen(initData(&fullscreen, false, "fullscreen", "Fullscreen mode"))
    ,cameraMode(initData(&cameraMode, "cameraMode", "Camera mode"))
    ,objectPickingMethod(initData(&objectPickingMethod, "objectPickingMethod","The method used to pick objects"))
{
    sofa::helper::OptionsGroup mode(2,"Perspective","Orthographic");
    cameraMode.setValue(mode);
    sofa::helper::OptionsGroup pickmethod(2,"Ray casting","Selection buffer");
    objectPickingMethod.setValue(pickmethod);
}

}

}

}
