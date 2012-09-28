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
#ifndef SOFA_COMPONENT_VISUALMODEL_INTERACTIVECAMERA_H
#define SOFA_COMPONENT_VISUALMODEL_INTERACTIVECAMERA_H

#include <sofa/component/visualmodel/BaseCamera.h>
#include <sofa/component/component.h>
namespace sofa
{

namespace component
{

namespace visualmodel
{

class SOFA_BASE_VISUAL_API InteractiveCamera : public BaseCamera
{
public:
    SOFA_CLASS(InteractiveCamera, BaseCamera);

    enum  { TRACKBALL_MODE, PAN_MODE, ZOOM_MODE, WHEEL_ZOOM_MODE, NONE_MODE };
    enum  { CAMERA_LOOKAT_PIVOT = 0, CAMERA_POSITION_PIVOT = 1, SCENE_CENTER_PIVOT = 2, WORLD_CENTER_PIVOT = 3};

    Data<double> p_zoomSpeed;
    Data<double> p_panSpeed;
    Data<int> p_pivot;
protected:
    InteractiveCamera();
    virtual ~InteractiveCamera();
public:
private:
    int currentMode;
    bool isMoving;
    int lastMousePosX, lastMousePosY;
    helper::gl::Trackball currentTrackball;

    void internalUpdate();
    void moveCamera(int x, int y);
    void manageEvent(core::objectmodel::Event* e);
    void processMouseEvent(core::objectmodel::MouseEvent* me);
    void processKeyPressedEvent(core::objectmodel::KeypressedEvent* kpe);
    void processKeyReleasedEvent(core::objectmodel::KeyreleasedEvent* kre);
};

} // namespace visualmodel

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_VISUALMODEL_INTERACTIVECAMERA_H
