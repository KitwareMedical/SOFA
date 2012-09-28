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
#include <sofa/component/visualmodel/InteractiveCamera.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace visualmodel
{

SOFA_DECL_CLASS(InteractiveCamera)

int InteractiveCameraClass = core::RegisterObject("InteractiveCamera")
        .add< InteractiveCamera >()
        .addAlias("Camera")
        ;


InteractiveCamera::InteractiveCamera()
    :p_zoomSpeed(initData(&p_zoomSpeed, (double) 250.0 , "zoomSpeed", "Zoom Speed"))
    ,p_panSpeed(initData(&p_panSpeed, (double) 0.1 , "panSpeed", "Pan Speed"))
    ,p_pivot(initData(&p_pivot, 0 , "pivot", "Pivot (0 => Camera lookAt, 1 => Camera position, 2 => Scene center, 3 => World center"))
    ,currentMode(InteractiveCamera::NONE_MODE)
    ,isMoving(false)
{
}

InteractiveCamera::~InteractiveCamera()
{
}

void InteractiveCamera::internalUpdate()
{
}

void InteractiveCamera::moveCamera(int x, int y)
{
    float x1, x2, y1, y2;
    Quat newQuat;
    const unsigned int widthViewport = p_widthViewport.getValue();
    const unsigned int heightViewport = p_heightViewport.getValue();

    //std::cout << "widthViewport: " << widthViewport << std::endl;
    //std::cout << "heightViewport: " << heightViewport << std::endl;
    //std::cout << "x: " << x << std::endl;
    //std::cout << "y: " << y << std::endl;

    if (isMoving)
    {
        if (currentMode == TRACKBALL_MODE)
        {
            x1 = (2.0f * widthViewport / 2.0f - widthViewport) / widthViewport;
            y1 = (heightViewport- 2.0f *heightViewport / 2.0f) /heightViewport;
            x2 = (2.0f * (x + (-lastMousePosX + widthViewport / 2.0f)) - widthViewport) / widthViewport;
            y2 = (heightViewport- 2.0f * (y + (-lastMousePosY +heightViewport / 2.0f))) /heightViewport;

            currentTrackball.ComputeQuaternion(x1, y1, x2, y2);
            //fetch rotation
            newQuat = currentTrackball.GetQuaternion();
            Vec3 pivot;
            switch (p_pivot.getValue())
            {
            case CAMERA_LOOKAT_PIVOT:
                pivot = getLookAt();
                break;
            case CAMERA_POSITION_PIVOT:
                pivot = getPosition();
                break;
            case WORLD_CENTER_PIVOT:
                pivot = Vec3(0.0, 0.0, 0.0);
                break;
            case SCENE_CENTER_PIVOT:
            default:
                pivot = sceneCenter;
                break;
            }
            //pivot = p_lookAt.getValue();
            //pivot = (Vec3(x,y, p_distance.getValue()));
            //std::cout << "Pivot : " <<  pivot << std::endl;
            //rotateCameraAroundPoint(newQuat, pivot);
            rotateWorldAroundPoint(newQuat, pivot, this->getOrientation());
        }
        else if (currentMode == ZOOM_MODE)
        {
            double zoomStep = p_zoomSpeed.getValue() *( 0.01*sceneRadius )/heightViewport;
            double zoomDistance = zoomStep * -(y - lastMousePosY);

            Vec3 trans(0.0, 0.0, zoomDistance);
            trans = cameraToWorldTransform(trans);
            translate(trans);
            //translateLookAt(trans);
            Vec3 newLookAt = cameraToWorldCoordinates(Vec3(0,0,-zoomStep));
            if (dot(getLookAt() - getPosition(), newLookAt - getPosition()) < 0)
                translateLookAt(newLookAt - getLookAt());
            getDistance(); // update distance between camera position and lookat
        }
        else if (currentMode == PAN_MODE)
        {
            Vec3 trans(lastMousePosX - x,  y-lastMousePosY, 0.0);
            trans = cameraToWorldTransform(trans)*p_panSpeed.getValue()*( 0.01*sceneRadius ) ;
            translate(trans);
            translateLookAt(trans);
        }
        //must call update afterwards

        lastMousePosX = x;
        lastMousePosY = y;
    }
    else if (currentMode == WHEEL_ZOOM_MODE)
    {
        double zoomStep = p_zoomSpeed.getValue() *( 0.01*sceneRadius )/heightViewport;
        double zoomDistance = zoomStep * -(y*0.5);

        Vec3 trans(0.0, 0.0, zoomDistance);
        trans = cameraToWorldTransform(trans);
        //serr << "old position = " << getPosition() << sendl;
        //serr << "old lookat   = " << getLookAt() << sendl;
        translate(trans);
        //serr << "new position = " << getPosition() << sendl;
        Vec3 newLookAt = cameraToWorldCoordinates(Vec3(0,0,-zoomStep));
        //serr << "new lookat   = " << newLookAt << sendl;
        //translateLookAt(trans);
        if (dot(getLookAt() - getPosition(), newLookAt - getPosition()) < 0)
        {
            //serr << "Moving target point" << sendl;
            translateLookAt(newLookAt - getLookAt());
        }
        getDistance(); // update distance between camera position and lookat

        currentMode = NONE_MODE;
    }

    computeZ();
}


void InteractiveCamera::manageEvent(core::objectmodel::Event* e)
{
    core::objectmodel::MouseEvent* me;
    core::objectmodel::KeypressedEvent* kpe;
    core::objectmodel::KeyreleasedEvent* kre;

    if(p_activated.getValue())
    {
        //Dispatch event
        if ((me = dynamic_cast<core::objectmodel::MouseEvent* > (e)))
            processMouseEvent(me);
        else if ((kpe = dynamic_cast<core::objectmodel::KeypressedEvent* > (e)))
            processKeyPressedEvent(kpe);
        else if ((kre = dynamic_cast<core::objectmodel::KeyreleasedEvent* > (e)))
            processKeyReleasedEvent(kre);

        internalUpdate();
    }
    else
    {
        isMoving = false;
        currentMode = NONE_MODE;
    }
}

void InteractiveCamera::processMouseEvent(core::objectmodel::MouseEvent* me)
{
    int posX = me->getPosX();
    int posY = me->getPosY();
    int wheelDelta = me->getWheelDelta();

    //Mouse Press
    if(me->getState() == core::objectmodel::MouseEvent::LeftPressed)
    {
        isMoving = true;
        currentMode = TRACKBALL_MODE;
        lastMousePosX = posX;
        lastMousePosY = posY;
    }
    else if(me->getState() == core::objectmodel::MouseEvent::RightPressed)
    {
        isMoving = true;
        currentMode = PAN_MODE;
        lastMousePosX = posX;
        lastMousePosY = posY;
    }
    else if(me->getState() == core::objectmodel::MouseEvent::MiddlePressed)
    {
        isMoving = true;
        currentMode = ZOOM_MODE;
        lastMousePosX = posX;
        lastMousePosY = posY;
    }
    else
        //Mouse Move
        if(me->getState() == core::objectmodel::MouseEvent::Move)
        {
            //isMoving = true;
        }
        else
            //Mouse Release
            if(me->getState() == core::objectmodel::MouseEvent::LeftReleased)
            {
                isMoving = false;
                currentMode = NONE_MODE;
            }
            else if(me->getState() == core::objectmodel::MouseEvent::RightReleased)
            {
                isMoving = false;
                currentMode = NONE_MODE;
            }
            else if(me->getState() == core::objectmodel::MouseEvent::MiddleReleased)
            {
                isMoving = false;
                currentMode = NONE_MODE;
            }
    //Mouse Wheel
    if(me->getState() == core::objectmodel::MouseEvent::Wheel)
    {
        isMoving = false;
        currentMode = WHEEL_ZOOM_MODE;
        posX = 0;
        posY = wheelDelta;
    }

    moveCamera(posX, posY);

    p_position.endEdit();

}

void InteractiveCamera::processKeyPressedEvent(core::objectmodel::KeypressedEvent* kpe)
{
    char keyPressed = kpe->getKey();

    switch(keyPressed)
    {
    case 'a':
    case 'A':
    {
        //glPushMatrix();
        //glLoadIdentity();
        //helper::gl::Axis(p_position.getValue(), p_orientation.getValue(), 10.0);
        //glPopMatrix();
        break;
    }
    default:
    {
        break;
    }
    }
}

void InteractiveCamera::processKeyReleasedEvent(core::objectmodel::KeyreleasedEvent* /* kre */)
{

}

} // namespace visualmodel

} // namespace component

} // namespace sofa
