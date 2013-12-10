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
#ifndef RECORDEDCAMERA_H
#define RECORDEDCAMERA_H

#include <sofa/component/visualmodel/BaseCamera.h>

namespace sofa
{

namespace component
{

namespace visualmodel
{

class SOFA_BASE_VISUAL_API RecordedCamera : public BaseCamera
{
public:
    SOFA_CLASS(RecordedCamera, BaseCamera);

    typedef BaseCamera::Vec3 Vec3;
    typedef BaseCamera::Quat Quat;
protected:
    RecordedCamera();
    virtual ~RecordedCamera() {}
public:
    virtual void init();

    virtual void reinit();

    virtual void reset();

    virtual void handleEvent(sofa::core::objectmodel::Event *);

    //virtual void rotateWorldAroundPoint(Quat &rotation, const Vec3 &point);

    enum  { TRACKBALL_MODE, PAN_MODE, ZOOM_MODE, WHEEL_ZOOM_MODE, NONE_MODE };
    enum  { SCENE_CENTER_PIVOT = 0, WORLD_CENTER_PIVOT = 1};

    Data<double> p_zoomSpeed;
    Data<double> p_panSpeed;
    Data<int> p_pivot;

    void draw(const core::visual::VisualParams* vparams);

private:
    int currentMode;
    bool isMoving;
    int lastMousePosX, lastMousePosY;
    helper::gl::Trackball currentTrackball;

    void moveCamera_rotation();
    void moveCamera_translation();
    void moveCamera_navigation();

    // Kepp functions for mouse interaction (TODO: removed them and allow interactive and recorded camera in same scene)
    void moveCamera_mouse(int x, int y);
    void manageEvent(core::objectmodel::Event* e);
    void processMouseEvent(core::objectmodel::MouseEvent* me);

    void configureRotation();
    void configureTranslation();
    void configureNavigation();
    void initializeViewUp();
    void drawRotation();

public:
    Data<SReal> m_startTime;
    Data<SReal> m_endTime;

    Data <bool> m_rotationMode;
    Data <bool> m_translationMode;
    Data <bool> m_navigationMode;
    Data <SReal> m_rotationSpeed;
    Data <Vec3> m_rotationCenter;
    Data <Vec3> m_rotationStartPoint;
    Data <Vec3> m_rotationLookAt;
    Data <Vec3> m_rotationAxis;
    Data <Vec3> m_cameraUp;

    Data <bool> p_drawRotation;
    Data <bool> p_drawTranslation;

    Data <sofa::helper::vector<Vec3> > m_translationPositions;
    Data <sofa::helper::vector<Quat> > m_translationOrientations;

protected:
    double m_nextStep;
    double m_angleStep;
    //double m_initAngle;
    //double m_radius;
    bool firstIterationforRotation;
    bool firstIterationforTranslation;
    bool firstIterationforNavigation;
  
    sofa::helper::vector <Vec3> m_rotationPoints;
};

} // namespace visualmodel

} // namespace component

} // namespace sofa

#endif // RecordedCamera_H
