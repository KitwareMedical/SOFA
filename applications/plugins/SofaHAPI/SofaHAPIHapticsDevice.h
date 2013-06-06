/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 4      *
*                (c) 2006-2009 MGH, INRIA, USTL, UJF, CNRS                    *
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
#ifndef SOFAHAPI_SOFAHAPIHAPTICSDEVICE_H
#define SOFAHAPI_SOFAHAPIHAPTICSDEVICE_H

#include "SofaHAPI.h"

//HAPI include
#include <H3DUtil/AutoRef.h>
#include <HAPI/HAPIHapticsDevice.h>
#include <HAPI/HapticSpring.h>

#include "SofaHAPIForceFeedbackEffect.h"

#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/behavior/BaseController.h>
#include <sofa/component/controller/Controller.h>
#include <sofa/component/controller/ForceFeedback.h>

namespace SofaHAPI
{

using sofa::helper::vector;
using sofa::defaulttype::Vec3d;
using sofa::defaulttype::Quat;
typedef sofa::defaulttype::SolidTypes<double>::Transform Transform;
using sofa::core::objectmodel::Data;
using sofa::core::objectmodel::BaseLink;
using sofa::core::objectmodel::MultiLink;
using sofa::core::objectmodel::KeypressedEvent;
using sofa::core::objectmodel::KeyreleasedEvent;
using sofa::core::behavior::MechanicalState;
using sofa::component::controller::Controller;
using sofa::component::controller::ForceFeedback;

/**
* HAPI Haptics Device
*/
class SofaHAPIHapticsDevice : public Controller
{
public:
    SOFA_CLASS(SofaHAPIHapticsDevice, Controller);
    Data<double> scale;
    Data<double> forceScale;
    Data<Vec3d> positionBase;
    Data<Quat> orientationBase;
    Data<Vec3d> positionTool;
    Data<Quat> orientationTool;
    Data<bool> permanent;
    Data<bool> toolSelector;
    Data<int> toolCount;
    Data<int> toolIndex;
    Data<double> toolTransitionSpringStiffness;
    Data<std::string> driverName;
    Data<bool> drawDevice;
    Data<float> drawHandleSize;
    Data<float> drawForceScale;

    //SofaHAPIHapticsDeviceData data;
    ForceFeedbackTransform data;

    SofaHAPIHapticsDevice();
    virtual ~SofaHAPIHapticsDevice();

    virtual void init();
    virtual void bwdInit();
    virtual void reset();
    void reinit();

    bool initDevice();
    void releaseDevice();

    void cleanup();
    virtual void draw(const sofa::core::visual::VisualParams* vparams);

    void setForceFeedbacks(vector<ForceFeedback*> ffs);

    void onKeyPressedEvent(KeypressedEvent *);
    void onKeyReleasedEvent(KeyreleasedEvent *);
    void onBeginAnimationStep(const double /*dt*/);
    void onEndAnimationStep(const double /*dt*/);

    void setDataValue();

protected:

    std::auto_ptr<HAPI::HAPIHapticsDevice> device;
    MultiLink<SofaHAPIHapticsDevice, SofaHAPIForceFeedbackEffect, BaseLink::FLAG_STRONGLINK> feedbackEffects;
    H3DUtil::AutoRef<HAPI::HapticSpring> transitionEffect;

    sofa::core::behavior::MechanicalState<Rigid3dTypes> *mState; ///< Controlled MechanicalState.

    bool isToolControlled;

    int fakeButtonState;
    int lastButtonState;
    Transform lastToolPosition;

    void setToolFeedback(int indice, bool enable = true, bool transfer = true);

    void sendHapticDeviceEvent();
};

} // namespace SofaHAPI

#endif // SOFAHAPI_SOFAHAPIHAPTICSDEVICE_H
