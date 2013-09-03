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
#ifndef SOFA_COMPONENT_CONTROLLER_OMNIEMU_H
#define SOFA_COMPONENT_CONTROLLER_OMNIEMU_H

//Sensable include
//#include <HD/hd.h>
//#include <HDU/hdu.h>
//#include <HDU/hduError.h>
//#include <HDU/hduVector.h>
#include <sofa/helper/LCPcalc.h>
#include <sofa/defaulttype/SolidTypes.h>

#include <sofa/core/behavior/BaseController.h>
#include <sofa/component/visualmodel/OglModel.h>
#include <sofa/component/controller/Controller.h>
#include <sofa/core/behavior/MechanicalState.h>

#include <sofa/helper/system/thread/CTime.h>
#include <pthread.h>

namespace sofa
{
namespace simulation { class Node; }

namespace component
{
namespace visualModel { class OglModel; }

namespace controller
{

class ForceFeedback;


using namespace sofa::defaulttype;
using namespace helper::system::thread;
using core::objectmodel::Data;

/** Holds data retrieved from HDAPI. */
typedef struct
{
    // changement unsigned int -> HDD
    //HHD id;
    int nupdates;
    int m_buttonState;					/* Has the device button has been pressed. */
    //hduVector3Dd m_devicePosition;	/* Current device coordinates. */
    //HDErrorInfo m_error;
    Vec3d pos;
    Quat quat;
    bool ready;
    bool stop;
} DeviceData;

typedef struct
{
    vector<ForceFeedback*> forceFeedbacks;
    //ForceFeedback* forceFeedback;
    // changement ajout
    int forceFeedbackIndice;
    simulation::Node *context;

    sofa::defaulttype::SolidTypes<double>::Transform endOmni_H_virtualTool;
    //Transform baseOmni_H_endOmni;
    sofa::defaulttype::SolidTypes<double>::Transform world_H_baseOmni;
    double forceScale;
    double scale;
    bool permanent_feedback;

    // API OMNI //
    DeviceData servoDeviceData;  // for the haptic loop
    DeviceData deviceData;		 // for the simulation loop

} OmniData;

/**
* Omni driver
*/
//changement NewOmni -> Omni
class OmniDriverEmu : public Controller
{

public:
    typedef Rigid3dTypes::Coord Coord;
    typedef Rigid3dTypes::VecCoord VecCoord;

    // changement NewOmni -> omni
    SOFA_CLASS(OmniDriverEmu, Controller);
    Data<double> forceScale;
    Data<double> scale;
    Data<Vec3d> positionBase;
    Data<Quat> orientationBase;
    Data<Vec3d> positionTool;
    Data<Quat> orientationTool;
    Data<bool> permanent;
    Data<bool> omniVisu;
    Data<int> simuFreq;
    Data<bool> simulateTranslation;
    Data<bool> toolSelector;
    Data<int> toolCount;

    OmniData	data;

    // changement NewOmni -> omni
    OmniDriverEmu();
    virtual ~OmniDriverEmu();

    virtual void init();
    virtual void bwdInit();
    virtual void reset();
    void reinit();

    int initDevice(OmniData& data);

    void cleanup();
    virtual void draw();

    //ajout
    void setForceFeedbacks(vector<ForceFeedback*> ffs);

    void onKeyPressedEvent(core::objectmodel::KeypressedEvent *);
    void onKeyReleasedEvent(core::objectmodel::KeyreleasedEvent *);

    void setDataValue();
    void reinitVisual();

    void setOmniSimThreadCreated(bool b) { omniSimThreadCreated = b;};

    bool afterFirstStep;
    SolidTypes<double>::Transform prevPosition;

    //neede for "omni simulation"
    CTime *thTimer;
    pthread_t hapSimuThread;
    double lastStep;
    bool executeAsynchro;
    Data<VecCoord> trajPts;
    Data<helper::vector<double> > trajTim;

    int getCurrentToolIndex() { return currentToolIndex;};
    void handleEvent(core::objectmodel::Event *);

private:

    void copyDeviceDataCallback(OmniData *pUserData);
    void stopCallback(OmniData *pUserData);
    sofa::component::visualmodel::OglModel::SPtr visu_base, visu_end;
    bool noDevice;

    bool moveOmniBase;
    Vec3d positionBase_buf;

    //ajout
    core::behavior::MechanicalState<Rigid3dTypes> *mState; ///< Controlled MechanicalState.

    bool omniSimThreadCreated;

    //ajout
    int currentToolIndex;
    //ajout
    bool isToolControlled;


};


} // namespace controller

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_CONTROLLER_OMNIEMU_H
