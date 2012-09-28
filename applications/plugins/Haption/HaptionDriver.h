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
#ifndef SOFA_COMPONENT_ODESOLVER_NEWOMNISOLVER_H
#define SOFA_COMPONENT_ODESOLVER_NEWOMNISOLVER_H

//Haption include
#include <sofa/helper/LCPcalc.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/helper/Quater.h>
#include <sofa/core/behavior/BaseController.h>
#include <sofa/component/visualmodel/OglModel.h>
#include <sofa/component/mapping/RigidMapping.h>
#include <sofa/component/controller/Controller.h>
#include <sofa/simulation/common/Node.h>
#include <sofa/simulation/common/Simulation.h>
#include <sofa/component/container/MechanicalObject.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/component/controller/MechanicalStateForceFeedback.h>
#include <sofa/component/controller/NullForceFeedbackT.h>
#include <sofa/simulation/common/AnimateBeginEvent.h>
#include <sofa/simulation/common/AnimateEndEvent.h>
#include <cstring>
#include <sofa/core/objectmodel/KeypressedEvent.h>
#include <sofa/core/objectmodel/KeyreleasedEvent.h>
//#include <sofa/core/objectmodel/MouseEvent.h>
#include <math.h>
#include <sofa/simulation/tree/GNode.h>
#include "virtuoseAPI.h"


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
using core::objectmodel::Data;

typedef struct
{
    simulation::Node *node;
    sofa::component::visualmodel::OglModel *visu;
    sofa::component::mapping::RigidMapping< Rigid3dTypes , ExtVec3fTypes  > *mapping;
} VisualComponent;

typedef struct
{
    VirtContext m_virtContext;
    MechanicalStateForceFeedback<Rigid3dTypes>* forceFeedback;
    float scale;
    float torqueScale;
    float forceScale;
} HaptionData;


class HaptionDriver : public Controller
{

public:
    SOFA_CLASS(HaptionDriver, Controller);
    typedef RigidTypes::VecCoord VecCoord;

    Data<double> scale;
    Data<bool> state_button;
    Data<bool> haptionVisu;
    Data<VecCoord> posBase;
    Data<double> torqueScale;
    Data<double> forceScale;
    Data< std::string > ip_haption;

    HaptionDriver();
    virtual ~HaptionDriver();

    virtual void init();
    virtual void reinit();
    virtual void bwdInit();
    virtual void reset();
    virtual void handleEvent(core::objectmodel::Event *);
    void onKeyPressedEvent(core::objectmodel::KeypressedEvent *);
    void onKeyReleasedEvent(core::objectmodel::KeyreleasedEvent *);
    void onAnimateBeginEvent();
    int initDevice(char* ip);
    void closeDevice();
    static void haptic_callback(VirtContext, void *);

    void setForceFeedback(MechanicalStateForceFeedback<Rigid3dTypes>* ff);

private:

    HaptionData myData;
    VirtIndexingType m_indexingMode;
    VirtCommandType m_typeCommand;
    float m_speedFactor;
    float m_forceFactor;
    float haptic_time_step;
    int connection_device;
    sofa::component::container::MechanicalObject<sofa::defaulttype::Rigid3dTypes> *rigidDOF;
    bool initCallback;
    simulation::Node *nodeHaptionVisual;
    sofa::component::container::MechanicalObject<sofa::defaulttype::Rigid3dTypes> *visualHaptionDOF;
    simulation::Node *nodeAxesVisual;
    sofa::component::container::MechanicalObject<sofa::defaulttype::Rigid3dTypes> *visualAxesDOF;
    VisualComponent visualNode[5];

    float oldScale;
    bool changeScale;
    bool visuAxes;
    bool modX,modY,modZ,modS;
    bool visuActif;
};

} // namespace controller

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_ODESOLVER_NEWOMNISOLVER_H
