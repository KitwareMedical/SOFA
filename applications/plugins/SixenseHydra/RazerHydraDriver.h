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

//Sixense includes
#include <sixense.h>
#include <sixense_math.hpp>
#ifdef WIN32
#include <sixense_utils/mouse_pointer.hpp>
#endif
#include <sixense_utils/derivatives.hpp>
#include <sixense_utils/button_states.hpp>
#include <sixense_utils/event_triggers.hpp>
#include <sixense_utils/controller_manager/controller_manager.hpp>


#include <sofa/core/behavior/BaseController.h>
#include <sofa/component/visualmodel/OglModel.h>
#include <sofa/component/controller/Controller.h>
#include <sofa/core/visual/VisualParams.h>

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/ObjectFactory.h>

#include <sofa/helper/system/gl.h>
#include <sofa/helper/gl/BasicShapes.h>
#include <sofa/helper/gl/glText.inl>

#include <deque>


namespace sofa
{
namespace simulation { class Node; }

namespace component
{

namespace controller
{

using namespace sofa::defaulttype;

class RazerHydraDriver : public Controller
{

public:
    SOFA_CLASS(RazerHydraDriver, Controller);	
	Data<double> scale;
	Data<Vec3d> positionBase;
    Data<Quat> orientationBase;
    Data<Vec3d> positionFirstTool, positionSecondTool;
    Data<Quat> orientationFirstTool, orientationSecondTool;
	Data< bool > triggerJustPressedFirstTool, triggerJustPressedSecondTool;
	Data< float > triggerValueFirstTool, triggerValueSecondTool;
	Data< bool > useBothTools;
	Data< bool > displayTools;

    RazerHydraDriver();
    virtual ~RazerHydraDriver();

    void init();
    void bwdInit();
    void reset();
    void reinit();

    void cleanup();
	void draw(const sofa::core::visual::VisualParams* vparams);

private:
	int first_controller_index, second_controller_index;
	sixenseAllControllerData acd;

    void handleEvent(core::objectmodel::Event *);
	void check_for_button_presses( sixenseAllControllerData *acd );
	void applyRotation (Vec3d* positionToRotate, Quat* orientationToRotate);

};

} // namespace controller

} // namespace component

} // namespace sofa
