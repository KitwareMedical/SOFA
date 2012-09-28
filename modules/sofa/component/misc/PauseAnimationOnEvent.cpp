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
#include <sofa/component/misc/PauseAnimationOnEvent.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/simulation/common/PauseEvent.h>
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace misc
{

PauseAnimationOnEvent::PauseAnimationOnEvent() : paused(false)
{
}


PauseAnimationOnEvent::~PauseAnimationOnEvent()
{

}

void PauseAnimationOnEvent::init()
{
    PauseAnimation::init();
    this->f_listening.setValue(true);
}

bool PauseAnimationOnEvent::isPaused()
{
    return paused;
}

void PauseAnimationOnEvent::handleEvent(sofa::core::objectmodel::Event* event)
{
    if (dynamic_cast<sofa::simulation::PauseEvent*>(event))
    {
        paused = true;
        pause();
    }
}

SOFA_DECL_CLASS(PauseAnimationOnEvent)

int PauseAnimationOnEventClass = core::RegisterObject("PauseAnimationOnEvent")
        .add< PauseAnimationOnEvent >();




} // namespace misc

} // namespace component

} // namespace sofa
