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
#ifndef SOFA_COMPONENT_COLLISION_CompliantAttachPerformer_CPP
#define SOFA_COMPONENT_COLLISION_CompliantAttachPerformer_CPP

#include "CompliantAttachPerformer.inl"
#include <sofa/defaulttype/Vec3Types.h>
#include <sofa/helper/Factory.inl>
#include <sofa/gui/PickHandler.h>
#include <sofa/component/collision/ComponentMouseInteraction.h>

namespace sofa
{
using namespace component::collision;

#ifdef WIN32
#ifndef SOFA_DOUBLE
#ifdef SOFA_DEV
helper::Creator<InteractionPerformer::InteractionPerformerFactory, CompliantAttachPerformer<defaulttype::Vec3fTypes> >  CompliantAttachPerformerVec3fClass("CompliantAttach",true);
#endif
#endif
#ifndef SOFA_FLOAT
#ifdef SOFA_DEV
helper::Creator<InteractionPerformer::InteractionPerformerFactory, CompliantAttachPerformer<defaulttype::Vec3dTypes> >  CompliantAttachPerformerVec3dClass("CompliantAttach",true);
#endif
#endif
#endif

namespace gui
{
//*******************************************************************************************
void CompliantAttachOperation::start()
{
    //Creation
    performer=component::collision::InteractionPerformer::InteractionPerformerFactory::getInstance()->createObject("CompliantAttach", pickHandle->getInteraction()->mouseInteractor.get());
    pickHandle->getInteraction()->mouseInteractor->addInteractionPerformer(performer);
    //Start
    performer->start();
}

void CompliantAttachOperation::execution()
{
    //do nothing
}

void CompliantAttachOperation::end()
{
    pickHandle->getInteraction()->mouseInteractor->removeInteractionPerformer(performer);
    delete performer; performer=0;
}

void CompliantAttachOperation::endOperation()
{
    pickHandle->getInteraction()->mouseInteractor->removeInteractionPerformer(performer);
}


}// gui


namespace component
{

namespace collision
{

#ifndef SOFA_DOUBLE
template class SOFA_Compliant_gui_API  CompliantAttachPerformer<defaulttype::Vec3fTypes>;
#endif
#ifndef SOFA_FLOAT
template class SOFA_Compliant_gui_API  CompliantAttachPerformer<defaulttype::Vec3dTypes>;
#endif


#ifndef SOFA_DOUBLE
helper::Creator<InteractionPerformer::InteractionPerformerFactory, CompliantAttachPerformer<defaulttype::Vec3fTypes> >  CompliantAttachPerformerVec3fClass("CompliantAttach",true);
#endif
#ifndef SOFA_FLOAT
helper::Creator<InteractionPerformer::InteractionPerformerFactory, CompliantAttachPerformer<defaulttype::Vec3dTypes> >  CompliantAttachPerformerVec3dClass("CompliantAttach",true);
#endif
}
}
}
#endif
