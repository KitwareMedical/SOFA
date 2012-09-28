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
#ifndef SOFA_COMPONENT_COLLISION_COMPONENTMOUSEINTERACTION_CPP
#define SOFA_COMPONENT_COLLISION_COMPONENTMOUSEINTERACTION_CPP


#include <sofa/component/collision/ComponentMouseInteraction.h>
#include <sofa/component/collision/ComponentMouseInteraction.inl>
#include <sofa/simulation/common/Node.h>
#include <sofa/simulation/common/DeleteVisitor.h>
#include <sofa/defaulttype/Vec3Types.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/component/component.h>
#include <sofa/helper/Factory.inl>


using namespace sofa::simulation;

namespace sofa
{

namespace component
{

namespace collision
{


ComponentMouseInteraction::ComponentMouseInteraction():
    nodeRayPick(NULL),
    mouseInSofa(NULL),
    mouseInteractor(NULL)
{
}

ComponentMouseInteraction::~ComponentMouseInteraction()
{
    if (nodeRayPick)
    {
        nodeRayPick->execute< simulation::DeleteVisitor >(sofa::core::ExecParams::defaultInstance());
        nodeRayPick.reset();
    }
}



void ComponentMouseInteraction::attach(Node* parentNode)
{
    if(parentNode)
    {
        if (!nodeRayPick)
        {
            nodeRayPick = parentNode->createChild("MouseInteraction");
            createInteractionComponents(parentNode,nodeRayPick.get());
            nodeRayPick->detachFromGraph();
        }
        parentNode->addChild(nodeRayPick);
    }
}

void ComponentMouseInteraction::detach()
{
    if (nodeRayPick)
        nodeRayPick->detachFromGraph();
}

void ComponentMouseInteraction::reset()
{
    if (mouseInteractor)
        mouseInteractor->cleanup();
}

#ifndef SOFA_DOUBLE
template class TComponentMouseInteraction<defaulttype::Vec3fTypes>;
template class TComponentMouseInteraction<defaulttype::Rigid3fTypes>;
#endif
#ifndef SOFA_FLOAT
template class TComponentMouseInteraction<defaulttype::Vec3dTypes>;
template class TComponentMouseInteraction<defaulttype::Rigid3dTypes>;
#endif

#ifndef SOFA_DOUBLE
helper::Creator<ComponentMouseInteraction::ComponentMouseInteractionFactory, TComponentMouseInteraction<defaulttype::Vec3fTypes> > ComponentMouseInteractionVec3fClass ("MouseSpringVec3f",true);
helper::Creator<ComponentMouseInteraction::ComponentMouseInteractionFactory, TComponentMouseInteraction<defaulttype::Rigid3fTypes> > ComponentMouseInteractionRigid3fClass ("MouseSpringRigid3f",true);
#endif
#ifndef SOFA_FLOAT
helper::Creator<ComponentMouseInteraction::ComponentMouseInteractionFactory, TComponentMouseInteraction<defaulttype::Vec3dTypes> > ComponentMouseInteractionVec3dClass ("MouseSpringVec3d",true);
helper::Creator<ComponentMouseInteraction::ComponentMouseInteractionFactory, TComponentMouseInteraction<defaulttype::Rigid3dTypes> > ComponentMouseInteractionRigid3dClass ("MouseSpringRigid3d",true);
#endif

}
}
namespace helper
{
template class SOFA_USER_INTERACTION_API Factory<std::string, component::collision::ComponentMouseInteraction, core::objectmodel::BaseContext*>;
}
}
#endif // SOFA_COMPONENT_COLLISION_COMPONENTMOUSEINTERACTION_CPP
