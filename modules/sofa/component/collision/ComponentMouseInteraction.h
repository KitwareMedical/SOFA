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
#ifndef SOFA_COMPONENT_COLLISION_COMPONENTMOUSEINTERACTION_H
#define SOFA_COMPONENT_COLLISION_COMPONENTMOUSEINTERACTION_H

#include <sofa/simulation/common/Simulation.h>
#include <sofa/simulation/common/DeleteVisitor.h>
#include <sofa/component/collision/MouseInteractor.h>
#include <sofa/component/mapping/IdentityMapping.h>
#include <sofa/component/component.h>
#include <sofa/core/Mapping.h>

namespace sofa
{

namespace simulation
{
class Node;
}
namespace component
{

namespace collision
{


class SOFA_USER_INTERACTION_API ComponentMouseInteraction
{
public:
    ComponentMouseInteraction();

    virtual ~ComponentMouseInteraction();

    virtual void createInteractionComponents(sofa::simulation::Node* parent,
            sofa::simulation::Node* current) = 0;

    void attach(simulation::Node* parentNode);

    void detach();

    void reset();

    virtual bool isCompatible( core::objectmodel::BaseContext *)const=0;

    typedef helper::Factory<std::string, ComponentMouseInteraction, core::objectmodel::BaseContext*> ComponentMouseInteractionFactory;

    template <class RealObject>
    static RealObject* create( RealObject*, core::objectmodel::BaseContext* /* context */)
    {
        return new RealObject;
    }

    //Components
    simulation::Node::SPtr nodeRayPick;
    sofa::core::behavior::BaseMechanicalState::SPtr mouseInSofa;
    sofa::component::collision::BaseMouseInteractor::SPtr mouseInteractor;
};



template <class DataTypes>
class TComponentMouseInteraction : public ComponentMouseInteraction
{
    typedef sofa::component::container::MechanicalObject< defaulttype::Vec3Types > MousePosition;
    typedef sofa::component::container::MechanicalObject< DataTypes > MouseContainer;
    typedef sofa::component::collision::MouseInteractor< DataTypes > Interactor;
    typedef sofa::component::mapping::IdentityMapping< defaulttype::Vec3Types, DataTypes > IdentityMechanicalMapping;
    typedef typename sofa::core::Mapping< defaulttype::Vec3Types, DataTypes >::SPtr MouseMapping;

public:


    void createInteractionComponents(sofa::simulation::Node* parent, sofa::simulation::Node* current);

    bool  isCompatible( core::objectmodel::BaseContext *context) const;
protected :
    MouseMapping mouseMapping;
};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_COLLISION_COMPONENTMOUSEINTERACTION_CPP)
#ifndef SOFA_DOUBLE
extern template class SOFA_USER_INTERACTION_API TComponentMouseInteraction<defaulttype::Vec3fTypes>;
extern template class SOFA_USER_INTERACTION_API TComponentMouseInteraction<defaulttype::Vec2fTypes>;
extern template class SOFA_USER_INTERACTION_API TComponentMouseInteraction<defaulttype::Rigid3fTypes>;
#endif
#ifndef SOFA_FLOAT
extern template class SOFA_USER_INTERACTION_API TComponentMouseInteraction<defaulttype::Vec2dTypes>;
extern template class SOFA_USER_INTERACTION_API TComponentMouseInteraction<defaulttype::Vec3dTypes>;
extern template class SOFA_USER_INTERACTION_API TComponentMouseInteraction<defaulttype::Rigid3dTypes>;
#endif
#endif
}
}

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_COLLISION_COMPONENTMOUSEINTERACTION_CPP)
namespace helper
{
extern template class SOFA_USER_INTERACTION_API Factory<std::string, component::collision::ComponentMouseInteraction, core::objectmodel::BaseContext*>;
}
#endif

}

#endif
