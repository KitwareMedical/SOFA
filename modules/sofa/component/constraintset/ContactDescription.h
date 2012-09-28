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
#ifndef SOFA_COMPONENT_CONSTRAINTSET_CONTACTDESCRIPTION_H
#define SOFA_COMPONENT_CONSTRAINTSET_CONTACTDESCRIPTION_H

#include <sofa/core/behavior/BaseLMConstraint.h>
#include <sofa/helper/vector.h>

namespace sofa
{

namespace component
{

namespace constraintset
{

enum ContactState {VANISHING, STICKING, SLIDING};

/// ContactDescription is a class providing precise information about the state of a contact
/// A contact can be in 3 different states: vanishing, sticking, sliding
/// with this states, we can solve precisely a contact, knowning its state, and the direction of the response
struct ContactDescription
{
    /// State of the contact
    /// @see ContactState
    ContactState state;
    /// When the contact is sliding, we have to know the direction of the reponse force.
    /// coeff is the linear combination of the normal, and two tangent constraint directions.
    SReal coeff[3];
};


/// Class handler to make the kink between constraint groups (a set of equations related to contact) and a description of the state of the contact
class ContactDescriptionHandler
{
    typedef  std::map< const core::behavior::ConstraintGroup*, ContactDescription> InternalData;
public:
    const ContactDescription& getContactDescription( const core::behavior::ConstraintGroup* contact) const
    {
        InternalData::const_iterator it = infos.find(contact);
        assert (it != infos.end());
        return it->second;
    };
    ContactDescription& getContactDescription(const core::behavior::ConstraintGroup* contact)
    {
        return infos[contact];
    };
protected:
    InternalData infos;
};

} // namespace constraintset

} // namespace component

} // namespace sofa

#endif
