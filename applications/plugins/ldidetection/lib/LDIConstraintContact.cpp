//
// C++ Interface: LDIConstraintContact
//
// Description: Contact used by the ContactManager to perform the collision response
//
//
// Author: Francois Faure, Sebastien Barbier, Jeremie Allard, Florent Falipou
//
// Licence: QPL, See LICENCE.txt file that comes with this distribution
// Copyright: See COPYING file that comes with this distribution
//
//
#include "LDIConstraintContact.inl"

namespace sofa
{

namespace component
{

namespace collision
{


SOFA_DECL_CLASS(LDIConstraintContact)

 Creator<Contact::Factory, LDIConstraintContact<TriangleModel, TriangleModel> > LDIConstraintContactClass("LDIConstraint",true);
 

} // namespace collision

} // namespace component

} // namespace sofa

