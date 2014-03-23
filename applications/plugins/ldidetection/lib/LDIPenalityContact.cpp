//
// C++ Interface: LDIPenalityContact
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
//#include <sofa/component/collision/LDIPenalityContact.inl>
#include "LDIPenalityContact.inl"

namespace sofa
{

namespace component
{

namespace collision
{


SOFA_DECL_CLASS(LDIPenalityContact)

 Creator<Contact::Factory, LDIPenalityContact<TriangleModel, TriangleModel> > LDIPenalityContactClass("LDI",true);
 

} // namespace collision

} // namespace component

} // namespace sofa

