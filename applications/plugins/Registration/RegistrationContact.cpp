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
#include "RegistrationContact.inl"
#include <sofa/component/collision/BarycentricContactMapper.inl>

#include <sofa/component/collision/IdentityContactMapper.inl>
#include <sofa/component/collision/RigidContactMapper.inl>
#include "RegistrationContactForceField.inl"

namespace sofa
{

namespace component
{

namespace collision
{

using namespace defaulttype;
using simulation::Node;

SOFA_DECL_CLASS(RegistrationContact)

Creator<Contact::Factory, RegistrationContact<SphereModel, SphereModel> > SphereSphereRegistrationContactClass("registration",true);
Creator<Contact::Factory, RegistrationContact<SphereModel, PointModel> > SpherePointRegistrationContactClass("registration",true);
//Creator<Contact::Factory, RegistrationContact<SphereTreeModel, SphereTreeModel> > SphereTreeSphereTreeRegistrationContactClass("registration", true);
//Creator<Contact::Factory, RegistrationContact<SphereTreeModel, TriangleModel> > SphereTreeTriangleContactClass("registration", true);
Creator<Contact::Factory, RegistrationContact<PointModel, PointModel> > PointPointRegistrationContactClass("registration",true);
Creator<Contact::Factory, RegistrationContact<LineModel, PointModel> > LinePointRegistrationContactClass("registration",true);
Creator<Contact::Factory, RegistrationContact<LineModel, LineModel> > LineLineRegistrationContactClass("registration",true);
Creator<Contact::Factory, RegistrationContact<LineModel, SphereModel> > LineSphereRegistrationContactClass("registration",true);
Creator<Contact::Factory, RegistrationContact<TriangleModel, SphereModel> > TriangleSphereRegistrationContactClass("registration",true);
Creator<Contact::Factory, RegistrationContact<TriangleModel, PointModel> > TrianglePointRegistrationContactClass("registration",true);
Creator<Contact::Factory, RegistrationContact<TriangleModel, LineModel> > TriangleLineRegistrationContactClass("registration",true);
Creator<Contact::Factory, RegistrationContact<TriangleModel, TriangleModel> > TriangleTriangleRegistrationContactClass("registration",true);
Creator<Contact::Factory, RegistrationContact<TetrahedronModel, SphereModel> > TetrahedronSphereRegistrationContactClass("registration",true);
Creator<Contact::Factory, RegistrationContact<TetrahedronModel, PointModel> > TetrahedronPointRegistrationContactClass("registration",true);
Creator<Contact::Factory, RegistrationContact<TetrahedronModel, LineModel> > TetrahedronLineRegistrationContactClass("registration",true);
Creator<Contact::Factory, RegistrationContact<TetrahedronModel, TriangleModel> > TetrahedronTriangleRegistrationContactClass("registration",true);
Creator<Contact::Factory, RegistrationContact<TetrahedronModel, TetrahedronModel> > TetrahedronTetrahedronRegistrationContactClass("registration",true);

Creator<Contact::Factory, RegistrationContact<RigidDistanceGridCollisionModel, RigidDistanceGridCollisionModel> > DistanceGridDistanceGridRegistrationContactClass("registration", true);
Creator<Contact::Factory, RegistrationContact<RigidDistanceGridCollisionModel, PointModel> > DistanceGridPointRegistrationContactClass("registration", true);
Creator<Contact::Factory, RegistrationContact<RigidDistanceGridCollisionModel, SphereModel> > DistanceGridSphereRegistrationContactClass("registration", true);
Creator<Contact::Factory, RegistrationContact<RigidDistanceGridCollisionModel, TriangleModel> > DistanceGridTriangleRegistrationContactClass("registration", true);

Creator<Contact::Factory, RegistrationContact<FFDDistanceGridCollisionModel, FFDDistanceGridCollisionModel> > FFDDistanceGridRegistrationContactClass("registration", true);
Creator<Contact::Factory, RegistrationContact<FFDDistanceGridCollisionModel, RigidDistanceGridCollisionModel> > FFDDistanceGridRigidDistanceGridRegistrationContactClass("registration", true);
Creator<Contact::Factory, RegistrationContact<FFDDistanceGridCollisionModel, PointModel> > FFDDistanceGridPoinRegistrationtContactClass("registration", true);
Creator<Contact::Factory, RegistrationContact<FFDDistanceGridCollisionModel, SphereModel> > FFDDistanceGridSphereRegistrationContactClass("registration", true);
Creator<Contact::Factory, RegistrationContact<FFDDistanceGridCollisionModel, TriangleModel> > FFDDistanceGridTriangleRegistrationContactClass("registration", true);


} // namespace collision

} // namespace component

} // namespace sofa

