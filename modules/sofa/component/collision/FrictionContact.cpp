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
#include <sofa/component/collision/FrictionContact.inl>

namespace sofa
{

namespace component
{

namespace collision
{

using namespace defaulttype;
using namespace sofa::helper;
using simulation::Node;

sofa::core::collision::DetectionOutput::ContactId Identifier::cpt=0;
std::list<sofa::core::collision::DetectionOutput::ContactId> Identifier::availableId;

SOFA_DECL_CLASS(FrictionContact)

Creator<Contact::Factory, FrictionContact<PointModel, PointModel> > PointPointFrictionContactClass("FrictionContact",true);
Creator<Contact::Factory, FrictionContact<LineModel, SphereModel> > LineSphereFrictionContactClass("FrictionContact",true);
Creator<Contact::Factory, FrictionContact<LineModel, PointModel> > LinePointFrictionContactClass("FrictionContact",true);
Creator<Contact::Factory, FrictionContact<LineModel, LineModel> > LineLineFrictionContactClass("FrictionContact",true);
Creator<Contact::Factory, FrictionContact<TriangleModel, SphereModel> > TriangleSphereFrictionContactClass("FrictionContact",true);
Creator<Contact::Factory, FrictionContact<TriangleModel, PointModel> > TrianglePointFrictionContactClass("FrictionContact",true);
Creator<Contact::Factory, FrictionContact<TriangleModel, LineModel> > TriangleLineFrictionContactClass("FrictionContact",true);
Creator<Contact::Factory, FrictionContact<TriangleModel, TriangleModel> > TriangleTriangleFrictionContactClass("FrictionContact",true);
Creator<Contact::Factory, FrictionContact<SphereModel, SphereModel> > SphereSphereFrictionContactClass("FrictionContact",true);
Creator<Contact::Factory, FrictionContact<SphereModel, PointModel> > SpherePointFrictionContactClass("FrictionContact",true);
Creator<Contact::Factory, FrictionContact<CapsuleModel, CapsuleModel> > CapsuleCapsuleFrictionContactClass("FrictionContact",true);
Creator<Contact::Factory, FrictionContact<CapsuleModel, TriangleModel> > CapsuleTriangleFrictionContactClass("FrictionContact",true);
Creator<Contact::Factory, FrictionContact<CapsuleModel, SphereModel> > CapsuleSphereFrictionContactClass("FrictionContact",true);
Creator<Contact::Factory, FrictionContact<OBBModel, OBBModel> > OBBOBBFrictionContactClass("FrictionContact",true);
Creator<Contact::Factory, FrictionContact<SphereModel, OBBModel> > SphereOBBFrictionContactClass("FrictionContact",true);
Creator<Contact::Factory, FrictionContact<CapsuleModel, OBBModel> > CapsuleOBBFrictionContactClass("FrictionContact",true);
Creator<Contact::Factory, FrictionContact<TriangleModel, OBBModel> > TriangleOBBFrictionContactClass("FrictionContact",true);
Creator<Contact::Factory, FrictionContact<RigidSphereModel, RigidSphereModel> > RigidSphereRigidSphereFrictionContactClass("FrictionContact",true);
Creator<Contact::Factory, FrictionContact<SphereModel, RigidSphereModel> > SphereRigidSphereFrictionContactClass("FrictionContact",true);
Creator<Contact::Factory, FrictionContact<LineModel, RigidSphereModel> > LineRigidSphereFrictionContactClass("FrictionContact",true);
Creator<Contact::Factory, FrictionContact<TriangleModel, RigidSphereModel> > TriangleRigidSphereFrictionContactClass("FrictionContact",true);
Creator<Contact::Factory, FrictionContact<RigidSphereModel, PointModel> > RigidSpherePointFrictionContactClass("FrictionContact",true);
Creator<Contact::Factory, FrictionContact<CapsuleModel, RigidSphereModel> > CapsuleRigidSphereFrictionContactClass("FrictionContact",true);
Creator<Contact::Factory, FrictionContact<RigidSphereModel, OBBModel> > RigidSphereOBBFrictionContactClass("FrictionContact",true);




} // namespace collision

} // namespace component

} // namespace sofa
