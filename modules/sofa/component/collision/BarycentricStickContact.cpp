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
#include <sofa/component/collision/BarycentricStickContact.inl>
#include <sofa/component/collision/BarycentricContactMapper.h>
#include <sofa/component/collision/IdentityContactMapper.h>
#include <sofa/component/collision/DistanceGridCollisionModel.h>

namespace sofa
{

namespace component
{

namespace collision
{

SOFA_DECL_CLASS(BarycentricStickContact)

Creator<Contact::Factory, BarycentricStickContact<SphereModel, SphereModel> > SphereSphereStickContactClass("stick",true);
Creator<Contact::Factory, BarycentricStickContact<SphereModel, PointModel> > SpherePointStickContactClass("stick",true);
Creator<Contact::Factory, BarycentricStickContact<PointModel, PointModel> > PointPointStickContactClass("stick",true);
Creator<Contact::Factory, BarycentricStickContact<LineModel, PointModel> > LinePointStickContactClass("stick",true);
Creator<Contact::Factory, BarycentricStickContact<LineModel, LineModel> > LineLineStickContactClass("stick",true);
Creator<Contact::Factory, BarycentricStickContact<LineModel, SphereModel> > LineSphereStickContactClass("stick",true);
Creator<Contact::Factory, BarycentricStickContact<TriangleModel, SphereModel> > TriangleSphereStickContactClass("stick",true);
Creator<Contact::Factory, BarycentricStickContact<TriangleModel, PointModel> > TrianglePointStickContactClass("stick",true);
Creator<Contact::Factory, BarycentricStickContact<TriangleModel, LineModel> > TriangleLineStickContactClass("stick",true);
Creator<Contact::Factory, BarycentricStickContact<TriangleModel, TriangleModel> > TriangleTriangleStickContactClass("stick",true);

Creator<Contact::Factory, BarycentricStickContact<RigidDistanceGridCollisionModel, RigidDistanceGridCollisionModel> > DistanceGridDistanceGridStickContactClass("stick", true);
Creator<Contact::Factory, BarycentricStickContact<RigidDistanceGridCollisionModel, PointModel> > DistanceGridPointStickContactClass("stick", true);
Creator<Contact::Factory, BarycentricStickContact<RigidDistanceGridCollisionModel, SphereModel> > DistanceGridSphereStickContactClass("stick", true);
Creator<Contact::Factory, BarycentricStickContact<RigidDistanceGridCollisionModel, TriangleModel> > DistanceGridTriangleStickContactClass("stick", true);

Creator<Contact::Factory, BarycentricStickContact<FFDDistanceGridCollisionModel, FFDDistanceGridCollisionModel> > FFDDistanceGridStickContactClass("stick", true);
Creator<Contact::Factory, BarycentricStickContact<FFDDistanceGridCollisionModel, RigidDistanceGridCollisionModel> > FFDDistanceGridRigidDistanceGridStickContactClass("stick", true);
Creator<Contact::Factory, BarycentricStickContact<FFDDistanceGridCollisionModel, PointModel> > FFDDistanceGridPointStickContactClass("stick", true);
Creator<Contact::Factory, BarycentricStickContact<FFDDistanceGridCollisionModel, SphereModel> > FFDDistanceGridSphereStickContactClass("stick", true);
Creator<Contact::Factory, BarycentricStickContact<FFDDistanceGridCollisionModel, TriangleModel> > FFDDistanceGridTriangleStickContactClass("stick", true);


} // namespace collision

} // namespace component

} // namespace sofa

