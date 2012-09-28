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
#ifndef SOFA_COMPONENT_COLLISION_BARYCENTRICPENALITYCONTACT_H
#define SOFA_COMPONENT_COLLISION_BARYCENTRICPENALITYCONTACT_H

#include <sofa/core/collision/Contact.h>
#include <sofa/core/collision/Intersection.h>
#include <sofa/component/mapping/BarycentricMapping.h>
//#ifdef SOFA_TEST_FRICTION
//#include <sofa/component/forcefield/PenalityContactFrictionForceField.h>
//#else
#include <sofa/component/interactionforcefield/PenalityContactForceField.h>
//#endif
#include <sofa/helper/Factory.h>

#include <sofa/component/collision/BarycentricContactMapper.h>
#include <sofa/component/collision/IdentityContactMapper.h>
#include <sofa/component/collision/RigidContactMapper.inl>


namespace sofa
{

namespace component
{

namespace collision
{

using namespace sofa::defaulttype;

template < class TCollisionModel1, class TCollisionModel2, class ResponseDataTypes = sofa::defaulttype::Vec3Types >
class BarycentricPenalityContact : public core::collision::Contact
{
public:
    SOFA_CLASS(SOFA_TEMPLATE3(BarycentricPenalityContact, TCollisionModel1, TCollisionModel2, ResponseDataTypes), core::collision::Contact);

    typedef TCollisionModel1 CollisionModel1;
    typedef TCollisionModel2 CollisionModel2;
    typedef core::collision::Intersection Intersection;
    typedef core::collision::DetectionOutputVector OutputVector;
    typedef core::collision::TDetectionOutputVector<CollisionModel1,CollisionModel2> TOutputVector;
    typedef ResponseDataTypes DataTypes1;
    typedef ResponseDataTypes DataTypes2;
    typedef core::behavior::MechanicalState<DataTypes1> MechanicalState1;
    typedef core::behavior::MechanicalState<DataTypes2> MechanicalState2;
    typedef typename CollisionModel1::Element CollisionElement1;
    typedef typename CollisionModel2::Element CollisionElement2;
//#ifdef SOFA_TEST_FRICTION
//	typedef forcefield::PenalityContactFrictionForceField<ResponseDataTypes> ResponseForceField;
//#else
    typedef interactionforcefield::PenalityContactForceField<ResponseDataTypes> ResponseForceField;
//#endif
protected:
    CollisionModel1* model1;
    CollisionModel2* model2;
    Intersection* intersectionMethod;

    ContactMapper<CollisionModel1,DataTypes1> mapper1;
    ContactMapper<CollisionModel2,DataTypes2> mapper2;

    typename ResponseForceField::SPtr ff;
    core::objectmodel::BaseContext* parent;

    typedef std::map<core::collision::DetectionOutput::ContactId,int> ContactIndexMap;
    /// Mapping of contactids to force element (+1, so that 0 means not active).
    /// This allows to ignore duplicate contacts, and preserve information associated with each contact point over time
    ContactIndexMap contactIndex;

    BarycentricPenalityContact(CollisionModel1* model1, CollisionModel2* model2, Intersection* intersectionMethod);
    ~BarycentricPenalityContact();

    void setInteractionTags(MechanicalState1* mstate1, MechanicalState2* mstate2);

public:
    void cleanup();

    std::pair<core::CollisionModel*,core::CollisionModel*> getCollisionModels() { return std::make_pair(model1,model2); }

    void setDetectionOutputs(OutputVector* outputs);

    void createResponse(core::objectmodel::BaseContext* group);

    void removeResponse();

    void draw(const core::visual::VisualParams* vparams);

};

} // namespace collision

} // namespace component

} // namespace sofa

#endif
