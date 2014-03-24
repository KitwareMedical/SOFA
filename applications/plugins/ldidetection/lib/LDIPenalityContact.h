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
#ifndef SOFA_COMPONENT_COLLISION_LDIPENALITYCONTACT_H
#define SOFA_COMPONENT_COLLISION_LDIPENALITYCONTACT_H

#include "LDIDetection.h"
#include <sofa/component/collision/BarycentricPenalityContact.h>
#include <sofa/core/collision/Contact.h>
#include <sofa/core/collision/Intersection.h>
#include <sofa/component/mapping/BarycentricMapping.h>
/* #include <sofa/component/forcefield/PenalityContactForceField.h> */

#include "LDIPenalityContactForceField.h"
#include <sofa/helper/Factory.h>
#include <sofa/component/collision/BarycentricContactMapper.h>


namespace sofa
{

  namespace component
  {

    namespace collision
    {

      using namespace sofa::defaulttype;

      template < class TCollisionModel1=TriangleModel, class TCollisionModel2=TriangleModel, class ResponseDataTypes = sofa::defaulttype::Vec3Types >
        class SOFA_LDIDETECTION_API LDIPenalityContact :  public core::collision::Contact
	{


public:
SOFA_CLASS(LDIPenalityContact,core::collision::Contact);
	typedef TCollisionModel1 CollisionModel1;
	typedef TCollisionModel2 CollisionModel2;
	typedef core::collision::Intersection Intersection;
        typedef core::collision::DetectionOutputVector OutputVector;
        typedef core::collision::TDetectionOutputVector<TriangleModel,TriangleModel> TOutputVector;
	typedef ResponseDataTypes DataTypes;
	//typedef core::componentmodel::behavior::MechanicalState<DataTypes> MechanicalState1;
	//typedef core::componentmodel::behavior::MechanicalState<DataTypes> MechanicalState2;
	typedef typename CollisionModel1::Element CollisionElement1;
	typedef typename CollisionModel2::Element CollisionElement2;
        //typedef forcefield::LDIPenalityContactForceField<ResponseDataTypes,ResponseDataTypes> ResponseForceField;
        typedef forcefield::BaseLDIPenalityContactForceField<ResponseDataTypes> ResponseForceField;
        typedef typename ResponseDataTypes::Coord::value_type Real;

protected:
	CollisionModel1* model1;
	CollisionModel2* model2;
	Intersection* intersectionMethod;

	sofa::helper::vector< Vector3 > dVdx[2];
	sofa::helper::vector< Vector3 > dVelocity[2];
	std::set< unsigned int > index_used[2];

    typename ResponseForceField::SPtr ff;
        core::behavior::BaseMechanicalState* mstate1;
        core::behavior::BaseMechanicalState* mstate2;
        core::objectmodel::BaseContext* parent;

	typedef std::map<core::collision::DetectionOutput::ContactId,int> ContactIndexMap;
	/// Mapping of contactids to force element (+1, so that 0 means not active).
	/// This allows to ignore duplicate contacts, and preserve information associated with each contact point over time
	 ContactIndexMap contactIndex;

    core::behavior::BaseMechanicalState* getBaseState(core::behavior::BaseMechanicalState* mstate0)
    {
      const bool printLog = this->f_printLog.getValue();
        while (mstate0)
        {
            core::BaseMapping* mapping;
            mstate0->getContext()->get(mapping);
            if (!mapping) break;
            if (mapping->isMechanical()) break;
            if (!mapping->getMechFrom()[0]) break;
            if (printLog) std::cout << "->"<<mapping->getClassName();
            mstate0 = mapping->getMechFrom()[0];
        }
        if (printLog)  std::cout << "="<<mstate0->getTemplateName() << std::endl;
        return mstate0;
    }

    template<class DataTypes1, class DataTypes2>
    bool tryCreateFF(core::behavior::BaseMechanicalState* mstate1, core::behavior::BaseMechanicalState* mstate2)
    {
        core::behavior::MechanicalState<DataTypes1>* m1 = dynamic_cast<core::behavior::MechanicalState<DataTypes1>*>(mstate1);
        core::behavior::MechanicalState<DataTypes2>* m2 = dynamic_cast<core::behavior::MechanicalState<DataTypes2>*>(mstate2);
        if (!m1 || !m2) return false;
        ff = sofa::core::objectmodel::New< forcefield::LDIPenalityContactForceField<DataTypes1,DataTypes2,Vec3Types> >(m1, m2);
//        std::cout << "LDIPenalityContactForceField<"<<DataTypes1::Name()<<","<<DataTypes2::Name()<<"> created."<<std::endl;
        this->mstate1 = m1;
        this->mstate2 = m2;
        return true;
    }
public:
	LDIPenalityContact(CollisionModel1* model1, CollisionModel2* model2, Intersection* intersectionMethod);
	~LDIPenalityContact();

	void cleanup();

	std::pair<core::CollisionModel*,core::CollisionModel*> getCollisionModels() { return std::make_pair(model1,model2); }

	void setDetectionOutputs(OutputVector* outputs);

	void createResponse(core::objectmodel::BaseContext* group);

	void removeResponse();

	// -- VisualModel interface
	void draw();
	void initTextures() { }
	void update() { }
	};
    } // namespace collision

  } // namespace component

} // namespace sofa

#endif