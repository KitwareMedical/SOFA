//
// C++ Interface: LDIConstraintContact
//
// Description: Interaction ForceField used to create a repulsion, and separate two colliding objects
//
//
// Author: Francois Faure, Sebastien Barbier, Jeremie Allard, Florent Falipou
//
// Licence: QPL, See LICENCE.txt file that comes with this distribution
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SOFA_COMPONENT_COLLISION_LDICONSTRAINTCONTACT_INL
#define SOFA_COMPONENT_COLLISION_LDICONSTRAINTCONTACT_INL

//#include <sofa/component/collision/LDIConstraintContact.h>
#include "LDIConstraintContact.h"
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/objectmodel/BaseContext.h>


#define DEBUG_LDICONTACTCONSTRAINT
#ifdef  DEBUG_LDICONTACTCONSTRAINT
#define DEBUG_LDICONTACTCONSTRAINT_OUT(c) \
  this->f_printLog.setValue(sofa::core::ExecParams::defaultInstance(),true); \
  c
#else
#define DEBUG_LDICONTACTCONSTRAINT_OUT(c)
#endif



namespace sofa
{

namespace component
{

namespace collision
{

using namespace sofa::defaulttype;
using namespace core::collision;
using simulation::Node;

template < class TCollisionModel1, class TCollisionModel2, class ResponseDataTypes >
LDIConstraintContact<TCollisionModel1,TCollisionModel2,ResponseDataTypes>::LDIConstraintContact(CollisionModel1* model1, CollisionModel2* model2, Intersection* intersectionMethod)
  : model1(model1), model2(model2), intersectionMethod(intersectionMethod), zeroVolume(NULL), mstate1(NULL), mstate2(NULL), parent(NULL)
{
}

template < class TCollisionModel1, class TCollisionModel2, class ResponseDataTypes >
LDIConstraintContact<TCollisionModel1,TCollisionModel2,ResponseDataTypes>::~LDIConstraintContact()
{
}


template < class TCollisionModel1, class TCollisionModel2, class ResponseDataTypes >
void LDIConstraintContact<TCollisionModel1,TCollisionModel2,ResponseDataTypes>::cleanup()
{
	if (zeroVolume!=NULL)
	{
          zeroVolume->resetConstraint();
          if (parent!=NULL) parent->removeObject(zeroVolume);


    }
}

template < class TCollisionModel1, class TCollisionModel2, class ResponseDataTypes >
void LDIConstraintContact<TCollisionModel1,TCollisionModel2,ResponseDataTypes>::setDetectionOutputs(OutputVector* o)
{  

  LDIDetectionOutputVector& outputs = *static_cast<LDIDetectionOutputVector*>(o);
  const bool printLog = this->f_printLog.getValue();

  int insize = outputs.size();

  if (insize == 0) { 
    if (zeroVolume) zeroVolume->setActive(false);
    return;
  }  

  if (zeroVolume==NULL )
  {
      if (printLog) this->sout << "mstate1=";
      core::behavior::BaseMechanicalState* mstate1 = getBaseState(static_cast< TriangleModel* >(this->model1)->getMechanicalState());
      if (printLog) this->sout << "mstate2=";
      core::behavior::BaseMechanicalState* mstate2 = getBaseState(static_cast< TriangleModel* >(this->model2)->getMechanicalState());
      if (printLog) this->sout << this->sendl;

      //TEMPORARY: need to modify ContactConstraint and more generally LMConstraint to take two templates
      core::behavior::MechanicalState<Vec3Types>* m1 = dynamic_cast<core::behavior::MechanicalState<Vec3Types>*>(mstate1);
      core::behavior::MechanicalState<Vec3Types>* m2 = dynamic_cast<core::behavior::MechanicalState<Vec3Types>*>(mstate2); 
      zeroVolume= sofa::core::objectmodel::New<ContactContraint>(m1,m2);
      zeroVolume->init();
      /////////////////////////////////////////////
}
  
  zeroVolume->setActive(true);

  

  // old index for each contact
  // >0 indicate preexisting contact
  // 0  indicate new contact
  // -1 indicate ignored duplicate contact
  std::vector<int> oldIndex(insize);
    
  int nbnew = 0;
  for (int i=0; i<(int)outputs.size(); i++ )
    {
      LDIDetectionOutput* o = &outputs.Vector[i];

      // find this contact in contactIndex, possibly creating a new entry initialized by 0
      int& index = contactIndex[o->id];
      if (index < 0) // duplicate contact
        {
	  int i2 = -1-index;
	  oldIndex[i] = oldIndex[i2];
	  oldIndex[i2] = -1;
        }
      else
        {
	  oldIndex[i] = index;
	  if (!index)
            {
	      ++nbnew;
	      if (printLog) std::cout << "LDIConstraintContact: New contact "<<o->id<<std::endl;
            }
        }
      index = -1-i; // save this index as a negative value in contactIndex map.
    }
  insize = outputs.size();
  // compute new index of each contact
  std::vector<int> newIndex(insize);
  // number of final contacts used in the response
  int size = 0;
  for (int i=0; i<insize; i++)
    {
      if (oldIndex[i] >= 0)
        {
	  ++size;
	  newIndex[i] = size;
        }
    }
    
  // update contactMap
  for (ContactIndexMap::iterator it = contactIndex.begin(), itend = contactIndex.end(); it != itend; )
    {
      int& index = it->second;
      if (index >= 0)
        {
	  if (printLog) std::cout << "LDIConstraintContact: Removed contact "<<it->first<<std::endl;
	  ContactIndexMap::iterator oldit = it;
	  ++it;
	  contactIndex.erase(oldit);
        }
      else
        {
	  index = newIndex[-1-index]; // write the final contact index
	  ++it;
        }
    }
  Vec<3,Real> noCollidingVolume = Vec<3,Real>(0,0,0);
  output.clear();
  
enum {R,G,B};


  for (int i=0; i<(int)outputs.size(); i++)
    {
      int index = oldIndex[i];
      if (index < 0) continue; // this contact is ignored
	
      LDIDetectionOutput* o = &outputs.Vector[i];
		      
      const int indexSubDiv = o->index;
      if (o->volume != NULL) noCollidingVolume[o->normal] += *o->volume;
          
      LDIOutput &out = output[indexSubDiv];


      const float alpha[2] = {(1-(o->point[0][0] + o->point[0][1])),(1-(o->point[1][0] + o->point[1][1]))};
      const float beta[2]  = {o->point[0][0],o->point[1][0]};
      const float gamma[2] = {o->point[0][1],o->point[1][1]};

      Vec<3,int> idxDOF[2];     
      
      
      float dScorrected = o->dS*o->subdivisionFactor;
      CollisionElement1 elem1(o->elem.first);
      CollisionElement2 elem2(o->elem.second);  
      idxDOF[0][R] = elem1.p1Index();idxDOF[1][R] = elem2.p1Index();
      idxDOF[0][G] = elem1.p2Index();idxDOF[1][G] = elem2.p2Index();
      idxDOF[0][B] = elem1.p3Index();idxDOF[1][B] = elem2.p3Index();
      
      int idxModel1, idxModel2;

      if (this->model1 == elem1.getCollisionModel())
        {
          idxModel1=0;
          idxModel2=1;
        }
      else
        {
          idxModel1=1;
          idxModel2=0;
        }

      {
        LDIOutput::SparseVec_dVdx &dVdx2=out.dVdx[idxModel2];
        const Vec<3,int> &dof=idxDOF[1];
        dVdx2[dof[R]][o->normal] += dScorrected*alpha[1];
        dVdx2[dof[G]][o->normal] += dScorrected*beta [1];
        dVdx2[dof[B]][o->normal] += dScorrected*gamma[1];
      }
      
      {
        LDIOutput::SparseVec_dVdx &dVdx1=out.dVdx[idxModel1];
        const Vec<3,int> &dof=idxDOF[0];
        dVdx1[dof[R]][o->normal] -= dScorrected*alpha[0];
        dVdx1[dof[G]][o->normal] -= dScorrected*beta [0];
        dVdx1[dof[B]][o->normal] -= dScorrected*gamma[0];
      }

      const sofa::helper::vector< sofa::defaulttype::Vector3 > pos[2] = {*(model1->getMechanicalState()->getX()),*(model2->getMechanicalState()->getX())};


      const Vector3 deltaX =
            (pos[idxModel1][idxDOF[0][R]]*alpha[0]+pos[idxModel2][idxDOF[1][R]]*alpha[1])+
            (pos[idxModel1][idxDOF[0][G]]*beta [0]+pos[idxModel2][idxDOF[1][G]]*beta [1])+
            (pos[idxModel1][idxDOF[0][B]]*gamma[0]+pos[idxModel2][idxDOF[1][B]]*gamma[1]);
          
      out.centralPoint += deltaX/2;


      ++out.contributions;
           
      out.volume[o->normal] += dScorrected*o->value;
      
    }


  LDIDetectionOutput* contact = &outputs.Vector[0];
  double frictionCoeff=(contact->elem.first.getContactFriction() + contact->elem.second.getContactFriction())*0.5;
  zeroVolume->setParameters(output, frictionCoeff); 
}

template < class TCollisionModel1, class TCollisionModel2, class ResponseDataTypes >
void LDIConstraintContact<TCollisionModel1,TCollisionModel2,ResponseDataTypes>::createResponse(core::objectmodel::BaseContext* group)
{
	if (zeroVolume!=NULL)
	{
		if (parent!=NULL)
		{
			parent->removeObject(this);
			parent->removeObject(zeroVolume);
		}
		parent = group;
		if (parent!=NULL)
		{
                    if (parent == this->model1->getContext() && mstate1)
                        parent = mstate1->getContext();
                    if (parent == this->model2->getContext() && mstate2)
                        parent = mstate2->getContext();
			//std::cout << "Attaching contact response to "<<parent->getName()<<std::endl;
			parent->addObject(this);
			parent->addObject(zeroVolume);
		}
	}
}

template < class TCollisionModel1, class TCollisionModel2, class ResponseDataTypes >
void LDIConstraintContact<TCollisionModel1,TCollisionModel2,ResponseDataTypes>::removeResponse()
{
	if (zeroVolume!=NULL)
	{
		if (parent!=NULL)
		{
			//std::cout << "Removing contact response from "<<parent->getName()<<std::endl;
			parent->removeObject(this);
			parent->removeObject(zeroVolume);
		}
		parent = NULL;
	}
}

template < class TCollisionModel1, class TCollisionModel2, class ResponseDataTypes >
void LDIConstraintContact<TCollisionModel1,TCollisionModel2,ResponseDataTypes>::draw()
{
}


} // namespace collision

} // namespace component

} // namespace sofa

#endif
