//
// C++ Interface: LDIPenalityContactForceField
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
#ifndef SOFA_COMPONENT_COLLISION_LDIPENALITYCONTACT_INL
#define SOFA_COMPONENT_COLLISION_LDIPENALITYCONTACT_INL

//#include <sofa/component/collision/LDIPenalityContact.h>
#include "LDIPenalityContact.h"
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/objectmodel/BaseContext.h>



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
LDIPenalityContact<TCollisionModel1,TCollisionModel2,ResponseDataTypes>::LDIPenalityContact(CollisionModel1* model1, CollisionModel2* model2, Intersection* intersectionMethod)
: model1(model1), model2(model2), intersectionMethod(intersectionMethod), ff(NULL), mstate1(NULL), mstate2(NULL), parent(NULL)
{
}

template < class TCollisionModel1, class TCollisionModel2, class ResponseDataTypes >
LDIPenalityContact<TCollisionModel1,TCollisionModel2,ResponseDataTypes>::~LDIPenalityContact()
{
}


template < class TCollisionModel1, class TCollisionModel2, class ResponseDataTypes >
void LDIPenalityContact<TCollisionModel1,TCollisionModel2,ResponseDataTypes>::cleanup()
{
  if (ff!=NULL)
    {
      ff->cleanup();
      if (parent!=NULL)
      {
          parent->removeObject(ff);
          parent->removeObject(this);
      }
    }
}

template < class TCollisionModel1, class TCollisionModel2, class ResponseDataTypes >
void LDIPenalityContact<TCollisionModel1,TCollisionModel2,ResponseDataTypes>::setDetectionOutputs(OutputVector* o)
{  
  LDIDetectionOutputVector& outputs = *static_cast<LDIDetectionOutputVector*>(o);
  const bool printLog = this->f_printLog.getValue();

  int insize = outputs.size();
  if (insize == 0) return;


  if (ff==NULL )
  {
      if (printLog) this->sout << "mstate1=";
      core::behavior::BaseMechanicalState* mstate1 = getBaseState(static_cast< TriangleModel* >(this->model1)->getMechanicalState());
      if (printLog) this->sout << "mstate2=";
      core::behavior::BaseMechanicalState* mstate2 = getBaseState(static_cast< TriangleModel* >(this->model2)->getMechanicalState());
      if (printLog) this->sout << this->sendl;
      
      tryCreateFF<Vec3Types,Vec3Types>(mstate1, mstate2);

      if (ff==NULL)
      {
          std::cerr << "ERROR: failed to create LDIPenalityContactForceField"<<std::endl;
          return;
      }
      ff->init();
}


  const sofa::helper::vector< sofa::defaulttype::Vector3 > velocity[2] = {*(model1->getMechanicalState()->getV()),*(model2->getMechanicalState()->getV())};
  const sofa::helper::vector< sofa::defaulttype::Vector3 > pos     [2] = {*(model1->getMechanicalState()->getX()),*(model2->getMechanicalState()->getX())};
 
  
  sofa::helper::vector< LDIDetectionOutput >::iterator it= outputs.Vector.begin();

  // old index for each contact
  // >0 indicate preexisting contact
  // 0  indicate new contact
  // -1 indicate ignored duplicate contact
  std::vector<int> oldIndex(insize);
    
  int nbnew = 0;
  for (int i=0; i<(int)outputs.size(); i++,it++)
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
	      if (printLog) std::cout << "LDIPenalityContact: New contact "<<o->id<<std::endl;
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
	  if (printLog) std::cout << "LDIPenalityContact: Removed contact "<<it->first<<std::endl;
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
  if (printLog) std::cout << "LDIPenalityContact: "<<insize<<" input contacts, "<<size<<" contacts used for response ("<<nbnew<<" new)."<<std::endl;

  ff->clear(size);

  unsigned int sizeDOF = this->model1->getMechanicalState()->getX()->size();
  if (dVdx[0].size() != sizeDOF)
    {
      dVdx[0].resize(sizeDOF) ;
      dVelocity[0].resize(sizeDOF) ;
    }
  sizeDOF = this->model2->getMechanicalState()->getX()->size();
  if (dVdx[1].size() != sizeDOF) 
    {
      dVdx[1].resize(sizeDOF) ;
      dVelocity[1].resize(sizeDOF) ;
    }

  Vec<3,Real> volume = Vec<3,Real>(0,0,0);
  Vec<3,Real> noCollidingVolume = Vec<3,Real>(0,0,0);
  
  LDIDetectionOutput* contact = &outputs.Vector[0];
  
enum {R,G,B};
	
  for (int i=0; i<(int)outputs.size(); i++)
    {
      int index = oldIndex[i];
      if (index < 0) continue; // this contact is ignored
	
      LDIDetectionOutput* o = &outputs.Vector[i];
	
      const float alpha[2] = {(1-(o->point[0][0] + o->point[0][1])),(1-(o->point[1][0] + o->point[1][1]))};
      const float beta[2]  = {o->point[0][0],o->point[1][0]};
      const float gamma[2] = {o->point[0][1],o->point[1][1]};

      Vec<3,int> idxDOF[2];     
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
	
      //Relative velocity at one pixel
      Vector3 deltaV = 
	(velocity[idxModel1][idxDOF[0][R]]*alpha[0]-velocity[idxModel2][idxDOF[1][R]]*alpha[1])+
	(velocity[idxModel1][idxDOF[0][G]]*beta [0]-velocity[idxModel2][idxDOF[1][G]]*beta [1])+
	(velocity[idxModel1][idxDOF[0][B]]*gamma[0]-velocity[idxModel2][idxDOF[1][B]]*gamma[1]);
      const Vector3 deltaX = 	
	(pos[idxModel1][idxDOF[0][R]]*alpha[0]-pos[idxModel2][idxDOF[1][R]]*alpha[1])+
	(pos[idxModel1][idxDOF[0][G]]*beta [0]-pos[idxModel2][idxDOF[1][G]]*beta [1])+
	(pos[idxModel1][idxDOF[0][B]]*gamma[0]-pos[idxModel2][idxDOF[1][B]]*gamma[1]);
	

      bool attractive = ( deltaX[o->normal]*o->dS<= 0);
      if (!attractive)
	{

	  if (o->volume != NULL) noCollidingVolume[o->normal] += *o->volume;

          const float dScorrected = o->dS*o->subdivisionFactor;

	  volume[o->normal] += o->value*dScorrected;
	  deltaV[o->normal] = 0;
	  const double dS_Viscosity = -dScorrected*o->viscosity;

          {
            helper::vector< Vector3 > &dVdx1     = dVdx[idxModel1];
            helper::vector< Vector3 > &dVelocity1= dVelocity[idxModel1];
            std::set< unsigned int >  &indexUsed = index_used[idxModel1];
            const Vec<3,int> &dof=idxDOF[0];

            dVdx1[dof[R]][o->normal] += dScorrected *alpha[0];	
            dVdx1[dof[G]][o->normal] += dScorrected *beta [0];	
            dVdx1[dof[B]][o->normal] += dScorrected *gamma[0];	

            dVelocity1[dof[R]]       += deltaV*dS_Viscosity*alpha[0];    
            dVelocity1[dof[G]]       += deltaV*dS_Viscosity*beta [0];    
            dVelocity1[dof[B]]       += deltaV*dS_Viscosity*gamma[0];  

            indexUsed.insert(dof[R]); 
            indexUsed.insert(dof[G]);
            indexUsed.insert(dof[B]);
          }

          {
            helper::vector< Vector3 > &dVdx2     = dVdx[idxModel2];
            helper::vector< Vector3 > &dVelocity2= dVelocity[idxModel2];
            std::set< unsigned int >  &indexUsed = index_used[idxModel2];
            const Vec<3,int> &dof=idxDOF[1];

            dVdx2[dof[R]][o->normal] -= dScorrected *alpha[1];	
            dVdx2[dof[G]][o->normal] -= dScorrected *beta [1];	
            dVdx2[dof[B]][o->normal] -= dScorrected *gamma[1];	

            dVelocity2[dof[R]]       -= deltaV*dS_Viscosity*alpha[1];    
            dVelocity2[dof[G]]       -= deltaV*dS_Viscosity*beta [1];    
            dVelocity2[dof[B]]       -= deltaV*dS_Viscosity*gamma[1];  

            indexUsed.insert(dof[R]); 
            indexUsed.insert(dof[G]);
            indexUsed.insert(dof[B]);          
          }
	
	}
    }
    
    const double stiffness = contact->K*(contact->elem.first.getContactStiffness() + contact->elem.second.getContactStiffness());
//     volume += noCollidingVolume;
  const Real V      = std::max(std::max(volume[0], volume[1]), volume[2]);
  //const Real GV      = std::max(std::max(gvolume[0], gvolume[1]), gvolume[2]);
  //std::cout << "V="<<V<<" GV="<<GV<<"\n";
  //std::cout<<"LDIPenalityContact<TCollisionModel1,TCollisionModel2,ResponseDataTypes>::setDetectionOutputs, contact->K = "<<contact->K<<", stiffness = "<<stiffness<<std::endl;

  ff->setParameters(this->model1, this->model2,
	            V, 
	            dVdx[0],dVelocity[0],index_used[0],dVdx[1],dVelocity[1], index_used[1],
		    stiffness);
}

template < class TCollisionModel1, class TCollisionModel2, class ResponseDataTypes >
void LDIPenalityContact<TCollisionModel1,TCollisionModel2,ResponseDataTypes>::createResponse(core::objectmodel::BaseContext* group)
{
	if (ff!=NULL)
	{
		if (parent!=NULL)
		{
			parent->removeObject(ff);
			parent->removeObject(this);
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
			parent->addObject(ff);
		}
	}
}

template < class TCollisionModel1, class TCollisionModel2, class ResponseDataTypes >
void LDIPenalityContact<TCollisionModel1,TCollisionModel2,ResponseDataTypes>::removeResponse()
{
  if (ff!=NULL)
    {
      //Clear the vector dVdx using the index of the previous used vertices
      for (std::set<unsigned int>::const_iterator it=index_used[0].begin(); it != index_used[0].end(); it++)
        {
          dVdx[0][(*it)] = Vector3();
          dVelocity[0][(*it)] = Vector3();
        }
  
  
      for (std::set<unsigned int>::const_iterator it=index_used[1].begin(); it != index_used[1].end(); it++)
        {
          dVdx[1][(*it)] = Vector3();
          dVelocity[1][(*it)] = Vector3();
        }
      
      //Vec<3,Real> gvolume = Vec<3,Real>(0,0,0);
      index_used[0].clear();
      index_used[1].clear();

      ff->clear();

      if (parent!=NULL)
        {
          //std::cout << "Removing contact response from "<<parent->getName()<<std::endl;
          parent->removeObject(this);
          parent->removeObject(ff);
        }
      parent = NULL;
    }
}

template < class TCollisionModel1, class TCollisionModel2, class ResponseDataTypes >
void LDIPenalityContact<TCollisionModel1,TCollisionModel2,ResponseDataTypes>::draw()
{
//	if (dynamic_cast<core::VisualModel*>(ff)!=NULL)
//		dynamic_cast<core::VisualModel*>(ff)->draw();
}


} // namespace collision

} // namespace component

} // namespace sofa

#endif
