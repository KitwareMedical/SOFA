/******************************************************************************
 *       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 3      *
 *                (c) 2006-2008 MGH, INRIA, USTL, UJF, CNRS                    *
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
#ifndef SOFA_COMPONENT_CONSTRAINT_CONTACTCONSTRAINT_INL
#define SOFA_COMPONENT_CONSTRAINT_CONTACTCONSTRAINT_INL

#include "ContactConstraint.h"
#include <sofa/helper/gl/Axis.h>
#include <sofa/core/visual/VisualParams.h>

//#define DEBUG_LDICONSTRAINTRESPONSE
#ifdef  DEBUG_LDICONSTRAINTRESPONSE
#define DEBUG_CONSTRAINTRESPONSE_OUT(c) \
  this->f_printLog.setValue(sofa::core::ExecParams::defaultInstance(),true); \
  c
#else
#define DEBUG_CONSTRAINTRESPONSE_OUT(c)
#endif




namespace sofa
{

  namespace component
  {

    namespace constraint
    {



#define M_PI 3.14159265358979323846
      //********************************************************
      //Functors: 

      template <typename DataTypes>
      struct CreateJacobian
      {
        typedef typename DataTypes::Deriv              Deriv;
        typedef typename DataTypes::MatrixDeriv MatrixDeriv;
        typedef typename DataTypes::MatrixDeriv::RowIterator MatrixDerivRowIterator;

        
        typedef LDIOutput::SparseVec_dVdx::value_type  Entry;
        typedef LDIOutput::SparseVec_dVdx::key_type    Index;
        typedef LDIOutput::SparseVec_dVdx::mapped_type Value;

        CreateJacobian(MatrixDerivRowIterator&_N, MatrixDerivRowIterator&_T1, MatrixDerivRowIterator&_T2,
          const Deriv &_n, const Deriv &_t1, const Deriv &_t2):N(_N), T1(_T1), T2(_T2), n(_n), t1(_t1), t2(_t2)
        {
        };

        void operator()( const Entry& v)
        {
          Index dof=v.first;
          const Value &value=v.second;
          const double Si=value*n;
          N.addCol(dof,  value);
          T1.addCol(dof, t1*Si);
          T2.addCol(dof, t2*Si);
        }

      protected:
        MatrixDerivRowIterator &N, &T1, &T2;
        const Deriv &n, &t1, &t2;
      };


      template<class DataTypes>
      void ContactConstraint<DataTypes>::setParameters(const std::map<unsigned int, LDIOutput > &_outputDetection, double _coeffFriction)
      {
        DEBUG_CONSTRAINTRESPONSE_OUT( this->sout << "=============" << this->sendl  );
        DEBUG_CONSTRAINTRESPONSE_OUT( this->sout << "setParameters" << this->sendl  );
        DEBUG_CONSTRAINTRESPONSE_OUT( this->sout << "=============" << this->sendl  );
        

        coeffFriction=_coeffFriction;
        if (coeffFriction < 0) coeffFriction=0;

        std::map<unsigned int, LDIOutput >::const_iterator it;
        for( it = _outputDetection.begin(); it!= _outputDetection.end(); ++it){
          const LDIOutput& outToTest = it->second;
          SReal v = std::max(outToTest.volume[0], std::max(outToTest.volume[1],outToTest.volume[2]) );
          LDIOutput out(outToTest);
          out.volumeValue = v;
          DEBUG_CONSTRAINTRESPONSE_OUT(this->sout << "Volume: " << out.volumeValue << this->sendl);
          out.centralPoint = out.centralPoint / out.contributions;
          DEBUG_CONSTRAINTRESPONSE_OUT( this->sout << "centralPoint: "<< out.centralPoint << this->sendl  );
          outputDetection.push_back(out);
        }



      }

      template <class DataTypes> 
      void ContactConstraint<DataTypes>::resetConstraint()
      {
        core::behavior::LMConstraint<DataTypes,DataTypes>::resetConstraint();
        clear();
      }

      template <class DataTypes>
      void ContactConstraint<DataTypes>::clear()
      {
        outputDetection.clear();
        constraintToContact.clear();
      }


      template<class DataTypes>
      void ContactConstraint<DataTypes>::buildConstraintMatrix(const core::ConstraintParams* /*cParams*/ /* PARAMS FIRST */, core::MultiMatrixDerivId cId, unsigned int &cIndex)
      {
        DEBUG_CONSTRAINTRESPONSE_OUT( this->sout << "=====================" << this->sendl  );
        DEBUG_CONSTRAINTRESPONSE_OUT( this->sout << "buildConstraintMatrix" << this->sendl  );
        DEBUG_CONSTRAINTRESPONSE_OUT( this->sout << "=====================" << this->sendl  );



        using namespace core;
        using namespace core::objectmodel;
        Data<MatrixDeriv>* dC1 = cId[this->constrainedObject1].write();
        helper::WriteAccessor<objectmodel::Data<MatrixDeriv> > c1 = *dC1;

        Data<MatrixDeriv>* dC2 = cId[this->constrainedObject2].write();
        helper::WriteAccessor<objectmodel::Data<MatrixDeriv> > c2 = *dC2;

        helper::vector<LDIOutput>::iterator it;
        
        for (it = outputDetection.begin(); it!= outputDetection.end(); ++it)
          {
            LDIOutput &out= *it;


            //-------------------------------------------------------------------------
            //Compute the Central points where the constraints will be applied
            //And the direction of the repulsion force
            //-------------------------------------------------------------------------

            defaulttype::Vec<3,SReal> repulsionForce;
            LDIOutput::SparseVec_dVdx::const_iterator entry_dVdx;
            for (entry_dVdx=out.dVdx[0].begin();entry_dVdx!=out.dVdx[0].end();++entry_dVdx)
              {
                repulsionForce -= entry_dVdx->second;
              }



            //-------------------------------------------------------------------------
            //Compute the orthonormal basis used to model the contact (n,t1,t2)
            //-------------------------------------------------------------------------
            //Compute normal direction
            out.n=repulsionForce;	
            out.n.normalize();    
            //Compute t1, and t2                                
            out.t1=Deriv(1,0,0);
            if (fabs(out.n*out.t1) >= 0.95) //t1 and direction are "almost" colinear
              out.t1=Deriv(0,1,0);
            
            out.t2=out.n.cross(out.t1);
            out.t1=out.t2.cross(out.n);      


            MatrixDerivRowIterator N1           = c1->writeLine(cIndex);
            MatrixDerivRowIterator N2           = c2->writeLine(cIndex);
            out.Nconstraint = cIndex++;
            MatrixDerivRowIterator T1_Friction1 = c1->writeLine(cIndex);
            MatrixDerivRowIterator T1_Friction2 = c2->writeLine(cIndex);
            out.T1constraint = cIndex++;
            MatrixDerivRowIterator T2_Friction1 = c1->writeLine(cIndex);
            MatrixDerivRowIterator T2_Friction2 = c2->writeLine(cIndex);
            out.T2constraint = cIndex++;
            //-------------------------------------------------------------------------
            //Write the equation L
            //-------------------------------------------------------------------------
          //Write Equations for Object1
            CreateJacobian<DataTypes> createJ1(N1,T1_Friction1,T2_Friction1,
                                               out.n, out.t1, out.t2);
            std::for_each(out.dVdx[0].begin(), out.dVdx[0].end(), createJ1);
            //Write Equations for Object2                
            CreateJacobian<DataTypes> createJ2(N2,T1_Friction2,T2_Friction2,
                                               out.n, out.t1, out.t2);
            std::for_each(out.dVdx[1].begin(), out.dVdx[1].end(), createJ2);


            //Using the Mask to reduce the passage through the mappings
            for (MatrixDerivColIterator it=N1.begin(); it!=N1.end();++it) this->constrainedObject1->forceMask.insertEntry(it.index());
            for (MatrixDerivColIterator it=N2.begin(); it!=N2.end();++it) this->constrainedObject2->forceMask.insertEntry(it.index());

          }



      }


      template<class DataTypes>
      void ContactConstraint<DataTypes>::writeConstraintEquations(unsigned int& lineNumber, core::MultiVecId id, ConstOrder order)
      {

        DEBUG_CONSTRAINTRESPONSE_OUT( this->sout << "========================" << this->sendl  );
        DEBUG_CONSTRAINTRESPONSE_OUT( this->sout << "writeConstraintEquations" << this->sendl  );
        DEBUG_CONSTRAINTRESPONSE_OUT( this->sout << "========================" << this->sendl  );
        using namespace core;
        using namespace core::objectmodel;
	      if (isSet)
	        {
            

            for (unsigned int c=0;c<outputDetection.size();++c)
              {
                LDIOutput& out = outputDetection[c];
                //------------------------------------------------------------
                //Compute Right Hand term                
                switch (order)
                {
                case core::ConstraintParams::VEL :
                  {
                   
                    ConstVecId v1 = id.getId(this->simulatedObject1);
                    ConstVecId v2 = id.getId(this->simulatedObject2);

                    core::behavior::ConstraintGroup *constraints = this->addGroupConstraint(order);                    
                    SReal correction_N_VEL  = this->simulatedObject1->getConstraintJacobianTimesVecDeriv(out.Nconstraint,v1);
                    correction_N_VEL       += this->simulatedObject2->getConstraintJacobianTimesVecDeriv(out.Nconstraint,v2);
                    constraints->addConstraint(lineNumber, out.Nconstraint, -correction_N_VEL);

                    if (coeffFriction > 0.0)
                    {
                      SReal correction_T1_VEL = this->simulatedObject1->getConstraintJacobianTimesVecDeriv(out.T1constraint,v1);
                      correction_T1_VEL      += this->simulatedObject2->getConstraintJacobianTimesVecDeriv(out.T1constraint,v2);
                      SReal correction_T2_VEL = this->simulatedObject1->getConstraintJacobianTimesVecDeriv(out.T2constraint,v1);
                      correction_T2_VEL      += this->simulatedObject2->getConstraintJacobianTimesVecDeriv(out.T2constraint,v2);

                      constraints->addConstraint(lineNumber, out.T1constraint, -correction_T1_VEL);
                      constraints->addConstraint(lineNumber, out.T2constraint, -correction_T2_VEL);
                    }

                    constraintToContact.insert(std::make_pair(constraints, &out) );

                    break;
                  }                     
                case  core::ConstraintParams::POS :
                  {

                    //--------------------------------------------------------------------------------------
                    //Adding the constraint to the VecConst
                    //Non Penetration Constraint: Direction along the normal of the contact
                    core::behavior::ConstraintGroup *constraints = this->addGroupConstraint(order);
                    constraints->addConstraint(lineNumber, out.Nconstraint, out.volumeValue);
                    constraintToContact.insert(std::make_pair(constraints, &out) );
                    break;
                  }
                default:
                  break;
                }
      
              }           
          }

      }
    
      template <class DataTypes> 
      void ContactConstraint<DataTypes>::LagrangeMultiplierEvaluation(const SReal* W, const SReal* c, SReal* Lambda, core::behavior::ConstraintGroup * group)
      {

        DEBUG_CONSTRAINTRESPONSE_OUT( this->sout << "============================" << this->sendl  );
        DEBUG_CONSTRAINTRESPONSE_OUT( this->sout << "LagrangeMultiplierEvaluation" << this->sendl  );
        DEBUG_CONSTRAINTRESPONSE_OUT( this->sout << "============================" << this->sendl  );

        const int numConstraintToProcess = group->getNumConstraint();
       

        switch (group->getOrder())
        {
        case core::ConstraintParams::VEL :
          {
            LDIOutput &out =  *(this->constraintToContact[group]);
           
            component::constraintset::ContactDescription &contact=this->getContactDescription(group);

            //The force cannot be attractive!
            if (Lambda[0] <= 0)
            {
              contact.state=component::constraintset::VANISHING;
              group->setActive(false);
              out.contactForce = Deriv();
              return;
            }

            if (numConstraintToProcess == 3)
            {
              const SReal normTangentForce=sqrt(pow(Lambda[1],2)+pow(Lambda[2],2));
              const SReal normNormalForce=fabs(Lambda[0]);

              //Test if we are outside the coulomb friction cone
              if ( normTangentForce > coeffFriction*normNormalForce)
              {
                contact.state=component::constraintset::SLIDING;
                //Force applied
                out.contactForce = out.n*Lambda[0]+out.t1*Lambda[1]+ out.t2*Lambda[2];
                //Outside: we project the force to the cone

                //directionCone <--> n' : unitary vector along the cone
                const SReal factor=coeffFriction*normNormalForce/normTangentForce;

                Deriv directionCone=out.n*Lambda[0]+(out.t1*Lambda[1]+ out.t2*Lambda[2])*factor;
                directionCone.normalize();
                

                contact.coeff[0]=out.n  *directionCone;
                contact.coeff[1]=out.t1 *directionCone;
                contact.coeff[2]=out.t2 *directionCone;

                const SReal value=W[0]*contact.coeff[0]+
                  W[1]*contact.coeff[1] +
                  W[2]*contact.coeff[2];


                if (value == 0)
                {
                  serr << "ERROR DIVISION BY ZERO AVOIDED: w=[" << W[0]  << "," << W[1] << "," << W[2]  << "] " << " DIRECTION CONE: " << directionCone << " BARY COEFF: " << contact.coeff[0] << ", " <<  contact.coeff[1] << ", " <<  contact.coeff[2] << std::endl;
                  group->setActive(false);
                  out.contactForce=Deriv();
                  return;
                }
                const SReal slidingLambda=c[0]/value;
                out.contactForce = directionCone*slidingLambda;
                //Then project the force to the border of the cone
                Lambda[0]=out.contactForce*out.n;
                //The force cannot be attractive!
                if (Lambda[0] <= 0)
                {
                  group->setActive(false);
                  out.contactForce=Deriv();
                  return;
                }
                Lambda[1]=out.contactForce*out.t1;
                Lambda[2]=out.contactForce*out.t2;

              }
              else contact.state=component::constraintset::STICKING;

              out.contactForce = out.n*Lambda[0]+out.t1*Lambda[1]+out.t2*Lambda[2];
            }
            else
              out.contactForce = out.n*Lambda[0];

            break;
          }
        case core::ConstraintParams::POS :
          {
            //The force cannot be attractive!
            if (Lambda[0] < 0)
            {                    
              group->setActive(false);
              return;
            }
          }
        default:{}
        }

        return;


      }



      template <class DataTypes>
      void ContactConstraint<DataTypes>::draw( const core::visual::VisualParams* vparams )
      {
         if ( !vparams->displayFlags().getShowInteractionForceFields() || !isSet ) return;
	 
            for (unsigned int i=0;i<outputDetection.size();++i)
              {
                LDIOutput& out = outputDetection[i];
                //Red: Contact Direction
                glColor4f(1.0f,0.0f,0.0f,1.0f);
                helper::gl::Axis::draw(out.centralPoint, out.centralPoint+out.n*3.0,(float)0.05);

                //Green: T1 Direction
                glColor4f(0.0f,1.0f,0.0f,1.0f);
                helper::gl::Axis::draw(out.centralPoint, out.centralPoint+out.t1*3.0,(float)0.05);

                //Blue: T2 Direction
                glColor4f(0.0f,0.0f,1.0f,1.0f);
                helper::gl::Axis::draw(out.centralPoint, out.centralPoint+out.t2*3.0,(float)0.05);

                //Direction of the Force
                Deriv f=out.contactForce; f.normalize(); //We normalize the response vector
                glColor4f(1.0f,0.0f,1.0f,1.0f);
                helper::gl::Axis::draw(out.centralPoint, out.centralPoint+f*3.0,(float)0.05);


              }

      }




    } // namespace constraint

  } // namespace component

} // namespace sofa

#endif
