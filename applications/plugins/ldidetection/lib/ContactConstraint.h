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
#ifndef SOFA_COMPONENT_CONSTRAINT_CONTACTCONSTRAINT_H
#define SOFA_COMPONENT_CONSTRAINT_CONTACTCONSTRAINT_H

#include "initldidetection.h"
#include <sofa/core/topology/BaseMeshTopology.h> 
#include <sofa/core/behavior/LMConstraint.h>
#include <sofa/helper/fixed_array.h>
#include <sofa/component/constraintset/ContactDescription.h>


namespace sofa
{

namespace component
{

namespace constraint
{

  struct LDIOutput
  {
    typedef std::map<unsigned int, defaulttype::Vec<3,SReal> > SparseVec_dVdx;
    LDIOutput()
    {
      centralPoint=defaulttype::Vec<3,SReal>();
      contributions = 0;
      volume.clear();
      dVdx[0].clear(); dVdx[1].clear();
    }
    defaulttype::Vec<3,SReal> volume;
    SparseVec_dVdx dVdx[2];



    SReal volumeValue;
    defaulttype::Vec<3,SReal> centralPoint;
    unsigned int contributions;
    //Index of the equations
    unsigned int Nconstraint;
    unsigned int T1constraint;
    unsigned int T2constraint;

    //Constrained Axis
    defaulttype::Vec<3,SReal> n,t1,t2;

    defaulttype::Vec<3,SReal> contactForce;

   
  };

  /// This class can be overridden if needed for additionnal storage within template specializations.
  template <class DataTypes>
  class SOFA_LDIDETECTION_API ContactConstraintInternalData
  {
  };




  /** Attach given particles to their initial positions.
  */
  template <class DataTypes>
  class ContactConstraint :  public core::behavior::LMConstraint<DataTypes,DataTypes>, public component::constraintset::ContactDescriptionHandler
  {
  public:
    SOFA_CLASS(SOFA_TEMPLATE(ContactConstraint,DataTypes),SOFA_TEMPLATE2(core::behavior::LMConstraint,DataTypes,DataTypes));
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename sofa::core::topology::BaseMeshTopology::SeqEdges SeqEdges;
    typedef typename sofa::core::topology::BaseMeshTopology::Edge Edge;

    typedef typename DataTypes::MatrixDeriv MatrixDeriv;
    typedef typename DataTypes::MatrixDeriv::RowIterator MatrixDerivRowIterator;
    typedef typename DataTypes::MatrixDeriv::ColIterator MatrixDerivColIterator;
    typedef typename DataTypes::MatrixDeriv::RowConstIterator MatrixDerivRowConstIterator;
    typedef typename DataTypes::MatrixDeriv::ColConstIterator MatrixDerivColConstIterator;

    typedef core::ConstraintParams::ConstOrder ConstOrder;


    typedef typename core::behavior::MechanicalState<DataTypes> MechanicalState;
    typedef typename core::VecId VecId;

  protected:
    ContactConstraintInternalData<DataTypes> data;
    friend class ContactConstraintInternalData<DataTypes>;

  public:
    ContactConstraint( MechanicalState *dof1, MechanicalState *dof2):
        core::behavior::LMConstraint<DataTypes,DataTypes>(dof1,dof2)
        {
          this->setName(dof1->getName() + "-" + dof2->getName());
        };
        ContactConstraint(){}

        ~ContactConstraint(){}; 

        // -- Constraint interface
        // -- Constraint interface

        void buildConstraintMatrix(const core::ConstraintParams* cParams /* PARAMS FIRST */, core::MultiMatrixDerivId cId, unsigned int &cIndex);
        void writeConstraintEquations(unsigned int& lineNumber, core::MultiVecId id, core::ConstraintParams::ConstOrder order);
        void LagrangeMultiplierEvaluation(const SReal* W, const SReal* c, SReal* Lambda, core::behavior::ConstraintGroup * group);


        void setParameters(const std::map<unsigned int, LDIOutput > &_outputDetection, double _coeffFriction);

        // -- virtual  
        void draw( const core::visual::VisualParams* vparams );
        void clear();
        void resetConstraint();
        // --

        void setActive(bool v){isSet=v;}

  protected:
    helper::vector<LDIOutput> outputDetection;
    std::map< core::behavior::ConstraintGroup*, LDIOutput*  > constraintToContact;

    double coeffFriction;
    bool isSet;
  };

} // namespace constraint

} // namespace component

} // namespace sofa


#endif
