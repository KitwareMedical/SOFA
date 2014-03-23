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

#ifndef SOFA_COMPONENT_FORCEFIELD_LDIPENALITYCONTACTFORCEFIELD_H
#define SOFA_COMPONENT_FORCEFIELD_LDIPENALITYCONTACTFORCEFIELD_H

#include <sofa/core/behavior/MixedInteractionForceField.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/Mapping.h>
#include <sofa/core/CollisionModel.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <vector>
#include <sofa/component/collision/TriangleModel.h>


namespace sofa
{

namespace component
{

namespace forcefield
{

  using namespace sofa::defaulttype;
  using namespace sofa::component::collision;

  template<class DataTypes>
  class BaseLDIPenalityContactForceField : public virtual core::objectmodel::BaseObject
  {
  public:
    SOFA_CLASS(BaseLDIPenalityContactForceField,core::objectmodel::BaseObject);
    typedef typename DataTypes::Real Real;
    virtual void clear(int reserve = 0) = 0;
    virtual void setParameters( TriangleModel* model1, TriangleModel* model2, const Real &_Volume, 
      sofa::helper::vector< Vector3>  &_dVdx1, sofa::helper::vector< Vector3>  &_dVelocity1, const std::set< unsigned int> &_index_used1,
      sofa::helper::vector< Vector3 > &_dVdx2, sofa::helper::vector< Vector3>  &_dVelocity2, const std::set< unsigned int> &_index_used2,
      const Real _K) = 0;
  };

  template<class DataTypes1, class DataTypes2, class ResponseDataTypes>
  class LDIPenalityContactForceField : public core::behavior::MixedInteractionForceField<DataTypes1, DataTypes2>, public BaseLDIPenalityContactForceField<ResponseDataTypes> //PenalityContactForceField<DataTypes>
  {
  public:
    SOFA_CLASS2(SOFA_TEMPLATE3(LDIPenalityContactForceField,DataTypes1,DataTypes2,ResponseDataTypes),SOFA_TEMPLATE2(core::behavior::MixedInteractionForceField,DataTypes1,DataTypes2),SOFA_TEMPLATE(BaseLDIPenalityContactForceField,ResponseDataTypes));
    typedef typename Vector3::value_type Real_Sofa;

    typedef typename core::behavior::MixedInteractionForceField<DataTypes1, DataTypes2> Inherit;
    typedef typename DataTypes1::VecCoord VecCoord1;
    typedef typename DataTypes1::VecDeriv VecDeriv1;
    typedef typename DataTypes1::Coord Coord1;
    typedef typename DataTypes1::Deriv Deriv1;
    typedef typename DataTypes1::MatrixDeriv MatrixDeriv1;
    typedef typename MatrixDeriv1::RowConstIterator MatrixDerivRowConstIterator1;
    typedef core::behavior::MechanicalState<DataTypes1> MechanicalState1;

    typedef typename DataTypes2::VecCoord VecCoord2;
    typedef typename DataTypes2::VecDeriv VecDeriv2;
    typedef typename DataTypes2::Coord Coord2;
    typedef typename DataTypes2::Deriv Deriv2;
    typedef typename DataTypes2::MatrixDeriv MatrixDeriv2;
    typedef typename MatrixDeriv2::RowConstIterator MatrixDerivRowConstIterator2;
    typedef core::behavior::MechanicalState<DataTypes2> MechanicalState2;

//    typedef core::objectmodel::Data<VecCoord1>    DataVecCoord1;
//    typedef core::objectmodel::Data<VecDeriv1>    DataVecDeriv1;
//    typedef core::objectmodel::Data<VecCoord2>    DataVecCoord2;
//    typedef core::objectmodel::Data<VecDeriv2>    DataVecDeriv2;

    typedef typename BaseLDIPenalityContactForceField<ResponseDataTypes>::Real Real;

    typedef core::objectmodel::Data<VecCoord1>    DataVecCoord1;
    typedef core::objectmodel::Data<VecDeriv1>    DataVecDeriv1;
    typedef core::objectmodel::Data<VecCoord2>    DataVecCoord2;
    typedef core::objectmodel::Data<VecDeriv2>    DataVecDeriv2;


  protected:
    typedef std::set< unsigned int>::const_iterator VertexIterator;
    Data< float > radiusForce;
    Real Volume;
    Real K;
    TriangleModel *model[2];
    sofa::helper::vector< Vector3>  dVelocity1;
    sofa::helper::vector< Vector3>  dVelocity2;
    sofa::helper::vector< Vector3>  dVdx1;
    sofa::helper::vector< Vector3>  dVdx2;
    std::set< unsigned int>   index_used1;
    std::set< unsigned int>   index_used2;
  
  public:
    LDIPenalityContactForceField(MechanicalState1* object1, MechanicalState2* object2)
      : Inherit(object1, object2)
      , radiusForce(sofa::core::objectmodel::Base::initData(&radiusForce, 0.01f,"radiusForce", "Radius of the displayed force"))
    {
    }

    LDIPenalityContactForceField(){	}



    void clear(int reserve = 0);
    void setParameters( TriangleModel* model1, TriangleModel* model2, const Real &_Volume, 
      sofa::helper::vector< Vector3>  &_dVdx1, sofa::helper::vector< Vector3>  &_dVelocity1, const std::set< unsigned int> &_index_used1,
      sofa::helper::vector< Vector3 > &_dVdx2, sofa::helper::vector< Vector3>  &_dVelocity2, const std::set< unsigned int> &_index_used2,
      const Real _K);

    void addForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, 
      DataVecDeriv1& data_f1, DataVecDeriv2& data_f2, 
      const DataVecCoord1& data_x1, const DataVecCoord2& data_x2, 
      const DataVecDeriv1& data_v1, const DataVecDeriv2& data_v2);

    void addDForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, 
      DataVecDeriv1& data_df1, DataVecDeriv2& data_df2, 
      const DataVecDeriv1& data_dx1, const DataVecDeriv2& data_dx2);

    double getPotentialEnergy(const core::MechanicalParams* mparams /* PARAMS FIRST */, 
      const DataVecCoord1& x1, const DataVecCoord2& x2) const ;

    virtual bool useMask(){return true;}
    void draw(const core::visual::VisualParams* vparams);
  };

} // namespace forcefield

} // namespace component

} // namespace sofa

#endif
