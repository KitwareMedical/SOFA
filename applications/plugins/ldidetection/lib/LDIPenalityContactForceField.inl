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
#ifndef SOFA_COMPONENT_FORCEFIELD_LDIPENALITYCONTACTFORCEFIELD_INL
#define SOFA_COMPONENT_FORCEFIELD_LDIPENALITYCONTACTFORCEFIELD_INL

#include <sofa/core/behavior/ForceField.inl>
// #include <sofa/component/forcefield/LDIPenalityContactForceField.h>
#include "LDIPenalityContactForceField.h"
#include <sofa/helper/system/config.h>
#include <assert.h>
// #include <GL/gl.h>
#include <sofa/helper/gl/template.h>
#include <iostream>
#include <sofa/helper/gl/Axis.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/helper/accessor.h>
#include <sofa/core/visual/VisualParams.h>

//#define DEBUG_LDIPENALTYCONTACTFORCEFIELD
#ifdef  DEBUG_LDIPENALTYCONTACTFORCEFIELD
#define DEBUG_LDIPENALTY_OUT(c) \
 this->f_printLog.setValue(sofa::core::ExecParams::defaultInstance(),true); \
  c
#else
#define DEBUG_LDIPENALTY_OUT(c)
#endif

namespace sofa
{

namespace component
{

namespace forcefield
{


  template<class DataTypes1,class DataTypes2, class ResponseDataTypes>
  void LDIPenalityContactForceField<DataTypes1,DataTypes2, ResponseDataTypes>::setParameters( TriangleModel* model1, TriangleModel* model2, const Real &_Volume, 
                                                                                             sofa::helper::vector< Vector3>  &_dVdx1, sofa::helper::vector< Vector3>  &_dVelocity1, const std::set< unsigned int> &_index_used1,
                                                                                             sofa::helper::vector< Vector3 > &_dVdx2, sofa::helper::vector< Vector3>  &_dVelocity2, const std::set< unsigned int> &_index_used2,
                                                                                             const Real _K)
  {
    using namespace core;
    using namespace core::objectmodel;
    model[0] = model1; model[1] = model2;
    Volume = _Volume;
    dVdx1 = _dVdx1; dVelocity1 = _dVelocity1; index_used1 = _index_used1;
    dVdx2 = _dVdx2; dVelocity2 = _dVelocity2; index_used2 = _index_used2;
    K = _K;
  }
  template<class DataTypes1,class DataTypes2, class ResponseDataTypes>
  void LDIPenalityContactForceField<DataTypes1,DataTypes2, ResponseDataTypes>::clear(int /*reserve*/)
  {
  }

  template<class DataTypes1,class DataTypes2, class ResponseDataTypes>
  void LDIPenalityContactForceField<DataTypes1,DataTypes2,ResponseDataTypes>::addForce(const core::MechanicalParams* /*mparams */, 
    DataVecDeriv1& data_f1, DataVecDeriv2& data_f2, 
    const DataVecCoord1& /*data_x1*/, const DataVecCoord2& /*data_x2*/, 
    const DataVecDeriv1& /*data_v1*/, const DataVecDeriv2& /*data_v2*/)
  {
    using namespace core;
    using namespace core::objectmodel;

    VecDeriv1& f1 = *data_f1.beginEdit();
    VecDeriv2& f2 = *data_f2.beginEdit();

    if (f1.size() == 0)
    {
      f1.resize(model[0]->getMechanicalState()->getX()->size());
    }
    if (f2.size() == 0)
    {
      f2.resize(model[1]->getMechanicalState()->getX()->size());
    }
    DEBUG_LDIPENALTY_OUT(this->sout << "Volume= " << Volume << this->sendl);
    const Real K_Volume = (-K*Volume);

   
    if( model[0]->isSimulated() )
    {
      for (VertexIterator it=index_used1.begin(); it!=index_used1.end();++it)
      {

        f1[*it] += (dVdx1[*it] * K_Volume); /*+ dVelocity1[*it]*/;
      }
    }
   
    if ( model[1]->isSimulated() )
    {
      for (VertexIterator it=index_used2.begin(); it!=index_used2.end();++it)
      {
        f2[*it] += (dVdx2[*it] * K_Volume) /*+ dVelocity2[*it]*/ ;
      }
    }
   

  
    data_f1.endEdit();
    data_f2.endEdit();
  }

  template<class DataTypes1,class DataTypes2, class ResponseDataTypes>
  void LDIPenalityContactForceField<DataTypes1,DataTypes2,ResponseDataTypes>::addDForce(const core::MechanicalParams* mparams , 
    DataVecDeriv1& data_df1, DataVecDeriv2& data_df2, 
    const DataVecDeriv1& data_dx1, const DataVecDeriv2& data_dx2)
  {
    VecDeriv1& df1 = *data_df1.beginEdit(mparams);
    VecDeriv2& df2 = *data_df2.beginEdit(mparams);
    const VecDeriv1& dx1 = data_dx1.getValue(mparams);
    const VecDeriv2& dx2 = data_dx2.getValue(mparams);
    double kfactor = mparams->kFactor();
    // Volume variation
    SReal deltaV = SReal();
 
    if( model[0]->isSimulated() )
    {
      for (VertexIterator it=index_used1.begin(); it!=index_used1.end() ;++it)
      {
        deltaV += dVdx1[(*it)] * dx1[(*it)];
        DEBUG_LDIPENALTY_OUT(this->sout << "dVdx1_" << *it << "= " << dVdx1[*it] << this->sendl);
        DEBUG_LDIPENALTY_OUT(this->sout << "dx1_" << *it << "= " << dx1[*it] << this->sendl);
      }
    }

    if( model[1]->isSimulated() )
    {
      for (VertexIterator it=index_used2.begin(); it!=index_used2.end() ;++it)
      {
        deltaV += dVdx2[(*it)] * dx2[(*it)];
        DEBUG_LDIPENALTY_OUT(this->sout << "dVdx2_" << *it << " = " << dVdx2[*it] << this->sendl);
        DEBUG_LDIPENALTY_OUT(this->sout << "dx2_" << *it << " = " << dx2[*it] << this->sendl);
      }
    }
  
    // Corresponding pressure variation
    SReal deltaP = -K*deltaV*kfactor;
    DEBUG_LDIPENALTY_OUT(this->sout << "deltaV=" << deltaV << this->sendl);
    DEBUG_LDIPENALTY_OUT(this->sout << "deltaP=" << deltaP << this->sendl);


    // Accumulate the resulting force variation

    if( model[0]->isSimulated() )
    {
      for (VertexIterator it=index_used1.begin(); it!=index_used1.end() ;++it)
      {
        df1[(*it)] += dVdx1[(*it)] * deltaP;
        DEBUG_LDIPENALTY_OUT(this->sout << "df1_" << *it << " = " << df1[*it] << this->sendl);
      }
    }
 
    if( model[1]->isSimulated() )
    {
      for (VertexIterator it=index_used2.begin(); it!=index_used2.end() ;++it)
      {
        df2[(*it)] += dVdx2[(*it)]*deltaP;
        DEBUG_LDIPENALTY_OUT(this->sout << "df2_" << *it << " = " << df2[*it] << this->sendl);
      }
    }
   

    data_df1.endEdit();
    data_df2.endEdit();

  }

 
  

  template<class DataTypes1,class DataTypes2, class ResponseDataTypes>
  void LDIPenalityContactForceField<DataTypes1,DataTypes2, ResponseDataTypes>::draw(const core::visual::VisualParams* vparams)
  {

    if (!((static_cast<core::behavior::BaseMechanicalState*>(this->mstate1) == static_cast<core::behavior::BaseMechanicalState*>(this->mstate2))?vparams->displayFlags().getShowForceFields():vparams->displayFlags().getShowInteractionForceFields())) return;

    const VecCoord1& p1 = *this->mstate1->getX();
    const VecCoord2& p2 = *this->mstate2->getX();
 
    
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_COLOR_MATERIAL);

    bool autoCollision = ((void *)this->mstate1 == (void *)this->mstate2);
    if (autoCollision)
      glColor4f(1.0f,1.0f,0.0f,1.0f);
    else
      glColor4f(1.0f,0.2f,0.2f,1.0f);

   
    Real K_Volume = (-K*Volume);


    for (VertexIterator it=index_used1.begin(); it!=index_used1.end() ;++it)
    {
      helper::gl::Axis::draw(p1[(*it)], p1[(*it)]+ dVdx1[(*it)]*K_Volume,radiusForce.getValue());
      // 	helper::gl::Axis::draw(p1[(*it)], p1[(*it)]+(dVdx[0][(*it)])*K_Volume,Volume*0.5);
    }

  
    for (VertexIterator it=index_used2.begin(); it!=index_used2.end() ;++it)
    {
      helper::gl::Axis::draw(p2[(*it)], p2[(*it)]+(dVdx2[(*it)])*K_Volume,radiusForce.getValue());
      // 	helper::gl::Axis::draw(p2[(*it)], p2[(*it)]+(dVdx[1][(*it)])*K_Volume,Volume*0.5);
    }

    glDisable(GL_LIGHTING);
    glDisable(GL_COLOR_MATERIAL);

  }

  template<class DataTypes1,class DataTypes2, class ResponseDataTypes>
  double LDIPenalityContactForceField<DataTypes1,DataTypes2,ResponseDataTypes>::getPotentialEnergy(const core::MechanicalParams* /*mparams */, const DataVecCoord1& /* x1*/ , const DataVecCoord2& /*x2*/ ) const
  {
    
    double e =  K * Volume * Volume * 0.5;
    return e;

  }


} // namespace forcefield

} // namespace component

} // namespace sofa

#endif
