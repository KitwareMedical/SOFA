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
*                               SOFA :: Plugins                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/

#ifndef SOFA_COMPONENT_CONSTRAINT_MyMappingPendulumInPlane_H
#define SOFA_COMPONENT_CONSTRAINT_MyMappingPendulumInPlane_H


#include <sofa/core/Mapping.h>
#include <sofa/component/component.h>
#include <sofa/helper/OptionsGroup.h>
#include <sofa/defaulttype/Vec.h>




namespace sofa
{

namespace component
{

namespace mapping
{
using helper::vector;
using defaulttype::Vec;

/** input: pendulum angle; output: coordinates of the endpoint of the pendulum
  */

template <class TIn, class TOut>
class MyMappingPendulumInPlane : public core::Mapping<TIn, TOut>
{
public:
    SOFA_CLASS( SOFA_TEMPLATE2(MyMappingPendulumInPlane,TIn,TOut), SOFA_TEMPLATE2(core::Mapping,TIn,TOut) );
    typedef core::Mapping<TIn, TOut> Inherit;
    typedef TIn In;
    typedef TOut Out;
    typedef typename In::Real InReal;
    typedef typename In::VecCoord VecInCoord;
    typedef typename In::Deriv InDeriv;
    typedef typename In::VecDeriv VecInDeriv;
    typedef typename In::MatrixDeriv MatrixInDeriv;
    typedef typename Out::Real OutReal;
    typedef typename Out::VecCoord VecOutCoord;
    typedef typename Out::Deriv OutDeriv;
    typedef typename Out::VecDeriv VecOutDeriv;
    typedef typename Out::MatrixDeriv MatrixOutDeriv;
protected:
    MyMappingPendulumInPlane();
    ~MyMappingPendulumInPlane();
public:
    Data<vector<OutReal> > f_length;

    virtual void init();
    virtual void draw(const core::visual::VisualParams*);

    virtual void apply(VecOutCoord& out, const VecInCoord& in);
    virtual void applyJ( VecOutDeriv& out, const VecInDeriv& in);
    virtual void applyJT( VecInDeriv& out, const VecOutDeriv& in);
    virtual void applyJT( MatrixInDeriv& out, const MatrixOutDeriv& in);
    virtual void applyDJT(const core::MechanicalParams* mparams /* PARAMS FIRST  = core::MechanicalParams::defaultInstance()*/, core::MultiVecDerivId parentForceChange, core::ConstMultiVecDerivId );


protected:
    typedef Vec<2,OutReal> Vec2;
    vector<Vec2> gap;

private:

};


}

}

}



#endif
