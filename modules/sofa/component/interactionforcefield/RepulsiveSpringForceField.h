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
#ifndef SOFA_COMPONENT_INTERACTIONFORCEFIELD_REPULSIVESPRINGFORCEFIELD_H
#define SOFA_COMPONENT_INTERACTIONFORCEFIELD_REPULSIVESPRINGFORCEFIELD_H

#include <sofa/component/interactionforcefield/StiffSpringForceField.h>
#include <sofa/core/MechanicalParams.h>

namespace sofa
{

namespace component
{

namespace interactionforcefield
{

using namespace sofa::core;



template<class DataTypes>
class RepulsiveSpringForceField : public sofa::component::interactionforcefield::StiffSpringForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(RepulsiveSpringForceField,DataTypes),
            SOFA_TEMPLATE(sofa::component::interactionforcefield::StiffSpringForceField,DataTypes));

    typedef sofa::component::interactionforcefield::StiffSpringForceField<DataTypes> Inherit;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename Coord::value_type Real;
    typedef typename Inherit::Mat Mat;
    typedef typename Inherit::Spring Spring;

    typedef core::objectmodel::Data<VecDeriv>    DataVecDeriv;
    typedef core::objectmodel::Data<VecCoord>    DataVecCoord;

    enum { N = Inherit::N };

protected:
    RepulsiveSpringForceField(core::behavior::MechanicalState<DataTypes>* object1, core::behavior::MechanicalState<DataTypes>* object2)
        : sofa::component::interactionforcefield::StiffSpringForceField<DataTypes>(object1, object2)
    {
    }

    RepulsiveSpringForceField()
    {
    }
public:
    virtual void addForce(const MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& data_f1, DataVecDeriv& data_f2, const DataVecCoord& data_x1, const DataVecCoord& data_x2, const DataVecDeriv& data_v1, const DataVecDeriv& data_v2 );
    ///SOFA_DEPRECATED_ForceField <<<virtual void addForce(VecDeriv& f1, VecDeriv& f2, const VecCoord& x1, const VecCoord& x2, const VecDeriv& v1, const VecDeriv& v2);

    virtual double getPotentialEnergy(const core::MechanicalParams* /* PARAMS FIRST */, const DataVecCoord&, const DataVecCoord& ) const ;
    ///SOFA_DEPRECATED_ForceField <<<virtual double getPotentialEnergy() const;
};

} // namespace interactionforcefield

} // namespace component

} // namespace sofa

#endif
