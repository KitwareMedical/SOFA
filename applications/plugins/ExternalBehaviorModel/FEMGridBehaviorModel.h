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
#ifndef SOFA_EXTERNALBEHAVIORMODEL_FEMGRIDBEHAVIORMODEL_H
#define SOFA_EXTERNALBEHAVIORMODEL_FEMGRIDBEHAVIORMODEL_H


// plugin includes
#include "initExternalBehaviorModel.h"

// SOFA includes
#include <sofa/component/misc/InteractingBehaviorModel.h>


// internal stuff, here it is using SOFA components that could be replaced by any library
#include <sofa/component/topology/MeshTopology.h>
#include <sofa/component/forcefield/HexahedronFEMForceFieldAndMass.h>
#include <sofa/simulation/tree/GNode.h>


namespace sofa
{

namespace externalBehaviorModel
{


/// Demo of how to implement a InteractingBehaviorModel to plug an external behavior model in a SOFA scene



/// In this simple example, the input SOFA dofs are an hexahedra (only 8 vertices) that are link by an elastic behavior model.
/// The behavior model is a fine FEM grid. All the "inside" dofs are invisible from the SOFA side, they are only visible inside this component.
/// All these dofs are independant, but on the SOFA side, the independant dofs are only the 8 hexahedron corners.

/// @warning Note that the internal behavior model is based on existing sofa components for simplicity,
/// but it could (should?) be built on any independant code or existing library

/// @warning Only a subset of the complete API is implemented, allowing only not assembled, implicit system solving.




// DataTypes describes the dof type
// here it is compiled for Vec3 corresponding to particles
template<class DataTypes>
class FEMGridBehaviorModel : public component::misc::InteractingBehaviorModel<DataTypes>
{

public:

    // SOFA black magic
    SOFA_CLASS( SOFA_TEMPLATE(FEMGridBehaviorModel, DataTypes), SOFA_TEMPLATE(component::misc::InteractingBehaviorModel, DataTypes) );

    typedef component::misc::InteractingBehaviorModel<DataTypes> Inherited;
    typedef typename Inherited::Dofs Dofs;


    // retrieve interesting types from DataTypes
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord; // encoding a dof position
    typedef typename DataTypes::Deriv Deriv; // encoding a dof variation of position or a force
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef Data<typename DataTypes::VecCoord> DataVecCoord;
    typedef Data<typename DataTypes::VecDeriv> DataVecDeriv;
    typedef typename defaulttype::Vec<3,Real> Vec3;




    /// @name Component API
    /// @{
    virtual void init(); /// call when initializing the simulation
    virtual void draw( const core::visual::VisualParams* vparams ); /// debug drawing
    virtual void handleEvent(sofa::core::objectmodel::Event *event);
    /// @}


    /// @name Internal force API
    /// @{
    /// f += K( x, v ) -> build the right part of the sytem (including gravity)
    virtual void addForce( const core::MechanicalParams* mparams, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v );
    /// df += K.dx -> call at each iteration by unassembled system solvers
    /// necessary to perform an implicit interaction
    virtual void addDForce( const core::MechanicalParams* mparams, DataVecDeriv& df, const DataVecDeriv& dx );
    /// @}


    /// @name Mass API
    /// @{
    /// f += factor M dx (for implicit solvers)
    virtual void addMDx( const core::MechanicalParams* mparams, DataVecDeriv& f, const DataVecDeriv& dx, double factor );
    /// @}



    // Data fields will automatically appear in qt-based GUI and in read/write XML scene files
    // in this simple example, there is only one young modulus for all elements and one similar mass to every particle
    // Note that Data must be initialized in the constructor, giving the default value, the field name and a description
    Data<Real> _youngModulus;
    Data<Real> _density;
    Data<unsigned> _subdivisions;




protected:

    FEMGridBehaviorModel();

    virtual ~FEMGridBehaviorModel() {}




    /// @name internal sofa stuff that could be replaced by any library
    /// @{
    typename Dofs::SPtr m_internalDofs;
    component::topology::RegularGridTopology::SPtr m_internalTopology;
    typename component::forcefield::HexahedronFEMForceFieldAndMass<DataTypes>::SPtr m_internalForceFieldAndMass;
    simulation::tree::GNode::SPtr m_internalNode;
    int mapExposedInternalIndices[8]; ///< identity mapping between exposed SOFA dofs and internal model dofs
    /// @}





}; // class FEMGridBehaviorModel




#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_EXTERNALBEHAVIORMODEL_FEMGRIDBEHAVIORMODEL_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_ExternalBehaviorModel_API FEMGridBehaviorModel<defaulttype::Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_ExternalBehaviorModel_API FEMGridBehaviorModel<defaulttype::Vec3fTypes>;
#endif
#endif



} // namespace externalBehaviorModel

} // namespace sofa



#endif // SOFA_EXTERNALBEHAVIORMODEL_FEMGRIDBEHAVIORMODEL_H
