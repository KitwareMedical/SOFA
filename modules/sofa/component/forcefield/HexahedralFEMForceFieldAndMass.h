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
#ifndef SOFA_COMPONENT_FORCEFIELD_HEXAHEDRALFEMFORCEFIELDANDMASS_H
#define SOFA_COMPONENT_FORCEFIELD_HEXAHEDRALFEMFORCEFIELDANDMASS_H


#include <sofa/component/forcefield/HexahedralFEMForceField.h>
#include <sofa/core/behavior/Mass.h>

#include <sofa/component/topology/TopologyData.h>

namespace sofa
{

namespace component
{

namespace forcefield
{

using namespace sofa::defaulttype;
using sofa::helper::vector;
using sofa::core::behavior::Mass;

/** Compute Finite Element forces based on hexahedral elements including continuum mass matrices
*/
template<class DataTypes>
class HexahedralFEMForceFieldAndMass : virtual public Mass<DataTypes>, virtual public HexahedralFEMForceField<DataTypes>
{
public:
    SOFA_CLASS2(SOFA_TEMPLATE(HexahedralFEMForceFieldAndMass,DataTypes), SOFA_TEMPLATE(Mass,DataTypes), SOFA_TEMPLATE(HexahedralFEMForceField,DataTypes));

    typedef HexahedralFEMForceField<DataTypes> HexahedralFEMForceFieldT;
    typedef Mass<DataTypes> MassT;

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef VecCoord Vector;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename Coord::value_type Real;
    typedef core::objectmodel::Data<VecDeriv> DataVecDeriv;
    typedef core::objectmodel::Data<VecCoord> DataVecCoord;
    typedef typename HexahedralFEMForceFieldT::Mat33 Mat33;
    typedef typename HexahedralFEMForceFieldT::Displacement Displacement;
    typedef typename HexahedralFEMForceFieldT::VecElement VecElement;
    typedef typename HexahedralFEMForceFieldT::VecElementStiffness VecElementMass;
    typedef typename HexahedralFEMForceFieldT::ElementStiffness ElementMass;
    typedef core::topology::BaseMeshTopology::index_type Index;
    typedef typename HexahedralFEMForceFieldT::HexahedronInformation HexahedronInformation;
    typedef typename HexahedralFEMForceFieldT::ElementStiffness ElementStiffness;
    typedef typename HexahedralFEMForceFieldT::Element Element;


protected:
    HexahedralFEMForceFieldAndMass();
public:
    virtual void init( );
    virtual void reinit( );

    // -- Mass interface
    virtual  void addMDx(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& f, const DataVecDeriv& dx, double factor);

    ///// WARNING this method only add diagonal elements in the given matrix !
    // virtual void addMToMatrix(defaulttype::BaseMatrix * matrix, double mFact, unsigned int &offset);
    virtual void addMToMatrix(const core::MechanicalParams* mparams /* PARAMS FIRST */, const sofa::core::behavior::MultiMatrixAccessor* matrix);

    ///// WARNING this method only add diagonal elements in the given matrix !
    // virtual void addKToMatrix(sofa::defaulttype::BaseMatrix *mat, SReal k, unsigned int &offset);
    virtual void addKToMatrix(const core::MechanicalParams* mparams /* PARAMS FIRST */, const sofa::core::behavior::MultiMatrixAccessor* matrix);

    ///// WARNING this method only add diagonal elements in the given matrix !
    // virtual void addMBKToMatrix(sofa::defaulttype::BaseMatrix * matrix,
    // double mFact, double bFact, double kFact, unsigned int &offset);
    virtual void addMBKToMatrix(const core::MechanicalParams* mparams /* PARAMS FIRST */, const sofa::core::behavior::MultiMatrixAccessor* matrix);

    virtual  void accFromF(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& a, const DataVecDeriv& f);

    virtual  void addForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v);

    virtual double getKineticEnergy(const core::MechanicalParams* /* mparams */ /* PARAMS FIRST */, const DataVecDeriv& /*v*/)  const ///< vMv/2 using dof->getV()
    {serr<<"HexahedralFEMForceFieldAndMass<DataTypes>::getKineticEnergy not yet implemented"<<sendl; return 0;}

    virtual void addDForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& df, const DataVecDeriv& dx);
    // virtual void addDForce(VecDeriv& df, const VecDeriv& dx);
    // virtual void addDForce(VecDeriv& df, const VecDeriv& dx, double kFactor, double);

    virtual void addGravityToV(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& d_v);

    virtual void draw(const core::visual::VisualParams* vparams);

    double getElementMass(unsigned int index);

    void setDensity(Real d) {_density.setValue( d );}
    Real getDensity() {return _density.getValue();}




protected:
    virtual void computeElementMasses( ); ///< compute the mass matrices
    Real integrateVolume( int signx, int signy, int signz, Real l0, Real l1, Real l2 );
    virtual void computeElementMass( ElementMass &Mass, Real& totalMass, const helper::fixed_array<Coord,8> &nodes); ///< compute the mass matrix of an element

    void computeParticleMasses();

    void computeLumpedMasses();

protected:
    //HFFHexahedronHandler* hexahedronHandler;

    Data<Real> _density;
    Data<bool> _useLumpedMass;

    HexahedronData<sofa::helper::vector<ElementMass> > _elementMasses; ///< mass matrices per element
    HexahedronData<sofa::helper::vector<Real> > _elementTotalMass; ///< total mass per element

    PointData<sofa::helper::vector<Real> > _particleMasses; ///< masses per particle in order to compute gravity
    PointData<sofa::helper::vector<Coord> > _lumpedMasses; ///< masses per particle computed by lumping mass matrices
};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_FORCEFIELD_HEXAHEDRALFEMFORCEFIELDANDMASS_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_SIMPLE_FEM_API HexahedralFEMForceFieldAndMass<defaulttype::Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_SIMPLE_FEM_API HexahedralFEMForceFieldAndMass<defaulttype::Vec3fTypes>;
#endif
#endif

} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FORCEFIELD_HEXAHEDRALFEMFORCEFIELDANDMASS_H
