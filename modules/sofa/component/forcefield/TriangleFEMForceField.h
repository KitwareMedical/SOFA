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
#ifndef SOFA_COMPONENT_FORCEFIELD_TRIANGLEFEMFORCEFIELD_H
#define SOFA_COMPONENT_FORCEFIELD_TRIANGLEFEMFORCEFIELD_H

#if !defined(__GNUC__) || (__GNUC__ > 3 || (_GNUC__ == 3 && __GNUC_MINOR__ > 3))
#pragma once
#endif

#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/Mat.h>


// corotational triangle from
// @InProceedings{NPF05,
//   author       = "Nesme, Matthieu and Payan, Yohan and Faure, Fran\c{c}ois",
//   title        = "Efficient, Physically Plausible Finite Elements",
//   booktitle    = "Eurographics (short papers)",
//   month        = "august",
//   year         = "2005",
//   editor       = "J. Dingliana and F. Ganovelli",
//   keywords     = "animation, physical model, elasticity, finite elements",
//   url          = "http://www-evasion.imag.fr/Publications/2005/NPF05"
// }



namespace sofa
{

namespace component
{

namespace forcefield
{

using namespace sofa::defaulttype;


/** Triangle FEM force field using the QR decomposition of the deformation gradient, inspired from http://www-evasion.imag.fr/Publications/2005/NPF05 , to handle large displacements.
  The material properties are uniform across the domain.
  Two methods are proposed, one for small displacements and one for large displacements.
  The method for small displacements has not been validated and we suspect that it is broke. Use it very carefully, and compare with the method for large displacements.
  */
template<class DataTypes>
class TriangleFEMForceField : public core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(TriangleFEMForceField, DataTypes), SOFA_TEMPLATE(core::behavior::ForceField, DataTypes));

    typedef core::behavior::ForceField<DataTypes> Inherited;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord    Coord   ;
    typedef typename DataTypes::Deriv    Deriv   ;
    typedef typename Coord::value_type   Real    ;

    typedef core::objectmodel::Data<VecCoord> DataVecCoord;
    typedef core::objectmodel::Data<VecDeriv> DataVecDeriv;

    typedef sofa::core::topology::BaseMeshTopology::index_type Index;
    typedef sofa::core::topology::BaseMeshTopology::Triangle Element;
    typedef sofa::core::topology::BaseMeshTopology::SeqTriangles VecElement;

    static const int SMALL = 1;										///< Symbol of small displacements triangle solver
    static const int LARGE = 0;										///< Symbol of large displacements triangle solver

protected:
    typedef Vec<6, Real> Displacement;								///< the displacement vector

    typedef Mat<3, 3, Real> MaterialStiffness;						///< the matrix of material stiffness
    typedef sofa::helper::vector<MaterialStiffness> VecMaterialStiffness;    ///< a vector of material stiffness matrices
    VecMaterialStiffness _materialsStiffnesses;						///< the material stiffness matrices vector

    typedef Mat<6, 3, Real> StrainDisplacement;						///< the strain-displacement matrix (the transpose, actually)
    typedef sofa::helper::vector<StrainDisplacement> VecStrainDisplacement;	///< a vector of strain-displacement matrices
    VecStrainDisplacement _strainDisplacements;						///< the strain-displacement matrices vector

    typedef Mat<3, 3, Real > Transformation;						///< matrix for rigid transformations like rotations

    /// Stiffness matrix ( = RJKJtRt  with K the Material stiffness matrix, J the strain-displacement matrix, and R the transformation matrix if any )
    typedef Mat<9, 9, Real> StiffnessMatrix;


    sofa::core::topology::BaseMeshTopology* _mesh;
    const VecElement *_indexedElements;
    Data< VecCoord > _initialPoints; ///< the intial positions of the points


    TriangleFEMForceField();
    virtual ~TriangleFEMForceField();



public:

    virtual void init();
    virtual void reinit();
    virtual void addForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v);
    virtual void addDForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& df, const DataVecDeriv& dx);

    void draw(const core::visual::VisualParams* vparams);

    int method;
    Data<std::string> f_method; ///< Choice of method: 0 for small, 1 for large displacements
    Data<Real> f_poisson;       ///< Poisson ratio of the material
    Data<Real> f_young;         ///< Young modulus of the material
    Data<Real> f_thickness;     ///< Thickness of the elements
//    Data<Real> f_damping;       ///< Damping coefficient of the material, currently unused
    Data<bool> f_planeStrain; ///< compute material stiffness corresponding to the plane strain assumption, or to the plane stress otherwise.

    Real getPoisson() { return f_poisson.getValue(); }
    void setPoisson(Real val) { f_poisson.setValue(val); }
    Real getYoung() { return f_young.getValue(); }
    void setYoung(Real val) { f_young.setValue(val); }
//    Real getDamping() { return f_damping.getValue(); }
//    void setDamping(Real val) { f_damping.setValue(val); }
    int  getMethod() { return method; }
    void setMethod(int val) { method = val; }

protected :

    /// f += Kx where K is the stiffness matrix and x a displacement
    virtual void applyStiffness( VecCoord& f, Real h, const VecCoord& x, const double &kFactor );
    void computeStrainDisplacement( StrainDisplacement &J, Coord a, Coord b, Coord c);
    void computeMaterialStiffnesses();
    void computeForce( Displacement &F, const Displacement &Depl, const MaterialStiffness &K, const StrainDisplacement &J );

    ////////////// small displacements method
    void initSmall();
    void accumulateForceSmall( VecCoord& f, const VecCoord & p, Index elementIndex, bool implicit = false );
//    void accumulateDampingSmall( VecCoord& f, Index elementIndex );
    void applyStiffnessSmall( VecCoord& f, Real h, const VecCoord& x, const double &kFactor );

    ////////////// large displacements method
    sofa::helper::vector< helper::fixed_array <Coord, 3> > _rotatedInitialElements;   ///< The initials positions in its frame
    sofa::helper::vector< Transformation > _rotations;
    void initLarge();
    void computeRotationLarge( Transformation &r, const VecCoord &p, const Index &a, const Index &b, const Index &c);
    void accumulateForceLarge( VecCoord& f, const VecCoord & p, Index elementIndex, bool implicit=false );
//    void accumulateDampingLarge( VecCoord& f, Index elementIndex );
    void applyStiffnessLarge( VecCoord& f, Real h, const VecCoord& x, const double &kFactor );

    //// stiffness matrix assembly
    void computeElementStiffnessMatrix( StiffnessMatrix& S, StiffnessMatrix& SR, const MaterialStiffness &K, const StrainDisplacement &J, const Transformation& Rot );
    void addKToMatrix(sofa::defaulttype::BaseMatrix *mat, SReal k, unsigned int &offset); // compute and add all the element stiffnesses to the global stiffness matrix

};

#ifndef SOFA_FLOAT
using sofa::defaulttype::Vec3dTypes;
#endif

#ifndef SOFA_DOUBLE
using sofa::defaulttype::Vec3fTypes;
#endif

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_FORCEFIELD_TRIANGLEFEMFORCEFIELD_CPP)

#ifndef SOFA_FLOAT
extern template class SOFA_SIMPLE_FEM_API TriangleFEMForceField<Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_SIMPLE_FEM_API TriangleFEMForceField<Vec3fTypes>;
#endif

#endif // defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_FORCEFIELD_TRIANGLEFEMFORCEFIELD_CPP)

} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FORCEFIELD_TRIANGLEFEMFORCEFIELD_H
