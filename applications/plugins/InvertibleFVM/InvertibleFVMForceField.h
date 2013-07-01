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
#ifndef SOFA_COMPONENT_FORCEFIELD_InvertibleFVMForceField_H
#define SOFA_COMPONENT_FORCEFIELD_InvertibleFVMForceField_H

#include "initPlugin.h"
#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/helper/vector.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/Mat.h>
//#include <sofa/component/component.h>
//#include <sofa/helper/OptionsGroup.h>



namespace sofa
{

namespace component
{

namespace forcefield
{

using namespace sofa::defaulttype;
using sofa::helper::vector;
using namespace sofa::core::topology;

template<class DataTypes>
class InvertibleFVMForceField;

/// This class can be overridden if needed for additionnal storage within template specializations.
template<class DataTypes>
class InvertibleFVMForceFieldInternalData
{
public:
};


/** Compute Finite Volume forces based on tetrahedral and hexahedral elements.
 * implementation of "invertible FEM..."
*/
template<class DataTypes>
class InvertibleFVMForceField : virtual public core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(InvertibleFVMForceField, DataTypes), SOFA_TEMPLATE(core::behavior::ForceField, DataTypes));

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::VecReal VecReal;
    typedef VecCoord Vector;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename Coord::value_type Real;

    typedef core::objectmodel::Data<VecDeriv>    DataVecDeriv;
    typedef core::objectmodel::Data<VecCoord>    DataVecCoord;

    typedef core::topology::BaseMeshTopology::index_type Index;
    typedef core::topology::BaseMeshTopology::Tetra Tetra;
    typedef core::topology::BaseMeshTopology::SeqTetrahedra VecTetra;



protected:

    /// @name Per tetrahedron data
    /// @{

    /// Displacement vector (deformation of the 4 corners of a tetrahedron)
    typedef VecNoInit<12, Real> Displacement;

    /// Rigid transformation (rotation) matrix
    typedef MatNoInit<3, 3, Real> Transformation;

    /// @}

    vector<Transformation> _rotationsU;
    vector<Transformation> _rotationsV;



    /* typedef std::pair<int,Real> Col_Value;
     typedef vector< Col_Value > CompressedValue;
     typedef vector< CompressedValue > CompressedMatrix;
     CompressedMatrix _stiffnesses;
     double m_potentialEnergy;*/



    core::topology::BaseMeshTopology* _mesh;
    const VecTetra *_indexedTetra;

    vector<Transformation> _initialTransformation;
    vector<Transformation> _initialRotation;

    vector<Transformation> _U;
    vector<Transformation> _V;
    vector<Vec<3,Coord> > _b;

    InvertibleFVMForceFieldInternalData<DataTypes> data;
    friend class InvertibleFVMForceFieldInternalData<DataTypes>;

public:


    Data< VecCoord > _initialPoints; ///< the intial positions of the points

    Data<Real> _poissonRatio;
    Data<VecReal > _youngModulus;
    Data<VecReal> _localStiffnessFactor;


    Data< bool > drawHeterogeneousTetra;
    Data< bool > drawAsEdges;

    Data< bool > _verbose;

    Real minYoung;
    Real maxYoung;
protected:
    InvertibleFVMForceField()
        : _mesh(NULL)
        , _indexedTetra(NULL)
        , _initialPoints(initData(&_initialPoints, "initialPoints", "Initial Position"))
        , _poissonRatio(initData(&_poissonRatio,(Real)0.45f,"poissonRatio","FEM Poisson Ratio [0,0.5["))
        , _youngModulus(initData(&_youngModulus,"youngModulus","FEM Young Modulus"))
        , _localStiffnessFactor(initData(&_localStiffnessFactor, "localStiffnessFactor","Allow specification of different stiffness per element. If there are N element and M values are specified, the youngModulus factor for element i would be localStiffnessFactor[i*M/N]"))
        , drawHeterogeneousTetra(initData(&drawHeterogeneousTetra,false,"drawHeterogeneousTetra","Draw Heterogeneous Tetra in different color"))
        , drawAsEdges(initData(&drawAsEdges,false,"drawAsEdges","Draw as edges instead of tetrahedra"))
        , _verbose(initData(&_verbose,false,"verbose","Print debug stuff"))
    {
        minYoung = 0.0;
        maxYoung = 0.0;
    }

    virtual ~InvertibleFVMForceField() {}

public:

    void setPoissonRatio(Real val) { this->_poissonRatio.setValue(val); }

    void setYoungModulus(Real val)
    {
        VecReal newY;
        newY.resize(1);
        newY[0] = val;
        _youngModulus.setValue(newY);
    }


    virtual void reset();
    virtual void init();
    virtual void reinit();

    virtual void addForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& d_f, const DataVecCoord& d_x, const DataVecDeriv& d_v);
    virtual void addDForce(const core::MechanicalParams*  /* PARAMS FIRST */, DataVecDeriv& , const DataVecDeriv& );

    virtual void addKToMatrix(sofa::defaulttype::BaseMatrix *m, SReal kFactor, unsigned int &offset);

    void draw(const core::visual::VisualParams* vparams);



};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_FORCEFIELD_InvertibleFVMForceField_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_InvertibleFVM_API InvertibleFVMForceField<defaulttype::Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_InvertibleFVM_API InvertibleFVMForceField<defaulttype::Vec3fTypes>;
#endif
#endif

} // namespace forcefield

} // namespace component

} // namespace sofa

#endif
