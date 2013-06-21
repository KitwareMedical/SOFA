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
#ifndef SOFA_COMPONENT_FORCEFIELD_TRIANGULARBIQUADRATICSPRINGSFORCEFIELD_H
#define SOFA_COMPONENT_FORCEFIELD_TRIANGULARBIQUADRATICSPRINGSFORCEFIELD_H

#if !defined(__GNUC__) || (__GNUC__ > 3 || (_GNUC__ == 3 && __GNUC_MINOR__ > 3))
#pragma once
#endif

#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/Mat.h>
#include <sofa/component/topology/TopologyData.h>


namespace sofa
{

namespace component
{


namespace forcefield
{
using namespace sofa::helper;
using namespace sofa::defaulttype;
using namespace sofa::component::topology;


template<class DataTypes>
class TriangularBiquadraticSpringsForceField : public core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(TriangularBiquadraticSpringsForceField, DataTypes), SOFA_TEMPLATE(core::behavior::ForceField, DataTypes));

    typedef core::behavior::ForceField<DataTypes> Inherited;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord    Coord   ;
    typedef typename DataTypes::Deriv    Deriv   ;
    typedef typename Coord::value_type   Real    ;

    typedef core::objectmodel::Data<VecCoord> DataVecCoord;
    typedef core::objectmodel::Data<VecDeriv> DataVecDeriv;

    class Mat3 : public fixed_array<Deriv,3>
    {
    public:
        Deriv operator*(const Deriv& v)
        {
            return Deriv((*this)[0]*v,(*this)[1]*v,(*this)[2]*v);
        }
        Deriv transposeMultiply(const Deriv& v)
        {
            return Deriv(v[0]*((*this)[0])[0]+v[1]*((*this)[1])[0]+v[2]*((*this)[2][0]),
                    v[0]*((*this)[0][1])+v[1]*((*this)[1][1])+v[2]*((*this)[2][1]),
                    v[0]*((*this)[0][2])+v[1]*((*this)[1][2])+v[2]*((*this)[2][2]));
        }
    };

protected:


    class EdgeRestInformation
    {
    public:
        Real  restSquareLength;	// the rest length
        Real  currentSquareLength; 	// the current edge length
        Real  deltaL2;  // the current unit direction
        Real stiffness;

        EdgeRestInformation()
        {
        }

        /// Output stream
        inline friend std::ostream& operator<< ( std::ostream& os, const EdgeRestInformation& /*eri*/ )
        {
            return os;
        }

        /// Input stream
        inline friend std::istream& operator>> ( std::istream& in, EdgeRestInformation& /*eri*/ )
        {
            return in;
        }
    };

    class TriangleRestInformation
    {
    public:
        Real  gamma[3];	// the angular stiffness
        Real stiffness[3]; // the elongation stiffness
        Mat3 DfDx[3]; /// the edge stiffness matrix

        Coord currentNormal;
        Coord lastValidNormal;
        Real area;
        Real restArea;
        Coord areaVector[3];
        Deriv dp[3];
        Real J;

        TriangleRestInformation()
        {
        }
        /// Output stream
        inline friend std::ostream& operator<< ( std::ostream& os, const TriangleRestInformation& /*tri*/ )
        {
            return os;
        }

        /// Input stream
        inline friend std::istream& operator>> ( std::istream& in, TriangleRestInformation& /*vec*/ )
        {
            return in;
        }
    };

    TriangleData<helper::vector<TriangleRestInformation> > triangleInfo;
    EdgeData<helper::vector<EdgeRestInformation> > edgeInfo;

    sofa::core::topology::BaseMeshTopology* _topology;
    Data < VecCoord >  _initialPoints;										///< the intial positions of the points

    bool updateMatrix;

    Data<Real> f_poissonRatio;
    Data<Real> f_youngModulus;
    Data<Real> f_dampingRatio;
    Data<bool> f_useAngularSprings; // whether angular springs should be included

    Data<bool> f_compressible; // whether the material is compressible or not
    /**** coefficient that controls how the material can cope with very compressible cases
    must be between 0 and 1 : if 0 then the deformation may diverge for large compression
    if 1 then the material can undergo large compression even inverse elements ***/
    Data<Real> f_stiffnessMatrixRegularizationWeight;

    Real lambda;  /// first Lam� coefficient
    Real mu;    /// second Lam� coefficient


    TriangularBiquadraticSpringsForceField();

    virtual ~TriangularBiquadraticSpringsForceField();
public:
    virtual void init();

    virtual void addForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& d_f, const DataVecCoord& d_x, const DataVecDeriv& d_v);
    virtual void addDForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& d_df, const DataVecDeriv& d_dx);

    virtual Real getLambda() const { return lambda;}
    virtual Real getMu() const { return mu;}

    void setYoungModulus(const double modulus)
    {
        f_youngModulus.setValue((Real)modulus);
    }
    void setPoissonRatio(const double ratio)
    {
        f_poissonRatio.setValue((Real)ratio);
    }

    void draw(const core::visual::VisualParams* vparams);
    /// compute lambda and mu based on the Young modulus and Poisson ratio
    void updateLameCoefficients();

    class TRBSEdgeHandler : public TopologyDataHandler<Edge,sofa::helper::vector<EdgeRestInformation> >
    {
    public:
        typedef typename TriangularBiquadraticSpringsForceField<DataTypes>::EdgeRestInformation EdgeRestInformation;

        TRBSEdgeHandler(TriangularBiquadraticSpringsForceField<DataTypes>* ff,
                EdgeData<sofa::helper::vector<EdgeRestInformation> >* data)
            :TopologyDataHandler<Edge,sofa::helper::vector<EdgeRestInformation> >(data)
            ,ff(ff)
        {
        }
        void applyCreateFunction(unsigned int, EdgeRestInformation &t, const Edge &,
                const sofa::helper::vector<unsigned int> &, const sofa::helper::vector<double> &);

    protected:
        TriangularBiquadraticSpringsForceField<DataTypes>* ff;
    };

    class TRBSTriangleHandler : public TopologyDataHandler<Triangle,sofa::helper::vector<TriangleRestInformation> >
    {
    public:
        typedef typename TriangularBiquadraticSpringsForceField<DataTypes>::TriangleRestInformation TriangleRestInformation;

        TRBSTriangleHandler(TriangularBiquadraticSpringsForceField<DataTypes>* ff,
                TriangleData<sofa::helper::vector<TriangleRestInformation> >* data)
            :TopologyDataHandler<Triangle,sofa::helper::vector<TriangleRestInformation> >(data)
            ,ff(ff)
        {
        }

        void applyCreateFunction(unsigned int, TriangleRestInformation &t,
                const Triangle &,
                const sofa::helper::vector<unsigned int> &,
                const sofa::helper::vector<double> &);
        void applyDestroyFunction(unsigned int, TriangleRestInformation &);


    protected:
        TriangularBiquadraticSpringsForceField<DataTypes>* ff;
    };



protected :
    TRBSEdgeHandler* edgeHandler;
    TRBSTriangleHandler* triangleHandler;

    EdgeData<helper::vector<EdgeRestInformation> > &getEdgeInfo() {return edgeInfo;}

};

#ifndef SOFA_FLOAT
using sofa::defaulttype::Vec3dTypes;
#endif
#ifndef SOFA_DOUBLE
using sofa::defaulttype::Vec3fTypes;
#endif

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_FORCEFIELD_TRIANGULARBIQUADRATICSPRINGSFORCEFIELD_CPP)

#ifndef SOFA_FLOAT
extern template class SOFA_DEFORMABLE_API TriangularBiquadraticSpringsForceField<Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_DEFORMABLE_API TriangularBiquadraticSpringsForceField<Vec3fTypes>;
#endif

#endif // defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_FORCEFIELD_TRIANGULARBIQUADRATICSPRINGSFORCEFIELD_CPP)

} //namespace forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FORCEFIELD_TRIANGULARBIQUADRATICSPRINGSFORCEFIELD_H
