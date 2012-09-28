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
#ifndef SOFA_COMPONENT_FORCEFIELD_HEXAHEDRALFEMFORCEFIELD_H
#define SOFA_COMPONENT_FORCEFIELD_HEXAHEDRALFEMFORCEFIELD_H

#include <sofa/core/behavior/ForceField.h>

#include <sofa/helper/vector.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/Mat.h>

#include <sofa/component/topology/HexahedronSetTopologyContainer.h>
#include <sofa/component/topology/TopologyData.h>

namespace sofa
{

namespace component
{

namespace forcefield
{

using namespace sofa::defaulttype;
using sofa::helper::vector;

using namespace sofa::component::topology;

/** Compute Finite Element forces based on hexahedral elements.
*
* Corotational hexahedron from
* @Article{NMPCPF05,
*   author       = "Nesme, Matthieu and Marchal, Maud and Promayon, Emmanuel and Chabanas, Matthieu and Payan, Yohan and Faure, Fran\c{c}ois",
*   title        = "Physically Realistic Interactive Simulation for Biological Soft Tissues",
*   journal      = "Recent Research Developments in Biomechanics",
*   volume       = "2",
*   year         = "2005",
*   keywords     = "surgical simulation physical animation truth cube",
*   url          = "http://www-evasion.imag.fr/Publications/2005/NMPCPF05"
* }
*
* indices ordering (same as in HexahedronSetTopology):
*
*     Y  7---------6
*     ^ /         /|
*     |/    Z    / |
*     3----^----2  |
*     |   /     |  |
*     |  4------|--5
*     | /       | /
*     |/        |/
*     0---------1-->X
*/
template<class DataTypes>
class HexahedralFEMForceField : virtual public core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(HexahedralFEMForceField, DataTypes), SOFA_TEMPLATE(core::behavior::ForceField, DataTypes));

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef VecCoord Vector;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename Coord::value_type Real;
    typedef core::objectmodel::Data<VecCoord> DataVecCoord;
    typedef core::objectmodel::Data<VecDeriv> DataVecDeriv;
    typedef helper::ReadAccessor< DataVecCoord > RDataRefVecCoord;
    typedef helper::WriteAccessor< DataVecDeriv > WDataRefVecDeriv;

    typedef core::topology::BaseMeshTopology::index_type Index;
    typedef core::topology::BaseMeshTopology::Hexa Element;
    typedef core::topology::BaseMeshTopology::SeqHexahedra VecElement;

    typedef Vec<24, Real> Displacement;		///< the displacement vector

    typedef Mat<6, 6, Real> MaterialStiffness;	///< the matrix of material stiffness
    typedef vector<MaterialStiffness> VecMaterialStiffness;  ///< a vector of material stiffness matrices
    typedef Mat<24, 24, Real> ElementMass;

    typedef Mat<24, 24, Real> ElementStiffness;
    typedef vector<ElementStiffness> VecElementStiffness;


    enum
    {
        LARGE = 0,   ///< Symbol of large displacements hexahedron solver
        POLAR = 1,   ///< Symbol of polar displacements hexahedron solver
    };

protected:




    typedef Mat<3, 3, Real> Mat33;
    typedef Mat33 Transformation; ///< matrix for rigid transformations like rotations

    typedef std::pair<int,Real> Col_Value;
    typedef vector< Col_Value > CompressedValue;
    typedef vector< CompressedValue > CompressedMatrix;

    /// the information stored for each hexahedron
    class HexahedronInformation
    {
    public:
        /// material stiffness matrices of each hexahedron
        MaterialStiffness materialMatrix;

        // large displacement method
        helper::fixed_array<Coord,8> rotatedInitialElements;

        Transformation rotation;
        ElementStiffness stiffness;

        HexahedronInformation() {}

        /// Output stream
        inline friend std::ostream& operator<< ( std::ostream& os, const HexahedronInformation& /*hi*/ )
        {
            return os;
        }

        /// Input stream
        inline friend std::istream& operator>> ( std::istream& in, HexahedronInformation& /*hi*/ )
        {
            return in;
        }
    };


    HexahedralFEMForceField();
    virtual ~HexahedralFEMForceField();
public:
    void setPoissonRatio(Real val) { this->f_poissonRatio.setValue(val); }

    void setYoungModulus(Real val) { this->f_youngModulus.setValue(val); }

    void setMethod(int val) { method = val; }

    virtual void init();
    virtual void reinit();

    virtual void addForce (const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v);

    virtual void addDForce (const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& df, const DataVecDeriv& dx);

    void addKToMatrix(const core::MechanicalParams* mparams /* PARAMS FIRST */, const sofa::core::behavior::MultiMatrixAccessor* matrix);

    void draw(const core::visual::VisualParams* vparams);

protected:

    virtual void computeElementStiffness( ElementStiffness &K, const MaterialStiffness &M, const Vec<8,Coord> &nodes);
    Mat33 integrateStiffness( int signx0, int signy0, int signz0, int signx1, int signy1, int signz1, const Real u, const Real v, const Real w, const Mat33& J_1  );

    /// compute the hookean material matrix
    void computeMaterialStiffness(MaterialStiffness &m, double youngModulus, double poissonRatio);

    void computeForce( Displacement &F, const Displacement &Depl, const ElementStiffness &K );


    ////////////// large displacements method
    void initLarge(const int i);
    void computeRotationLarge( Transformation &r, Coord &edgex, Coord &edgey);
    virtual void accumulateForceLarge( WDataRefVecDeriv& f, RDataRefVecCoord& p, const int i);

    ////////////// polar decomposition method
    void initPolar(const int i);
    void computeRotationPolar( Transformation &r, Vec<8,Coord> &nodes);
    virtual void accumulateForcePolar( WDataRefVecDeriv& f, RDataRefVecCoord & p, const int i);

public:
    int method;
    Data<std::string> f_method; ///< the computation method of the displacements
    Data<Real> f_poissonRatio;
    Data<Real> f_youngModulus;
    Data<bool> f_drawing;
    /// container that stotes all requires information for each hexahedron
    HexahedronData<sofa::helper::vector<HexahedronInformation> > hexahedronInfo;

    class HFFHexahedronHandler : public TopologyDataHandler<Hexahedron,sofa::helper::vector<HexahedronInformation> >
    {
    public:
        typedef typename HexahedralFEMForceField<DataTypes>::HexahedronInformation HexahedronInformation;

        HFFHexahedronHandler(HexahedralFEMForceField<DataTypes>* ff, HexahedronData<sofa::helper::vector<HexahedronInformation> >* data )
            :TopologyDataHandler<Hexahedron,sofa::helper::vector<HexahedronInformation> >(data)
            ,ff(ff)
        {
        }

        void applyCreateFunction(unsigned int, HexahedronInformation &t, const Hexahedron &,
                const sofa::helper::vector<unsigned int> &, const sofa::helper::vector<double> &);
    protected:
        HexahedralFEMForceField<DataTypes>* ff;
    };



protected:
    HFFHexahedronHandler* hexahedronHandler;

    HexahedronSetTopologyContainer* _topology;

    Mat<8,3,int> _coef; ///< coef of each vertices to compute the strain stress matrix
};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_FORCEFIELD_HEXAHEDRALFEMFORCEFIELD_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_SIMPLE_FEM_API HexahedralFEMForceField<defaulttype::Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_SIMPLE_FEM_API HexahedralFEMForceField<defaulttype::Vec3fTypes>;
#endif
#endif

} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FORCEFIELD_HEXAHEDRALFEMFORCEFIELD_H
