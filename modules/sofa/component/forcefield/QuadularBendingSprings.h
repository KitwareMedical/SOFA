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
#ifndef SOFA_COMPONENT_FORCEFIELD_QUADULARBENDINGSPRINGS_H
#define SOFA_COMPONENT_FORCEFIELD_QUADULARBENDINGSPRINGS_H

#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/Mat.h>
#include <sofa/component/topology/TopologyData.h>

#include <sofa/defaulttype/Mat.h>
#include <sofa/helper/fixed_array.h>

#include <map>
#include <set>

#include <sofa/component/component.h>

namespace sofa
{

namespace component
{

namespace forcefield
{

using namespace sofa::helper;
using namespace sofa::defaulttype;
using namespace sofa::component::topology;

/**
Bending springs added between vertices of quads sharing a common edge.
The springs connect the vertices not belonging to the common edge. It compresses when the surface bends along the common edge.
*/
template<class DataTypes>
class QuadularBendingSprings : public core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(QuadularBendingSprings,DataTypes), SOFA_TEMPLATE(core::behavior::ForceField,DataTypes));

    typedef core::behavior::ForceField<DataTypes> Inherited;
    //typedef typename DataTypes::Real Real;
    typedef typename DataTypes::VecCoord VecCoord;

    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename Coord::value_type Real;

    typedef core::objectmodel::Data<VecDeriv>    DataVecDeriv;
    typedef core::objectmodel::Data<VecCoord>    DataVecCoord;

    enum { N=DataTypes::spatial_dimensions };
    typedef defaulttype::Mat<N,N,Real> Mat;

    QuadularBendingSprings();

    ~QuadularBendingSprings();

protected:

    class EdgeInformation
    {
    public:
        Mat DfDx; /// the edge stiffness matrix

        int     m1, m2;  /// the two extremities of the first spring: masses m1 and m2
        int     m3, m4;  /// the two extremities of the second spring: masses m3 and m4

        double  ks;      /// spring stiffness (initialized to the default value)
        double  kd;      /// damping factor (initialized to the default value)

        double  restlength1; /// rest length of the first spring
        double  restlength2; /// rest length of the second spring

        bool is_activated;

        bool is_initialized;

        EdgeInformation(int m1=0, int m2=0, int m3=0, int m4=0, double restlength1=0.0, double restlength2=0.0, bool is_activated=false, bool is_initialized=false)
            : m1(m1), m2(m2), m3(m3), m4(m4), restlength1(restlength1), restlength2(restlength2), is_activated(is_activated), is_initialized(is_initialized)
        {
        }

        /// Output stream
        inline friend std::ostream& operator<< ( std::ostream& os, const EdgeInformation& /*ei*/ )
        {
            return os;
        }

        /// Input stream
        inline friend std::istream& operator>> ( std::istream& in, EdgeInformation& /*ei*/ )
        {
            return in;
        }
    };



public:
    /// Searches quad topology and creates the bending springs
    virtual void init();

    virtual void addForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& d_f, const DataVecCoord& d_x, const DataVecDeriv& d_v);
    virtual void addDForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& d_df, const DataVecDeriv& d_dx);

    virtual double getKs() const { return f_ks.getValue();}
    virtual double getKd() const { return f_kd.getValue();}

    void setKs(const double ks)
    {
        f_ks.setValue((double)ks);
    }
    void setKd(const double kd)
    {
        f_kd.setValue((double)kd);
    }

    // -- VisualModel interface
    void draw(const core::visual::VisualParams* vparams);
    void initTextures() { }
    void update() { }

    EdgeData<sofa::helper::vector<EdgeInformation> > &getEdgeInfo() {return edgeInfo;}


    class EdgeBSHandler : public TopologyDataHandler<Edge,sofa::helper::vector<EdgeInformation> >
    {
    public:
        typedef typename QuadularBendingSprings<DataTypes>::EdgeInformation EdgeInformation;

        EdgeBSHandler(QuadularBendingSprings<DataTypes>* ff, EdgeData<sofa::helper::vector<EdgeInformation> >* data )
            :TopologyDataHandler<Edge,sofa::helper::vector<EdgeInformation> >(data)
            ,ff(ff)
        {
        }

        void applyCreateFunction(unsigned int edgeIndex, EdgeInformation& ei,
                const Edge &,
                const sofa::helper::vector< unsigned int > &,
                const sofa::helper::vector< double > &);

        void applyQuadCreation(const sofa::helper::vector<unsigned int> & quadAdded,
                const sofa::helper::vector<Edge> &,
                const sofa::helper::vector<sofa::helper::vector<unsigned int> > &,
                const sofa::helper::vector<sofa::helper::vector<double> > &);

        void applyQuadDestruction(const sofa::helper::vector<unsigned int> & quadRemoved);

        void applyPointDestruction(const sofa::helper::vector<unsigned int> &pointIndices);

        void applyPointRenumbering(const sofa::helper::vector<unsigned int> &pointToRenumber);

    protected:
        QuadularBendingSprings<DataTypes>* ff;
    };



protected:
    bool updateMatrix;

    Data<double> f_ks;
    Data<double> f_kd;

    double m_potentialEnergy;

    EdgeData<sofa::helper::vector<EdgeInformation> > edgeInfo;

    /// Pointer to the current topology
    sofa::core::topology::BaseMeshTopology* _topology;

    /// Handler for subset Data
    EdgeBSHandler* edgeHandler;
};

using sofa::defaulttype::Vec3dTypes;
using sofa::defaulttype::Vec3fTypes;

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_FORCEFIELD_QUADULARBENDINGSPRINGS_CPP)

#ifndef SOFA_FLOAT
extern template class SOFA_DEFORMABLE_API QuadularBendingSprings<Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_DEFORMABLE_API QuadularBendingSprings<Vec3fTypes>;
#endif


#endif // defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_FORCEFIELD_QUADULARBENDINGSPRINGS_CPP)

} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FORCEFIELD_QUADULARBENDINGSPRINGS_H
