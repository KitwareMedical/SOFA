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
#ifndef SOFA_COMPONENT_INTERACTIONFORCEFIELD_DISTANCEGRIDFORCEFIELD_H
#define SOFA_COMPONENT_INTERACTIONFORCEFIELD_DISTANCEGRIDFORCEFIELD_H

#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/component/container/DistanceGrid.h>
#include <sofa/component/component.h>

namespace sofa
{

namespace component
{

namespace forcefield
{

/// This class can be overridden if needed for additionnal storage within template specializations.
template<class DataTypes>
class DistanceGridForceFieldInternalData
{
public:
};

template<class DataTypes>
class DistanceGridForceField : public core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(DistanceGridForceField, DataTypes), SOFA_TEMPLATE(core::behavior::ForceField, DataTypes));

    typedef core::behavior::ForceField<DataTypes> Inherit;
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef Data<VecCoord> DataVecCoord;
    typedef Data<VecDeriv> DataVecDeriv;
    typedef container::DistanceGrid DistanceGrid;

protected:
    DistanceGrid* grid;

    class Contact
    {
    public:
        int index;
        Coord normal;
        Real fact;
        Contact( int index=0, Coord normal=Coord(),Real fact=Real(0))
            : index(index),normal(normal),fact(fact)
        {
        }

        inline friend std::istream& operator >> ( std::istream& in, Contact& c )
        {
            in>>c.index>>c.normal>>c.fact;
            return in;
        }

        inline friend std::ostream& operator << ( std::ostream& out, const Contact& c )
        {
            out << c.index << " " << c.normal << " " << c.fact ;
            return out;
        }

    };

    sofa::helper::vector<bool> pOnBorder;
    Data<sofa::helper::vector<Contact> > contacts;


    class TContact
    {
    public:
        helper::fixed_array<unsigned int,3> index;
        Coord normal,B,C;
        Real fact;

        inline friend std::istream& operator >> ( std::istream& in, TContact& c )
        {
            in>>c.index>>c.normal>>c.B>>c.C>>c.fact;
            return in;
        }

        inline friend std::ostream& operator << ( std::ostream& out, const TContact& c )
        {
            out << c.index << " " << c.normal << " " << c.B << " " << c.C << " " << c.fact ;
            return out;
        }

    };

    Data<sofa::helper::vector<TContact> > tcontacts;

    class VContact
    {
    public:
        helper::fixed_array<unsigned int,4> index;
        Coord A,B,C;
        Real fact;

        inline friend std::istream& operator >> ( std::istream& in, VContact& c )
        {
            in>>c.index>>c.A>>c.B>>c.C>>c.fact;
            return in;
        }

        inline friend std::ostream& operator << ( std::ostream& out, const VContact& c )
        {
            out << c.index << " " << c.A << " " << c.B << " " << c.C << " " << c.fact ;
            return out;
        }

    };

    Data<sofa::helper::vector<VContact> > vcontacts;

    DistanceGridForceFieldInternalData<DataTypes> data;

public:

    // Input data parameters
    sofa::core::objectmodel::DataFileName fileDistanceGrid;
    Data< double > scale;
    Data< helper::fixed_array<DistanceGrid::Coord,2> > box;
    Data< int > nx;
    Data< int > ny;
    Data< int > nz;

    Data<Real> stiffnessIn, stiffnessOut;
    Data<Real> damping;
    Data<Real> maxDist;
    Data<Real> minArea;
    Data<Real> stiffnessArea;
    Data<Real> minVolume;
    Data<Real> stiffnessVolume;
    bool flipNormals;

    Data<defaulttype::Vec3f> color;
    Data<bool> bDraw;
    Data<Real> drawSize;

    /// optional range of local DOF indices. Any computation involving only indices outside of this range are discarded (useful for parallelization using mesh partitionning)
    Data< defaulttype::Vec<2,int> > localRange;
protected:
    DistanceGridForceField()
        : grid(NULL)
        , fileDistanceGrid( initData( &fileDistanceGrid, "fileDistanceGrid", "load distance grid from specified file"))
        , scale( initData( &scale, 1.0, "scale", "scaling factor for input file"))
        , box( initData( &box, "box", "Field bounding box defined by xmin,ymin,zmin, xmax,ymax,zmax") )
        , nx( initData( &nx, 64, "nx", "number of values on X axis") )
        , ny( initData( &ny, 64, "ny", "number of values on Y axis") )
        , nz( initData( &nz, 64, "nz", "number of values on Z axis") )
        , stiffnessIn(initData(&stiffnessIn, (Real)500, "stiffnessIn", "force stiffness when inside of the object"))
        , stiffnessOut(initData(&stiffnessOut, (Real)0, "stiffnessOut", "force stiffness when outside of the object"))
        , damping(initData(&damping, (Real)0.01, "damping", "force damping coefficient"))
        , maxDist(initData(&maxDist, (Real)1.0, "maxdist", "max distance of the surface after which no more force is applied"))
        , minArea(initData(&minArea, (Real)0, "minArea", "minimal area for each triangle, as seen from the direction of the local surface (i.e. a flipped triangle will have a negative area)"))
        , stiffnessArea(initData(&stiffnessArea, (Real)100, "stiffnessArea", "force stiffness if a triangle have an area less than minArea"))
        , minVolume(initData(&minVolume, (Real)0, "minVolume", "minimal volume for each tetrahedron (a flipped triangle will have a negative volume)"))
        , stiffnessVolume(initData(&stiffnessVolume, (Real)0, "stiffnessVolume", "force stiffness if a tetrahedron have an volume less than minVolume"))
        , color(initData(&color, defaulttype::Vec3f(0.0f,.5f,.2f), "color", "display color"))
        , bDraw(initData(&bDraw, false, "draw", "enable/disable drawing of distancegrid"))
        , drawSize(initData(&drawSize, (Real)10.0f, "drawSize", "display size if draw is enabled"))
        , localRange( initData(&localRange, defaulttype::Vec<2,int>(-1,-1), "localRange", "optional range of local DOF indices. Any computation involving only indices outside of this range are discarded (useful for parallelization using mesh partitionning)" ) )
    {
        this->addAlias(&stiffnessIn,"stiffness");
        this->addAlias(&stiffnessOut,"stiffness");
        this->addAlias(&fileDistanceGrid,"filename");
    }
public:
    virtual void init();

    void setMState(  core::behavior::MechanicalState<DataTypes>* mstate ) { this->mstate = mstate; }

    void setStiffness(Real stiffIn, Real stiffOut)
    {
        stiffnessIn.setValue( stiffIn );
        stiffnessOut.setValue( stiffOut );
    }

    void setDamping(Real damp)
    {
        damping.setValue( damp );
    }

    virtual void addForce(const sofa::core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecDeriv &  dataF, const DataVecCoord &  dataX , const DataVecDeriv & dataV ) ;
    virtual void addDForce(const sofa::core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv&   datadF , const DataVecDeriv&   datadX ) ;
    virtual void addKToMatrix(const sofa::core::MechanicalParams* mparams /* PARAMS FIRST */, const sofa::core::behavior::MultiMatrixAccessor* matrix) ;

    void draw(const core::visual::VisualParams* vparams);
    void drawDistanceGrid(const core::visual::VisualParams*, float size=0.0f);


};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_INTERACTIONFORCEFIELD_DISTANCEGRIDFORCEFIELD_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_VOLUMETRIC_DATA_API DistanceGridForceField<defaulttype::Vec3dTypes>;
//extern template class SOFA_VOLUMETRIC_DATA_API DistanceGridForceField<defaulttype::Vec2dTypes>;
//extern template class SOFA_VOLUMETRIC_DATA_API DistanceGridForceField<defaulttype::Vec1dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_VOLUMETRIC_DATA_API DistanceGridForceField<defaulttype::Vec3fTypes>;
//extern template class SOFA_VOLUMETRIC_DATA_API DistanceGridForceField<defaulttype::Vec2fTypes>;
//extern template class SOFA_VOLUMETRIC_DATA_API DistanceGridForceField<defaulttype::Vec1fTypes>;
#endif
#endif

} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_INTERACTIONFORCEFIELD_DISTANCEGRIDFORCEFIELD_H
