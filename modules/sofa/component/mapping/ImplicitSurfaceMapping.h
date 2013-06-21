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
#ifndef SOFA_COMPONENT_MAPPING_IMPLICITSURFACEMAPPING_H
#define SOFA_COMPONENT_MAPPING_IMPLICITSURFACEMAPPING_H

#include <sofa/core/Mapping.h>
#include <sofa/component/topology/MeshTopology.h>
#include <sofa/helper/MarchingCubeUtility.h>
#include <sofa/defaulttype/VecTypes.h>

namespace sofa
{

namespace component
{

namespace mapping
{

template <class In, class Out>
class ImplicitSurfaceMapping : public core::Mapping<In, Out>, public topology::MeshTopology
{
public:
    SOFA_CLASS2(SOFA_TEMPLATE2(ImplicitSurfaceMapping, In, Out), SOFA_TEMPLATE2(core::Mapping, In, Out), topology::MeshTopology);

    typedef core::Mapping<In, Out> Inherit;
    typedef typename Out::VecCoord OutVecCoord;
    typedef typename Out::VecDeriv OutVecDeriv;
    typedef typename Out::Coord OutCoord;
    typedef typename Out::Deriv OutDeriv;
    typedef typename OutCoord::value_type OutReal;
    typedef typename In::VecCoord InVecCoord;
    typedef typename In::VecDeriv InVecDeriv;
    typedef typename In::Coord InCoord;
    typedef typename In::Deriv InDeriv;
    typedef typename InCoord::value_type InReal;
protected:
    ImplicitSurfaceMapping()
        : Inherit(),
          mStep(initData(&mStep,0.5,"step","Step")),
          mRadius(initData(&mRadius,2.0,"radius","Radius")),
          mIsoValue(initData(&mIsoValue,0.5,"isoValue","Iso Value")),
          mGridMin(initData(&mGridMin,InCoord(-100,-100,-100),"min","Grid Min")),
          mGridMax(initData(&mGridMax,InCoord(100,100,100),"max","Grid Max"))
    {
    }

    virtual ~ImplicitSurfaceMapping()
    {
    }
public:
    virtual void init();

    virtual void parse(core::objectmodel::BaseObjectDescription* arg);

    double getStep() const { return mStep.getValue(); }
    void setStep(double val) { mStep.setValue(val); }

    double getRadius() const { return mRadius.getValue(); }
    void setRadius(double val) { mRadius.setValue(val); }

    double getIsoValue() const { return mIsoValue.getValue(); }
    void setIsoValue(double val) { mIsoValue.setValue(val); }

    const InCoord& getGridMin() const { return mGridMin.getValue(); }
    void setGridMin(const InCoord& val) { mGridMin.setValue(val); }
    void setGridMin(double x, double y, double z) { mGridMin.setValue( InCoord((InReal)x,(InReal)y,(InReal)z)); }

    const InCoord& getGridMax() const { return mGridMax.getValue(); }
    void setGridMax(const InCoord& val) { mGridMax.setValue(val); }
    void setGridMax(double x, double y, double z) { mGridMax.setValue( InCoord((InReal)x,(InReal)y,(InReal)z)); }

    void apply(const core::MechanicalParams *mparams /* PARAMS FIRST */, Data<OutVecCoord>& out, const Data<InVecCoord>& in);
    void applyJ(const core::MechanicalParams *mparams /* PARAMS FIRST */, Data<OutVecDeriv>& out, const Data<InVecDeriv>& in);

    //void applyJT( InVecDeriv& out, const OutVecDeriv& in );
protected:
    Data <double > mStep;
    Data <double > mRadius;
    Data <double > mIsoValue;

    Data< InCoord > mGridMin;
    Data< InCoord > mGridMax;

    // Marching cube data

    /// For each cube, store the vertex indices on each 3 first edges, and the data value
    struct CubeData
    {
        int p[3];
        OutReal data;
        inline friend std::istream& operator >> ( std::istream& in, CubeData& c)
        {
            in >> c.p[0] >> c.p[1] >> c.p[2] >> c.data;

            return in;
        }

        inline friend std::ostream& operator << ( std::ostream& out, const CubeData& c)
        {
            out << c.p[0] << " " << c.p[1] << " " << c.p[2] << " " << c.data ;
            return out;
        }
    };

    Data < sofa::helper::vector<CubeData> > planes;
    typename sofa::helper::vector<CubeData>::iterator P0; /// Pointer to first plane
    typename sofa::helper::vector<CubeData>::iterator P1; /// Pointer to second plane

    void newPlane();

    template<int C>
    int addPoint(OutVecCoord& out, int x,int y,int z, OutReal v0, OutReal v1, OutReal iso)
    {
        int p = out.size();
        OutCoord pos = OutCoord((OutReal)x,(OutReal)y,(OutReal)z);
        pos[C] -= (iso-v0)/(v1-v0);
        out.resize(p+1);
        out[p] = pos * mStep.getValue();
        return p;
    }

    int addFace(int p1, int p2, int p3, int nbp)
    {
        if ((unsigned)p1<(unsigned)nbp &&
            (unsigned)p2<(unsigned)nbp &&
            (unsigned)p3<(unsigned)nbp)
        {
            SeqTriangles& triangles = *seqTriangles.beginEdit();
            int f = triangles.size();
            triangles.push_back(Triangle(p1, p3, p2));
            seqTriangles.endEdit();
            return f;
        }
        else
        {
            serr << "ERROR: Invalid face "<<p1<<" "<<p2<<" "<<p3<<sendl;
            return -1;
        }
    }

};

// MARCHING CUBE TABLES
// ( table copied from an article of Paul Bourke )
// based on code by Cory Gene Bloyd

/* Convention:

       Z
       ^
       |
       4----4----5
      /|        /|
     7 |       5 |
    /  8      /  9
   7---+6----6   |
   |   |     |   |
   |   0----0+---1--> X
   11  /     10  /
   | 3       | 1
   |/        |/
   3----2----2
   /
  /
|_
Y

*/

#ifndef SOFA_FLOAT
using sofa::defaulttype::Vec3dTypes;
#endif

#ifndef SOFA_DOUBLE
using sofa::defaulttype::Vec3fTypes;
using sofa::defaulttype::ExtVec3fTypes;
#endif

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_MAPPING_IMPLICITSURFACEMAPPING_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_VOLUMETRIC_DATA_API ImplicitSurfaceMapping< Vec3dTypes, Vec3dTypes >;
extern template class SOFA_VOLUMETRIC_DATA_API ImplicitSurfaceMapping< Vec3dTypes, ExtVec3fTypes >;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_VOLUMETRIC_DATA_API ImplicitSurfaceMapping< Vec3fTypes, Vec3fTypes >;
extern template class SOFA_VOLUMETRIC_DATA_API ImplicitSurfaceMapping< Vec3fTypes, ExtVec3fTypes >;
#endif

#ifndef SOFA_FLOAT
#ifndef SOFA_DOUBLE
extern template class SOFA_VOLUMETRIC_DATA_API ImplicitSurfaceMapping< Vec3dTypes, Vec3fTypes >;
extern template class SOFA_VOLUMETRIC_DATA_API ImplicitSurfaceMapping< Vec3fTypes, Vec3dTypes >;
#endif
#endif
#endif


} // namespace mapping

} // namespace component

} // namespace sofa

#endif
