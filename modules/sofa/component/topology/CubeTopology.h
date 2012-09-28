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
#ifndef SOFA_COMPONENT_TOPOLOGY_CUBETOPOLOGY_H
#define SOFA_COMPONENT_TOPOLOGY_CUBETOPOLOGY_H

#include <sofa/component/topology/MeshTopology.h>
#include <sofa/defaulttype/Vec.h>

namespace sofa
{

namespace component
{

namespace topology
{

using namespace sofa::defaulttype;

class SOFA_BASE_TOPOLOGY_API CubeTopology : public MeshTopology
{
public:
    SOFA_CLASS(CubeTopology,MeshTopology);
protected:
    CubeTopology(int nx, int ny, int nz);
    CubeTopology();
public:
    void setSize(int nx, int ny, int nz);

    void parse(core::objectmodel::BaseObjectDescription* arg);

    int getNx() const { return nx.getValue(); }
    int getNy() const { return ny.getValue(); }
    int getNz() const { return nz.getValue(); }

    void setNx(int n) { nx.setValue(n); setSize(); }
    void setNy(int n) { ny.setValue(n); setSize(); }
    void setNz(int n) { nz.setValue(n); setSize(); }

    virtual void init();
    virtual void reinit();

    //int getNbQuads();
    //Quad getQuad(int i);
    //Quad getQuad(int x, int y, int z);

    enum Plane { PLANE_UNKNOWN=0,
            PLANE_X0,
            PLANE_X1,
            PLANE_Y0,
            PLANE_Y1,
            PLANE_Z0,
            PLANE_Z1
               };

    int point(int x, int y, int z, Plane p = PLANE_UNKNOWN) const; // { return x+nx.getValue()*(y+ny.getValue()*z); }

    void setP0(const Vector3& val) { p0 = val; }
    void setDx(const Vector3& val) { dx = val; inv_dx2 = 1/(dx*dx); }
    void setDy(const Vector3& val) { dy = val; inv_dy2 = 1/(dy*dy); }
    void setDz(const Vector3& val) { dz = val; inv_dz2 = 1/(dz*dz); }

    void setPos(SReal xmin, SReal xmax, SReal ymin, SReal ymax, SReal zmin, SReal zmax);

    const Vector3& getP0() const { return p0; }
    const Vector3& getDx() const { return dx; }
    const Vector3& getDy() const { return dy; }
    const Vector3& getDz() const { return dz; }

    Vector3   getMin() const { return min.getValue();}
    Vector3   getMax() const { return max.getValue();}

    Vector3 getPoint(int i) const;
    virtual Vector3 getPoint(int x, int y, int z) const;
    bool hasPos()  const { return true; }
    double getPX(int i)  const { return getPoint(i)[0]; }
    double getPY(int i) const { return getPoint(i)[1]; }
    double getPZ(int i) const { return getPoint(i)[2]; }

    void setSplitNormals(bool b) {splitNormals.setValue(b);}

protected:
    Data<int> nx;
    Data<int> ny;
    Data<int> nz;
    Data<bool> internalPoints;
    Data<bool> splitNormals;

    Data< Vector3 > min, max;
    /// Position of point 0
    Vector3 p0;
    /// Distance between points in the grid. Must be perpendicular to each other
    Vector3 dx,dy,dz;
    SReal inv_dx2, inv_dy2, inv_dz2;

    virtual void setSize();
    void updateEdges();
    void updateQuads();
    //void updateHexahedra();
};

} // namespace topology

} // namespace component

} // namespace sofa

#endif
