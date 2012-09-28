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
#ifndef SOFA_COMPONENT_COLLISION_TRIANGLEOCTREE_H
#define SOFA_COMPONENT_COLLISION_TRIANGLEOCTREE_H

//#include <sofa/component/collision/TriangleOctreeModel.h>

#include <sofa/core/CollisionModel.h>
#include <sofa/component/container/MechanicalObject.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/defaulttype/Vec3Types.h>
/*THIS STATIC CUBE SIZE MUST BE CHANGE, it represents the size of the occtree cube*/
#define CUBE_SIZE 800
#define bb_max(a,b) (((a)>(b))?(a):(b))
#define bb_max3(a,b,c) bb_max((bb_max(a,b)),c)

#define bb_min(a,b) (((a)<(b))?(a):(b))
#define bb_min3(a,b,c) bb_min(bb_min(a,b),c)

namespace sofa
{

namespace component
{

namespace collision
{

using namespace sofa::defaulttype;

class TriangleOctree;

class SOFA_MESH_COLLISION_API TriangleOctreeRoot
{
public:
    typedef sofa::core::topology::BaseMeshTopology::SeqTriangles SeqTriangles;
    typedef sofa::core::topology::BaseMeshTopology::Triangle Tri;
    typedef sofa::defaulttype::Vec3Types::VecCoord VecCoord;
    typedef sofa::defaulttype::Vec3Types::Coord Coord;
    /// the triangles used as input to construct the octree
    const SeqTriangles* octreeTriangles;
    /// the positions of vertices used as input to construct the octree
    const VecCoord* octreePos;
    /// the first node of the octree
    TriangleOctree* octreeRoot;
    /// the size of the octree cube
    int cubeSize;

    TriangleOctreeRoot();
    ~TriangleOctreeRoot();

    void buildOctree();
    void buildOctree(const sofa::core::topology::BaseMeshTopology::SeqTriangles* triangles, const sofa::defaulttype::Vec3Types::VecCoord* pos)
    {
        this->octreeTriangles = triangles;
        this->octreePos = pos;
        buildOctree();
    }

protected:
    /// used to add a triangle  to the octree
    int fillOctree (int t, int d = 0, Vector3 v = Vector3 (0, 0, 0));
    /// used to compute the Bounding Box for each triangle
    void calcTriangleAABB(int t, double* bb, double& size);
};

class SOFA_MESH_COLLISION_API TriangleOctree
{
public:
    class traceResult
    {
    public:
        traceResult():tid(-1),t(0),u(0),v(0) {}
        int tid;
        double t,u,v;
        bool operator == (const traceResult& r) const { return tid == r.tid && t == r.t && u == r.u && v == r.v; }
        bool operator != (const traceResult& r) const { return tid != r.tid || t != r.t || u != r.u || v != r.v; }
        bool operator < (const traceResult& r) const { return t < r.t; }
        bool operator <= (const traceResult& r) const { return t <= r.t; }
        bool operator > (const traceResult& r) const { return t > r.t; }
        bool operator >= (const traceResult& r) const { return t >= r.t; }
    };

    double x, y, z;
    bool visited;

    double size;
    bool val;
    bool is_leaf;
    bool internal;
    TriangleOctreeRoot *tm;
    vector < int >objects;
    TriangleOctree *childVec[8];

    ~TriangleOctree ();
    /*the default cube has xmin=-CUBE_SIZE xmax=CUBE_SIZE, ymin=-CUBE_SIZE, ymax=CUBE_SIZE, zmin=-CUBE_SIZE,zmax=CUBE_SIZE*/
    TriangleOctree (TriangleOctreeRoot * _tm, double _x = (double)-CUBE_SIZE, double _y = (double)-CUBE_SIZE, double _z =(double) -CUBE_SIZE, double _size = 2 * CUBE_SIZE)
        : x (_x), y (_y), z (_z), size (_size)
        , tm(_tm)
    {
        is_leaf = true;
        internal = false;
        for (int i = 0; i < 8; i++)
            childVec[i] = NULL;
    }

    void draw (const core::visual::VisualParams* vparams);

    /// Find the nearest triangle intersecting the given ray, or -1 of not found
    int trace (Vector3 origin, Vector3 direction, traceResult &result);

    /// Find all triangles intersecting the given ray
    void traceAll (Vector3 origin, Vector3 direction, vector<traceResult>& results);

    /// Find all triangles intersecting the given ray
    void traceAllCandidates(Vector3 origin, Vector3 direction, std::set<int>& results);

    /// Find all triangles intersecting the given ray
    void bboxAllCandidates(Vector3 bbmin, Vector3 bbmax, std::set<int>& results);

    friend class TriangleOctreeRoot;

protected:
    int trace (const Vector3 & origin, const Vector3 & direction,
            double tx0, double ty0, double tz0, double tx1, double ty1,
            double tz1, unsigned int a, unsigned int b,Vector3 &origin1,Vector3 &direction1, traceResult &result);

    template<class Res>
    void traceAllStart (Vector3 origin, Vector3 direction, Res& results);

    template<class Res>
    void traceAll (const Vector3 & origin, const Vector3 & direction,
            double tx0, double ty0, double tz0, double tx1, double ty1,
            double tz1, unsigned int a, unsigned int b,Vector3 &origin1,Vector3 &direction1, Res& results);

    template<class Res>
    void bbAll (const Vector3 & bbmin, const Vector3 & bbmax, Res& results);

    int nearestTriangle (int minIndex, const Vector3 & origin,
            const Vector3 & direction,traceResult &result);

    void allTriangles (const Vector3 & origin,
            const Vector3 & direction, vector<traceResult>& results);

    void allTriangles (const Vector3 & origin,
            const Vector3 & direction, std::set<int>& results);

    void bbAllTriangles (const Vector3 & bbmin,
            const Vector3 & bbmax, std::set<int>& results);

    void insert (double _x, double _y, double _z, double _inc, int t);

};

} // namespace collision

} // namespace component

} // namespace sofa

#endif
