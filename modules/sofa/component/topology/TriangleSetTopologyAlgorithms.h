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
#ifndef SOFA_COMPONENT_TOPOLOGY_TRIANGLESETTOPOLOGYALGORITHMS_H
#define SOFA_COMPONENT_TOPOLOGY_TRIANGLESETTOPOLOGYALGORITHMS_H

#include <sofa/component/topology/EdgeSetTopologyAlgorithms.h>

namespace sofa
{
namespace component
{
namespace topology
{
class TriangleSetTopologyContainer;

class TriangleSetTopologyModifier;

template < class DataTypes >
class TriangleSetGeometryAlgorithms;

using core::topology::BaseMeshTopology;
typedef BaseMeshTopology::TriangleID TriangleID;
typedef BaseMeshTopology::Triangle Triangle;
typedef BaseMeshTopology::SeqTriangles SeqTriangles;
typedef BaseMeshTopology::TrianglesAroundVertex TrianglesAroundVertex;
typedef BaseMeshTopology::TrianglesAroundEdge TrianglesAroundEdge;
typedef BaseMeshTopology::EdgesInTriangle EdgesInTriangle;

/**
* A class that performs topology algorithms on an TriangleSet.
*/
template < class DataTypes >
class TriangleSetTopologyAlgorithms : public EdgeSetTopologyAlgorithms<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(TriangleSetTopologyAlgorithms,DataTypes), SOFA_TEMPLATE(EdgeSetTopologyAlgorithms,DataTypes));

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecCoord VecDeriv;
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
protected:
    TriangleSetTopologyAlgorithms()
        : EdgeSetTopologyAlgorithms<DataTypes>()
        , m_listTriRemove( initData(&m_listTriRemove,  "RemoveTrianglesByIndex", "Debug : Remove a triangle or a list of triangles by using their indices (only while animate)."))
        , m_listTriAdd( initData(&m_listTriAdd,  "addTrianglesByIndex", "Debug : Add a triangle or a list of triangles by using their indices (only while animate)."))
    {
    }

    virtual ~TriangleSetTopologyAlgorithms() {}
public:
    virtual void init();

    virtual void reinit();

    /** \brief  Moves and fixes the two closest points of two triangles to their median point
     */
    bool Suture2Points(unsigned int ind_ta, unsigned int ind_tb, unsigned int &ind1, unsigned int &ind2);

    /** \brief Removes triangles along the list of points (ind_edge,coord) intersected by the vector from point a to point b and the triangular mesh
     */
    void RemoveAlongTrianglesList(const sofa::defaulttype::Vec<3,double>& a,
            const sofa::defaulttype::Vec<3,double>& b,
            const unsigned int ind_ta, const unsigned int ind_tb);

    /** \brief Incises along the list of points (ind_edge,coord) intersected by the sequence of input segments (list of input points) and the triangular mesh
     */
    void InciseAlongLinesList(const sofa::helper::vector< sofa::defaulttype::Vec<3,double> >& input_points,
            const sofa::helper::vector< unsigned int > &input_triangles);

    /** \brief Duplicates the given edge. Only works if at least one of its points is adjacent to a border.
     * @returns the number of newly created points, or -1 if the incision failed.
     */
    virtual int InciseAlongEdge(unsigned int edge, int* createdPoints = NULL);


    /** \brief Split triangles to create edges along a path given as a the list of existing edges and triangles crossed by it.
     * Each end of the path is given either by an existing point or a point inside the first/last triangle. If the first/last triangle is (TriangleID)-1, it means that to path crosses the boundary of the surface.
     * @returns the indice of the end point, or -1 if the incision failed.
     */
    virtual int SplitAlongPath(unsigned int pa, Coord& a, unsigned int pb, Coord& b,
            sofa::helper::vector< sofa::core::topology::TopologyObjectType>& topoPath_list,
            sofa::helper::vector<unsigned int>& indices_list,
            sofa::helper::vector< sofa::defaulttype::Vec<3, double> >& coords_list,
            sofa::helper::vector<EdgeID>& new_edges, double epsilonSnapPath = 0.0, double epsilonSnapBorder = 0.0);



    /* void SnapAlongPath (sofa::helper::vector<TriangleID>& triangles_list, sofa::helper::vector<EdgeID>& edges_list,
      sofa::helper::vector<double>& coords_list, sofa::helper::vector<double>& points2Snap);*/

    void SnapAlongPath (sofa::helper::vector< sofa::core::topology::TopologyObjectType>& topoPath_list,
            sofa::helper::vector<unsigned int>& indices_list, sofa::helper::vector< sofa::defaulttype::Vec<3, double> >& coords_list,
            sofa::helper::vector< sofa::helper::vector<double> >& points2Snap,
            double epsilonSnapPath);

    void SnapBorderPath (unsigned int pa, Coord& a, unsigned int pb, Coord& b,
            sofa::helper::vector< sofa::core::topology::TopologyObjectType>& topoPath_list,
            sofa::helper::vector<unsigned int>& indices_list,
            sofa::helper::vector< sofa::defaulttype::Vec<3, double> >& coords_list,
            sofa::helper::vector< sofa::helper::vector<double> >& points2Snap,
            double epsilonSnapBorder);



    /** \brief Duplicates the given edges. Only works if at least the first or last point is adjacent to a border.
     * @returns true if the incision succeeded.
     */
    virtual bool InciseAlongEdgeList(const sofa::helper::vector<unsigned int>& edges, sofa::helper::vector<unsigned int>& new_points, sofa::helper::vector<unsigned int>& end_points, bool& reachBorder);

    unsigned int getOtherPointInTriangle(const Triangle& t, unsigned int p1, unsigned int p2) const
    {
        if (t[0] != p1 && t[0] != p2) return t[0];
        else if (t[1] != p1 && t[1] != p2) return t[1];
        else return t[2];
    }


protected:
    Data< sofa::helper::vector< unsigned int> > m_listTriRemove;
    Data< sofa::helper::vector< Triangle> > m_listTriAdd;

private:
    TriangleSetTopologyContainer*				m_container;
    TriangleSetTopologyModifier*				m_modifier;
    TriangleSetGeometryAlgorithms< DataTypes >*		m_geometryAlgorithms;


};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_TOPOLOGY_TRIANGLESETTOPOLOGYALGORITHMS_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_BASE_TOPOLOGY_API TriangleSetTopologyAlgorithms<defaulttype::Vec3dTypes>;
extern template class SOFA_BASE_TOPOLOGY_API TriangleSetTopologyAlgorithms<defaulttype::Vec2dTypes>;
extern template class SOFA_BASE_TOPOLOGY_API TriangleSetTopologyAlgorithms<defaulttype::Vec1dTypes>;
//extern template class SOFA_BASE_TOPOLOGY_API TriangleSetTopologyAlgorithms<defaulttype::Rigid3dTypes>;
//extern template class SOFA_BASE_TOPOLOGY_API TriangleSetTopologyAlgorithms<defaulttype::Rigid2dTypes>;
#endif

#ifndef SOFA_DOUBLE
extern template class SOFA_BASE_TOPOLOGY_API TriangleSetTopologyAlgorithms<defaulttype::Vec3fTypes>;
extern template class SOFA_BASE_TOPOLOGY_API TriangleSetTopologyAlgorithms<defaulttype::Vec2fTypes>;
extern template class SOFA_BASE_TOPOLOGY_API TriangleSetTopologyAlgorithms<defaulttype::Vec1fTypes>;
//extern template class SOFA_BASE_TOPOLOGY_API TriangleSetTopologyAlgorithms<defaulttype::Rigid3fTypes>;
//extern template class SOFA_BASE_TOPOLOGY_API TriangleSetTopologyAlgorithms<defaulttype::Rigid2fTypes>;
#endif
#endif

} // namespace topology

} // namespace component

} // namespace sofa

#endif
