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
#ifndef SOFA_COMPONENT_TOPOLOGY_MESH2POINTTOPOLOGICALMAPPING_H
#define SOFA_COMPONENT_TOPOLOGY_MESH2POINTTOPOLOGICALMAPPING_H

#include <sofa/core/topology/TopologicalMapping.h>
#include <sofa/component/topology/PointSetTopologyModifier.h>

#include <sofa/defaulttype/Vec.h>
#include <map>
#include <set>

#include <sofa/core/BaseMapping.h>
#include <sofa/component/topology/TopologyData.h>


namespace sofa
{
namespace component
{
namespace topology
{
using namespace sofa::defaulttype;
using namespace sofa::component::topology;
using namespace sofa::core::topology;
using namespace sofa::core;

/**
 * This class, called Mesh2PointTopologicalMapping, is a specific implementation of the interface TopologicalMapping where :
 *
 * INPUT TOPOLOGY = any MeshTopology
 * OUTPUT TOPOLOGY = A PointSetTopologie, as the boundary of the INPUT TOPOLOGY
 *
 * Each primitive in the input Topology will be mapped to a point in the output topology computed from a parameter vector (pointBaryCoords, edgeBaryCoords, triangleBaryCoords, quadBaryCoords, tetraBaryCoords, hexaBaryCoords)
 *
 * Mesh2PointTopologicalMapping class is templated by the pair (INPUT TOPOLOGY, OUTPUT TOPOLOGY)
 *
*/

class SOFA_TOPOLOGY_MAPPING_API Mesh2PointTopologicalMapping : public TopologicalMapping
{
public:
    SOFA_CLASS(Mesh2PointTopologicalMapping,TopologicalMapping);
protected:
    /** \brief Constructor.
     *
     */
    Mesh2PointTopologicalMapping ();

    /** \brief Destructor.
     *
         * Does nothing.
         */
    virtual ~Mesh2PointTopologicalMapping() {};
public:
    /** \brief Initializes the target BaseTopology from the source BaseTopology.
     */
    virtual void init();

    /// Method called at each topological changes propagation which comes from the INPUT topology to adapt the OUTPUT topology :
    virtual void updateTopologicalMappingTopDown();

    virtual unsigned int getGlobIndex(unsigned int ind)
    {
        if(ind<pointSource.size())
        {
            return pointSource[ind].second;
        }
        else
        {
            return 0;
        }
    }

    virtual unsigned int getFromIndex(unsigned int ind)
    {
        return ind;
    }

    enum Element
    {
        POINT = 0,
        EDGE,
        TRIANGLE,
        QUAD,
        TETRA,
        HEXA,
        NB_ELEMENTS
    };

    const vector< vector<int> >& getPointsMappedFromPoint() const { return pointsMappedFrom[POINT]; }
    const vector< vector<int> >& getPointsMappedFromEdge() const { return pointsMappedFrom[EDGE]; }
    const vector< vector<int> >& getPointsMappedFromTriangle() const { return pointsMappedFrom[TRIANGLE]; }
    const vector< vector<int> >& getPointsMappedFromQuad() const { return pointsMappedFrom[QUAD]; }
    const vector< vector<int> >& getPointsMappedFromTetra() const { return pointsMappedFrom[TETRA]; }
    const vector< vector<int> >& getPointsMappedFromHexa() const { return pointsMappedFrom[HEXA]; }

    const vector< Vec3d >& getPointBaryCoords() const { return pointBaryCoords.getValue(); }
    const vector< Vec3d >& getEdgeBaryCoords() const { return edgeBaryCoords.getValue(); }
    const vector< Vec3d >& getTriangleBaryCoords() const { return triangleBaryCoords.getValue(); }
    const vector< Vec3d >& getQuadBaryCoords() const { return quadBaryCoords.getValue(); }
    const vector< Vec3d >& getTetraBaryCoords() const { return tetraBaryCoords.getValue(); }
    const vector< Vec3d >& getHexaBaryCoords() const { return hexaBaryCoords.getValue(); }

    const vector< std::pair<Element,int> >& getPointSource() const { return pointSource;}

protected:

    Data< vector< Vec3d > > pointBaryCoords; ///< Coordinates for the points of the output topology created from the points of the input topology
    Data< vector< Vec3d > > edgeBaryCoords; ///< Coordinates for the points of the output topology created from the edges of the input topology
    Data< vector< Vec3d > > triangleBaryCoords; ///< Coordinates for the points of the output topology created from the triangles of the input topology
    Data< vector< Vec3d > > quadBaryCoords; ///< Coordinates for the points of the output topology created from the quads of the input topology
    Data< vector< Vec3d > > tetraBaryCoords; ///< Coordinates for the points of the output topology created from the tetra of the input topology
    Data< vector< Vec3d > > hexaBaryCoords; ///< Coordinates for the points of the output topology created from the hexa of the input topology

    Data< bool > copyEdges; ///< Activate mapping of input edges into the output topology (requires at least one item in pointBaryCoords)
    Data< bool > copyTriangles; ///< Activate mapping of input triangles into the output topology (requires at least one item in pointBaryCoords)

    helper::fixed_array< vector< vector<int> >, NB_ELEMENTS > pointsMappedFrom; ///< Points mapped from the differents elements (see the enum Element declared before)

    vector< std::pair<Element,int> > pointSource; ///< Correspondance between the points mapped and the elements from which are mapped

    std::set<unsigned int> pointsToRemove;

    size_t addInputPoint(unsigned int i, PointSetTopologyModifier* toPointMod=NULL); ///< Returns the number of points added inside the output topology. 
    void addInputEdge(unsigned int i, PointSetTopologyModifier* toPointMod=NULL);
    void addInputTriangle(unsigned int i, PointSetTopologyModifier* toPointMod=NULL);

    void swapInput(Element elem, int i1, int i2);
    void removeInput(Element elem, const sofa::helper::vector<unsigned int>& tab );
    void renumberInput(Element elem, const sofa::helper::vector<unsigned int>& index );

    void swapOutputPoints(int i1, int i2, bool removeLast = false);
    void removeOutputPoints( const sofa::helper::vector<unsigned int>& tab );

protected:
    bool internalCheck(const char* step, const helper::fixed_array <int, NB_ELEMENTS >& nbInputRemoved);
    
    bool internalCheck(const char* step)
    {
        helper::fixed_array <int, NB_ELEMENTS > nbInputRemoved;
        nbInputRemoved.assign(0);
        return internalCheck(step, nbInputRemoved);
    }
    bool initDone;
};

} // namespace topology
} // namespace component
} // namespace sofa

#endif // SOFA_COMPONENT_TOPOLOGY_MESH2POINTTOPOLOGICALMAPPING_H
