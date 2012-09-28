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
#ifndef SOFA_COMPONENT_TOPOLOGY_QUADSETTOPOLOGYCONTAINER_H
#define SOFA_COMPONENT_TOPOLOGY_QUADSETTOPOLOGYCONTAINER_H

#include <sofa/component/topology/EdgeSetTopologyContainer.h>

namespace sofa
{
namespace component
{
namespace topology
{
class QuadSetTopologyModifier;

using core::topology::BaseMeshTopology;

typedef BaseMeshTopology::PointID			PointID;
typedef BaseMeshTopology::EdgeID			EdgeID;
typedef BaseMeshTopology::QuadID			QuadID;
typedef BaseMeshTopology::Edge				Edge;
typedef BaseMeshTopology::Quad				Quad;
typedef BaseMeshTopology::SeqQuads			SeqQuads;
typedef BaseMeshTopology::EdgesInQuad			EdgesInQuad;
typedef BaseMeshTopology::QuadsAroundVertex		QuadsAroundVertex;
typedef BaseMeshTopology::QuadsAroundEdge		QuadsAroundEdge;
typedef sofa::helper::vector<QuadID>                  VecQuadID;

/** Object that stores a set of quads and provides access
to each quad and its edges and vertices */
class SOFA_BASE_TOPOLOGY_API QuadSetTopologyContainer : public EdgeSetTopologyContainer
{
    friend class QuadSetTopologyModifier;

public:
    SOFA_CLASS(QuadSetTopologyContainer,EdgeSetTopologyContainer);
protected:
    QuadSetTopologyContainer();

    virtual ~QuadSetTopologyContainer() {}
public:
    virtual void init();


    /// Procedural creation methods
    /// @{
    virtual void clear();
    virtual void addQuad( int a, int b, int c, int d );
    /// @}



    /// BaseMeshTopology API
    /// @{

    /** \brief Returns the quad array.
     *
     */
    virtual const SeqQuads& getQuads()
    {
        return getQuadArray();
    }

    /** \brief Returns a reference to the Data of quads array container. */
    Data< sofa::helper::vector<Quad> >& getQuadDataArray() {return d_quad;}

    /** \brief Returns the quad corresponding to the QuadID i.
     *
     * @param ID of a Quad.
     * @return The corresponding Quad.
     */
    virtual const Quad getQuad(QuadID i);


    /** Returns the indices of a quad given four vertex indices.
     *
     * @param the four vertex indices.
     * @return the ID of the corresponding quad.
     * @return -1 if none
     */
    virtual int getQuadIndex(PointID v1, PointID v2, PointID v3, PointID v4);


    /** \brief Returns the set of edges adjacent to a given quad.
     *
     * @param ID of a quad.
     * @return EdgesInQuad list composing the input quad.
    */
    virtual const EdgesInQuad& getEdgesInQuad(QuadID i) ;


    /** \brief Returns the set of quads adjacent to a given vertex.
     *
     * @param ID of a vertex.
     * @return QuadsAroundVertex list around the input vertex.
     */
    virtual const QuadsAroundVertex& getQuadsAroundVertex(PointID i);


    /** \brief Returns the set of quads adjacent to a given edge.
     *
     * @param ID of an edge.
     * @return QuadsAroundEdge list around the input edge.
     */
    virtual const QuadsAroundEdge& getQuadsAroundEdge(EdgeID i);


    /** \brief Returns the index (either 0, 1, 2, 3) of the vertex whose global index is vertexIndex.
     *
     * @param Ref to a quad.
     * @param Id of a vertex.
     * @return the position of this vertex in the quad (i.e. either 0, 1, 2, 3).
     * @return -1 if none.
     */
    virtual int getVertexIndexInQuad(const Quad &t, PointID vertexIndex) const;


    /** \brief Returns the index (either 0, 1, 2, 3) of the edge whose global index is edgeIndex.
     *
     * @param Ref to an EdgesInQuad.
     * @param Id of an edge .
     * @return the position of this edge in the quad (i.e. either 0, 1, 2, 3).
     * @return -1 if none.
     */
    virtual int getEdgeIndexInQuad(const EdgesInQuad &t, EdgeID edheIndex) const;

    /// @}



    /// Dynamic Topology API
    /// @{


    /** \brief Checks if the topology is coherent
    *
    * Check if the shell arrays are coherent
    * @see m_quad
    * @see m_edgesInQuad
    * @see m_quadsAroundVertex
    * @see m_quadsAroundEdge
    */
    virtual bool checkTopology() const;


    /// Get information about connexity of the mesh
    /// @{
    /** \brief Checks if the topology has only one connected component
      *
      * @return true if only one connected component
      */
    virtual bool checkConnexity();

    /// Returns the number of connected component.
    virtual unsigned int getNumberOfConnectedComponent();

    /// Returns the set of element indices connected to an input one (i.e. which can be reached by topological links)
    virtual const VecQuadID getConnectedElement(QuadID elem);

    /// Returns the set of element indices adjacent to a given element (i.e. sharing a link)
    virtual const VecQuadID getElementAroundElement(QuadID elem);
    /// Returns the set of element indices adjacent to a given list of elements (i.e. sharing a link)
    virtual const VecQuadID getElementAroundElements(VecQuadID elems);
    /// @}


    /** \brief Returns the number of quads in this topology.
     * The difference to getNbQuads() is that this method does not generate the quad array if it does not exist.
     */
    unsigned int getNumberOfQuads() const;


    /** \brief Returns the number of topological element of the current topology.
     * This function avoids to know which topological container is in used.
     */
    virtual unsigned int getNumberOfElements() const;


    /** \brief Returns the Quad array. */
    const sofa::helper::vector<Quad> &getQuadArray();


    /** \brief Returns the EdgesInQuadArray array (i.e. provide the 4 edge indices for each quad) */
    const sofa::helper::vector< EdgesInQuad > &getEdgesInQuadArray() ;


    /** \brief Returns the QuadsAroundVertex array (i.e. provide the quad indices adjacent to each vertex). */
    const sofa::helper::vector< QuadsAroundVertex > &getQuadsAroundVertexArray();


    /** \brief Returns the QuadsAroundEdge array (i.e. provide the quad indices adjacent to each edge). */
    const sofa::helper::vector< QuadsAroundEdge > &getQuadsAroundEdgeArray() ;


    bool hasQuads() const;

    bool hasEdgesInQuad() const;

    bool hasQuadsAroundVertex() const;

    bool hasQuadsAroundEdge() const;

    /// @}


protected:

    /** \brief Creates the QuadSet array.
     *
     * This function must be implemented by derived classes to create a list of quads from a set of hexahedra for instance.
     */
    virtual void createQuadSetArray();


    /** \brief Creates the EdgeSet array.
     *
     * Create the set of edges when needed.
     */
    virtual void createEdgeSetArray();


    /** \brief Creates the array of edge indices for each quad.
     *
     * This function is only called if the EdgesInQuad array is required.
     * m_edgesInQuad[i] contains the 4 indices of the 4 edges composing the ith quad.
     */
    virtual void createEdgesInQuadArray();


    /** \brief Creates the QuadsAroundVertex Array.
     *
     * This function is only called if the QuadsAroundVertex array is required.
     * m_quadsAroundVertex[i] contains the indices of all quads adjacent to the ith vertex.
     */
    virtual void createQuadsAroundVertexArray();


    /** \brief Creates the quadsAroundEdge Array.
     *
     * This function is only called if the QuadsAroundVertex array is required.
     * m_quadsAroundEdge[i] contains the indices of all quads adjacent to the ith edge
     */
    virtual void createQuadsAroundEdgeArray();


    void clearQuads();

    void clearEdgesInQuad();

    void clearQuadsAroundVertex();

    void clearQuadsAroundEdge();


protected:

    /** \brief Returns a non-const list of quad indices around a given DOF for subsequent modification.
     *
     * @return QuadsAroundVertex lists in non-const.
     * @see getQuadsAroundVertex()
     */
    virtual QuadsAroundVertex& getQuadsAroundVertexForModification(const PointID vertexIndex);


    /** \brief Returns a non-const list of quad indices around a given edge for subsequent modification.
     *
     * @return QuadsAroundEdge lists in non-const.
     * @see getQuadsAroundEdge()
     */
    virtual QuadsAroundEdge& getQuadsAroundEdgeForModification(const EdgeID edgeIndex);



    /// \brief Function creating the data graph linked to d_quad
    virtual void updateTopologyEngineGraph();


    /// Use a specific boolean @see m_quadTopologyDirty in order to know if topology Data is dirty or not.
    /// Set/Get function access to this boolean
    void setQuadTopologyToDirty() {m_quadTopologyDirty = true;}
    void cleanQuadTopologyFromDirty() {m_quadTopologyDirty = false;}
    const bool& isQuadTopologyDirty() {return m_quadTopologyDirty;}

protected:

    /// provides the set of quads.
    Data< sofa::helper::vector<Quad> > d_quad;

    /// provides the 4 edges in each quad.
    sofa::helper::vector<EdgesInQuad> m_edgesInQuad;

    /// for each vertex provides the set of quads adjacent to that vertex.
    sofa::helper::vector< QuadsAroundVertex > m_quadsAroundVertex;

    /// for each edge provides the set of quads adjacent to that edge.
    sofa::helper::vector< QuadsAroundEdge > m_quadsAroundEdge;


    /// Boolean used to know if the topology Data of this container is dirty
    bool m_quadTopologyDirty;

    /// List of engines related to this specific container
    sofa::helper::list <sofa::core::topology::TopologyEngine *> m_enginesList;

    /// \brief variables used to display the graph of Data/DataEngines linked to this Data array.
    sofa::helper::vector < sofa::helper::vector <std::string> > m_dataGraph;
    sofa::helper::vector < sofa::helper::vector <std::string> > m_enginesGraph;
};

} // namespace topology

} // namespace component

} // namespace sofa

#endif
