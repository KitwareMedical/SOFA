/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *
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
*                              SOFA :: Framework                              *
*                                                                             *
* Authors: The SOFA Team (see Authors.txt)                                    *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_CORE_TOPOLOGY_BASETOPOLOGY_H
#define SOFA_CORE_TOPOLOGY_BASETOPOLOGY_H

#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/topology/BaseTopologyObject.h>
#include <sofa/core/VecId.h>

#include <sofa/helper/list.h>


namespace sofa
{

namespace core
{

namespace topology
{
using core::topology::BaseMeshTopology;

// forward declarations:

/// Provides high-level topology algorithms (e.g. CutAlongPlane, DecimateTopology, etc).
class TopologyAlgorithms;

/// Provides some geometric functions (e.g. ComputeTriangleNormal, ComputeShell, etc).
class GeometryAlgorithms;

/// Provides low-level topology methods (e.g. AddPoint, RemoveEdge, etc).
class TopologyModifier;

/// Contains the actual topology data and give acces to it.
class TopologyContainer;

/// Translates topology events (TopologyChange objects) from a topology so that they apply on another one.
class TopologicalMapping;

/// Allow topological handle events
class TopologyEngine;


/** A class that contains a set of high-level (user frisendly) methods that perform topological changes */
class SOFA_CORE_API TopologyAlgorithms : public sofa::core::topology::BaseTopologyObject
{
public:
    SOFA_CLASS(TopologyAlgorithms, BaseTopologyObject);

protected:
    /** \brief Constructor.
    *
    */
    TopologyAlgorithms()
    {}


    /// Destructor
    virtual ~TopologyAlgorithms()
    {}
public:
    virtual void init();

protected:
    /** \brief Adds a TopologyChange object to the list of the topology this object describes.
    */
    void addTopologyChange(const TopologyChange *topologyChange);

protected:
    /// Contains the actual topology data and give acces to it (nature of these data heavily depends on the kind of topology).
    TopologyContainer *m_topologyContainer;
};

/** A class that contains a set of methods that describes the geometry of the object */
class SOFA_CORE_API GeometryAlgorithms : public sofa::core::topology::BaseTopologyObject
{
public:
    SOFA_CLASS(GeometryAlgorithms, BaseTopologyObject);

protected:
    /** \brief Constructor.
    *
    */
    GeometryAlgorithms()
    {}


    /// Destructor
    virtual ~GeometryAlgorithms()
    {}
public:
    virtual void init();

    /** \brief Called by the MechanicalObject state change callback to initialize added
    * points according to the topology (topology element & local coordinates)
    *
    * \param ancestorElems are the ancestors topology info used in the points modifications
    */
    virtual void initPointsAdded(const helper::vector< unsigned int > &indices, const helper::vector< AncestorElem > &ancestorElems
        , const helper::vector< core::VecCoordId >& coordVecs, const helper::vector< core::VecDerivId >& derivVecs );
};

/** A class that contains a set of low-level methods that perform topological changes */
class SOFA_CORE_API TopologyModifier : public sofa::core::topology::BaseTopologyObject
{
public:
    SOFA_CLASS(TopologyModifier, BaseTopologyObject);

protected:
    /** \brief Constructor.
    *
    */
    TopologyModifier()
        : m_topologyContainer(NULL)
    { }


    /// Destructor
    virtual ~TopologyModifier()
    { }
public:
    virtual void init();

    /** \brief Called by a topology to warn the Mechanical Object component that points have been added or will be removed.
    *
    * StateChangeList should contain all TopologyChange objects corresponding to vertex changes in this topology
    * that just happened (in the case of creation) or are about to happen (in the case of destruction) since
    * last call to propagateTopologicalChanges.
    *
    * @sa firstChange()
    * @sa lastChange()
    */
    virtual void propagateStateChanges();

    /** \brief Called by a topology to warn specific topologies linked to it that TopologyChange objects happened.
    *
    * ChangeList should contain all TopologyChange objects corresponding to changes in this topology
    * that just happened (in the case of creation) or are about to happen (in the case of destruction) since
    * last call to propagateTopologicalChanges.
    *
    * @sa firstChange()
    * @sa lastChange()
    */
    virtual void propagateTopologicalChanges();

    /** \notify the end for the current sequence of topological change events.
    */
    virtual void notifyEndingEvent();

    /** \brief Generic method to remove a list of items.
    */
    virtual void removeItems(sofa::helper::vector< unsigned int >& /*items*/);

protected:
    /** \brief Adds a TopologyChange object to the list of the topology this object describes.
    */
    void addTopologyChange(const TopologyChange *topologyChange);

    /** \brief Adds a StateChange object to the list of the topology this object describes.
    */
    void addStateChange(const TopologyChange *topologyChange);

protected:
    /// Contains the actual topology data and give acces to it (nature of these data heavily depends on the kind of topology).
    TopologyContainer *m_topologyContainer;
};



/** A class that contains a description of the topology (set of edges, triangles, adjacency information, ...) */
class SOFA_CORE_API TopologyContainer : public sofa::core::topology::BaseTopologyObject,
    public core::topology::BaseMeshTopology
{
public:
    SOFA_CLASS2(TopologyContainer, BaseTopologyObject, BaseMeshTopology);

protected:

    /** \brief Constructor.
    *
    */
    TopologyContainer()
    {}


    /// Destructor
    virtual ~TopologyContainer();
public:
    virtual void init();

    /// BaseMeshTopology API
    /// @{
    virtual const SeqEdges& getEdges()         { static SeqEdges     empty; return empty; }
    virtual const SeqTriangles& getTriangles() { static SeqTriangles empty; return empty; }
    virtual const SeqQuads& getQuads()         { static SeqQuads     empty; return empty; }
    virtual const SeqTetrahedra& getTetrahedra()       { static SeqTetrahedra    empty; return empty; }
    virtual const SeqHexahedra& getHexahedra()         { static SeqHexahedra     empty; return empty; }

    /** \brief Get the current revision of this mesh.
    *
    * This can be used to detect changes, however topological changes event should be used whenever possible.
    */
    virtual int getRevision() const { return m_changeList.getCounter(); }

    /// @}

    /// TopologyChange interactions
    /// @{
    const sofa::helper::list<const TopologyChange *> &getChangeList() const { return m_changeList.getValue(); }

    const sofa::helper::list<const TopologyChange *> &getStateChangeList() const { return m_stateChangeList.getValue(); }

    const Data <sofa::helper::list<const TopologyChange *> > &getDataChangeList() const { return m_changeList; }

    const Data <sofa::helper::list<const TopologyChange *> > &getDataStateChangeList() const { return m_stateChangeList; }

    /** \brief Adds a TopologyChange to the list.
    *
    * Needed by topologies linked to this one to know what happened and what to do to take it into account.
    *
    */
    virtual void addTopologyChange(const TopologyChange *topologyChange);

    /** \brief Adds a StateChange to the list.
    *
    * Needed by topologies linked to this one to know what happened and what to do to take it into account.
    *
    */
    virtual void addStateChange(const TopologyChange *topologyChange);

    /** \brief Provides an iterator on the first element in the list of TopologyChange objects.
     */
    sofa::helper::list<const TopologyChange *>::const_iterator beginChange() const;

    /** \brief Provides an iterator on the last element in the list of TopologyChange objects.
     */
    sofa::helper::list<const TopologyChange *>::const_iterator endChange() const;

    /** \brief Provides an iterator on the first element in the list of StateChange objects.
     */
    sofa::helper::list<const TopologyChange *>::const_iterator beginStateChange() const;

    /** \brief Provides an iterator on the last element in the list of StateChange objects.
     */
    sofa::helper::list<const TopologyChange *>::const_iterator endStateChange() const;


    /** \brief Free each Topology changes in the list and remove them from the list
    *
    */
    virtual void resetTopologyChangeList();

    /** \brief Free each State changes in the list and remove them from the list
    *
    */
    virtual void resetStateChangeList();

    ///@}

    /// TopologyEngine interactions
    ///@{
    const sofa::helper::list<TopologyEngine *> &getTopologyEngineList() const { return m_topologyEngineList; }

    /** \brief Adds a TopologyEngine to the list.
    */
    virtual void addTopologyEngine(TopologyEngine* _topologyEngine);


    /** \brief Provides an iterator on the first element in the list of TopologyEngine objects.
     */
    sofa::helper::list<TopologyEngine *>::const_iterator beginTopologyEngine() const;

    /** \brief Provides an iterator on the last element in the list of TopologyEngine objects.
     */
    sofa::helper::list<TopologyEngine *>::const_iterator endTopologyEngine() const;

    /** \brief Free each Topology changes in the list and remove them from the list
    *
    */
    void resetTopologyEngineList();

    ///@}


protected:

    virtual void updateTopologyEngineGraph() {}

    /// Array of topology modifications that have already occured (addition) or will occur next (deletion).
    Data <sofa::helper::list<const TopologyChange *> >m_changeList;

    /// Array of state modifications that have already occured (addition) or will occur next (deletion).
    Data <sofa::helper::list<const TopologyChange *> >m_stateChangeList;

    /// List of topology engines which will interact on all topological Data.
    sofa::helper::list<TopologyEngine *> m_topologyEngineList;
};




} // namespace topology

} // namespace core

} // namespace sofa

#endif // SOFA_CORE_BASICTOPOLOGY_H
