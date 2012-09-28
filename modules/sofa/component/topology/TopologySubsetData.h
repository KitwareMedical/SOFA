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
#ifndef SOFA_COMPONENT_TOPOLOGY_TOPOLOGYSUBSETDATA_H
#define SOFA_COMPONENT_TOPOLOGY_TOPOLOGYSUBSETDATA_H

#include <sofa/helper/vector.h>
#include <sofa/component/component.h>

#include <sofa/core/topology/BaseTopologyData.h>
#include <sofa/component/topology/TopologyEngine.h>
#include <sofa/component/topology/TopologySubsetDataHandler.h>

namespace sofa
{

namespace component
{

namespace topology
{

////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////   Generic Topology Data Implementation   /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

/** \brief A class for storing Edge related data. Automatically manages topology changes.
*
* This class is a wrapper of class helper::vector that is made to take care transparently of all topology changes that might
* happen (non exhaustive list: Edges added, removed, fused, renumbered).
*/
template< class TopologyElementType, class VecT>
class TopologySubsetDataImpl : public sofa::core::topology::BaseTopologyData<VecT>
{

public:
    typedef VecT container_type;
    typedef typename container_type::value_type value_type;

    /// size_type
    typedef typename container_type::size_type size_type;
    /// reference to a value (read-write)
    typedef typename container_type::reference reference;
    /// const reference to a value (read only)
    typedef typename container_type::const_reference const_reference;
    /// const iterator
    typedef typename container_type::const_iterator const_iterator;


    /// Constructor
    TopologySubsetDataImpl( const typename sofa::core::topology::BaseTopologyData< VecT >::InitData& data)
        : sofa::core::topology::BaseTopologyData< VecT >(data),
          m_topologicalEngine(NULL),
          m_topology(NULL),
          m_topologyHandler(NULL)
    {}

    virtual ~TopologySubsetDataImpl();


    /** Public functions to handle topological engine creation */
    /// To create topological engine link to this Data. Pointer to current topology is needed.
    virtual void createTopologicalEngine(sofa::core::topology::BaseMeshTopology* _topology, sofa::core::topology::TopologyHandler* _topologyHandler);

    /** Public functions to handle topological engine creation */
    /// To create topological engine link to this Data. Pointer to current topology is needed.
    virtual void createTopologicalEngine(sofa::core::topology::BaseMeshTopology* _topology);

    /// Allow to add additionnal dependencies to others Data.
    void addInputData(sofa::core::objectmodel::BaseData* _data);

    /// Function to link the topological Data with the engine and the current topology. And init everything.
    /// This function should be used at the end of the all declaration link to this Data while using it in a component.
    void registerTopologicalData();


    value_type& operator[](int i)
    {
        container_type& data = *(this->beginEdit());
        value_type& result = data[i];
        this->endEdit();
        return result;
    }


    /// Link Data to topology arrays
    void linkToPointDataArray();
    void linkToEdgeDataArray();
    void linkToTriangleDataArray();
    void linkToQuadDataArray();
    void linkToTetrahedronDataArray();
    void linkToHexahedronDataArray();


protected:
    virtual void linkToElementDataArray() {}

    virtual void createTopologyHandler() {}

    typename sofa::component::topology::TopologyEngineImpl<VecT>::SPtr m_topologicalEngine;
    sofa::core::topology::BaseMeshTopology* m_topology;
    sofa::component::topology::TopologySubsetDataHandler<TopologyElementType,VecT>* m_topologyHandler;
};




////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////   Point Topology Data Implementation   /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

template< class VecT >
class PointSubsetData : public TopologySubsetDataImpl<Point, VecT>
{
public:
    typedef typename TopologySubsetDataImpl<Point, VecT>::container_type container_type;
    typedef typename TopologySubsetDataImpl<Point, VecT>::value_type value_type;

    PointSubsetData( const typename sofa::core::topology::BaseTopologyData< VecT >::InitData& data)
        : TopologySubsetDataImpl<Point, VecT>(data)
    {}

protected:
    void linkToElementDataArray() {this->linkToPointDataArray();}
};




////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////   Edge Topology Data Implementation   /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

template< class VecT >
class EdgeSubsetData : public TopologySubsetDataImpl<Edge, VecT>
{
public:
    typedef typename TopologySubsetDataImpl<Edge, VecT>::container_type container_type;
    typedef typename TopologySubsetDataImpl<Edge, VecT>::value_type value_type;

    EdgeSubsetData( const typename sofa::core::topology::BaseTopologyData< VecT >::InitData& data)
        : TopologySubsetDataImpl<Edge, VecT>(data)
    {}

protected:
    void linkToElementDataArray() {this->linkToEdgeDataArray();}

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   Triangle Topology Data Implementation   ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

template< class VecT >
class TriangleSubsetData : public TopologySubsetDataImpl<Triangle, VecT>
{
public:
    typedef typename TopologySubsetDataImpl<Triangle, VecT>::container_type container_type;
    typedef typename TopologySubsetDataImpl<Triangle, VecT>::value_type value_type;

    TriangleSubsetData( const typename sofa::core::topology::BaseTopologyData< VecT >::InitData& data)
        : TopologySubsetDataImpl<Triangle, VecT>(data)
    {}

protected:
    void linkToElementDataArray() {this->linkToTriangleDataArray();}

};



////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////   Quad Topology Data Implementation   /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

template< class VecT >
class QuadSubsetData : public TopologySubsetDataImpl<Quad, VecT>
{
public:
    typedef typename TopologySubsetDataImpl<Quad, VecT>::container_type container_type;
    typedef typename TopologySubsetDataImpl<Quad, VecT>::value_type value_type;

    QuadSubsetData( const typename sofa::core::topology::BaseTopologyData< VecT >::InitData& data)
        : TopologySubsetDataImpl<Quad, VecT>(data)
    {}

protected:
    void linkToElementDataArray() {this->linkToQuadDataArray();}

};



////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////   Tetrahedron Topology Data Implementation   /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

template< class VecT >
class TetrahedronSubsetData : public TopologySubsetDataImpl<Tetrahedron, VecT>
{
public:
    typedef typename TopologySubsetDataImpl<Tetrahedron, VecT>::container_type container_type;
    typedef typename TopologySubsetDataImpl<Tetrahedron, VecT>::value_type value_type;

    TetrahedronSubsetData( const typename sofa::core::topology::BaseTopologyData< VecT >::InitData& data)
        : TopologySubsetDataImpl<Tetrahedron, VecT>(data)
    {}

protected:
    void linkToElementDataArray() {this->linkToTetrahedronDataArray();}

};




////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   Hexahedron Topology Data Implementation   /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

template< class VecT >
class HexahedronSubsetData : public TopologySubsetDataImpl<Hexahedron, VecT>
{
public:
    typedef typename TopologySubsetDataImpl<Hexahedron, VecT>::container_type container_type;
    typedef typename TopologySubsetDataImpl<Hexahedron, VecT>::value_type value_type;

    HexahedronSubsetData( const typename sofa::core::topology::BaseTopologyData< VecT >::InitData& data)
        : TopologySubsetDataImpl<Hexahedron, VecT>(data)
    {}

protected:
    void linkToElementDataArray() {this->linkToHexahedronDataArray();}

};





} // namespace topology

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_TOPOLOGY_TOPOLOGYSUBSETDATA_H
