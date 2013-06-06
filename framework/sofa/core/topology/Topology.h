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
#ifndef SOFA_CORE_TOPOLOGY_TOPOLOGY_H
#define SOFA_CORE_TOPOLOGY_TOPOLOGY_H

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/helper/list.h>
#include <sofa/core/DataEngine.h>

#include <sofa/helper/fixed_array.h>
#include <iostream>
#include <sofa/helper/vector.h>
#include <stdlib.h>
#include <string>


namespace sofa
{

namespace core
{

namespace topology
{

/// The enumeration used to give unique identifiers to Topological objects.
enum TopologyObjectType
{
    POINT,
    EDGE,
    TRIANGLE,
    QUAD,
    TETRAHEDRON,
    HEXAHEDRON
};



class SOFA_CORE_API Topology : public virtual core::objectmodel::BaseObject
{
public:
    /// Topology global typedefs
    //typedef int index_type;
    typedef unsigned int index_type;
    enum { InvalidID = (unsigned)-1 };
    typedef index_type	        	    PointID;
    typedef index_type          		    EdgeID;
    typedef index_type                          TriangleID;
    typedef index_type                 	    QuadID;
    typedef index_type	                    TetraID;
    typedef index_type	                    TetrahedronID;
    typedef index_type	                    HexaID;
    typedef index_type	                    HexahedronID;


    typedef sofa::helper::vector<index_type>                  SetIndex;
    typedef sofa::helper::vector<index_type>                  SetIndices;

    typedef PointID                             Point;
    // in the following types, we use wrapper classes to have different types for each element, otherwise Quad and Tetrahedron would be the same
    class Edge : public sofa::helper::fixed_array<PointID,2>
    {
    public:
        Edge() {}
        Edge(PointID a, PointID b) : sofa::helper::fixed_array<PointID,2>(a,b) {}
    };
    class Triangle : public sofa::helper::fixed_array<PointID,3>
    {
    public:
        Triangle() {}
        Triangle(PointID a, PointID b, PointID c) : sofa::helper::fixed_array<PointID,3>(a,b,c) {}
    };
    class Quad : public sofa::helper::fixed_array<PointID,4>
    {
    public:
        Quad() {}
        Quad(PointID a, PointID b, PointID c, PointID d) : sofa::helper::fixed_array<PointID,4>(a,b,c,d) {}
    };
    class Tetrahedron : public sofa::helper::fixed_array<PointID,4>
    {
    public:
        Tetrahedron() {}
        Tetrahedron(PointID a, PointID b, PointID c, PointID d) : sofa::helper::fixed_array<PointID,4>(a,b,c,d) {}
    };
    typedef Tetrahedron                         Tetra;
    class Hexahedron : public sofa::helper::fixed_array<PointID,8>
    {
    public:
        Hexahedron() {}
        Hexahedron(PointID a, PointID b, PointID c, PointID d, PointID e, PointID f, PointID g, PointID h) : sofa::helper::fixed_array<PointID,8>(a,b,c,d,e,f,g,h) {}
    };
    typedef Hexahedron                          Hexa;

    SOFA_CLASS(Topology, core::objectmodel::BaseObject);
protected:
    Topology():BaseObject() {}
    virtual ~Topology()
    {}
public:
    // Access to embedded position information (in case the topology is a regular grid for instance)
    // This is not very clean and is quit slow but it should only be used during initialization

    virtual bool hasPos() const { return false; }
    virtual int getNbPoints() const { return 0; }
    virtual void setNbPoints(int /*n*/) {}
    virtual double getPX(int /*i*/) const { return 0.0; }
    virtual double getPY(int /*i*/) const { return 0.0; }
    virtual double getPZ(int /*i*/) const { return 0.0; }
};

template<class TopologyElement>
struct TopologyElementInfo;

template<>
struct TopologyElementInfo<Topology::Point>
{
    static TopologyObjectType type() { return POINT; }
    static const char* name() { return "Point"; }
};

template<>
struct TopologyElementInfo<Topology::Edge>
{
    static TopologyObjectType type() { return EDGE; }
    static const char* name() { return "Edge"; }
};

template<>
struct TopologyElementInfo<Topology::Triangle>
{
    static TopologyObjectType type() { return TRIANGLE; }
    static const char* name() { return "Triangle"; }
};

template<>
struct TopologyElementInfo<Topology::Quad>
{
    static TopologyObjectType type() { return QUAD; }
    static const char* name() { return "Quad"; }
};

template<>
struct TopologyElementInfo<Topology::Tetrahedron>
{
    static TopologyObjectType type() { return TETRAHEDRON; }
    static const char* name() { return "Tetrahedron"; }
};

template<>
struct TopologyElementInfo<Topology::Hexahedron>
{
    static TopologyObjectType type() { return HEXAHEDRON; }
    static const char* name() { return "Hexahedron"; }
};

} // namespace topology

} // namespace core

} // namespace sofa

// Specialization of the defaulttype::DataTypeInfo type traits template

namespace sofa
{

namespace defaulttype
{

template<>
struct DataTypeInfo< sofa::core::topology::Topology::Edge > : public FixedArrayTypeInfo<sofa::helper::fixed_array<unsigned int,2> >
{
    static std::string name() { return "Edge"; }
};

template<>
struct DataTypeInfo< sofa::core::topology::Topology::Triangle > : public FixedArrayTypeInfo<sofa::helper::fixed_array<unsigned int,3> >
{
    static std::string name() { return "Triangle"; }
};

template<>
struct DataTypeInfo< sofa::core::topology::Topology::Quad > : public FixedArrayTypeInfo<sofa::helper::fixed_array<unsigned int,4> >
{
    static std::string name() { return "Quad"; }
};

template<>
struct DataTypeInfo< sofa::core::topology::Topology::Tetrahedron > : public FixedArrayTypeInfo<sofa::helper::fixed_array<unsigned int,4> >
{
    static std::string name() { return "Tetrahedron"; }
};

template<>
struct DataTypeInfo< sofa::core::topology::Topology::Hexahedron > : public FixedArrayTypeInfo<sofa::helper::fixed_array<unsigned int,8> >
{
    static std::string name() { return "Hexahedron"; }
};

} // namespace defaulttype

} // namespace sofa

#endif
