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
#ifndef SOFA_COMPONENT_TOPOLOGY_BASETOPOLOGYDATA_H
#define SOFA_COMPONENT_TOPOLOGY_BASETOPOLOGYDATA_H

#include <sofa/core/topology/Topology.h>

namespace sofa
{

namespace core
{

namespace topology
{



/** A class that define topological Data general methods

      DEPRECATED

*/
template < class T = void* >
class BaseTopologyData : public sofa::core::objectmodel::Data <T>
{
public:
    //SOFA_CLASS(SOFA_TEMPLATE2(BaseTopologyData,T,VecT), SOFA_TEMPLATE(sofa::core::objectmodel::Data, T));

    class InitData : public sofa::core::objectmodel::BaseData::BaseInitData
    {
    public:
        InitData() : value(T()) {}
        InitData(const T& v) : value(v) {}
        InitData(const sofa::core::objectmodel::BaseData::BaseInitData& i) : sofa::core::objectmodel::BaseData::BaseInitData(i), value(T()) {}

        T value;
    };

    /** Constructor
            this constructor should be used through the initData() methods
         */
    explicit BaseTopologyData(const sofa::core::objectmodel::BaseData::BaseInitData& init)
        : Data<T>(init)
    {
    }

    /** Constructor
            this constructor should be used through the initData() methods
         */
    explicit BaseTopologyData(const InitData& init)
        : Data<T>(init)
    {
    }


    /** Constructor
        \param helpMsg help on the field
         */
    BaseTopologyData( const char* helpMsg=0, bool isDisplayed=true, bool isReadOnly=false)
        : Data<T>(helpMsg, isDisplayed, isReadOnly)
    {

    }

    /** Constructor
        \param value default value
        \param helpMsg help on the field
         */
    BaseTopologyData( const T& value, const char* helpMsg=0, bool isDisplayed=true, bool isReadOnly=false)
        : Data<T>(helpMsg, isDisplayed, isReadOnly)
    {
    }


    // Generic methods to apply changes on the Data
    //{
    /// Apply adding points elements.
    virtual void applyCreatePointFunction(const sofa::helper::vector<unsigned int>& ) {}
    /// Apply removing points elements.
    virtual void applyDestroyPointFunction(const sofa::helper::vector<unsigned int>& ) {}

    /// Apply adding edges elements.
    virtual void applyCreateEdgeFunction(const sofa::helper::vector<unsigned int>& ) {}
    /// Apply removing edges elements.
    virtual void applyDestroyEdgeFunction(const sofa::helper::vector<unsigned int>& ) {}

    /// Apply adding triangles elements.
    virtual void applyCreateTriangleFunction(const sofa::helper::vector<unsigned int>& ) {}
    /// Apply removing triangles elements.
    virtual void applyDestroyTriangleFunction(const sofa::helper::vector<unsigned int>& ) {}

    /// Apply adding quads elements.
    virtual void applyCreateQuadFunction(const sofa::helper::vector<unsigned int>& ) {}
    /// Apply removing quads elements.
    virtual void applyDestroyQuadFunction(const sofa::helper::vector<unsigned int>& ) {}

    /// Apply adding tetrahedra elements.
    virtual void applyCreateTetrahedronFunction(const sofa::helper::vector<unsigned int>& ) {}
    /// Apply removing tetrahedra elements.
    virtual void applyDestroyTetrahedronFunction(const sofa::helper::vector<unsigned int>& ) {}

    /// Apply adding hexahedra elements.
    virtual void applyCreateHexahedronFunction(const sofa::helper::vector<unsigned int>& ) {}
    /// Apply removing hexahedra elements.
    virtual void applyDestroyHexahedronFunction(const sofa::helper::vector<unsigned int>& ) {}
    //}

    /// Add some values. Values are added at the end of the vector.
    virtual void add(unsigned int ,
            const sofa::helper::vector< sofa::helper::vector< unsigned int > >& ,
            const sofa::helper::vector< sofa::helper::vector< double > >& ) {}

    /// Temporary Hack: find a way to have a generic description of topological element:
    /// add Edge
    virtual void add( unsigned int ,
            const sofa::helper::vector< sofa::helper::fixed_array<unsigned int,2> >& ,
            const sofa::helper::vector< sofa::helper::vector< unsigned int > > &,
            const sofa::helper::vector< sofa::helper::vector< double > >& ) {}

    /// add Triangle
    virtual void add( unsigned int ,
            const sofa::helper::vector< sofa::helper::fixed_array<unsigned int,3> >& ,
            const sofa::helper::vector< sofa::helper::vector< unsigned int > > &,
            const sofa::helper::vector< sofa::helper::vector< double > >& ) {}

    /// add Quad & Tetrahedron
    virtual void add( unsigned int ,
            const sofa::helper::vector< sofa::helper::fixed_array<unsigned int,4> >& ,
            const sofa::helper::vector< sofa::helper::vector< unsigned int > > &,
            const sofa::helper::vector< sofa::helper::vector< double > >& ) {}

    /// add Hexahedron
    virtual void add( unsigned int ,
            const sofa::helper::vector< sofa::helper::fixed_array<unsigned int,8> >& ,
            const sofa::helper::vector< sofa::helper::vector< unsigned int > > &,
            const sofa::helper::vector< sofa::helper::vector< double > >& ) {}

    /// Remove the values corresponding to the points removed.
    virtual void remove( const sofa::helper::vector<unsigned int>& ) {}

    /// Swaps values at indices i1 and i2.
    virtual void swap( unsigned int , unsigned int ) {}

    /// Reorder the values.
    virtual void renumber( const sofa::helper::vector<unsigned int>& ) {}

    /// Move a list of points
    virtual void move( const sofa::helper::vector<unsigned int>& ,
            const sofa::helper::vector< sofa::helper::vector< unsigned int > >& ,
            const sofa::helper::vector< sofa::helper::vector< double > >& ) {}



};


} // namespace topology

} // namespace component

} // namespace sofa

#endif //SOFA_COMPONENT_TOPOLOGY_BASETOPOLOGYDATA_H
