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
#ifndef SOFA_COMPONENT_COLLISION_TRIANGLELOCALMINDISTANCEFILTER_H
#define SOFA_COMPONENT_COLLISION_TRIANGLELOCALMINDISTANCEFILTER_H

#include <sofa/component/collision/LocalMinDistanceFilter.h>
#include <sofa/component/collision/LineLocalMinDistanceFilter.h>
#include <sofa/component/collision/PointLocalMinDistanceFilter.h>
#include <sofa/component/collision/TriangleModel.h>
#include <sofa/component/topology/TopologyData.h>

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/helper/AdvancedTimer.h>

namespace sofa
{

namespace component
{

namespace collision
{
using sofa::helper::AdvancedTimer;

/**
 * @brief LocalMinDistance cone information class for a Triangle collision primitive.
 */
class TriangleInfo : public InfoFilter //< Triangle >
{
    typedef sofa::core::topology::BaseMeshTopology::Triangle Triangle;
public:
    /**
     * @brief Default constructor.
     */
    TriangleInfo(LocalMinDistanceFilter *lmdFilters)
        : InfoFilter(lmdFilters)
    {

    }

    /**
     * @brief Empty constructor. Required by TriangleData<>.
     */
    TriangleInfo()
        : InfoFilter(NULL)
    {

    }

    /**
     * @brief Default destructor.
     */
    virtual ~TriangleInfo() {};

    /**
     * @brief Returns the validity of a detected contact according to this TriangleInfo.
     */
    virtual bool validate(const unsigned int /*p*/, const defaulttype::Vector3 & /*PQ*/);
    /**
     * @brief Output stream.
     */
    inline friend std::ostream& operator<< ( std::ostream& os, const TriangleInfo& /*ti*/ )
    {
        return os;
    }

    /**
     * @brief Input stream.
     */
    inline friend std::istream& operator>> ( std::istream& in, TriangleInfo& /*ti*/ )
    {
        return in;
    }

    /**
     * @brief Computes the region of interest cone of the Triangle primitive.
     */
    //virtual void buildFilter(const Triangle & /*t*/);
    virtual void buildFilter(unsigned int /*t*/);

protected:


    Vector3 m_normal; ///< Stored normal of the triangle.
};



/**
 * @brief
 */
class SOFA_MESH_COLLISION_API TriangleLocalMinDistanceFilter : public LocalMinDistanceFilter
{
public:
    SOFA_CLASS(TriangleLocalMinDistanceFilter, LocalMinDistanceFilter);

protected:
    TriangleLocalMinDistanceFilter();
    virtual ~TriangleLocalMinDistanceFilter();

public:

    /**
     * @brief Scene graph initialization method.
     */
    void init();

    /**
     * @brief Handle topological changes.
     */
    void handleTopologyChange();

    /**
     * @name These methods check the validity of a found intersection.
     * According to the local configuration around the found intersected primitive, we build a "Region Of Interest" geometric cone.
     * Pertinent intersections have to belong to this cone, others are not taking into account anymore.
     * If the filtration cone is unbuilt or invalid, these methods launch the build or update.
     */
    //@{

    /**
     * @brief Point Collision Primitive validation method.
     */
    bool validPoint(const int pointIndex, const defaulttype::Vector3 &PQ);

    /**
     * @brief Line Collision Primitive validation method.
     */
    bool validLine(const int lineIndex, const defaulttype::Vector3 &PQ);

    /**
     * @brief Triangle Collision Primitive validation method.
     */
    bool validTriangle(const int triangleIndex, const defaulttype::Vector3 &PQ);

    //@}

    /**
     * @brief New Points creations callback.
     */
    class PointInfoHandler : public topology::TopologyDataHandler<topology::Point, helper::vector<PointInfo> >
    {
    public:
        PointInfoHandler(TriangleLocalMinDistanceFilter* _f, topology::PointData<helper::vector<PointInfo> >* _data) : topology::TopologyDataHandler<topology::Point, helper::vector<PointInfo> >(_data), f(_f) {}

        void applyCreateFunction(unsigned int pointIndex, PointInfo& m, const sofa::helper::vector< unsigned int > &,
                const sofa::helper::vector< double > &);
    protected:
        TriangleLocalMinDistanceFilter* f;
    };

    /**
     * @brief New Edges creations callback.
     */
    class LineInfoHandler : public topology::TopologyDataHandler<topology::Edge, helper::vector<LineInfo> >
    {
    public:
        LineInfoHandler(TriangleLocalMinDistanceFilter* _f, topology::EdgeData<helper::vector<LineInfo> >* _data) : topology::TopologyDataHandler<topology::Edge, helper::vector<LineInfo> >(_data), f(_f) {}

        void applyCreateFunction(unsigned int edgeIndex, LineInfo& m, const topology::Edge&, const sofa::helper::vector< unsigned int > &,
                const sofa::helper::vector< double > &);
    protected:
        TriangleLocalMinDistanceFilter* f;
    };

    /**
     * @brief New Triangles creations callback.
     */
    class TriangleInfoHandler : public topology::TopologyDataHandler<topology::Triangle, helper::vector<TriangleInfo> >
    {
    public:
        TriangleInfoHandler(TriangleLocalMinDistanceFilter* _f, topology::TriangleData<helper::vector<TriangleInfo> >* _data) : topology::TopologyDataHandler<topology::Triangle, helper::vector<TriangleInfo> >(_data), f(_f) {}

        void applyCreateFunction(unsigned int triangleIndex, TriangleInfo& m, const topology::Triangle&, const sofa::helper::vector< unsigned int > &,
                const sofa::helper::vector< double > &);
    protected:
        TriangleLocalMinDistanceFilter* f;
    };

private:
    topology::PointData< sofa::helper::vector<PointInfo> > m_pointInfo;
    topology::EdgeData< sofa::helper::vector<LineInfo> > m_lineInfo;
    topology::TriangleData< sofa::helper::vector<TriangleInfo> > m_triangleInfo;

    PointInfoHandler* pointInfoHandler;
    LineInfoHandler* lineInfoHandler;
    TriangleInfoHandler* triangleInfoHandler;

    core::topology::BaseMeshTopology *bmt;
};


} // namespace collision

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_COLLISION_TRIANGLELOCALMINDISTANCEFILTER_H
