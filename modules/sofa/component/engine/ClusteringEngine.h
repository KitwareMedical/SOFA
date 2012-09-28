/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 4      *
*                (c) 2006-2009 MGH, INRIA, USTL, UJF, CNRS                    *
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
#ifndef SOFA_COMPONENT_ENGINE_CLUSTERING_H
#define SOFA_COMPONENT_ENGINE_CLUSTERING_H

#if !defined(__GNUC__) || (__GNUC__ > 3 || (_GNUC__ == 3 && __GNUC_MINOR__ > 3))
#pragma once
#endif

#include <sofa/core/DataEngine.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/VecId.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/component/topology/TriangleSetTopologyContainer.h>
#include <sofa/component/topology/TriangleSetGeometryAlgorithms.h>
#include <sofa/component/topology/TriangleSetGeometryAlgorithms.inl>
#include <sofa/defaulttype/Vec.h>
#include <sofa/helper/SVector.h>

#include <sofa/component/component.h>

namespace sofa
{

namespace component
{

namespace engine
{

using namespace core::behavior;
using namespace core::topology;
using namespace core::objectmodel;

/**
 * This class groups points into overlapping clusters according to a user defined number of clusters and radius
 * It takes input positions (and optionally a meshtopology if geodesic distances are prefered)
 * ouput clusters can then be fed to
 *     - shape matching engine
 *     - blendSkinningMapping
 *
 */

template <class DataTypes>
class ClusteringEngine : public core::DataEngine
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(ClusteringEngine,DataTypes),core::DataEngine);
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename Coord::value_type Real;

    typedef BaseMeshTopology::PointID ID;
    typedef helper::vector<ID> VI;
    typedef helper::vector<VI> VVI;

    typedef helper::vector<Real> VD;
    typedef helper::vector<VD> VVD;

    typedef defaulttype::Vec<2,unsigned int> indicesType;

public:

    ClusteringEngine();

    virtual ~ClusteringEngine() {}

    void init();
    void update();

    void draw(const core::visual::VisualParams* vparams);

    Data<bool> useTopo;
    //Data<unsigned int> maxIter;

    Data<Real> radius;
    Data<Real> fixedRadius;
    Data<int> number;
    Data< VecCoord > fixedPosition;  ///< input (non mechanical particle reference position)
    Data< VecCoord > position; ///< input (reference mstate position)
    Data< VVI > cluster;       ///< result

    sofa::core::objectmodel::DataFileName input_filename;
    sofa::core::objectmodel::DataFileName output_filename;

    virtual std::string getTemplateName() const    { return templateName(this);    }
    static std::string templateName(const ClusteringEngine<DataTypes>* = NULL) {   return DataTypes::Name(); }

private:
    MechanicalState<DataTypes>* mstate;
    BaseMeshTopology* topo;

    // recursively add to cluster[i] neighbors of lastNeighbors if dist<radius or if in voronoi(i)+one ring
    void AddNeighborhoodFromNeighborhood(VI& lastN, const unsigned int i, const VI& voronoi);

    // recursively add farthest point from already selected points
    void farthestPointSampling(VI& indices,VI& voronoi);

    // voronoi from a set of points -> returns voronoi and distances
    void Voronoi(const VI& indices , VD& distances, VI& voronoi);

    // dijkstra from a set of points -> returns voronoi and distances = Voronoi function with geodesic distances
    void dijkstra(const VI& indices , VD& distances, VI& voronoi);

    // relax farthest point sampling using LLoyd (k-means) algorithm
    void LLoyd();

    // IO
    bool load();
    std::string loadedFilename;
    bool save();
};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_ENGINE_CLUSTERING_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_ENGINE_API ClusteringEngine<defaulttype::Vec3dTypes>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
extern template class SOFA_ENGINE_API ClusteringEngine<defaulttype::Vec3fTypes>;
#endif //SOFA_DOUBLE
#endif

} // namespace engine

} // namespace component

} // namespace sofa

#endif
