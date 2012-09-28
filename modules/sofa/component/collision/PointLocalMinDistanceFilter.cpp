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

#include <sofa/component/collision/PointLocalMinDistanceFilter.h>
#include <sofa/core/visual/VisualParams.h>

#include <sofa/component/collision/LineModel.h>
#include <sofa/component/topology/TopologyData.inl>

#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/topology/Topology.h>

#include <sofa/simulation/common/Node.h>

#include <sofa/core/ObjectFactory.h>


namespace sofa
{

namespace component
{

namespace collision
{


void PointInfo::buildFilter(unsigned int p_index)
{
    using sofa::simulation::Node;
    using sofa::helper::vector;
    using sofa::core::topology::BaseMeshTopology;


    bool debug=false;
    if((int)p_index==-1)
        debug=true;
    //std::cout<<"buildFilter for point"<<p_index;
    m_noLineModel = false;



    // get the positions:
    const sofa::helper::vector<sofa::defaulttype::Vector3>& x = *this->position_filtering;
    const Vector3 &pt = x[p_index];

    //std::cout<<"  pt"<<pt<<std::endl;

    // get the topology
    BaseMeshTopology* bmt = this->base_mesh_topology;
    const vector< unsigned int >& edgesAroundVertex = bmt->getEdgesAroundVertex(p_index);
    const vector< unsigned int >& trianglesAroundVertex = bmt->getTrianglesAroundVertex(p_index);

    if(edgesAroundVertex.size() ==0)
    {
        std::cerr<<"WARNING no topology defined: no filtering"<<std::endl;
        std::cout<<"Mesh Topology found :"<<bmt->getName()<<std::endl;
        m_noLineModel = true;
        setValid();
        return;
    }


    // compute the normal (nMean) of the point : IS IT STORED ANYWHERE ELSE ?
    // 1. using triangle around the point
    vector< unsigned int >::const_iterator triIt = trianglesAroundVertex.begin();
    vector< unsigned int >::const_iterator triItEnd = trianglesAroundVertex.end();
    Vector3 nMean;
    while (triIt != triItEnd)
    {
        const BaseMeshTopology::Triangle& triangle = bmt->getTriangle(*triIt);
        Vector3 nCur = (x[triangle[1]] - x[triangle[0]]).cross(x[triangle[2]] - x[triangle[0]]);
        nCur.normalize();
        nMean += nCur;
        ++triIt;
    }

    // 2. if no triangle around the point: compute an other normal using edges
    if (trianglesAroundVertex.empty())
    {
        if(debug)
            std::cout<<" trianglesAroundVertex.empty !"<<std::endl;
        vector< unsigned int >::const_iterator edgeIt = edgesAroundVertex.begin();
        vector< unsigned int >::const_iterator edgeItEnd = edgesAroundVertex.end();

        while (edgeIt != edgeItEnd)
        {
            const BaseMeshTopology::Edge& edge = bmt->getEdge(*edgeIt);

            Vector3 l = (pt - x[edge[0]]) + (pt - x[edge[1]]);
            l.normalize();
            nMean += l;
            ++edgeIt;
        }
    }

    // 3. test to verify the normal value and normalize it
    if (nMean.norm() > 1e-20)
        nMean.normalize();
    else
    {
        std::cerr << "WARNING PointInfo m_nMean is null" << std::endl;
    }

    if (debug)
        std::cout<<"  nMean ="<<nMean<<std::endl;


    // Build the set of unit vector that are normal to the planes that defines the cone
    // for each plane, we can "extend" the cone: allow for a larger cone
    vector< unsigned int >::const_iterator edgeIt = edgesAroundVertex.begin();
    vector< unsigned int >::const_iterator edgeItEnd = edgesAroundVertex.end();

    m_computedData.clear();
    while (edgeIt != edgeItEnd)
    {
        const BaseMeshTopology::Edge& edge = bmt->getEdge(*edgeIt);

        Vector3 l = (pt - x[edge[0]]) + (pt - x[edge[1]]);
        l.normalize();



        double computedAngleCone = dot(nMean , l) * m_lmdFilters->getConeExtension();

        if (computedAngleCone < 0)
            computedAngleCone = 0.0;

        computedAngleCone += m_lmdFilters->getConeMinAngle();
        //std::cout<<"  add filtration with l="<<l<<"    and angle="<<computedAngleCone<<std::endl;
        m_computedData.push_back(std::make_pair(l, computedAngleCone));
        ++edgeIt;

        if (debug)
            std::cout<<"  l ="<<l<<"computedAngleCone ="<< computedAngleCone<<"  for edge ["<<edge[0]<<"-"<<edge[1]<<"]"<<std::endl;

    }


    setValid();
}



bool PointInfo::validate(const unsigned int p, const defaulttype::Vector3 &PQ)
{

    bool debug=false;
    if ((int)p==-1)
        debug=true;

    if (isValid())
    {
        if(debug)
            std::cout<<"Point "<<p<<" is valid"<<std::endl;

        if (m_noLineModel)
        {
            std::cout<<"Warning : No Line Model"<<std::endl;
            return true;
        }

        TDataContainer::const_iterator it = m_computedData.begin();
        TDataContainer::const_iterator itEnd = m_computedData.end();

        while (it != itEnd)
        {
            if(debug)
                std::cout<<" test avec direction : "<<it->first <<"   dot(it->first , PQ)="<<dot(it->first , PQ)<<"    (-it->second * PQ.norm()) ="<<(-it->second * PQ.norm())<<std::endl;
            if (dot(it->first , PQ) <= (-it->second * PQ.norm()))
                return false;

            ++it;
        }

        return true;
    }
    else
    {
        if(debug)
            std::cout<<"Point "<<p<<" is not valid ------------ build"<<std::endl;
        buildFilter(p);
        return validate(p, PQ);
    }
}


PointLocalMinDistanceFilter::PointLocalMinDistanceFilter()
    : m_pointInfo(initData(&m_pointInfo, "pointInfo", "point filter data"))
    , pointInfoHandler(NULL)
    , bmt(NULL)
{
}

PointLocalMinDistanceFilter::~PointLocalMinDistanceFilter()
{
    if (pointInfoHandler) delete pointInfoHandler;
}

void PointLocalMinDistanceFilter::init()
{
    bmt = getContext()->getMeshTopology();

    if (bmt != 0)
    {
        helper::vector< PointInfo >& pInfo = *(m_pointInfo.beginEdit());
        pInfo.resize(bmt->getNbPoints());
        m_pointInfo.endEdit();

        pointInfoHandler = new PointInfoHandler(this,&m_pointInfo);
        m_pointInfo.createTopologicalEngine(bmt, pointInfoHandler);
        m_pointInfo.registerTopologicalData();
    }
    if(this->isRigid())
    {
        // Precomputation of the filters in the rigid case
        //points:
        helper::vector< PointInfo >& pInfo = *(m_pointInfo.beginEdit());
        for(unsigned int p=0; p<pInfo.size(); p++)
        {
            pInfo[p].buildFilter(p);

        }
        m_pointInfo.endEdit();

    }

}



void PointLocalMinDistanceFilter::handleTopologyChange()
{
    if(this->isRigid())
    {
        serr<<"WARNING: filters optimization needed for topological change on rigid collision model"<<sendl;
        this->invalidate(); // all the filters will be recomputed, not only those involved in the topological change
    }

    //core::topology::BaseMeshTopology *bmt = getContext()->getMeshTopology();

    //assert(bmt != 0);

    //std::list< const core::topology::TopologyChange * >::const_iterator itBegin = bmt->beginChange();
    //std::list< const core::topology::TopologyChange * >::const_iterator itEnd = bmt->endChange();

    //m_pointInfo.handleTopologyEvents(itBegin, itEnd);
}



void PointLocalMinDistanceFilter::PointInfoHandler::applyCreateFunction(unsigned int /*pointIndex*/, PointInfo &pInfo, const sofa::helper::vector< unsigned int > &, const sofa::helper::vector< double >&)
{

    std::cout<<" LMDFilterPointCreationFunction is called"<<std::endl;
    const PointLocalMinDistanceFilter *pLMDFilter = this->f;
    pInfo.setLMDFilters(pLMDFilter);

    sofa::core::topology::BaseMeshTopology * bmt = pLMDFilter->bmt; //getContext()->getMeshTopology();
    pInfo.setBaseMeshTopology(bmt);
    /////// TODO : template de la classe
    component::container::MechanicalObject<Vec3Types>*  mstateVec3d= dynamic_cast<component::container::MechanicalObject<Vec3Types>*>(pLMDFilter->getContext()->getMechanicalState());
    if(pLMDFilter->isRigid())
    {
        /////// TODO : template de la classe
        if(mstateVec3d != NULL)
        {
            pInfo.setPositionFiltering(mstateVec3d->getX0());
        }

    }
    else
    {
        /////// TODO : template de la classe
        if(mstateVec3d != NULL)
        {
            pInfo.setPositionFiltering(mstateVec3d->getX());
        }
    }

}



SOFA_DECL_CLASS(PointLocalMinDistanceFilter)

int PointLocalMinDistanceFilterClass = core::RegisterObject("This class manages Point collision models cones filters computations and updates.")
        .add< PointLocalMinDistanceFilter >()
        ;

} // namespace collision

} // namespace component

} // namespace sofa
