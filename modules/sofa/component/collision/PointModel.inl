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

#ifndef SOFA_COMPONENT_COLLISION_POINTMODEL_INL
#define SOFA_COMPONENT_COLLISION_POINTMODEL_INL

#include <sofa/helper/system/config.h>
#include <sofa/helper/proximity.h>
#include <sofa/defaulttype/Mat.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/core/collision/Intersection.inl>
#include <iostream>
#include <algorithm>




#include <sofa/component/collision/PointModel.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/component/collision/PointLocalMinDistanceFilter.h>
#include <sofa/component/collision/CubeModel.h>
#include <sofa/core/ObjectFactory.h>
#include <vector>
#include <sofa/helper/system/gl.h>
#include <sofa/helper/gl/template.h>
#include <sofa/core/collision/Intersection.inl>

#include <sofa/core/topology/BaseMeshTopology.h>

#include <sofa/simulation/common/Simulation.h>

namespace sofa
{

namespace component
{

namespace collision
{

using namespace sofa::defaulttype;
using namespace sofa::core::collision;
using namespace helper;

template<class DataTypes>
TPointModel<DataTypes>::TPointModel()
    : bothSide(initData(&bothSide, false, "bothSide", "activate collision on both side of the point model (when surface normals are defined on these points)") )
    , mstate(NULL)
    , computeNormals( initData(&computeNormals, false, "computeNormals", "activate computation of normal vectors (required for some collision detection algorithms)") )
    , PointActiverPath(initData(&PointActiverPath,"PointActiverPath", "path of a component PointActiver that activate or deactivate collision point during execution") )
    , m_lmdFilter( NULL )
    , m_displayFreePosition(initData(&m_displayFreePosition, false, "displayFreePosition", "Display Collision Model Points free position(in green)") )
{
    enum_type = POINT_TYPE;
}

template<class DataTypes>
void TPointModel<DataTypes>::resize(int size)
{
    this->core::CollisionModel::resize(size);
}

template<class DataTypes>
void TPointModel<DataTypes>::init()
{
    this->CollisionModel::init();
    mstate = dynamic_cast< core::behavior::MechanicalState<DataTypes>* > (getContext()->getMechanicalState());

    if (mstate==NULL)
    {
        serr<<"ERROR: PointModel requires a Vec3 Mechanical Model" << sendl;
        return;
    }

    simulation::Node* node = dynamic_cast< simulation::Node* >(this->getContext());
    if (node != 0)
    {
        m_lmdFilter = node->getNodeObject< PointLocalMinDistanceFilter >();
    }

    const int npoints = mstate->getX()->size();
    resize(npoints);
    if (computeNormals.getValue()) updateNormals();


    const std::string path = PointActiverPath.getValue();

    if (path.size()==0)
    {

        myActiver = new PointActiver();
        sout<<"path = "<<path<<" no Point Activer found for PointModel "<<this->getName()<<sendl;
    }
    else
    {

        core::objectmodel::BaseObject *activer=NULL;
        this->getContext()->get(activer ,path  );

        if (activer != NULL)
            sout<<" Activer named"<<activer->getName()<<" found"<<sendl;
        else
            serr<<"wrong path for PointActiver"<<sendl;


        myActiver = dynamic_cast<PointActiver *> (activer);



        if (myActiver==NULL)
        {
            myActiver = new PointActiver();


            serr<<"no dynamic cast possible for Point Activer for PointModel "<< this->getName() <<sendl;
        }
        else
        {
            sout<<"PointActiver named"<<activer->getName()<<" found !! for PointModel "<< this->getName() <<sendl;
        }
    }
}

template<class DataTypes>
void TPointModel<DataTypes>::draw(const core::visual::VisualParams* ,int index)
{
#ifndef SOFA_NO_OPENGL
    TPoint<DataTypes> p(this,index);
    if (!p.activated())
        return;
    glBegin(GL_POINTS);
    helper::gl::glVertexT(p.p());
    glEnd();
    if ((unsigned)index < normals.size())
    {
        glBegin(GL_LINES);
        helper::gl::glVertexT(p.p());
        helper::gl::glVertexT(p.p()+normals[index]*0.1f);
        glEnd();
    }
#endif /* SOFA_NO_OPENGL */
}

template<class DataTypes>
void TPointModel<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    if (vparams->displayFlags().getShowCollisionModels())
    {
        if (vparams->displayFlags().getShowWireFrame())
            vparams->drawTool()->setPolygonMode(0,true);

        // Check topological modifications
        const int npoints = mstate->getX()->size();
        if (npoints != size)
        {
            resize(npoints);
        }

        std::vector< Vector3 > pointsP;
        std::vector< Vector3 > pointsL;
        for (int i = 0; i < size; i++)
        {
            TPoint<DataTypes> p(this,i);
            if (p.activated())
            {
                pointsP.push_back(p.p());
                if ((unsigned)i < normals.size())
                {
                    pointsL.push_back(p.p());
                    pointsL.push_back(p.p()+normals[i]*0.1f);
                }
            }
        }

        vparams->drawTool()->drawPoints(pointsP, 3, Vec<4,float>(getColor4f()));
        vparams->drawTool()->drawLines(pointsL, 1, Vec<4,float>(getColor4f()));

        if (m_displayFreePosition.getValue())
        {
            std::vector< Vector3 > pointsPFree;

            for (int i = 0; i < size; i++)
            {
                TPoint<DataTypes> p(this,i);
                if (p.activated())
                {
                    pointsPFree.push_back(p.pFree());
                }
            }

            vparams->drawTool()->drawPoints(pointsPFree, 3, Vec<4,float>(0.0f,1.0f,0.2f,1.0f));
        }

        if (vparams->displayFlags().getShowWireFrame())
            vparams->drawTool()->setPolygonMode(0,false);
    }

    if (getPrevious()!=NULL && vparams->displayFlags().getShowBoundingCollisionModels())
        getPrevious()->draw(vparams);
}

template<class DataTypes>
bool TPointModel<DataTypes>::canCollideWithElement(int index, CollisionModel* model2, int index2)
{

    if (!this->bSelfCollision.getValue()) return true; // we need to perform this verification process only for the selfcollision case.
    if (this->getContext() != model2->getContext()) return true;


    bool debug= false;


    if (model2 == this)
    {

        if (index<=index2) // to avoid to have two times the same auto-collision we only consider the case when index > index2
            return false;

        sofa::core::topology::BaseMeshTopology* topology = this->getMeshTopology();




        // in the neighborhood, if we find a point in common, we cancel the collision
        const helper::vector <unsigned int>& verticesAroundVertex1 =topology->getVerticesAroundVertex(index);
        const helper::vector <unsigned int>& verticesAroundVertex2 =topology->getVerticesAroundVertex(index2);

        for (unsigned int i1=0; i1<verticesAroundVertex1.size(); i1++)
        {

            unsigned int v1 = verticesAroundVertex1[i1];

            for (unsigned int i2=0; i2<verticesAroundVertex2.size(); i2++)
            {

                if (debug)
                    std::cout<<"v1 = "<<v1 <<"  verticesAroundVertex2[i2]"<<verticesAroundVertex2[i2]<<std::endl;
                if (v1==verticesAroundVertex2[i2] || v1==(unsigned int)index2 || index == (int)verticesAroundVertex2[i2])
                {
                    if(debug)
                        std::cout<<" return false"<<std::endl;
                    return false;
                }
            }
        }
        if(debug)
            std::cout<<" return true"<<std::endl;
        return true;
    }
    else
        return model2->canCollideWithElement(index2, this, index);
}

template<class DataTypes>
void TPointModel<DataTypes>::computeBoundingTree(int maxDepth)
{
    CubeModel* cubeModel = createPrevious<CubeModel>();
    const int npoints = mstate->getX()->size();
    bool updated = false;
    if (npoints != size)
    {
        resize(npoints);
        updated = true;
    }
    if (updated) cubeModel->resize(0);
    if (!isMoving() && !cubeModel->empty() && !updated) return; // No need to recompute BBox if immobile

    if (computeNormals.getValue()) updateNormals();

    cubeModel->resize(size);
    if (!empty())
    {
        //VecCoord& x = *mstate->getX();
        const SReal distance = this->proximity.getValue();
        for (int i=0; i<size; i++)
        {
            TPoint<DataTypes> p(this,i);
            const Vector3& pt = p.p();
            cubeModel->setParentOf(i, pt - Vector3(distance,distance,distance), pt + Vector3(distance,distance,distance));
        }
        cubeModel->computeBoundingTree(maxDepth);
    }

    if (m_lmdFilter != 0)
    {
        m_lmdFilter->invalidate();
    }
}

template<class DataTypes>
void TPointModel<DataTypes>::computeContinuousBoundingTree(double dt, int maxDepth)
{
    CubeModel* cubeModel = createPrevious<CubeModel>();
    const int npoints = mstate->getX()->size();
    bool updated = false;
    if (npoints != size)
    {
        resize(npoints);
        updated = true;
    }
    if (!isMoving() && !cubeModel->empty() && !updated) return; // No need to recompute BBox if immobile

    if (computeNormals.getValue()) updateNormals();

    Vector3 minElem, maxElem;

    cubeModel->resize(size);
    if (!empty())
    {
        //VecCoord& x = *mstate->getX();
        //VecDeriv& v = *mstate->getV();
        const SReal distance = (SReal)this->proximity.getValue();
        for (int i=0; i<size; i++)
        {
            TPoint<DataTypes> p(this,i);
            const Vector3& pt = p.p();
            const Vector3 ptv = pt + p.v()*dt;

            for (int c = 0; c < 3; c++)
            {
                minElem[c] = pt[c];
                maxElem[c] = pt[c];
                if (ptv[c] > maxElem[c]) maxElem[c] = ptv[c];
                else if (ptv[c] < minElem[c]) minElem[c] = ptv[c];
                minElem[c] -= distance;
                maxElem[c] += distance;
            }
            cubeModel->setParentOf(i, minElem, maxElem);
        }
        cubeModel->computeBoundingTree(maxDepth);
    }
}

template<class DataTypes>
void TPointModel<DataTypes>::updateNormals()
{
    const VecCoord& x = *this->mstate->getX();
    int n = x.size();
    normals.resize(n);
    for (int i=0; i<n; ++i)
    {
        normals[i].clear();
    }
    core::topology::BaseMeshTopology* mesh = getContext()->getMeshTopology();
    if (mesh->getNbTetrahedra()+mesh->getNbHexahedra() > 0)
    {
        if (mesh->getNbTetrahedra()>0)
        {
            const core::topology::BaseMeshTopology::SeqTetrahedra &elems = mesh->getTetrahedra();
            for (unsigned int i=0; i < elems.size(); ++i)
            {
                const core::topology::BaseMeshTopology::Tetra &e = elems[i];
                const Coord& p1 = x[e[0]];
                const Coord& p2 = x[e[1]];
                const Coord& p3 = x[e[2]];
                const Coord& p4 = x[e[3]];
                Coord& n1 = normals[e[0]];
                Coord& n2 = normals[e[1]];
                Coord& n3 = normals[e[2]];
                Coord& n4 = normals[e[3]];
                Coord n;
                n = cross(p3-p1,p2-p1); n.normalize();
                n1 += n;
                n2 += n;
                n3 += n;
                n = cross(p4-p1,p3-p1); n.normalize();
                n1 += n;
                n3 += n;
                n4 += n;
                n = cross(p2-p1,p4-p1); n.normalize();
                n1 += n;
                n4 += n;
                n2 += n;
                n = cross(p3-p2,p4-p2); n.normalize();
                n2 += n;
                n4 += n;
                n3 += n;
            }
        }
        /// @TODO Hexahedra
    }
    else if (mesh->getNbTriangles()+mesh->getNbQuads() > 0)
    {
        if (mesh->getNbTriangles()>0)
        {
            const core::topology::BaseMeshTopology::SeqTriangles &elems = mesh->getTriangles();
            for (unsigned int i=0; i < elems.size(); ++i)
            {
                const core::topology::BaseMeshTopology::Triangle &e = elems[i];
                const Coord& p1 = x[e[0]];
                const Coord& p2 = x[e[1]];
                const Coord& p3 = x[e[2]];
                Coord& n1 = normals[e[0]];
                Coord& n2 = normals[e[1]];
                Coord& n3 = normals[e[2]];
                Coord n;
                n = cross(p2-p1,p3-p1); n.normalize();
                n1 += n;
                n2 += n;
                n3 += n;
            }
        }
        if (mesh->getNbQuads()>0)
        {
            const core::topology::BaseMeshTopology::SeqQuads &elems = mesh->getQuads();
            for (unsigned int i=0; i < elems.size(); ++i)
            {
                const core::topology::BaseMeshTopology::Quad &e = elems[i];
                const Coord& p1 = x[e[0]];
                const Coord& p2 = x[e[1]];
                const Coord& p3 = x[e[2]];
                const Coord& p4 = x[e[3]];
                Coord& n1 = normals[e[0]];
                Coord& n2 = normals[e[1]];
                Coord& n3 = normals[e[2]];
                Coord& n4 = normals[e[3]];
                Coord n;
                n = cross(p3-p1,p4-p2); n.normalize();
                n1 += n;
                n2 += n;
                n3 += n;
                n4 += n;
            }
        }
    }
    for (int i=0; i<n; ++i)
    {
        SReal l = normals[i].norm();
        if (l > 1.0e-3)
            normals[i] *= 1/l;
        else
            normals[i].clear();
    }
}

template<class DataTypes>
bool TPoint<DataTypes>::testLMD(const Vector3 &PQ, double &coneFactor, double &coneExtension)
{

    Vector3 pt = p();

    sofa::core::topology::BaseMeshTopology* mesh = this->model->getMeshTopology();
    helper::vector<Vector3> x = (*this->model->mstate->getX());

    const helper::vector <unsigned int>& trianglesAroundVertex = mesh->getTrianglesAroundVertex(this->index);
    const helper::vector <unsigned int>& edgesAroundVertex = mesh->getEdgesAroundVertex(this->index);


    Vector3 nMean;

    for (unsigned int i=0; i<trianglesAroundVertex.size(); i++)
    {
        unsigned int t = trianglesAroundVertex[i];
        const fixed_array<unsigned int,3>& ptr = mesh->getTriangle(t);
        Vector3 nCur = (x[ptr[1]]-x[ptr[0]]).cross(x[ptr[2]]-x[ptr[0]]);
        nCur.normalize();
        nMean += nCur;
    }

    if (trianglesAroundVertex.size()==0)
    {
        for (unsigned int i=0; i<edgesAroundVertex.size(); i++)
        {
            unsigned int e = edgesAroundVertex[i];
            const fixed_array<unsigned int,2>& ped = mesh->getEdge(e);
            Vector3 l = (pt - x[ped[0]]) + (pt - x[ped[1]]);
            l.normalize();
            nMean += l;
        }
    }

    if (nMean.norm()> 0.0000000001)
        nMean.normalize();


    for (unsigned int i=0; i<edgesAroundVertex.size(); i++)
    {
        unsigned int e = edgesAroundVertex[i];
        const fixed_array<unsigned int,2>& ped = mesh->getEdge(e);
        Vector3 l = (pt - x[ped[0]]) + (pt - x[ped[1]]);
        l.normalize();
        double computedAngleCone = dot(nMean , l) * coneFactor;
        if (computedAngleCone<0)
            computedAngleCone=0.0;
        computedAngleCone+=coneExtension;
        if (dot(l , PQ) < -computedAngleCone*PQ.norm())
            return false;
    }
    return true;


}

template<class DataTypes>
PointLocalMinDistanceFilter *TPointModel<DataTypes>::getFilter() const
{
    return m_lmdFilter;
}

template<class DataTypes>
void TPointModel<DataTypes>::setFilter(PointLocalMinDistanceFilter *lmdFilter)
{
    m_lmdFilter = lmdFilter;
}


} // namespace collision

} // namespace component

} // namespace sofa

#endif
