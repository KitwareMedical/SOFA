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
#ifndef SOFA_COMPONENT_ENGINE_SUBSETTOPOLOGY_INL
#define SOFA_COMPONENT_ENGINE_SUBSETTOPOLOGY_INL

#if !defined(__GNUC__) || (__GNUC__ > 3 || (_GNUC__ == 3 && __GNUC_MINOR__ > 3))
#pragma once
#endif

#include <sofa/component/engine/SubsetTopology.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/gl/template.h>
#include <sofa/helper/gl/BasicShapes.h>
#include <sofa/defaulttype/BoundingBox.h>
#include <limits>

namespace sofa
{

namespace component
{

namespace engine
{

using namespace sofa::helper;
using namespace sofa::defaulttype;
using namespace core::objectmodel;
using namespace core::topology;

template <class DataTypes>
SubsetTopology<DataTypes>::SubsetTopology()
    : boxes( initData(&boxes, "box", "Box defined by xmin,ymin,zmin, xmax,ymax,zmax") )
    , centers( initData(&centers, "centers", "Center(s) of the sphere(s)") )
    , radii( initData(&radii, "radii", "Radius(i) of the sphere(s)") )
    , direction( initData(&direction, "direction", "Edge direction(if edgeAngle > 0)") )
    , normal( initData(&normal, "normal", "Normal direction of the triangles (if triAngle > 0)") )
    , edgeAngle( initData(&edgeAngle, (Real)0, "edgeAngle", "Max angle between the direction of the selected edges and the specified direction") )
    , triAngle( initData(&triAngle, (Real)0, "triAngle", "Max angle between the normal of the selected triangle and the specified normal direction") )
    , f_X0( initData (&f_X0, "rest_position", "Rest position coordinates of the degrees of freedom") )
    , f_edges(initData (&f_edges, "edges", "Edge Topology") )
    , f_triangles(initData (&f_triangles, "triangles", "Triangle Topology") )
    , f_tetrahedra(initData (&f_tetrahedra, "tetrahedra", "Tetrahedron Topology") )
    , f_indices( initData(&f_indices,"indices","Indices of the points contained in the ROI") )
    , f_edgeIndices( initData(&f_edgeIndices,"edgeIndices","Indices of the edges contained in the ROI") )
    , f_triangleIndices( initData(&f_triangleIndices,"triangleIndices","Indices of the triangles contained in the ROI") )
    , f_tetrahedronIndices( initData(&f_tetrahedronIndices,"tetrahedronIndices","Indices of the tetrahedra contained in the ROI") )
    , f_pointsInROI( initData(&f_pointsInROI,"pointsInROI","Points contained in the ROI") )
    , f_pointsOutROI( initData(&f_pointsOutROI,"pointsOutROI","Points out of the ROI") )
    , f_edgesInROI( initData(&f_edgesInROI,"edgesInROI","Edges contained in the ROI") )
    , f_edgesOutROI( initData(&f_edgesOutROI,"edgesOutROI","Edges out of the ROI") )
    , f_trianglesInROI( initData(&f_trianglesInROI,"trianglesInROI","Triangles contained in the ROI") )
    , f_trianglesOutROI( initData(&f_trianglesOutROI,"trianglesOutROI","Triangles out of the ROI") )
    , f_tetrahedraInROI( initData(&f_tetrahedraInROI,"tetrahedraInROI","Tetrahedra contained in the ROI") )
    , f_tetrahedraOutROI( initData(&f_tetrahedraOutROI,"tetrahedraOutROI","Tetrahedra out of the ROI") )
    , f_nbrborder( initData(&f_nbrborder,(unsigned int)0,"nbrborder","If localIndices option is activated, will give the number of vertices on the border of the ROI (being the n first points of each output Topology). ") )
    , p_localIndices( initData(&p_localIndices,false,"localIndices","If true, will compute local dof indices in topological elements") )
    , p_drawROI( initData(&p_drawROI,false,"drawROI","Draw ROI") )
    , p_drawPoints( initData(&p_drawPoints,false,"drawPoints","Draw Points") )
    , p_drawEdges( initData(&p_drawEdges,false,"drawEdges","Draw Edges") )
    , p_drawTriangles( initData(&p_drawTriangles,false,"drawTriangle","Draw Triangles") )
    , p_drawTetrahedra( initData(&p_drawTetrahedra,false,"drawTetrahedra","Draw Tetrahedra") )
    , _drawSize( initData(&_drawSize,0.0,"drawSize","rendering size for box and topological elements") )
{
    boxes.beginEdit()->push_back(Vec6(0,0,0,1,1,1));
    boxes.endEdit();
    addAlias(&f_X0,"rest_position");

    f_indices.beginEdit()->push_back(0);
    f_indices.endEdit();
    typeROI = BOX;
}

template <class DataTypes>
void SubsetTopology<DataTypes>::init()
{
    if (centers.isSet())
        typeROI = SPHERE;

    if (!f_X0.isSet())
    {
        MechanicalState<DataTypes>* mstate;
        this->getContext()->get(mstate);
        if (mstate)
        {
            BaseData* parent = mstate->findField("position");
            if (parent)
            {
                f_X0.setParent(parent);
                f_X0.setReadOnly(true);
            }
        }
        else
        {
            core::loader::MeshLoader* loader = NULL;
            this->getContext()->get(loader);
            if (loader)
            {
                BaseData* parent = loader->findField("rest_position");
                if (parent)
                {
                    f_X0.setParent(parent);
                    f_X0.setReadOnly(true);
                }
            }
        }
    }
    if (!f_edges.isSet() || !f_triangles.isSet() || !f_tetrahedra.isSet())
    {
        BaseMeshTopology* topology;
        this->getContext()->get(topology);
        if (topology)
        {
            if (!f_edges.isSet())
            {
                BaseData* eparent = topology->findField("edges");
                if (eparent)
                {
                    f_edges.setParent(eparent);
                    f_edges.setReadOnly(true);
                }
            }
            if (!f_triangles.isSet())
            {
                BaseData* tparent = topology->findField("triangles");
                if (tparent)
                {
                    f_triangles.setParent(tparent);
                    f_triangles.setReadOnly(true);
                }
            }
            if (!f_tetrahedra.isSet())
            {
                BaseData* tparent = topology->findField("tetrahedra");
                if (tparent)
                {
                    f_tetrahedra.setParent(tparent);
                    f_tetrahedra.setReadOnly(true);
                }
            }
        }
    }

    addInput(&centers);
    addInput(&radii);
    addInput(&direction);
    addInput(&normal);
    addInput(&edgeAngle);
    addInput(&triAngle);

    addInput(&f_X0);
    addInput(&f_edges);
    addInput(&f_triangles);
    addInput(&f_tetrahedra);

    addOutput(&f_indices);
    addOutput(&f_edgeIndices);
    addOutput(&f_triangleIndices);
    addOutput(&f_tetrahedronIndices);
    addOutput(&f_pointsInROI);
    addOutput(&f_pointsOutROI);
    addOutput(&f_edgesInROI);
    addOutput(&f_edgesOutROI);
    addOutput(&f_trianglesInROI);
    addOutput(&f_trianglesOutROI);
    addOutput(&f_tetrahedraInROI);
    addOutput(&f_tetrahedraOutROI);
    setDirtyValue();
}

template <class DataTypes>
void SubsetTopology<DataTypes>::reinit()
{
    update();
}

template <class DataTypes>
bool SubsetTopology<DataTypes>::isPointInROI(const CPos &p, unsigned int idROI)
{
    if (typeROI == 0)
    {
        const Vec6& b = boxes.getValue()[idROI];
        return ( p[0] >= b[0] && p[0] <= b[3] && p[1] >= b[1] && p[1] <= b[4] && p[2] >= b[2] && p[2] <= b[5] );
    }
    else
    {
        const Vec3& c = centers.getValue()[idROI];
        const Real& r = radii.getValue()[idROI];

        if((p-c).norm() > r)
            return false;
        else
            return true;
    }
}

template <class DataTypes>
bool SubsetTopology<DataTypes>::isPointInROI(const PointID &pid, unsigned int idROI)
{
    const VecCoord* x0 = &f_X0.getValue();
    CPos p =  DataTypes::getCPos((*x0)[pid]);
    return ( isPointInROI(p,idROI) );
}

template <class DataTypes>
bool SubsetTopology<DataTypes>::isEdgeInROI(const Edge &e, unsigned int idROI)
{
    const VecCoord* x0 = &f_X0.getValue();
    for (unsigned int i=0; i<2; ++i)
    {
        CPos p =  DataTypes::getCPos((*x0)[e[i]]);
        if (!isPointInROI(p, idROI))
            return false;
    }
    return true;
}

template <class DataTypes>
bool SubsetTopology<DataTypes>::isTriangleInROI(const Triangle &t, unsigned int idROI)
{
    const VecCoord* x0 = &f_X0.getValue();
    for (unsigned int i=0; i<3; ++i)
    {
        CPos p =  DataTypes::getCPos((*x0)[t[i]]);
        if (!isPointInROI(p, idROI))
            return false;
    }
    return true;
}

template <class DataTypes>
bool SubsetTopology<DataTypes>::isTetrahedronInROI(const Tetra &t, unsigned int idROI)
{
    const VecCoord* x0 = &f_X0.getValue();
    for (unsigned int i=0; i<4; ++i)
    {
        CPos p =  DataTypes::getCPos((*x0)[t[i]]);
        if (!isPointInROI(p, idROI))
            return false;
    }
    return true;
}

template <class DataTypes>
void SubsetTopology<DataTypes>::findVertexOnBorder(const Triangle &t, unsigned int idROI)
{
    const VecCoord* x0 = &f_X0.getValue();

    sofa::helper::vector<unsigned int> onborder;
    bool findIn = false;
    bool findOut = false;

    for (unsigned int i = 0; i<3; ++i)
    {
        CPos p =  DataTypes::getCPos((*x0)[t[i]]);
        if (isPointInROI(p, idROI))
        {
            findIn = true;
            onborder.push_back(t[i]);
        }
        else
            findOut = true;
    }

    if (findIn && findOut) // triangle on the border
    {
        bool find;
        for (unsigned int j=0; j<onborder.size(); ++j)
        {
            find = false;
            for (unsigned int i = 0; i<listOnBorder.size(); ++i)
                if (onborder[j] == listOnBorder[i])
                {
                    find = true;
                    break;
                }

            if (!find)
                listOnBorder.push_back(onborder[j]);
        }
    }
}



template <class DataTypes>
void SubsetTopology<DataTypes>::findVertexOnBorder(const Tetra &t, unsigned int idROI)
{
    const VecCoord* x0 = &f_X0.getValue();

    sofa::helper::vector<unsigned int> onborder;
    bool findIn = false;
    bool findOut = false;

    for (unsigned int i = 0; i<4; ++i)
    {
        CPos p =  DataTypes::getCPos((*x0)[t[i]]);
        if (isPointInROI(p, idROI))
        {
            findIn = true;
            onborder.push_back(t[i]);
        }
        else
            findOut = true;
    }

    if (findIn && findOut) // tetrahedron on the border
    {
        bool find;
        for (unsigned int j=0; j<onborder.size(); ++j)
        {
            find = false;
            for (unsigned int i = 0; i<listOnBorder.size(); ++i)
                if (onborder[j] == listOnBorder[i])
                {
                    find = true;
                    break;
                }

            if (!find)
                listOnBorder.push_back(onborder[j]);
        }
    }
}


template <class DataTypes>
void SubsetTopology<DataTypes>::update()
{
    cleanDirty();
    unsigned int ROInum = 0;
    const helper::vector<Vec3>& cen = (centers.getValue());
    const helper::vector<Real>& rad = (radii.getValue());

    helper::vector<Vec6>& vb = *(boxes.beginEdit());
    for (unsigned int bi=0; bi<vb.size(); ++bi)
    {
        if (vb[bi][0] > vb[bi][3]) std::swap(vb[bi][0],vb[bi][3]);
        if (vb[bi][1] > vb[bi][4]) std::swap(vb[bi][1],vb[bi][4]);
        if (vb[bi][2] > vb[bi][5]) std::swap(vb[bi][2],vb[bi][5]);
    }

    ROInum = vb.size();
    boxes.endEdit();

    if (typeROI == SPHERE)
    {
        if (cen.empty())
            return;

        if (cen.size() != rad.size())
        {
            serr << "WARNING: number of sphere centers and radius doesn't match." <<sendl;
            return;
        }
        ROInum = cen.size();
    }

    // Write accessor for topological element indices in ROI
    SetIndex& indices = *(f_indices.beginEdit());
    SetIndex& edgeIndices = *(f_edgeIndices.beginEdit());
    SetIndex& triangleIndices = *(f_triangleIndices.beginEdit());
    SetIndex& tetrahedronIndices = *f_tetrahedronIndices.beginEdit();


    // Read accessor for input topology
    helper::ReadAccessor< Data<helper::vector<Edge> > > edges = f_edges;
    helper::ReadAccessor< Data<helper::vector<Triangle> > > triangles = f_triangles;
    helper::ReadAccessor< Data<helper::vector<Tetra> > > tetrahedra = f_tetrahedra;

    // Write accessor for toplogical element in ROI
    helper::WriteAccessor< Data<VecCoord > > pointsInROI = f_pointsInROI;
    helper::WriteAccessor< Data<VecCoord > > pointsOutROI = f_pointsOutROI;
    helper::WriteAccessor< Data<helper::vector<Edge> > > edgesInROI = f_edgesInROI;
    helper::WriteAccessor< Data<helper::vector<Edge> > > edgesOutROI = f_edgesOutROI;
    helper::WriteAccessor< Data<helper::vector<Triangle> > > trianglesInROI = f_trianglesInROI;
    helper::WriteAccessor< Data<helper::vector<Triangle> > > trianglesOutROI = f_trianglesOutROI;
    helper::WriteAccessor< Data<helper::vector<Tetra> > > tetrahedraInROI = f_tetrahedraInROI;
    helper::WriteAccessor< Data<helper::vector<Tetra> > > tetrahedraOutROI = f_tetrahedraOutROI;

    // Clear lists
    indices.clear();
    edgeIndices.clear();
    triangleIndices.clear();
    tetrahedronIndices.clear();

    pointsInROI.clear();
    edgesInROI.clear();
    trianglesInROI.clear();
    tetrahedraInROI.clear();
    pointsOutROI.clear();
    edgesOutROI.clear();
    trianglesOutROI.clear();
    tetrahedraOutROI.clear();

    const VecCoord* x0 = &f_X0.getValue();

    const bool local = p_localIndices.getValue();
    unsigned int cpt_in = 0, cpt_out = 0, cpt_border = 0;
    unsigned int& nbrBorder = *f_nbrborder.beginEdit();
    nbrBorder = 0;

    if (local)
    {
        if ((f_tetrahedra.getValue().empty()) )
        {
            for(unsigned int i=0 ; i<triangles.size() ; i++)
            {
                Triangle t = triangles[i];
                this->findVertexOnBorder(t, 0);
            }
        }
        else
        {
            for(unsigned int i=0 ; i<tetrahedra.size() ; i++)
            {
                Tetra t = tetrahedra[i];
                this->findVertexOnBorder(t, 0);
            }
        }

        localIndices.resize(x0->size());
        nbrBorder = listOnBorder.size();
        pointsInROI.resize(nbrBorder);
        pointsOutROI.resize(nbrBorder);
        indices.resize(nbrBorder);

        cpt_in = cpt_out = nbrBorder;

        // reverse sort listOnBorder
        std::sort(listOnBorder.begin(), listOnBorder.end());
        std::reverse (listOnBorder.begin(), listOnBorder.end());
    }

    //Points
    for( unsigned i=0; i<x0->size(); ++i )
    {
        // local reorder
        if (local)
        {
            bool find = false;
            for (int j = listOnBorder.size()-1; j>=0; j--)
            {
                if (listOnBorder[j] == i)
                {
                    localIndices[i] = cpt_border;
                    indices[cpt_border] = i;
                    pointsInROI[cpt_border] = (*x0)[i];
                    pointsOutROI[cpt_border] = (*x0)[i];
                    cpt_border++;
                    find = true;
                    break;
                }
            }

            if (find)
            {
                listOnBorder.pop_back();
                continue;
            }
        }

        bool inside = false;
        for (unsigned int bi=0; bi<ROInum; ++bi)
        {
            if (isPointInROI(i, bi))
            {
                indices.push_back(i);
                pointsInROI.push_back((*x0)[i]);
                inside = true;

                if (local)
                {
                    localIndices[i] = cpt_in;
                    cpt_in++;
                }

                break;
            }
        }

        if (!inside)
        {
            pointsOutROI.push_back((*x0)[i]);
            if (local)
            {
                localIndices[i] = cpt_out;
                cpt_out++;
            }
        }
    }

    //Edges
    for(unsigned int i=0 ; i<edges.size() ; i++)
    {
        bool inside = false;
        Edge e = edges[i];
        for (unsigned int bi=0; bi<ROInum; ++bi)
        {
            if (isEdgeInROI(e, bi))
            {
                if (local) { e[0] = localIndices[e[0]]; e[1] = localIndices[e[1]]; }
                edgeIndices.push_back(i);
                edgesInROI.push_back(e);
                inside = true;
                break;
            }
        }

        if (!inside)
        {
            if (local) { e[0] = localIndices[e[0]]; e[1] = localIndices[e[1]]; }
            edgesOutROI.push_back(e);
        }
    }

    //Triangles
    for(unsigned int i=0 ; i<triangles.size() ; i++)
    {
        bool inside = false;
        Triangle t = triangles[i];
        for (unsigned int bi=0; bi<ROInum; ++bi)
        {
            if (isTriangleInROI(t, bi))
            {
                if (local) { t[0] = localIndices[t[0]]; t[1] = localIndices[t[1]]; t[2] = localIndices[t[2]];}
                triangleIndices.push_back(i);
                trianglesInROI.push_back(t);
                inside = true;
                break;
            }
        }

        if (!inside)
        {
            if (local) { t[0] = localIndices[t[0]]; t[1] = localIndices[t[1]]; t[2] = localIndices[t[2]];}
            trianglesOutROI.push_back(t);
        }
    }

    //Tetrahedra
    for(unsigned int i=0 ; i<tetrahedra.size() ; i++)
    {
        bool inside = false;
        Tetra t = tetrahedra[i];
        for (unsigned int bi=0; bi<ROInum; ++bi)
        {
            if (isTetrahedronInROI(t, bi))
            {
                if (local) { t[0] = localIndices[t[0]]; t[1] = localIndices[t[1]]; t[2] = localIndices[t[2]]; t[3] = localIndices[t[3]];}
                tetrahedronIndices.push_back(i);
                tetrahedraInROI.push_back(t);
                inside = true;
                break;
            }
        }

        if (!inside)
        {
            if (local) { t[0] = localIndices[t[0]]; t[1] = localIndices[t[1]]; t[2] = localIndices[t[2]]; t[3] = localIndices[t[3]];}
            tetrahedraOutROI.push_back(t);
        }
    }

    f_indices.endEdit();
    f_edgeIndices.endEdit();
    f_triangleIndices.endEdit();
    f_tetrahedronIndices.endEdit();
    f_nbrborder.endEdit();
}

template <class DataTypes>
void SubsetTopology<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    if (!vparams->displayFlags().getShowBehaviorModels())
        return;

    const VecCoord* x0 = &f_X0.getValue();
    glColor3f(0.0, 1.0, 1.0);
    if( p_drawROI.getValue())
    {
        ///draw the boxes
        glBegin(GL_LINES);
        const helper::vector<Vec6>& vb=boxes.getValue();
        for (unsigned int bi=0; bi<vb.size(); ++bi)
        {
            const Vec6& b=vb[bi];
            const Real& Xmin=b[0];
            const Real& Xmax=b[3];
            const Real& Ymin=b[1];
            const Real& Ymax=b[4];
            const Real& Zmin=b[2];
            const Real& Zmax=b[5];
            glVertex3d(Xmin,Ymin,Zmin);
            glVertex3d(Xmin,Ymin,Zmax);
            glVertex3d(Xmin,Ymin,Zmin);
            glVertex3d(Xmax,Ymin,Zmin);
            glVertex3d(Xmin,Ymin,Zmin);
            glVertex3d(Xmin,Ymax,Zmin);
            glVertex3d(Xmin,Ymax,Zmin);
            glVertex3d(Xmax,Ymax,Zmin);
            glVertex3d(Xmin,Ymax,Zmin);
            glVertex3d(Xmin,Ymax,Zmax);
            glVertex3d(Xmin,Ymax,Zmax);
            glVertex3d(Xmin,Ymin,Zmax);
            glVertex3d(Xmin,Ymin,Zmax);
            glVertex3d(Xmax,Ymin,Zmax);
            glVertex3d(Xmax,Ymin,Zmax);
            glVertex3d(Xmax,Ymax,Zmax);
            glVertex3d(Xmax,Ymin,Zmax);
            glVertex3d(Xmax,Ymin,Zmin);
            glVertex3d(Xmin,Ymax,Zmax);
            glVertex3d(Xmax,Ymax,Zmax);
            glVertex3d(Xmax,Ymax,Zmin);
            glVertex3d(Xmax,Ymin,Zmin);
            glVertex3d(Xmax,Ymax,Zmin);
            glVertex3d(Xmax,Ymax,Zmax);
        }
        glEnd();
    }
    if( p_drawPoints.getValue())
    {
        ///draw points in boxes
        glBegin(GL_POINTS);
        glPointSize(5.0);
        helper::ReadAccessor< Data<VecCoord > > pointsInROI = f_pointsInROI;
        for (unsigned int i=0; i<pointsInROI.size() ; ++i)
        {
            CPos p = DataTypes::getCPos(pointsInROI[i]);
            helper::gl::glVertexT(p);
        }
        glEnd();
    }
    if( p_drawEdges.getValue())
    {
        ///draw edges in boxes
        glBegin(GL_LINES);
        helper::ReadAccessor< Data<helper::vector<Edge> > > edgesInROI = f_edgesInROI;
        for (unsigned int i=0; i<edgesInROI.size() ; ++i)
        {
            Edge e = edgesInROI[i];
            for (unsigned int j=0 ; j<2 ; j++)
            {
                CPos p = DataTypes::getCPos((*x0)[e[j]]);
                helper::gl::glVertexT(p);
            }
        }
        glEnd();
    }
    if( p_drawTriangles.getValue())
    {
        ///draw triangles in boxes
        glBegin(GL_TRIANGLES);
        helper::ReadAccessor< Data<helper::vector<Triangle> > > trianglesInROI = f_trianglesInROI;
        for (unsigned int i=0; i<trianglesInROI.size() ; ++i)
        {
            Triangle t = trianglesInROI[i];
            for (unsigned int j=0 ; j<3 ; j++)
            {
                CPos p = DataTypes::getCPos((*x0)[t[j]]);
                helper::gl::glVertexT(p);
            }
        }
        glEnd();
    }

    if( p_drawTetrahedra.getValue())
    {
        ///draw tetrahedra in boxes
        glBegin(GL_LINES);
        helper::ReadAccessor< Data<helper::vector<Tetra> > > tetrahedraInROI = f_tetrahedraInROI;
        for (unsigned int i=0; i<tetrahedraInROI.size() ; ++i)
        {
            Tetra t = tetrahedraInROI[i];
            for (unsigned int j=0 ; j<4 ; j++)
            {
                CPos p = DataTypes::getCPos((*x0)[t[j]]);
                helper::gl::glVertexT(p);
                p = DataTypes::getCPos((*x0)[t[(j+1)%4]]);
                helper::gl::glVertexT(p);
            }

            CPos p = DataTypes::getCPos((*x0)[t[0]]);
            helper::gl::glVertexT(p);
            p = DataTypes::getCPos((*x0)[t[2]]);
            helper::gl::glVertexT(p);
            p = DataTypes::getCPos((*x0)[t[1]]);
            helper::gl::glVertexT(p);
            p = DataTypes::getCPos((*x0)[t[3]]);
            helper::gl::glVertexT(p);
        }
        glEnd();
    }
}

template <class DataTypes>
void SubsetTopology<DataTypes>::computeBBox(const core::ExecParams*  params )
{
    const helper::vector<Vec6>& vb=boxes.getValue();
    const Real max_real = std::numeric_limits<Real>::max();
    const Real min_real = std::numeric_limits<Real>::min();
    Real maxBBox[3] = {min_real,min_real,min_real};
    Real minBBox[3] = {max_real,max_real,max_real};

    for (unsigned int bi=0; bi<vb.size(); ++bi)
    {
        const Vec6& b=vb[bi];
        if (b[0] < minBBox[0]) minBBox[0] = b[0];
        if (b[1] < minBBox[1]) minBBox[1] = b[1];
        if (b[2] < minBBox[2]) minBBox[2] = b[2];
        if (b[3] > maxBBox[0]) maxBBox[0] = b[3];
        if (b[4] > maxBBox[1]) maxBBox[1] = b[4];
        if (b[5] > maxBBox[2]) maxBBox[2] = b[5];
    }
    this->f_bbox.setValue(params,sofa::defaulttype::TBoundingBox<Real>(minBBox,maxBBox));
}

} // namespace engine

} // namespace component

} // namespace sofa

#endif
