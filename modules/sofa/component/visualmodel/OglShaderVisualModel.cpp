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
/*
 * OglShaderVisualModel.cpp
 *
 *  Created on: 9 févr. 2009
 *      Author: froy
 */

#include <sofa/component/visualmodel/OglShaderVisualModel.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/ObjectFactory.h>

#include <sofa/core/topology/TopologyChange.h>
#include <sofa/component/topology/PointSetGeometryAlgorithms.h>
#include <sofa/simulation/common/Node.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/component/visualmodel/OglAttribute.inl>

namespace sofa
{

namespace component
{

namespace visualmodel
{

using namespace sofa::core::topology;
using namespace sofa::core::behavior;

SOFA_DECL_CLASS(OglShaderVisualModel)

int OglShaderVisualModelClass = core::RegisterObject("Visual model for OpenGL display using Glew extensions")
        .add< OglShaderVisualModel >()
        ;

OglShaderVisualModel::OglShaderVisualModel()
    : shader(NULL)
    , restPosition_lastUpdate(-1)
    , vrestpositions(NULL)
    , vrestnormals(NULL)
    , modelMatrixUniform(NULL)
{
}

OglShaderVisualModel::~OglShaderVisualModel()
{
}

void OglShaderVisualModel::pushTransformMatrix(float* matrix)
{
    OglModel::pushTransformMatrix(matrix);

    helper::vector<float> tempModelMatrixValue;

    for ( unsigned int i = 0; i < 16; i++ )
        tempModelMatrixValue.push_back(matrix[i]);
    if (modelMatrixUniform)
        modelMatrixUniform->setValue(tempModelMatrixValue);

    if (shader)
    {
        shader->stop();
        if (modelMatrixUniform)
            modelMatrixUniform->pushValue();
        shader->start();
    }
}

void OglShaderVisualModel::popTransformMatrix()
{
    OglModel::popTransformMatrix();
    /*
    	if (shader)
    		shader->stop();
    */
}

void OglShaderVisualModel::init()
{
    if (!shader)
        OglModel::init();

    sofa::core::objectmodel::BaseContext* context = this->getContext();
    shader = context->core::objectmodel::BaseContext::get<OglShader>();

    if( shader)
    {
        if (!vrestpositions)
            context->get(vrestpositions, "restPosition");
        if (!vrestnormals)
            context->get(vrestnormals, "restNormal");
        if (!modelMatrixUniform)
            context->get(modelMatrixUniform, "modelMatrix");
        //add restPosition as Attribute
        if (!vrestpositions)
        {
            vrestpositions = new OglFloat3Attribute;
            vrestpositions->setName("restPosition");
            this->getContext()->addObject(vrestpositions);
            vrestpositions->setID( std::string("restPosition"));
            vrestpositions->setIndexShader(0);
            vrestpositions->init();
        }

        if (!vrestnormals)
        {
            vrestnormals = new OglFloat3Attribute;
            vrestnormals->setName("restNormal");
            this->getContext()->addObject(vrestnormals);
            vrestnormals->setID( std::string("restNormal"));
            vrestnormals->setIndexShader(0);
            vrestnormals->init();
        }
        /*
            ResizableExtVector<Coord>& vrestpos = * ( vrestpositions.beginEdit() );
            const ResizableExtVector<Coord>& vertices = m_vertices.getValue();
            vrestpos.resize (vertices.size() );
            for ( unsigned int i = 0; i < vertices.size(); i++ )
            {
              vrestpos[i] = vertices[i];
            }
            vrestpositions.endEdit();

            //add restNormal as Attribute
            vrestnormals.setContext( this->getContext());
        	vrestnormals.setID( std::string("restNormal") );
        	vrestnormals.setIndexShader( 0);
            vrestnormals.init();

            computeRestNormals();
        */
        computeRestPositions();
//
//    //add Model Matrix as Uniform
        if (!modelMatrixUniform)
        {
            modelMatrixUniform = new OglMatrix4Variable;
            modelMatrixUniform->setName("modelMatrix");
            this->getContext()->addObject(modelMatrixUniform);
            modelMatrixUniform->setID( std::string("modelMatrix") );
            modelMatrixUniform->setIndexShader( 0);
            modelMatrixUniform->init();
        }
    }
}


void OglShaderVisualModel::initVisual()
{
    OglModel::initVisual();
}

void OglShaderVisualModel::updateVisual()
{
    OglModel::updateVisual();
    computeRestPositions();
}

void OglShaderVisualModel::computeRestPositions()
{
    if (!vrestpositions) return;
    int counter = m_restPositions.getCounter();
    if (counter == restPosition_lastUpdate)
        return;
    restPosition_lastUpdate = counter;

    helper::ReadAccessor< Data<ResizableExtVector<Coord> > > positions = (m_positions.getValue().size() != m_restPositions.getValue().size()) ? m_positions : m_restPositions;
    ResizableExtVector<Coord>& vrestpos = * ( vrestpositions->beginEdit() );
    vrestpos.resize ( positions.size() );

    for ( unsigned int i = 0; i < positions.size(); i++ )
    {
        vrestpos[i] = positions[i];
    }

    vrestpositions->endEdit();
    computeRestNormals();
}

void OglShaderVisualModel::handleTopologyChange()
{
    //init();
    OglModel::handleTopologyChange();

    if (m_topology && shader)
    {
        bool update=false;
        std::list<const TopologyChange *>::const_iterator itBegin=m_topology->beginChange();
        std::list<const TopologyChange *>::const_iterator itEnd=m_topology->endChange();

        while( itBegin != itEnd )
        {
            core::topology::TopologyChangeType changeType = (*itBegin)->getChangeType();
            if ((changeType==core::topology::TRIANGLESREMOVED) ||
                (changeType==core::topology::TRIANGLESADDED) ||
                (changeType==core::topology::QUADSADDED) ||
                (changeType==core::topology::QUADSREMOVED))
                update=true;
            itBegin++;
        }
        if (update)
        {
            computeRestPositions();
            computeRestNormals();
//            std::cerr<< "OglShaderVisualModel - Updating Rest Normals"<<std::endl;
        }
    }

    //TODO// update the rest position when inserting a point. Then, call computeNormals() to update the attributes.
    // Not done here because we don't have the rest position of the model.
    // For the moment, the only class using dynamic topology is HexaToTriangleTopologicalMapping which update itself the attributes...
}

void OglShaderVisualModel::bwdDraw(core::visual::VisualParams* /*vp*/)
{
}

void OglShaderVisualModel::fwdDraw(core::visual::VisualParams* /*vp*/)
{
}

void OglShaderVisualModel::computeRestNormals()
{
    if (!vrestpositions || !vrestnormals) return;
    const ResizableExtVector<Coord>& vrestpos = vrestpositions->getValue();
    const ResizableExtVector<Triangle>& triangles = m_triangles.getValue();
    const ResizableExtVector<Quad>& quads = m_quads.getValue();
    ResizableExtVector<Coord>& restNormals = * ( vrestnormals->beginEdit() );
    restNormals.resize(vrestpos.size());
    for (unsigned int i = 0; i < restNormals.size(); i++)
    {
        restNormals[i].clear();
    }
    for (unsigned int i = 0; i < triangles.size() ; i++)
    {
        if (triangles[i][0] >= vrestpos.size() || triangles[i][1] >= vrestpos.size() || triangles[i][2] >= vrestpos.size())
            continue;
        const Coord  v1 = vrestpos[triangles[i][0]];
        const Coord  v2 = vrestpos[triangles[i][1]];
        const Coord  v3 = vrestpos[triangles[i][2]];
        Coord n = cross(v2-v1, v3-v1);

        n.normalize();
        restNormals[triangles[i][0]] += n;
        restNormals[triangles[i][1]] += n;
        restNormals[triangles[i][2]] += n;
    }
    for (unsigned int i = 0; i < quads.size() ; i++)
    {
        if (quads[i][0] >= vrestpos.size() || quads[i][1] >= vrestpos.size() || quads[i][2] >= vrestpos.size() || quads[i][3] >= vrestpos.size())
            continue;
        const Coord & v1 = vrestpos[quads[i][0]];
        const Coord & v2 = vrestpos[quads[i][1]];
        const Coord & v3 = vrestpos[quads[i][2]];
        const Coord & v4 = vrestpos[quads[i][3]];
        Coord n1 = cross(v2-v1, v4-v1);
        Coord n2 = cross(v3-v2, v1-v2);
        Coord n3 = cross(v4-v3, v2-v3);
        Coord n4 = cross(v1-v4, v3-v4);
        n1.normalize(); n2.normalize(); n3.normalize(); n4.normalize();
        restNormals[quads[i][0]] += n1;
        restNormals[quads[i][1]] += n2;
        restNormals[quads[i][2]] += n3;
        restNormals[quads[i][3]] += n4;
    }
    for (unsigned int i = 0; i < restNormals.size(); i++)
    {
        restNormals[i].normalize();
    }
    vrestnormals->endEdit();
}


} //namespace visualmodel

} //namespace component

} //namespace sofa
