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
 * GlslModel.h
 *
 *  Created on: 9 févr. 2009
 *      Author: froy
 */

#ifndef OGLSHADERVISUALMODEL_H_
#define OGLSHADERVISUALMODEL_H_

#include <sofa/component/visualmodel/OglModel.h>
#include <sofa/component/visualmodel/OglShader.h>
#include <sofa/component/visualmodel/OglAttribute.h>
#include <sofa/component/visualmodel/OglVariable.h>

namespace sofa
{

namespace component
{

namespace visualmodel
{


class SOFA_OPENGL_VISUAL_API OglShaderVisualModel : public OglModel
{
public:
    SOFA_CLASS(OglShaderVisualModel, OglModel);

protected:

    typedef ExtVec3fTypes::Coord Coord;
    typedef ExtVec3fTypes::VecCoord VecCoord;

    GLuint abo;
    OglShader* shader;
    int restPosition_lastUpdate;
public:
    // These attributes are public due to dynamic topologies updates.
    OglFloat3Attribute* vrestpositions;
    OglFloat3Attribute* vrestnormals;

    OglMatrix4Variable* modelMatrixUniform;
protected:
    OglShaderVisualModel();
    virtual ~OglShaderVisualModel();
public:
    void init();
    void initVisual();

    void updateVisual();

    //void putRestPositions(const Vec3fTypes::VecCoord& positions);

    virtual void bwdDraw(core::visual::VisualParams*);
    virtual void fwdDraw(core::visual::VisualParams*);

    // handle topological changes
    virtual void handleTopologyChange();
    void computeRestPositions();
    void computeRestNormals();

private:
    virtual void pushTransformMatrix(float* matrix);
    virtual void popTransformMatrix();


};

} //namespace visualmodel

} //namespace component

} //namespace sofa

#endif /* OGLSHADERVISUALMODEL_H_ */
