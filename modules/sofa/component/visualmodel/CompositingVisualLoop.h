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
* Authors: The SOFA Team (see Authors.txt)                                    *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
/*
 * CompositingVisualLoop.h
 *
 *  Created on: 16 janv. 2012
 *      Author: Jeremy Ringard
 */

#ifndef SOFA_SIMULATION_COMPOSITINGVISUALLOOP_H
#define SOFA_SIMULATION_COMPOSITINGVISUALLOOP_H

#include <sofa/component/component.h>
#include <sofa/simulation/common/DefaultVisualManagerLoop.h>
#include <sofa/core/visual/VisualParams.h>

#ifdef SOFA_HAVE_GLEW
#include <sofa/component/visualmodel/OglShader.h>
#include <sofa/helper/gl/FrameBufferObject.h>
#include <sofa/component/visualmodel/VisualManagerPass.h>
#endif

#include <sofa/core/objectmodel/DataFileName.h>
#include <sofa/core/objectmodel/Event.h>

using namespace sofa::core::objectmodel;
using namespace sofa::core::behavior;
using namespace sofa::simulation;

namespace sofa
{

namespace component
{

namespace visualmodel
{

/**
 *  \Compositing visual loop: render multiple passes and composite them into one single rendered frame
 */

class SOFA_OPENGL_VISUAL_API CompositingVisualLoop : public simulation::DefaultVisualManagerLoop
{
public:
    SOFA_CLASS(CompositingVisualLoop,simulation::DefaultVisualManagerLoop);

    ///Files where vertex shader is defined
    sofa::core::objectmodel::DataFileName vertFilename;
    ///Files where fragment shader is defined
    sofa::core::objectmodel::DataFileName fragFilename;

private:

    void traceFullScreenQuad();
    void defaultRendering(sofa::core::visual::VisualParams* vparams);

protected:
    CompositingVisualLoop(simulation::Node* gnode = NULL);

    virtual ~CompositingVisualLoop();

public:

    virtual void init();
    virtual void initVisual();
    virtual void drawStep(sofa::core::visual::VisualParams* vparams);
};

} // namespace visualmodel

} // namespace component

} //sofa
#endif  /* SOFA_SIMULATION_COMPOSITINGVISUALLOOP_H */
