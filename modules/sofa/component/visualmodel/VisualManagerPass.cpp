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
 * VisualManagerPass.cpp
 *
 *  Created on: 16 janv. 2012
 *      Author: Jeremy Ringard
 */

#include <sofa/component/visualmodel/VisualManagerPass.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/simulation/common/VisualVisitor.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/system/FileRepository.h>

#include <sofa/core/objectmodel/KeypressedEvent.h>
#include <sofa/core/objectmodel/KeyreleasedEvent.h>
#include <sofa/core/objectmodel/Data.h>

namespace sofa
{

namespace component
{

namespace visualmodel
{

using namespace helper::gl;
using namespace simulation;
using namespace core::visual;

SOFA_DECL_CLASS(VisualManagerPass)
//Register LightManager in the Object Factory
int VisualManagerPassClass = core::RegisterObject("VisualManagerPass")
        .add< VisualManagerPass >()
        ;

VisualManagerPass::VisualManagerPass()
    : factor(initData(&factor, (float)1.0, "factor","set the resolution factor for the output pass. default value:1.0")),
      renderToScreen(initData(&renderToScreen, "renderToScreen", "if true, this pass will be displayed on screen (only one renderPass in the scene must be defined as renderToScreen)")),
      outputName(initData(&outputName, "outputName","name the output texture"))
{
    if(factor.getValue()==0.0)
    {
        cerr<<this->getName()<<":\"factor\" attribute shall not be null. Using 1.0 instead..."<<endl;
        factor.setValue(1.0);
    }

    prerendered=false;
    //listen by default, in order to get the keys to activate shadows
    if(!f_listening.isSet())
        f_listening.setValue(true);
}

std::string VisualManagerPass::getOutputName()
{
    if(outputName.getValue().empty())
        return this->getName();
    else
        return outputName.getValue();
}

VisualManagerPass::~VisualManagerPass()
{}

bool VisualManagerPass::checkMultipass(sofa::core::objectmodel::BaseContext* con)
{
    sofa::component::visualmodel::CompositingVisualLoop* isMultipass=NULL;
    isMultipass= con->core::objectmodel::BaseContext::get<sofa::component::visualmodel::CompositingVisualLoop>();
    return (isMultipass!=NULL);
}

void VisualManagerPass::init()
{
    sofa::core::objectmodel::BaseContext* context = this->getContext();
    multiPassEnabled=checkMultipass(context);
    fbo = new FrameBufferObject(true, true, true);
}

/* herited from VisualModel */
void VisualManagerPass::initVisual()
{
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);

    passWidth = viewport[2]*(GLint)factor.getValue();
    passHeight = viewport[3]*(GLint)factor.getValue();

    fbo->init(passWidth, passHeight);
}

void VisualManagerPass::fwdDraw(core::visual::VisualParams* )
{
}

void VisualManagerPass::bwdDraw(core::visual::VisualParams* )
{
}

void VisualManagerPass::draw(const core::visual::VisualParams* )
{
}
/***************************/

void VisualManagerPass::preDrawScene(VisualParams* vp)
{
    if(renderToScreen.getValue() || (!multiPassEnabled))
        return;

    //const VisualParams::Viewport& viewport = vp->viewport();
    fbo->setSize(passWidth, passHeight);
    fbo->start();

    glViewport(0,0,passWidth,passHeight);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_BUFFER_BIT);

    //render opaque meshes
    vp->pass() = sofa::core::visual::VisualParams::Std;
    VisualDrawVisitor act ( vp );
    act.setTags(this->getTags());
    act.execute ( getContext() );
    //render transparent meshes
    vp->pass() = sofa::core::visual::VisualParams::Transparent;
    VisualDrawVisitor act2 ( vp );
    act2.setTags(this->getTags());
    act2.execute ( getContext() );

    fbo->stop();
    prerendered=true;
}

bool VisualManagerPass::drawScene(VisualParams* vp)
{
    if(!multiPassEnabled)
        return false;

    if(renderToScreen.getValue())
    {
        glViewport(0,0,passWidth,passHeight);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
        glEnable(GL_DEPTH_BUFFER_BIT);

        //render opaque meshes
        vp->pass() = sofa::core::visual::VisualParams::Std;
        VisualDrawVisitor act ( vp );
        act.setTags(this->getTags());
        act.execute ( getContext() );

        //render transparent meshes
        vp->pass() = sofa::core::visual::VisualParams::Transparent;
        VisualDrawVisitor act2 ( vp );
        act2.setTags(this->getTags());
        act2.execute ( getContext() );

        return true;
    }
    else
        return false;
}

void VisualManagerPass::postDrawScene(VisualParams* /*vp*/)
{
    prerendered=false;
}


//keyboard event management. Not sure what I'm gonna do with that for the moment, but I'm quite sure it should be usefull in the future
void VisualManagerPass::handleEvent(sofa::core::objectmodel::Event* /*event*/)
{
//   if (sofa::core::objectmodel::KeypressedEvent* ev = dynamic_cast<sofa::core::objectmodel::KeypressedEvent*>(event))
//   {
//     switch(ev->getKey())
//     {
//     case 'P':
//     break;
//     }
//   }
}

bool VisualManagerPass::hasFilledFbo()
{
    return prerendered;
}


}//namespace visualmodel

}//namespace component

}//namespace sofa
