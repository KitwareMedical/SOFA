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
*                               SOFA :: Plugins                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef OGRESOFAVIEWER_H
#define OGRESOFAVIEWER_H

#include <sofa/gui/qt/viewer/VisualModelPolicy.h>
#include <sofa/gui/qt/viewer/SofaViewer.h>
#include "DrawToolOGRE.h"
namespace sofa
{
namespace gui
{
namespace qt
{
namespace viewer
{

class OgreVisualModelPolicy : public VisualModelPolicy
{
protected:
    sofa::core::ObjectFactory::ClassEntry* classVisualModel;
    sofa::core::ObjectFactory::ClassEntry* classOglModel;
    sofa::core::visual::DrawToolOGRE drawToolOGRE;
public:
    void load()
    {
        // Replace OpenGL visual models with OgreVisualModel
        sofa::core::ObjectFactory::AddAlias("OglModel", "OgreVisualModel", true, &classOglModel);
        sofa::core::ObjectFactory::AddAlias("VisualModel", "OgreVisualModel", true, &classVisualModel);
        vparams->drawTool() = &drawToolOGRE;
//    vparams->setSupported(sofa::core::visual::API_OGRE); // ?
    }

    void unload()
    {
        sofa::core::ObjectFactory::ResetAlias("OglModel", classOglModel);
        sofa::core::ObjectFactory::ResetAlias("VisualModel", classVisualModel);
        vparams->drawTool() = NULL;
    }

};

typedef sofa::gui::qt::viewer::CustomPolicySofaViewer< OgreVisualModelPolicy > OgreSofaViewer;

}
}
}
}


#endif // OGRESOFAVIEWER_H
