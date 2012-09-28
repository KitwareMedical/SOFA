/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU General Public License as published by the Free  *
* Software Foundation; either version 2 of the License, or (at your option)   *
* any later version.                                                          *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    *
* more details.                                                               *
*                                                                             *
* You should have received a copy of the GNU General Public License along     *
* with this program; if not, write to the Free Software Foundation, Inc., 51  *
* Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.                   *
*******************************************************************************
*                            SOFA :: Applications                             *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_GUI_QT_VIEWER_VISUALMODELPOLICY_H
#define SOFA_GUI_QT_VIEWER_VISUALMODELPOLICY_H

#include <sofa/gui/qt/SofaGUIQt.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/visual/DrawToolGL.h>

namespace sofa
{
namespace gui
{
namespace qt
{
namespace viewer
{

class SOFA_SOFAGUIQT_API VisualModelPolicy
{
public:
    VisualModelPolicy(core::visual::VisualParams* vparams = core::visual::VisualParams::defaultInstance())
        :vparams(vparams)
    {}
    virtual ~VisualModelPolicy() {}
    virtual void load() = 0;
    virtual void unload() = 0;
protected:
    sofa::core::visual::VisualParams* vparams;

};


class OglModelPolicy : public VisualModelPolicy
{
protected:
    sofa::core::ObjectFactory::ClassEntry* classVisualModel;
    sofa::core::visual::DrawToolGL drawTool;
public:
    void load()
    {
        // Replace generic visual models with OglModel
        sofa::core::ObjectFactory::AddAlias("VisualModel", "OglModel", true,
                &classVisualModel);
        vparams->drawTool() = &drawTool;
        vparams->setSupported(sofa::core::visual::API_OpenGL);
    }
    void unload()
    {
        sofa::core::ObjectFactory::ResetAlias("VisualModel", classVisualModel);
        vparams->drawTool() = NULL;
    }
};

}
}
}
}



#endif // SOFA_GUI_QT_VIEWER_VISUALMODELPOLICY_H
