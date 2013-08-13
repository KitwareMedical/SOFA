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
#ifndef LABELIMAGETOOLBOX_H
#define LABELIMAGETOOLBOX_H 

#include <QObject>

#include "initImage.h"
#include "ImageTypes.h"
#include "sofa/defaulttype/defaulttype.h"
#include "sofa/defaulttype/VecTypes.h"
#include <sofa/core/DataEngine.h>
#include <sofa/component/component.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/visual/VisualParams.h>

#include <sofa/core/objectmodel/Event.h>
#include <sofa/simulation/common/AnimateEndEvent.h>



//#include "labelimagetoolboxaction.h"



namespace sofa
{
namespace gui
{
namespace qt
{
class LabelImageToolBoxAction;
}
}
}

namespace sofa
{

namespace component
{

namespace engine
{

using helper::vector;
using defaulttype::Vec;
using defaulttype::Mat;
using namespace cimg_library;

/**
 * This class coorespond to a label visualized by imagetoolbox
 */

class LabelImageToolBox : public core::DataEngine
{
public:
    typedef core::DataEngine Inherited;
    SOFA_CLASS(LabelImageToolBox,Inherited);
    
    
    
    
    Data< bool > d_islinkedToToolBox;
    Data< sofa::defaulttype::Vec4d > d_color;

//    virtual std::string getTemplateName() const    { return templateName(this);    }
//    static std::string templateName(const LabelImageToolBox* = NULL) { return ImageTypes::Name();    }

    LabelImageToolBox();

    virtual void init()
    {
        //addInput(&image);
        //addOutput(&triangles);
        setDirtyValue();
    }

    virtual void reinit() { update(); }

protected:

    unsigned int time;

    virtual void update()
    {
        cleanDirty();

    }

    void handleEvent(sofa::core::objectmodel::Event */*event*/)
    {
    }

    virtual void draw(const core::visual::VisualParams* /*vparams*/)
    {
    }

public:
    
    virtual sofa::gui::qt::LabelImageToolBoxAction* createTBAction(QWidget* /*parent*/=NULL )=0;

};


} // namespace engine

} // namespace component

} // namespace sofa


#endif // LABELIMAGETOOLBOX_H
