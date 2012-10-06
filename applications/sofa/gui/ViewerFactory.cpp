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
#define SOFA_GUI_VIEWERFACTORY_CPP

#include <sofa/helper/Factory.h>
#include <sofa/helper/Factory.inl>
#include "ViewerFactory.h"

namespace sofa
{
namespace gui
{
}
}

namespace sofa
{
namespace helper
{

template class SOFA_SOFAGUI_API Factory< std::string, sofa::gui::BaseViewer, sofa::gui::BaseViewerArgument& >;

SofaViewerFactory*  SofaViewerFactory::getInstance()
{
    static SofaViewerFactory instance;
    return &instance;
}

const char* SofaViewerFactory::getViewerName(Key key)
{

    Creator* creator;
    std::multimap<Key, Creator*>::iterator it = this->registry.lower_bound(key);
    std::multimap<Key, Creator*>::iterator end = this->registry.upper_bound(key);
    while (it != end)
    {
        creator = (*it).second;
        const char* viewerName = creator->viewerName();
        if(viewerName != NULL )
        {
            return viewerName;
        }
        ++it;
    }
    //	std::cerr<<"Object type "<<key<<" creation failed."<<std::endl;
    return NULL;
}

const char* SofaViewerFactory::getAcceleratedViewerName(Key key)
{

    Creator* creator;
    std::multimap<Key, Creator*>::iterator it = this->registry.lower_bound(key);
    std::multimap<Key, Creator*>::iterator end = this->registry.upper_bound(key);
    while (it != end)
    {
        creator = (*it).second;
        const char* acceleratedName = creator->acceleratedName();
        if(acceleratedName != NULL )
        {
            return acceleratedName;
        }
        ++it;
    }
    //	std::cerr<<"Object type "<<key<<" creation failed."<<std::endl;
    return NULL;

}
}
}
