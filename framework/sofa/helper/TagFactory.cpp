/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *
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
*                              SOFA :: Framework                              *
*                                                                             *
* Authors: The SOFA Team (see Authors.txt)                                    *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#include <sofa/helper/TagFactory.h>
#include <iostream>
namespace sofa
{

namespace helper
{

TagFactory::TagFactory()
{
    tagsList.push_back(std::string("0")); // ID 0 == "0" or empty string
    // Add standard tags
    tagsList.push_back(std::string("Visual"));
}

/*SOFA_HELPER_API*/ unsigned int TagFactory::getID(std::string name)
{
    if (name.empty()) return 0;
    TagFactory * tagfac = TagFactory::getInstance();
    std::vector<std::string>::iterator it = tagfac->tagsList.begin();
    unsigned int i=0;

    while(it != tagfac->tagsList.end() && (*it)!= name)
    {
        it++;
        i++;
    }

    if (it!=tagfac->tagsList.end())
        return i;
    else
    {
        std::cout << "TagFactory: creating new tag "<<i<<": "<<name<<std::endl;
        tagfac->tagsList.push_back(name);
        return i;
    }
}

/*SOFA_HELPER_API*/ std::string TagFactory::getName(unsigned int id)
{
    if( id < getInstance()->tagsList.size() )
        return getInstance()->tagsList[id];
    else
        return "";
}

/*SOFA_HELPER_API*/ TagFactory* TagFactory::getInstance()
{
    static TagFactory instance;
    return &instance;
}


} // namespace helper

} // namespace sofa

