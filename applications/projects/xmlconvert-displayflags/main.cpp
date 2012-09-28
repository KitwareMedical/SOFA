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
#include <tinyxml.h>
#include "xmlvisitor.h"
#include <sofa/core/visual/DisplayFlags.h>

TiXmlDocument* loadFromFile(const char *filename)
{

#ifndef WIN32
    // Reset local settings to make sure that floating-point values are interpreted correctly
    setlocale(LC_ALL,"C");
    setlocale(LC_NUMERIC,"C");
#endif
    //
    // this initialize the library and check potential ABI mismatches
    // between the version it was compiled for and the actual shared
    // library used.
    //

    TiXmlDocument* doc = new TiXmlDocument; // the resulting document tree

    // xmlSubstituteEntitiesDefault(1);

    if (!(doc->LoadFile(filename)))
    {
        std::cerr << "Failed to open " << filename << "\n" << doc->ErrorDesc() << " at line " << doc->ErrorRow() << " row " << doc->ErrorCol() << std::endl;
        delete doc;
        return NULL;
    }
    return doc;
}

using namespace sofa::core::visual;
using namespace sofa::xml;

int main(int argc, char** argv)
{
    if( argc <  2 ) return -1;
    TiXmlDocument* doc = loadFromFile(argv[1]);
    if(doc )
    {
        DiscoverNodes v_nodes;
        DiscoverDisplayFlagsVisitor v_displayFlags;
        doc->Accept(&v_nodes);
        doc->Accept(&v_displayFlags);

        for(unsigned int i=0; i < v_nodes.leaves.size(); ++i )
        {
            createVisualStyleVisitor(v_nodes.leaves[i],v_displayFlags.map_displayFlags);
        }

        for( unsigned int i=0; i < v_nodes.nodes.size(); ++i)
        {
            removeShowAttributes(v_nodes.nodes[i]);
        }

        if( v_nodes.rootNode() )
        {
            std::map<const TiXmlElement*,sofa::core::visual::DisplayFlags*>::iterator it_root;
            it_root = v_displayFlags.map_displayFlags.find(v_nodes.rootNode());
            if( it_root != v_displayFlags.map_displayFlags.end() )
            {
                DisplayFlags* root_flags = it_root->second;
                convert_false_to_neutral(*root_flags);
                TiXmlElement* visualStyle;
                if( ! root_flags->isNeutral() )
                {
                    if( (visualStyle = v_nodes.rootNode()->FirstChildElement("VisualStyle") ) == NULL )
                    {
                        TiXmlElement* visualstyle = new TiXmlElement("VisualStyle");
                        std::ostringstream oss;
                        oss << *root_flags;
                        visualstyle->SetAttribute("displayFlags",oss.str());
                        TiXmlElement* first_child = v_nodes.rootNode()->FirstChildElement();
                        if(first_child) v_nodes.rootNode()->InsertBeforeChild(first_child,*visualstyle);
                        else v_nodes.rootNode()->LinkEndChild(visualstyle);
                    }
                    else
                    {
                        std::ostringstream oss;
                        oss << *root_flags;
                        visualStyle->SetAttribute("displayFlags",oss.str());
                    }
                }
                else
                {
                    if( (visualStyle = v_nodes.rootNode()->FirstChildElement("VisualStyle")) != NULL )
                    {
                        v_nodes.rootNode()->RemoveChild(visualStyle);
                    }
                }
            }
        }

        doc->Print();
        std::cout.flush();
        delete doc;
        return 0;
    }
    else
    {
        return -1;
    }


}


