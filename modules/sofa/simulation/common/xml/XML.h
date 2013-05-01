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
#ifndef SOFA_SIMULATION_COMMON_XML_XML_H
#define SOFA_SIMULATION_COMMON_XML_XML_H

#include <sofa/simulation/common/common.h>
#include <sofa/simulation/common/xml/Element.h>

#ifdef SOFA_XML_PARSER_TINYXML
#ifndef TIXML_USE_STL
#define TIXML_USE_STL
#endif
#include <tinyxml.h>
#endif
#ifdef SOFA_XML_PARSER_LIBXML
#include <libxml/parser.h>
#include <libxml/tree.h>
#endif


namespace sofa
{

namespace simulation
{

namespace xml
{

#ifdef SOFA_XML_PARSER_TINYXML
SOFA_SIMULATION_COMMON_API BaseElement* processXMLLoading(const char *filename, const TiXmlDocument &doc);
#endif
#ifdef SOFA_XML_PARSER_LIBXML
SOFA_SIMULATION_COMMON_API BaseElement* processXMLLoading(const char *filename, const xmlDocPtr &doc);
#endif

SOFA_SIMULATION_COMMON_API BaseElement* loadFromFile(const char *filename);

SOFA_SIMULATION_COMMON_API BaseElement* loadFromMemory(const char *filename, const char *data, unsigned int size );


SOFA_SIMULATION_COMMON_API bool save(const char *filename, BaseElement* root);

extern int SOFA_SIMULATION_COMMON_API numDefault;

} // namespace xml

} // namespace simulation

} // namespace sofa

#endif
