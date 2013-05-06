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
#include "SceneLoaderXML.h"

#include <sofa/simulation/common/xml/NodeElement.h>
#include <sofa/simulation/common/FindByTypeVisitor.h>

namespace sofa
{

namespace simulation
{

// register the loader in the factory
static SceneLoader* loaderXML = SceneLoaderFactory::getInstance()->addEntry(new SceneLoaderXML());



bool SceneLoaderXML::canLoadFileExtension(const char *extension)
{
    std::string ext = extension;
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    return (ext=="xml" || ext=="scn");
}

/// get the file type description
std::string SceneLoaderXML::getFileTypeDesc()
{
    return "Scenes";
}

/// get the list of file extensions
void SceneLoaderXML::getExtensionList(ExtensionList* list)
{
    list->clear();
    list->push_back("xml");
    list->push_back("scn");
}

sofa::simulation::Node::SPtr SceneLoaderXML::load(const char *filename)
{
    sofa::simulation::Node::SPtr root;

    if (!canLoadFileName(filename))
        return 0;

    //::sofa::simulation::init();
    // 				std::cerr << "Loading simulation XML file "<<filename<<std::endl;
    xml::BaseElement* xml = xml::loadFromFile ( filename );

    root = processXML(xml, filename);

    // 				std::cout << "load done."<<std::endl;
    delete xml;

    return root;
}

/// Load a scene from a file
Node::SPtr SceneLoaderXML::processXML(xml::BaseElement* xml, const char *filename)
{
    if ( xml==NULL )
    {
        return NULL;
    }
    sofa::core::ExecParams* params = sofa::core::ExecParams::defaultInstance();

    // We go the the current file's directory so that all relative path are correct
    helper::system::SetDirectory chdir ( filename );

#if !defined(WIN32) && !defined(PS3)
    // Reset local settings to make sure that floating-point values are interpreted correctly
    setlocale(LC_ALL,"C");
    setlocale(LC_NUMERIC,"C");
#endif

    // 				std::cout << "Initializing objects"<<std::endl;
    sofa::simulation::xml::NodeElement* nodeElt = dynamic_cast<sofa::simulation::xml::NodeElement *>(xml);
    if( nodeElt==NULL )
    {
        std::cerr << "LOAD ERROR: XML Root Node is not an Element."<<std::endl;
        std::exit(1);
    }
    else if( !(nodeElt->init()) )
    {
        std::cerr << "LOAD ERROR: Node initialization failed."<<std::endl;
    }

    Node::SPtr root = dynamic_cast<Node*> ( xml->getObject() );
    if ( root == NULL )
    {
        std::cerr << "LOAD ERROR: Objects initialization failed."<<std::endl;
        return NULL;
    }

    // 				std::cout << "Initializing simulation "<<sRoot->getName() <<std::endl;

    // Find the Simulation component in the scene
    FindByTypeVisitor<Simulation> findSimu(params);
    findSimu.execute(root.get());
    if( !findSimu.found.empty() )
        setSimulation( findSimu.found[0] );

    return root;
}

/// Load from a string in memory
Node::SPtr SceneLoaderXML::loadFromMemory ( const char *filename, const char *data, unsigned int size )
{
    //::sofa::simulation::init();
    // 				std::cerr << "Loading simulation XML file "<<filename<<std::endl;
    xml::BaseElement* xml = xml::loadFromMemory (filename, data, size );

    Node::SPtr root = processXML(xml, filename);

    // 				std::cout << "load done."<<std::endl;
    delete xml;

    return root;
}


} // namespace simulation

} // namespace sofa

