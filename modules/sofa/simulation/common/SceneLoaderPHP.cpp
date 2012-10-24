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
#include "SceneLoaderPHP.h"

#include <sofa/simulation/common/SceneLoaderXML.h>
#include <sofa/helper/system/PipeProcess.h>
#include <sofa/simulation/common/xml/NodeElement.h>

namespace sofa
{

namespace simulation
{

// register the loader in the factory
static SceneLoader* loaderPHP = SceneLoaderFactory::getInstance()->addEntry(new SceneLoaderPHP());




bool SceneLoaderPHP::canLoadFileExtension(const char *extension)
{
    std::string ext = extension;
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    return (ext=="php" || ext=="pscn");
}

/// get the file type description
std::string SceneLoaderPHP::getFileTypeDesc()
{
    return "Php Scenes";
}

/// get the list of file extensions
void SceneLoaderPHP::getExtensionList(ExtensionList* list)
{
    list->clear();
    list->push_back("pscn");
//    list->push_back("php");
}


sofa::simulation::Node::SPtr SceneLoaderPHP::load(const char *filename)
{
    sofa::simulation::Node::SPtr root;

    if (!canLoadFileName(filename))
        return 0;

    std::string out="",error="";
    std::vector<std::string> args;


    //TODO : replace when PipeProcess will get file as stdin
    //at the moment, the filename is given as an argument
    args.push_back(std::string("-f" + std::string(filename)));
    //args.push_back("-w");
    std::string newFilename="";
    //std::string newFilename=filename;

    helper::system::FileRepository fp("PATH", ".");
#ifdef WIN32
    std::string command = "php.exe";
#else
    std::string command = "php";
#endif
    if (!fp.findFile(command,""))
    {
        std::cerr << "Simulation : Error : php not found in your PATH environment" << std::endl;
        return NULL;
    }

    sofa::helper::system::PipeProcess::executeProcess(command.c_str(), args,  newFilename, out, error);

    if(error != "")
    {
        std::cerr << "Simulation : load : "<< error << std::endl;
        if (out == "")
            return NULL;
    }
    root = loadFromMemory(filename, out.c_str(), out.size());

    return root;
}

/// Load from a string in memory
Node::SPtr SceneLoaderPHP::loadFromMemory ( const char *filename, const char *data, unsigned int size )
{
    //::sofa::simulation::init();
    // 				std::cerr << "Loading simulation XML file "<<filename<<std::endl;
    xml::BaseElement* xml = xml::loadFromMemory (filename, data, size );

    Node::SPtr root = SceneLoaderXML::processXML(xml, filename);

    // 				std::cout << "load done."<<std::endl;
    delete xml;

    return root;
}


} // namespace simulation

} // namespace sofa

