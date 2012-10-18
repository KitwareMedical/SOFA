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
#ifndef SOFA_SIMULATION_SCENELOADERFACTORY_H
#define SOFA_SIMULATION_SCENELOADERFACTORY_H

#include <sofa/core/core.h>
#include <sofa/simulation/common/Node.h>
#include <sofa/helper/system/SetDirectory.h>


namespace sofa
{

namespace simulation
{

/**
 *  \brief Main class used to register scene file loaders
 *
 *  It uses the Factory design pattern, where each class is registered in a map,
 *  and dynamically retrieved given the type name.
 *
 */
class SOFA_SIMULATION_COMMON_API SceneLoaderFactory
{

public:
    class SceneLoader;

    typedef std::vector<SceneLoader*> SceneLoaderList;

    /// Abstract interface of a scene loader
    class SceneLoader
    {
    public:
        typedef std::vector<std::string> ExtensionList;

        /// Pre-loading check
        virtual bool canLoadFileName(const char *filename)
        {
            std::string ext = sofa::helper::system::SetDirectory::GetExtension(filename);
            return canLoadFileExtension(ext.c_str());
        }

        virtual bool canLoadFileExtension(const char *extension) = 0;

        /// load the file
        virtual sofa::simulation::Node::SPtr load(const char *filename) = 0;

        /// get the file type description
        virtual std::string getFileTypeDesc() = 0;

        /// get the list of file extensions
        virtual void getExtensionList(ExtensionList* list) = 0;

    };

    /// Get the ObjectFactory singleton instance
    static SceneLoaderFactory* getInstance();

protected:

    /// Main class registry
    SceneLoaderList registry;

public:
    /// Get an entry given a file extension
    SceneLoader* getEntryFileExtension(std::string extension);

    /// Get an entry given a file name
    SceneLoader* getEntryFileName(std::string filename);

    /// Add a scene loader
    SceneLoader* addEntry(SceneLoader *loader);

    /// Get the list of loaders
    SceneLoaderList* getEntries() {return &registry;}

};

} // namespace simulation

} // namespace sofa


#endif // SOFA_SIMULATION_SCENELOADERFACTORY_H

