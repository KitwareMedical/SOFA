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

#include "GUIManager.h"
#include "BaseGUI.h"
#include <sofa/component/init.h>
#include <sofa/simulation/common/xml/initXml.h>

namespace sofa
{

namespace gui
{

/*STATIC FIELD DEFINITIONS */
BaseGUI* GUIManager::currentGUI = NULL;
std::list<GUIManager::GUICreator> GUIManager::guiCreators;
std::vector<std::string> GUIManager::guiOptions;
const char* GUIManager::valid_guiname = NULL;


BaseGUI* GUIManager::getGUI()
{
    return currentGUI;
}

void GUIManager::AddGUIOption(const char* option)
{
    guiOptions.push_back(option);
}

const std::string &GUIManager::GetCurrentGUIName() {return currentGUI->GetGUIName();};

int GUIManager::RegisterGUI(const char* name, CreateGUIFn* creator, InitGUIFn* init, int priority)
{
	if(guiCreators.size())
	{
		std::list<GUICreator>::iterator it = guiCreators.begin();
		std::list<GUICreator>::iterator itend = guiCreators.end();
		while (it != itend && strcmp(name, it->name))
			++it;
		if (it != itend)
		{
			std::cerr << "ERROR(GUIManager): GUI "<<name<<" duplicate registration."<<std::endl;
			return 1;
		}
	}

    GUICreator entry;
    entry.name = name;
    entry.creator = creator;
    entry.init = init;
    entry.priority = priority;
    guiCreators.push_back(entry);
    return 0;
}

std::vector<std::string> GUIManager::ListSupportedGUI()
{
    std::vector<std::string> names;
    for(std::list<GUICreator>::iterator it = guiCreators.begin(), itend = guiCreators.end(); it != itend; ++it)
    {
        names.push_back(it->name);
    }
    return names;
}

std::string GUIManager::ListSupportedGUI(char separator)
{
    std::string names;
    bool first = true;
    for(std::list<GUICreator>::iterator it =guiCreators.begin(), itend =guiCreators.end(); it != itend; ++it)
    {
        if (!first) names += separator; else first = false;
        names += it->name;
    }
    return names;
}

const char* GUIManager::GetValidGUIName()
{
    const char* name;
    std::string lastGuiFilename = "config/lastUsedGUI.ini";
    if (guiCreators.empty())
    {
        std::cerr << "ERROR(SofaGUI): No GUI registered."<<std::endl;
        return NULL;
    }
    else
    {
        //Check the config file for the last used GUI type
        if(sofa::helper::system::DataRepository.findFile(lastGuiFilename))
        {
            std::string configPath = sofa::helper::system::DataRepository.getFile(lastGuiFilename);
            std::string lastGuiName;
            std::ifstream lastGuiStream(configPath.c_str());
            std::getline(lastGuiStream,lastGuiName);
            lastGuiStream.close();

            const char* lastGuiNameChar = lastGuiName.c_str();

            // const char* lastGuiNameChar = "qt";
            std::list<GUICreator>::iterator it1 = guiCreators.begin();
            std::list<GUICreator>::iterator itend1 = guiCreators.end();
            while(++it1 != itend1)
            {
                if( strcmp(lastGuiNameChar, it1->name) == 0 )
                {
                    return it1->name;
                }
            }
            std::cerr << "WARNING(SofaGUI): Previously used GUI not registered. Using default GUI." << std::endl;
        }

        std::list<GUICreator>::iterator it =guiCreators.begin();
        std::list<GUICreator>::iterator itend =guiCreators.end();
        name = it->name;
        int prio = it->priority;
        while (++it != itend)
        {
            if (it->priority > prio)
            {
                name = it->name;
                prio = it->priority;
            }
        }
    }
    return name;
}

GUIManager::GUICreator* GUIManager::GetGUICreator(const char* name)
{
    if (!name) name = GetValidGUIName();
    std::list<GUICreator>::iterator it =guiCreators.begin();
    std::list<GUICreator>::iterator itend =guiCreators.end();
    while (it != itend && strcmp(name, it->name))
        ++it;
    if (it == itend)
    {
        std::cerr << "ERROR(GUIManager): GUI "<<name<<" creation failed."<<std::endl;
        std::cerr << "Available GUIs:" << ListSupportedGUI(' ') << std::endl;
        return NULL;
    }
    else
        return &(*it);
}

int GUIManager::Init(const char* argv0, const char* name /* = "" */)
{
    BaseGUI::SetProgramName(argv0);
    sofa::component::init();
    sofa::simulation::xml::initXml();
    GUICreator* creator;

    if (currentGUI)
        return 0; // already initialized

    if (guiCreators.empty())
    {
        std::cerr << "ERROR(SofaGUI): No GUI registered."<<std::endl;
        return 1;
    }

    if( strcmp(name,"") == 0 || name == NULL)
    {
        name = GetValidGUIName(); // get the default gui name
    }
    creator = GetGUICreator(name);
    if(!creator)
    {
        return 1;
    }
    valid_guiname = name; // at this point we must have a valid name for the gui.

    if (creator->init)
        return (*creator->init)(valid_guiname, guiOptions);
    else
        return 0;
}


int GUIManager::createGUI(sofa::simulation::Node::SPtr groot, const char* filename)
{
    if (!currentGUI)
    {
        GUICreator* creator = GetGUICreator(valid_guiname);
        if (!creator)
        {
            return 1;
        }
        currentGUI = (*creator->creator)(valid_guiname, guiOptions, groot, filename);
        if (!currentGUI)
        {
            std::cerr << "ERROR(GUIManager): GUI "<<valid_guiname<<" creation failed."<<std::endl;
            return 1;
        }
        //Save this GUI type as the last used GUI
        std::string lastGUIfileName;
        std::string path = sofa::helper::system::DataRepository.getFirstPath();
        lastGUIfileName = path.append("/config/lastUsedGUI.ini");

        std::ofstream out(lastGUIfileName.c_str(),std::ios::out);
        out << valid_guiname << std::endl;
        out.close();
    }
    return 0;
}

void GUIManager::closeGUI()
{
    if(currentGUI) currentGUI->closeGUI();
}

void GUIManager::Redraw()
{
    if (currentGUI) currentGUI->redraw();
}

sofa::simulation::Node* GUIManager::CurrentSimulation()
{
    if (currentGUI)
        return currentGUI->currentSimulation();
    else
        return NULL;
}

void GUIManager::SetScene(sofa::simulation::Node::SPtr groot, const char* filename /*=NULL*/, bool temporaryFile /*=false*/ )
{
    if (currentGUI)
    {
        currentGUI->setScene(groot,filename,temporaryFile);
        currentGUI->configureGUI(groot);
    }

}

int GUIManager::MainLoop(sofa::simulation::Node::SPtr groot, const char* filename)
{
    int ret = 0;
    if (!currentGUI)
    {
        createGUI(groot, filename);
    }
    ret = currentGUI->mainLoop();
    if (ret)
    {
        std::cerr << "ERROR(SofaGUI): GUI "<<currentGUI->GetGUIName()<<" main loop failed (code "<<ret<<")."<<std::endl;
        return ret;
    }
    return ret;
}
void GUIManager::SetDimension(int  width , int  height )
{
    if (currentGUI) currentGUI->setViewerResolution(width,height);
}
void GUIManager::SetFullScreen()
{
    if (currentGUI) currentGUI->setFullScreen();
}



}
// namespace gui

}
// namespace sofa
