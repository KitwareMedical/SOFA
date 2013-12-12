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
*                               SOFA :: Plugins                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#include "initSceneCreator.h"

namespace sofa
{

namespace component
{

	//Here are just several convenient functions to help user to know what contains the plugin

	extern "C" {
                SOFA_SceneCreator_API void initExternalModule();
                SOFA_SceneCreator_API const char* getModuleName();
                SOFA_SceneCreator_API const char* getModuleVersion();
                SOFA_SceneCreator_API const char* getModuleLicense();
                SOFA_SceneCreator_API const char* getModuleDescription();
                SOFA_SceneCreator_API const char* getModuleComponentList();
	}
	
	void initExternalModule()
	{
		static bool first = true;
		if (first)
		{
			first = false;
		}
	}

	const char* getModuleName()
	{
	  return "SceneCreator";
	}

	const char* getModuleVersion()
	{
		return "0.2";
	}

	const char* getModuleLicense()
	{
		return "LGPL";
	}


	const char* getModuleDescription()
	{
		return "TODO: replace this with the description of your plugin";
	}

	const char* getModuleComponentList()
	{
	  /// string containing the names of the classes provided by the plugin
	  return "";
	  //return "MyMappingPendulumInPlane, MyBehaviorModel, MyProjectiveConstraintSet";
	}



} 

} 

/// Use the SOFA_LINK_CLASS macro for each class, to enable linking on all platforms
//SOFA_LINK_CLASS(MyMappingPendulumInPlane)
//SOFA_LINK_CLASS(MyBehaviorModel)
//SOFA_LINK_CLASS(MyProjectiveConstraintSet)

