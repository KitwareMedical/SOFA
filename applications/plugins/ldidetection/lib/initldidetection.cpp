/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 3      *
*                (c) 2006-2008 MGH, INRIA, USTL, UJF, CNRS                    *
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
#include "initldidetection.h"



namespace sofa
{

namespace component
{

	extern "C" {
		SOFA_LDIDETECTION_API void initExternalModule();
		SOFA_LDIDETECTION_API const char* getModuleName();
                SOFA_LDIDETECTION_API const char* getModuleVersion();
                SOFA_LDIDETECTION_API const char* getModuleLicense();
		SOFA_LDIDETECTION_API const char* getModuleDescription();
		SOFA_LDIDETECTION_API const char* getModuleComponentList();
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
		return "LDI Detection";
	}

        const char* getModuleVersion()
        {
                return "0.5";
        }

        const char* getModuleDescription()
	{
		return "Layered Depth peeling Contact detection, based on GPU";
	}

	const char* getModuleLicense()
	{
		return "QPL";
	}

	const char* getModuleComponentList()
	{
		return "LayeredDepthImagesPipeline, LDIMasterSolver, ContactConstraint, LDIConstraintContact, LDIPenalityContact, LDIPenalityContactForceField, RelaxationSolver";
	}



} 

} 

SOFA_LINK_CLASS(LayeredDepthImagesPipeline);
#ifdef SOFA_HAVE_EIGEN2
SOFA_LINK_CLASS(ContactConstraint)
SOFA_LINK_CLASS(LDIConstraintContact)
#endif
SOFA_LINK_CLASS(LDIDetection)
SOFA_LINK_CLASS(LDIPenalityContact)
SOFA_LINK_CLASS(LDIPenalityContactForceField)

