/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 4      *
*                (c) 2006-2009 MGH, INRIA, USTL, UJF, CNRS                    *
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
#include "initMultiThreading.h"


namespace sofa
{

namespace component
{

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
		return "MultiThreading";
	}

	const char* getModuleVersion()
	{
		return "1.0";
	}

	const char* getModuleLicense()
	{
		return "LGPL";
	}

	const char* getModuleDescription()
	{
		return "MultiThreading SOFA Framework";
	}

	const char* getModuleComponentList()
	{
		return "DataExchange, AnimationLoopParallelScheduler ";
	}





} // namespace component


} // namespace sofa

#ifdef SOFA_HAVE_BOOST
	SOFA_LINK_CLASS(AnimationLoopParallelScheduler)
#endif
SOFA_LINK_CLASS(DataExchange)

