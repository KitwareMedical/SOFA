#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/system/PluginManager.h>
#include <sofa/component/init.h>
#include <sofa/helper/vector.h>

#include <sofa/helper/ArgumentParser.h>
#include <sofa/core/objectmodel/BaseClass.h>
#include <sofa/helper/system/SetDirectory.h>

#include <iostream>
#include <algorithm>

namespace
{
	const std::string authors(
		"\
		/******************************************************************************\n\
		*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *\n\
		*                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *\n\
		*                                                                             *\n\
		* This library is free software; you can redistribute it and/or modify it     *\n\
		* under the terms of the GNU Lesser General Public License as published by    *\n\
		* the Free Software Foundation; either version 2.1 of the License, or (at     *\n\
		* your option) any later version.                                             *\n\
		*                                                                             *\n\
		* This library is distributed in the hope that it will be useful, but WITHOUT *\n\
		* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *\n\
		* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *\n\
		* for more details.                                                           *\n\
		*                                                                             *\n\
		* You should have received a copy of the GNU Lesser General Public License    *\n\
		* along with this library; if not, write to the Free Software Foundation,     *\n\
		* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *\n\
		*******************************************************************************\n\
		*                               SOFA :: Modules                               *\n\
		*                                                                             *\n\
		* Authors: The SOFA Team and external contributors (see Authors.txt)          *\n\
		*                                                                             *\n\
		* Contact information: contact@sofa-framework.org                             *\n\
		******************************************************************************/\n\
		"		);
	const std::string defaultIncludes(
		"\
		//Default files containing the declaration of the vector type\n\
		#include <sofa/defaulttype/VecTypes.h>\n\
		#include <sofa/defaulttype/RigidTypes.h>\n\
		#include <sofa/defaulttype/Mat.h>\n\n\
		\n\
		#ifdef SOFA_GPU_CUDA\n\
		#include <sofa/gpu/cuda/CudaTypes.h>\n\
		#endif\n\
		#ifdef SOFA_GPU_OPENCL\n\
		#include <sofa/gpu/opencl/OpenCLTypes.h>\n\
		#endif\n");
}

using namespace sofa;
using namespace sofa::core;
using namespace sofa::helper::system;

int main(int argc, char** argv)
{

	std::vector<std::string> plugins;
	std::vector<std::string> targets;
	sofa::helper::parse("This application generate header files for each SOFA binary.\n\
						Each header file contains includes to the SOFA components belonging to the binary.\n\
						A glorious Sofa.h header is also created which contains includes to all the other headers generated\n\
						by this application.\n\
						Headers are created in $SOFA_DIR/modules/sofa")
						.option(&plugins,'l',"load","load given plugins")
						.option(&targets,'t',"targets","generate typedefs for the given targets only")
						(argc,argv);

	sofa::component::init();
	for (unsigned int i=0; i<plugins.size(); i++)
		sofa::helper::system::PluginManager::getInstance().loadPlugin(plugins[i]);

	sofa::helper::system::PluginManager::getInstance().initRecentlyOpened();
	sofa::helper::system::PluginManager::getInstance().init();


	helper::vector<ObjectFactory::ClassEntry*> registry;
	helper::vector<ObjectFactory::ClassEntry*>::iterator it;

	//std::string outputPath = SetDirectory::GetRelativeFromProcess("../modules/sofa/");
	std::string outputPath = "/modules/sofa/";

	std::string fullPath = SetDirectory::GetParentDir(SetDirectory::GetParentDir(SetDirectory::GetProcessFullPath("").c_str()).c_str());
	int absoluteBasePathSize = fullPath.length();

	std::cout << "absoluteBasePathSize = " << absoluteBasePathSize << std::endl;
	std::cout << "outputPath[" << outputPath.length() << "] = " << outputPath << std::endl;

	// retrieve creators

	typedef std::map< std::string, helper::vector< ObjectFactory::CreatorList::value_type >  > TargetCreatorMap;
	typedef std::map< std::string, helper::vector< ObjectFactory::ClassEntry* >  > TargetClassEntryMap;
	TargetCreatorMap targetCreatorMap;
	TargetClassEntryMap targetClassEntryMap;

	sofa::core::ObjectFactory::getInstance()->getAllEntries(registry);

	for( it = registry.begin(); it != registry.end(); ++it)
	{
		ObjectFactory::ClassEntry* entry = *it;
		ObjectFactory::CreatorList::iterator itc;
		itc = entry->creatorList.begin();
		if( itc != entry->creatorList.end() )
		{
			targetClassEntryMap[itc->second->getTarget()].push_back(entry);
		}

		for(itc = entry->creatorList.begin(); itc != entry->creatorList.end(); ++itc )
		{
			ObjectFactory::Creator* creator = itc->second;

			if(!std::string(creator->getTarget()).empty())
				targetCreatorMap[creator->getTarget()].push_back(*itc);
		}
	}

	std::ostringstream timestamp;
	timestamp << "// File automatically generated by \"SofaTypedef\" -- " << std::string(__DATE__) << " at " << std::string(__TIME__) << std::endl;


	// generate target headers

	for( TargetCreatorMap::iterator it_target = targetCreatorMap.begin(); it_target != targetCreatorMap.end(); ++it_target)
	{
		std::ostringstream target_h_location,target_h_guarding_block;
		target_h_location << fullPath << outputPath << it_target->first << ".h";

		std::ofstream target_h(target_h_location.str().c_str());
		{
			std::string target_name(it_target->first);
			std::transform(target_name.begin(),target_name.end(),target_name.begin(),::toupper);
			target_h_guarding_block << "#ifndef " << target_name << "_H" << std::endl;
			target_h_guarding_block << "#define " << target_name << "_H" << std::endl;
		}

		target_h << authors;
		target_h << timestamp.str();
		target_h << target_h_guarding_block.str();
		target_h << defaultIncludes;

		helper::vector< ObjectFactory::CreatorList::value_type > targetCreatorList = it_target->second;
		helper::vector< ObjectFactory::CreatorList::value_type >::iterator it_creatorlist;
		std::pair<TargetClassEntryMap::iterator,TargetClassEntryMap::iterator> range = targetClassEntryMap.equal_range(it_target->first);
		TargetClassEntryMap::iterator it_vecEntry;
		for(it_vecEntry = range.first; it_vecEntry != range.second; ++it_vecEntry)
		{
			helper::vector< ObjectFactory::ClassEntry*> &  vecEntry = it_vecEntry->second;
			helper::vector< ObjectFactory::ClassEntry*>::iterator it_entry;
			for(it_entry = vecEntry.begin(); it_entry != vecEntry.end(); ++it_entry)
			{
				ObjectFactory::ClassEntry* entry = *it_entry;
				if( entry->creatorList.begin() != entry->creatorList.end() )
				{
					ObjectFactory::Creator* creator = entry->creatorList.begin()->second;
					std::cout << "getHeaderFileLocation["<< std::string(creator->getHeaderFileLocation()).length() << "] = " << creator->getHeaderFileLocation() << std::endl;
					std::string includePath = std::string(creator->getHeaderFileLocation());
#if defined(WIN32)
					includePath = std::string("sofa/")+includePath.substr(absoluteBasePathSize + outputPath.length());
#else
					std::size_t pos = includePath.find(outputPath);
					if (pos != std::string::npos)
						includePath = std::string("sofa/") + includePath.substr(pos+outputPath.length());
					else
						includePath = std::string("sofa/component/") + includePath;
#endif
					std::cout << " => " << includePath << std::endl;
					target_h << "#include <" << includePath << ">" << std::endl;
				}
			}
		}

		target_h << std::endl;

		for(it_creatorlist = targetCreatorList.begin(); it_creatorlist != targetCreatorList.end(); ++it_creatorlist)
		{
			ObjectFactory::Creator* creator = it_creatorlist->second;
			std::string templateName = creator->getClass()->templateName;

			{
				std::string::iterator itEnd;

				itEnd = std::remove(templateName.begin(),templateName.end(),'[');
				templateName.erase(itEnd, templateName.end());

				itEnd = std::remove(templateName.begin(),templateName.end(),']');
				templateName.erase(itEnd, templateName.end());
			}

			std::replace(templateName.begin(),templateName.end(),',','_');
			std::replace(templateName.begin(),templateName.end(),' ','_');
			std::replace(templateName.begin(),templateName.end(),'<','_');
			std::replace(templateName.begin(),templateName.end(),'>','_');

			// trim right underscores
			{
				for(unsigned size = templateName.size(); size != 0; --size)
				{
					if('_' != templateName[size - 1])
					{
						if(templateName.size() != size)
							templateName.resize(size);

						break;
					}
				}
			}

			const std::type_info& type = creator->type();

			size_t curPos = 0;
			size_t oldPos = 0;
			std::string namespaceName = creator->getClass()->namespaceName;
			size_t bracketCount = 0;
			helper::vector<std::string> namespaces;
			while( ( curPos = namespaceName.find(std::string("::"),oldPos) ) != std::string::npos)
			{
				std::string currentNamespace = namespaceName.substr(oldPos,curPos-oldPos);
				namespaces.push_back(currentNamespace);
				oldPos = curPos+2;
				++bracketCount;
			}
			namespaces.push_back(namespaceName.substr(oldPos));
			helper::vector<std::string>::iterator it_namespace;
			for(it_namespace = namespaces.begin(); it_namespace != namespaces.end(); ++it_namespace)
			{
				target_h << "namespace " << *it_namespace << std::endl;
				target_h << "{" << std::endl;
			}

			if(!templateName.empty())
			{
				target_h << "typedef " <<  sofa::core::objectmodel::BaseClass::decodeFullName(type)
					<< " " << creator->getClass()->className << "_" << templateName << ";";
			}
			else
			{
				target_h << "typedef " <<  sofa::core::objectmodel::BaseClass::decodeFullName(type)
					<< " " << creator->getClass()->className << ";";
			}
			target_h << std::endl;
			for( size_t i=0; i<namespaces.size(); ++i)
			{
				target_h << "}" << std::endl;
			}

			target_h << std::endl;
		}
		target_h << "#endif " << std::endl;
		target_h.close();
	}

	// generate sofa.h header

	std::string sofa_name("sofa");
	std::ostringstream sofa_h_location,sofa_h_guarding_block;
	sofa_h_location << fullPath << outputPath << sofa_name << ".h";
	std::ofstream sofa_h(sofa_h_location.str().c_str());
	{
		std::transform(sofa_name.begin(),sofa_name.end(),sofa_name.begin(),::toupper);
		sofa_h_guarding_block << "#ifndef " << sofa_name << "_H" << std::endl;
		sofa_h_guarding_block << "#define " << sofa_name << "_H" << std::endl;
	}

	sofa_h << authors << std::endl;
	sofa_h << timestamp.str();
	sofa_h << sofa_h_guarding_block.str();

	for( TargetCreatorMap::iterator it_target = targetCreatorMap.begin(); it_target != targetCreatorMap.end(); ++it_target)
	{
		sofa_h << "#include <sofa/" << it_target->first << ".h>" << std::endl;
	}

	sofa_h << std::endl << "#endif" << std::endl;
	sofa_h.close();
}
