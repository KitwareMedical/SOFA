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

#include "SofaLibrary.h"
#include <sofa/core/ObjectFactory.h>


//Automatically create and destroy all the components available: easy way to verify the default constructor and destructor
//#define TEST_CREATION_COMPONENT

namespace sofa
{

namespace core
{
void SofaLibrary::build( const std::vector< std::string >& examples)
{
    exampleFiles=examples;
    //-----------------------------------------------------------------------
    //Read the content of the Object Factory
    //-----------------------------------------------------------------------
    std::vector< ClassEntry* > entries;
    sofa::core::ObjectFactory::getInstance()->getAllEntries(entries);
    //Set of categories found in the Object Factory
    std::set< std::string > mainCategories;
    //Data containing all the entries for a given category
    std::multimap< std::string, ClassEntry* > inventory;

    for (unsigned int i=0; i<entries.size(); ++i)
    {
#ifdef      TEST_CREATION_COMPONENT
        {
            sofa::core::objectmodel::BaseObject::SPtr object;
            std::cerr << "Creating " << entries[i]->className << std::endl;
            if (entries[i]->creatorMap.find(entries[i]->defaultTemplate) != entries[i]->creatorMap.end())
            {
                object = entries[i]->creatorMap.find(entries[i]->defaultTemplate)->second->createInstance(NULL, NULL);
            }
            else
            {
                object = entries[i]->creatorList.begin()->second->createInstance(NULL, NULL);
            }
            std::cerr << "Deleting " << entries[i]->className << std::endl;
            object.reset();
            std::cerr << "Ok for " << entries[i]->className << std::endl;
        }
#endif

        //Insert Template specification
        std::set< std::string >::iterator it;
        for (it = entries[i]->baseClasses.begin(); it != entries[i]->baseClasses.end(); ++it)
        {
            mainCategories.insert((*it));
            inventory.insert(std::make_pair((*it), entries[i]));
        }
        //If no inheritance was found for the given component, we store it in a default category
        if (entries[i]->baseClasses.empty())
        {
            mainCategories.insert("_Miscellaneous");
            inventory.insert(std::make_pair("_Miscellaneous", entries[i]));
        }
    }

    //-----------------------------------------------------------------------
    //Using the inventory, Add each component to the Sofa Library
    //-----------------------------------------------------------------------
    std::set< std::string >::iterator itCategory;
    typedef std::multimap< std::string, ClassEntry* >::iterator IteratorInventory;


    //We add the components category by category
    for (itCategory = mainCategories.begin(); itCategory != mainCategories.end(); ++itCategory)
    {
        const std::string& categoryName = *itCategory;
        IteratorInventory itComponent;

        std::pair< IteratorInventory,IteratorInventory > rangeCategory;
        rangeCategory = inventory.equal_range(categoryName);
        const unsigned int numComponentInCategory = inventory.count(categoryName);



        CategoryLibrary *category = createCategory(categoryName,numComponentInCategory);

        //Process all the component of the current category, and add them to the group
        for (itComponent=rangeCategory.first; itComponent != rangeCategory.second; ++itComponent)
        {
            ClassEntry *entry = itComponent->second;
            const std::string &componentName=entry->className;

            //Special Case of Mass Component: they are also considered as forcefield. We remove their occurence of the force field category group
            if (categoryName == "ForceField")
            {
                std::set< std::string >::iterator inheritanceClass;
                bool needToRemove=false;
                for (inheritanceClass = entry->baseClasses.begin(); inheritanceClass != entry->baseClasses.end(); ++inheritanceClass)
                {
                    if (*inheritanceClass == "Mass") { needToRemove=true; break;};
                }
                if (needToRemove) continue;
            }
            //Special Case of TopologyObject: they are also considered as Topology. We remove their occurence of the topology category group
            else if (categoryName == "Topology")
            {
                std::set< std::string >::iterator inheritanceClass;
                bool needToRemove=false;
                for (inheritanceClass = entry->baseClasses.begin(); inheritanceClass != entry->baseClasses.end(); ++inheritanceClass)
                {
                    if (*inheritanceClass == "TopologyObject") { needToRemove=true; break;};
                }
                if (needToRemove) continue;
            }

            //Add the component to the category
            category->addComponent(componentName, entry, exampleFiles);
        }
        category->endConstruction();
        addCategory(category);
    }
    computeNumComponents();
}

void SofaLibrary::computeNumComponents()
{
    numComponents=0;
    for (unsigned int cat=0; cat<categories.size(); ++cat)
    {
        numComponents += categories[cat]->getNumComponents();
    }

}

void SofaLibrary::addCategory(CategoryLibrary *category)
{
    categories.push_back(category);
}


std::string SofaLibrary::getComponentDescription( const std::string &componentName ) const
{
    const ComponentLibrary *component = getComponent(componentName);
    if (component) return component->getDescription();
    else return "";
}

const CategoryLibrary *SofaLibrary::getCategory( const std::string &categoryName) const
{
    for (VecCategoryIterator it=categories.begin(); it != categories.end(); ++it)
    {
        if ((*it)->getName().find(categoryName) != std::string::npos)
            return *it;
    }
    return NULL;
}

const ComponentLibrary *SofaLibrary::getComponent( const std::string &componentName ) const
{
    //Look into all the categories
    for (unsigned int cat=0; cat<categories.size(); ++cat)
    {
        //For each category, look at all the components if one has the name wanted
        const std::vector< ComponentLibrary* > &components = categories[cat]->getComponents();
        for (unsigned int comp=0; comp<components.size(); ++comp)
        {
            if (componentName == components[comp]->getName()) return components[comp];
        }
    }
    return NULL;
}

void SofaLibrary::clear()
{
    for (unsigned int i=0; i<categories.size(); ++i)
    {
        delete categories[i];
    }
    categories.clear();
}

}
}
