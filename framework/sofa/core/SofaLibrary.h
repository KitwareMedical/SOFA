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
#ifndef SOFA_SOFALIBRARY_H
#define SOFA_SOFALIBRARY_H

#include "CategoryLibrary.h"


namespace sofa
{

namespace core
{


/**
 *  \brief An Generic Library
 *
 *  It reads the content of the Object Factory and builds a library of components sorted inside categories.
 *  This Interface is used for the Modeler mainly.
 *
 */
class SOFA_CORE_API SofaLibrary
{
public:
    typedef std::vector< CategoryLibrary* > VecCategory;
    typedef VecCategory::const_iterator VecCategoryIterator;

public:
    virtual ~SofaLibrary() {};

    virtual void build(const std::vector< std::string >& examples=std::vector< std::string >());
    virtual void clear();

    std::string getComponentDescription( const std::string &componentName) const;

    const VecCategory& getCategories() const {return categories;};

    const CategoryLibrary  *getCategory(  const std::string &categoryName ) const;
    const ComponentLibrary *getComponent( const std::string &componentName) const;
    unsigned int getNumComponents() const {return numComponents;}

protected:
    virtual CategoryLibrary *createCategory(const std::string &category , unsigned int/*  numComponent */) {return new CategoryLibrary(category);};
    virtual void addCategory(CategoryLibrary *);
    void computeNumComponents();

    VecCategory categories;
    std::vector< std::string > exampleFiles;
    int numComponents;

};

}
}

#endif
