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
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_COMPONENT_ENGINE_GROUPFILTERYOUNGMODULUS_H
#define SOFA_COMPONENT_ENGINE_GROUPFILTERYOUNGMODULUS_H

#if !defined(__GNUC__) || (__GNUC__ > 3 || (_GNUC__ == 3 && __GNUC_MINOR__ > 3))
#pragma once
#endif

#include <sofa/defaulttype/Vec.h>
#include <sofa/core/DataEngine.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/loader/PrimitiveGroup.h>

namespace sofa
{

namespace component
{

namespace engine
{

/**
 * This class returns a vector of Young modulus, according a list of groups
 */

template <class DataTypes>
class GroupFilterYoungModulus : public core::DataEngine
{
public:
    SOFA_CLASS(GroupFilterYoungModulus,core::DataEngine);

    typedef typename DataTypes::Real Real;

protected:

    GroupFilterYoungModulus();
    ~GroupFilterYoungModulus() {}
public:
    void init();
    void reinit();
    void update();

    //Input
    Data<helper::vector<sofa::core::loader::PrimitiveGroup > > f_groups;
    Data<helper::vector<unsigned int> > f_primitives; //not mandatory
    Data<helper::vector<int > > f_elementsGroup;
    //Output
    Data<helper::vector<Real> > f_youngModulus;
    //Parameters
    Data<std::string> p_mapGroupModulus;
    Data<Real> p_defaultModulus;

    virtual std::string getTemplateName() const
    {
        return templateName(this);
    }

    static std::string templateName(const GroupFilterYoungModulus<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }

};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_ENGINE_GROUPFILTERYOUNGMODULUS_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_ENGINE_API GroupFilterYoungModulus<defaulttype::Vec3dTypes>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
extern template class SOFA_ENGINE_API GroupFilterYoungModulus<defaulttype::Vec3fTypes>;
#endif //SOFA_DOUBLE
#endif

} // namespace engine

} // namespace component

} // namespace sofa

#endif
