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
#ifndef SOFA_HELPER_FNDISPATCHER_INL
#define SOFA_HELPER_FNDISPATCHER_INL

#include <sofa/helper/FnDispatcher.h>
#include <sofa/helper/Factory.h> // for gettypename()
#include <iostream>
#include <string>

namespace sofa
{

namespace helper
{

// template <class BaseClass, typename ResulT>
// BasicDispatcher<BaseClass, ResulT>::~BasicDispatcher()
// {
// }

template <class BaseClass, typename ResulT>
ResulT BasicDispatcher<BaseClass, ResulT>::defaultFn(BaseClass& arg1, BaseClass& arg2)
{
    std::cerr << "ERROR DISPATCH ("
            << gettypename(typeid(arg1)) << ", "
            << gettypename(typeid(arg2)) << ")\n";
    this->ignore(typeid(arg1), typeid(arg2));
    return ResulT();
}

template <class BaseClass, typename ResulT>
ResulT BasicDispatcher<BaseClass, ResulT>::ignoreFn(BaseClass& /*arg1*/, BaseClass& /*arg2*/)
{
    return ResulT();
}

template <class BaseClass, typename ResulT>
ResulT BasicDispatcher<BaseClass, ResulT>::go(BaseClass &arg1,BaseClass &arg2)
{
    typename MapType::iterator itt = this->callBackMap.find(KeyType(TypeInfo(typeid(arg1)),TypeInfo(typeid(arg2))));
    if (itt==callBackMap.end())
        return defaultFn(arg1,arg2);
    else
        return (itt->second)(arg1,arg2);
}

/// Return true if a pair of argument correspond to a callback function (different than ignoreFn)
template <class BaseClass, typename ResulT>
bool BasicDispatcher<BaseClass, ResulT>::isSupported(BaseClass &arg1, BaseClass &arg2)
{
    typename MapType::iterator itt = this->callBackMap.find(KeyType(TypeInfo(typeid(arg1)),TypeInfo(typeid(arg2))));
    if (itt==callBackMap.end())
        return false;
    else
        return itt->second != ignoreFn;
}

template <class BaseClass, typename ResulT>
SingletonFnDispatcher<BaseClass, ResulT>::SingletonFnDispatcher()
{
}

template <class BaseClass, typename ResulT>
SingletonFnDispatcher<BaseClass, ResulT>* SingletonFnDispatcher<BaseClass, ResulT>::getInstance()
{
    static SingletonFnDispatcher<BaseClass, ResulT> instance;
    return &instance;
}

} // namespace helper

} // namespace sofa


#endif
