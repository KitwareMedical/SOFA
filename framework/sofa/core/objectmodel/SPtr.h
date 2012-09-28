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
#ifndef SOFA_CORE_OBJECTMODEL_SPTR_H
#define SOFA_CORE_OBJECTMODEL_SPTR_H

#include <sofa/helper/system/config.h>
#include <sofa/helper/vector.h>
#include <sofa/core/core.h>

namespace sofa
{

namespace core
{

namespace objectmodel
{

/**
 *  \brief new operator for classes with smart pointers (such as all components deriving from Base)
 *
 *  This class should be used as :
 *     MyT::SPtr p = sofa::core::objectmodel::New<MyT>(myargs);
 *  instead of :
 *     MyT* p = new MyT(myargs);
 *
 *  The use of this New operator and SPtr pointers insures that all created objects are :
 *    - destroyed (no leak),
 *    - only once (no double desctructions),
 *    - and only after the last reference to them are erased (no invalid pointers).
 *
 */
template<class T>
class New : public T::SPtr
{
    typedef typename T::SPtr SPtr;
public:
    New() : SPtr(new T) {}
    template <class A1>
    New(A1 a1) : SPtr(new T(a1)) {}
    template <class A1, class A2>
    New(A1 a1, A2 a2) : SPtr(new T(a1,a2)) {}
    template <class A1, class A2, class A3>
    New(A1 a1, A2 a2, A3 a3) : SPtr(new T(a1,a2,a3)) {}
    template <class A1, class A2, class A3, class A4>
    New(A1 a1, A2 a2, A3 a3, A4 a4) : SPtr(new T(a1,a2,a3,a4)) {}
};

/// dynamic_cast operator for SPtr
template<class T>
class SPtr_dynamic_cast : public T::SPtr
{
public:
    template<class UPtr>
    SPtr_dynamic_cast(UPtr p) : T::SPtr(dynamic_cast<T*>(p.get())) {}
};

/// static_cast operator for SPtr
template<class T>
class SPtr_static_cast : public T::SPtr
{
public:
    template<class UPtr>
    SPtr_static_cast(UPtr p) : T::SPtr(static_cast<T*>(p.get())) {}
};

/// const_cast operator for SPtr
template<class T>
class SPtr_const_cast : public T::SPtr
{
public:
    template<class UPtr>
    SPtr_const_cast(UPtr p) : T::SPtr(const_cast<T*>(p.get())) {}
};

} // namespace objectmodel

} // namespace core

} // namespace sofa



#endif

