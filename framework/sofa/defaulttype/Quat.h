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
#ifndef SOFA_DEFAULTTYPE_QUAT_H
#define SOFA_DEFAULTTYPE_QUAT_H

#include <sofa/helper/Quater.h>
#include <sofa/defaulttype/DataTypeInfo.h>

namespace sofa
{

namespace defaulttype
{
typedef helper::Quater<double> Quatd; ///< alias
typedef helper::Quater<float>  Quatf; ///< alias
#ifdef SOFA_FLOAT
typedef Quatf Quat; ///< alias
#else
typedef Quatd Quat; ///< alias
#endif
typedef Quat Quaternion; ///< alias

// Specialization of the defaulttype::DataTypeInfo type traits template

template<class T>
struct DataTypeInfo< sofa::helper::Quater<T> > : public FixedArrayTypeInfo< sofa::helper::Quater<T> >
{
    static std::string name() { std::ostringstream o; o << "Quater<" << DataTypeName<T>::name() << ">"; return o.str(); }
};

// The next line hides all those methods from the doxygen documentation
/// \cond TEMPLATE_OVERRIDES

template<> struct DataTypeName<defaulttype::Quatf> { static const char* name() { return "Quatf"; } };
template<> struct DataTypeName<defaulttype::Quatd> { static const char* name() { return "Quatd"; } };

/// \endcond

} // namespace defaulttype

} // namespace sofa

#endif

