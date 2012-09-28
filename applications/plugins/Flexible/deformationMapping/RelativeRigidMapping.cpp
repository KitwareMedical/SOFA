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
#define SOFA_COMPONENT_MAPPING_RelativeRigidMapping_CPP

#include "../initFlexible.h"
#include "../deformationMapping/RelativeRigidMapping.inl"
#include <sofa/core/ObjectFactory.h>

#include <sofa/helper/vector.h>
#include <sofa/defaulttype/DataTypeInfo.h>

namespace sofa
{

namespace component
{

namespace mapping
{

SOFA_DECL_CLASS(RelativeRigidMapping)

using namespace defaulttype;

// Register in the Factory
int RelativeRigidMappingClass = core::RegisterObject("Computes relative rigid configurations")
#ifndef SOFA_FLOAT
        .add< RelativeRigidMapping< Rigid3dTypes, Rigid3dTypes > >()
#endif
#ifndef SOFA_DOUBLE
        .add< RelativeRigidMapping< Rigid3fTypes, Rigid3fTypes > >()
#endif
        ;

#ifndef SOFA_FLOAT
template class SOFA_Flexible_API RelativeRigidMapping< Rigid3dTypes, Rigid3dTypes >;
#endif

#ifndef SOFA_DOUBLE
template class SOFA_Flexible_API RelativeRigidMapping< Rigid3fTypes, Rigid3fTypes >;
#endif




} // namespace mapping

} // namespace component

} // namespace sofa

