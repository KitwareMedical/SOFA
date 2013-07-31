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
*                               SOFA :: Plugins                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#define SOFA_COMPONENT_MAPPING_RelativeMAPPING_CPP


#include "RelativeStrainMapping.h"
#include <sofa/core/ObjectFactory.h>

#include "../types/DeformationGradientTypes.h"
#include "../types/StrainTypes.h"

namespace sofa
{
namespace component
{
namespace mapping
{

SOFA_DECL_CLASS(RelativeStrainMapping);

using namespace defaulttype;

// Register in the Factory
int RelativeStrainMappingClass = core::RegisterObject("Map a total strain to an elastic strain + offset")

        .add< RelativeStrainMapping< E331Types > >(true)
        .add< RelativeStrainMapping< E321Types > >()
        .add< RelativeStrainMapping< E332Types > >()
        .add< RelativeStrainMapping< E333Types > >()

//        .add< RelativeStrainMapping< D331Types > >()
//        .add< RelativeStrainMapping< D321Types > >()
//        .add< RelativeStrainMapping< D332Types > >()
        ;

template class SOFA_Flexible_API RelativeStrainMapping< E331Types >;
template class SOFA_Flexible_API RelativeStrainMapping< E321Types >;
template class SOFA_Flexible_API RelativeStrainMapping< E332Types >;
template class SOFA_Flexible_API RelativeStrainMapping< E333Types >;

//template class SOFA_Flexible_API RelativeStrainMapping< D331Types >;
//template class SOFA_Flexible_API RelativeStrainMapping< D321Types >;
//template class SOFA_Flexible_API RelativeStrainMapping< D332Types >;

} // namespace mapping
} // namespace component
} // namespace sofa

