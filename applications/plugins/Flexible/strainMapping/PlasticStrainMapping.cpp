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
#define SOFA_COMPONENT_MAPPING_CorotationalStrainMAPPING_CPP


#include "PlasticStrainMapping.h"
#include <sofa/core/ObjectFactory.h>

#include "../types/DeformationGradientTypes.h"
#include "../types/StrainTypes.h"

namespace sofa
{
namespace component
{
namespace mapping
{

SOFA_DECL_CLASS(PlasticStrainMapping);

using namespace defaulttype;

// Register in the Factory
int PlasticStrainMappingClass = core::RegisterObject("Map a total strain to an elastic strain + a plastic strain.")

        .add< PlasticStrainMapping< E331Types > >(true)
        .add< PlasticStrainMapping< E321Types > >()
        .add< PlasticStrainMapping< E332Types > >()
        .add< PlasticStrainMapping< E333Types > >()

        .add< PlasticStrainMapping< D331Types > >()
        .add< PlasticStrainMapping< D321Types > >()
        .add< PlasticStrainMapping< D332Types > >()
        ;

template class SOFA_Flexible_API PlasticStrainMapping< E331Types >;
template class SOFA_Flexible_API PlasticStrainMapping< E321Types >;
template class SOFA_Flexible_API PlasticStrainMapping< E332Types >;
template class SOFA_Flexible_API PlasticStrainMapping< E333Types >;

template class SOFA_Flexible_API PlasticStrainMapping< D331Types >;
template class SOFA_Flexible_API PlasticStrainMapping< D321Types >;
template class SOFA_Flexible_API PlasticStrainMapping< D332Types >;

} // namespace mapping
} // namespace component
} // namespace sofa

