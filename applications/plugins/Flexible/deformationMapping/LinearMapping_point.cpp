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
#define SOFA_COMPONENT_MAPPING_LINEARMAPPING_point_CPP

#include "../initFlexible.h"
#include "../deformationMapping/LinearMapping.h"
#include <sofa/core/ObjectFactory.h>

#include <sofa/defaulttype/VecTypes.h>
#include "../types/DeformationGradientTypes.h"
#include "../types/AffineTypes.h"

namespace sofa
{
namespace component
{
namespace mapping
{

SOFA_DECL_CLASS(LinearMapping_point);

using namespace defaulttype;

// Register in the Factory
int LinearMappingClass_point = core::RegisterObject("Map child positions as a linear combination of parents.")

        .add< LinearMapping< Vec3Types, Vec3Types > >(true)
        .add< LinearMapping< Vec3Types, ExtVec3fTypes > >()
        .add< LinearMapping< Vec3Types, F331Types > >()
        .add< LinearMapping< Vec3Types, F332Types > >()
        .add< LinearMapping< Vec3Types, F321Types > >()
        .add< LinearMapping< Vec3Types, F311Types > >()
        .add< LinearMapping< Vec2Types, Vec2Types > >()
        .add< LinearMapping< Vec2Types, F221Types > >()
        .add< LinearMapping< Vec3Types, Affine3Types > >()
        ;

template class SOFA_Flexible_API LinearMapping< Vec3Types, Vec3Types >;
template class SOFA_Flexible_API LinearMapping< Vec3Types, ExtVec3fTypes >;
template class SOFA_Flexible_API LinearMapping< Vec3Types, F331Types >;
template class SOFA_Flexible_API LinearMapping< Vec3Types, F332Types >;
template class SOFA_Flexible_API LinearMapping< Vec3Types, F321Types >;
template class SOFA_Flexible_API LinearMapping< Vec3Types, F311Types >;
template class SOFA_Flexible_API LinearMapping< Vec2Types, Vec2Types >;
template class SOFA_Flexible_API LinearMapping< Vec2Types, F221Types >;
template class SOFA_Flexible_API LinearMapping< Vec3Types, Affine3Types >;

} // namespace mapping
} // namespace component
} // namespace sofa

