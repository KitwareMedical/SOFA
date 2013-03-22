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
#define SOFA_COMPONENT_MAPPING_MLSMAPPING_point_CPP

#include "../initFlexible.h"
#include "../deformationMapping/MLSMapping.h"
#include <sofa/core/ObjectFactory.h>

#include <sofa/defaulttype/VecTypes.h>
#include "../types/DeformationGradientTypes.h"

namespace sofa
{
namespace component
{
namespace mapping
{

SOFA_DECL_CLASS(MLSMapping_point);

using namespace defaulttype;

// Register in the Factory
int MLSMappingClass_point = core::RegisterObject("Map child positions using moving least squares.")

        .add< MLSMapping< Vec3Types, Vec3Types > >(true)
//        .add< MLSMapping< Vec3Types, ExtVec3fTypes > >()
//        .add< MLSMapping< Vec3Types, F331Types > >()
//        .add< MLSMapping< Vec3Types, F332Types > >()
//        .add< MLSMapping< Vec3Types, F321Types > >()
//        .add< MLSMapping< Vec3Types, F311Types > >()
        ;

template class SOFA_Flexible_API MLSMapping< Vec3Types, Vec3Types >;
//template class SOFA_Flexible_API MLSMapping< Vec3Types, ExtVec3fTypes >;
//template class SOFA_Flexible_API MLSMapping< Vec3Types, F331Types >;
//template class SOFA_Flexible_API MLSMapping< Vec3Types, F332Types >;
//template class SOFA_Flexible_API MLSMapping< Vec3Types, F321Types >;
//template class SOFA_Flexible_API MLSMapping< Vec3Types, F311Types >;

} // namespace mapping
} // namespace component
} // namespace sofa

