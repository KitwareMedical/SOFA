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
#define SOFA_COMPONENT_TOPOLOGY_EDGESETGEOMETRYALGORITHMS_CPP
#include <sofa/component/topology/EdgeSetGeometryAlgorithms.h>
#include <sofa/component/topology/EdgeSetGeometryAlgorithms.inl>

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/ObjectFactory.h>
namespace sofa
{

namespace component
{

namespace topology
{
using namespace sofa::defaulttype;
SOFA_DECL_CLASS(EdgeSetGeometryAlgorithms)
int EdgeSetGeometryAlgorithmsClass = core::RegisterObject("Edge set geometry algorithms")

#ifdef SOFA_FLOAT
        .add< EdgeSetGeometryAlgorithms<sofa::defaulttype::Vec3fTypes> >(true) // default template
#else
        .add< EdgeSetGeometryAlgorithms<sofa::defaulttype::Vec3dTypes> >(true) // default template
#ifndef SOFA_DOUBLE
        .add< EdgeSetGeometryAlgorithms<sofa::defaulttype::Vec3fTypes> >() // default template
#endif
#endif
#ifndef SOFA_FLOAT
        .add< EdgeSetGeometryAlgorithms<Vec2dTypes> >()
        .add< EdgeSetGeometryAlgorithms<Vec1dTypes> >()
        .add< EdgeSetGeometryAlgorithms<Rigid3dTypes> >()
        .add< EdgeSetGeometryAlgorithms<Rigid2dTypes> >()
#endif
#ifndef SOFA_DOUBLE
        .add< EdgeSetGeometryAlgorithms<Vec2fTypes> >()
        .add< EdgeSetGeometryAlgorithms<Vec1fTypes> >()
        .add< EdgeSetGeometryAlgorithms<Rigid3fTypes> >()
        .add< EdgeSetGeometryAlgorithms<Rigid2fTypes> >()
#endif
        ;

#ifndef SOFA_FLOAT
template class SOFA_BASE_TOPOLOGY_API EdgeSetGeometryAlgorithms<sofa::defaulttype::Vec3dTypes>;
template class SOFA_BASE_TOPOLOGY_API EdgeSetGeometryAlgorithms<Vec2dTypes>;
template class SOFA_BASE_TOPOLOGY_API EdgeSetGeometryAlgorithms<Vec1dTypes>;
template class SOFA_BASE_TOPOLOGY_API EdgeSetGeometryAlgorithms<Rigid3dTypes>;
template class SOFA_BASE_TOPOLOGY_API EdgeSetGeometryAlgorithms<Rigid2dTypes>;
#endif

#ifndef SOFA_DOUBLE
template class SOFA_BASE_TOPOLOGY_API EdgeSetGeometryAlgorithms<sofa::defaulttype::Vec3fTypes>;
template class SOFA_BASE_TOPOLOGY_API EdgeSetGeometryAlgorithms<Vec2fTypes>;
template class SOFA_BASE_TOPOLOGY_API EdgeSetGeometryAlgorithms<Vec1fTypes>;
template class SOFA_BASE_TOPOLOGY_API EdgeSetGeometryAlgorithms<Rigid3fTypes>;
template class SOFA_BASE_TOPOLOGY_API EdgeSetGeometryAlgorithms<Rigid2fTypes>;
#endif

} // namespace topology

} // namespace component

} // namespace sofa

