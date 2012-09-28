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
#define SOFA_COMPONENT_ENGINE_MERGEPOINTS_CPP
#include <sofa/component/engine/MergePoints.inl>
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace engine
{

SOFA_DECL_CLASS(MergePoints)

int MergePointsClass = core::RegisterObject("Merge 2 cordinate vectors")
#ifdef SOFA_FLOAT
        .add< MergePoints<defaulttype::Vec3fTypes> >(true) // default template
#else
        .add< MergePoints<defaulttype::Vec3dTypes> >(true) // default template
#ifndef SOFA_DOUBLE
        .add< MergePoints<defaulttype::Vec3fTypes> >()
#endif
#endif
#ifndef SOFA_FLOAT
        .add< MergePoints<defaulttype::Vec1dTypes> >()
        .add< MergePoints<defaulttype::Vec2dTypes> >()
        .add< MergePoints<defaulttype::Rigid2dTypes> >()
        .add< MergePoints<defaulttype::Rigid3dTypes> >()
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
        .add< MergePoints<defaulttype::Vec1fTypes> >()
        .add< MergePoints<defaulttype::Vec2fTypes> >()
        .add< MergePoints<defaulttype::Rigid2fTypes> >()
        .add< MergePoints<defaulttype::Rigid3fTypes> >()
#endif //SOFA_DOUBLE
        ;

#ifndef SOFA_FLOAT
template class SOFA_ENGINE_API MergePoints<defaulttype::Vec1dTypes>;
template class SOFA_ENGINE_API MergePoints<defaulttype::Vec2dTypes>;
template class SOFA_ENGINE_API MergePoints<defaulttype::Vec3dTypes>;
template class SOFA_ENGINE_API MergePoints<defaulttype::Rigid2dTypes>;
template class SOFA_ENGINE_API MergePoints<defaulttype::Rigid3dTypes>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
template class SOFA_ENGINE_API MergePoints<defaulttype::Vec1fTypes>;
template class SOFA_ENGINE_API MergePoints<defaulttype::Vec2fTypes>;
template class SOFA_ENGINE_API MergePoints<defaulttype::Vec3fTypes>;
template class SOFA_ENGINE_API MergePoints<defaulttype::Rigid2fTypes>;
template class SOFA_ENGINE_API MergePoints<defaulttype::Rigid3fTypes>;
#endif //SOFA_DOUBLE


} // namespace constraint

} // namespace component

} // namespace sofa

