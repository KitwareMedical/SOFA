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
#define SOFA_COMPONENT_ENGINE_MERGEVECTORS_CPP
#include <sofa/component/engine/MergeVectors.inl>
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace engine
{

SOFA_DECL_CLASS(MergeVectors)

int MergeVectorsClass = core::RegisterObject("Apply a merge operation to combine several inputs")
#if defined(SOFA_DOUBLE)
    .add< MergeVectors< helper::vector<double> > >(true)
#elif defined(SOFA_FLOAT)
    .add< MergeVectors< helper::vector<float> > >(true)
#else
    .add< MergeVectors< helper::vector<double> > >(true)
    .add< MergeVectors< helper::vector<float> > >()
#endif
    .add< MergeVectors< helper::vector<int> > >()
    .add< MergeVectors< helper::vector<bool> > >()
    .add< MergeVectors< helper::vector<std::string> > >()
#ifndef SOFA_FLOAT
    .add< MergeVectors< helper::vector<defaulttype::Vec2d> > >()
    .add< MergeVectors< helper::vector<defaulttype::Vec3d> > >()
    .add< MergeVectors< defaulttype::Rigid2dTypes::VecCoord > >()
    .add< MergeVectors< defaulttype::Rigid2dTypes::VecDeriv > >()
    .add< MergeVectors< defaulttype::Rigid3dTypes::VecCoord > >()
    .add< MergeVectors< defaulttype::Rigid3dTypes::VecDeriv > >()
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
    .add< MergeVectors< helper::vector<defaulttype::Vec2f> > >()
    .add< MergeVectors< helper::vector<defaulttype::Vec3f> > >()
    .add< MergeVectors< defaulttype::Rigid2fTypes::VecCoord > >()
    .add< MergeVectors< defaulttype::Rigid2fTypes::VecDeriv > >()
    .add< MergeVectors< defaulttype::Rigid3fTypes::VecCoord > >()
    .add< MergeVectors< defaulttype::Rigid3fTypes::VecDeriv > >()
#endif //SOFA_DOUBLE
        ;

template class SOFA_ENGINE_API MergeVectors< helper::vector<int> >;
template class SOFA_ENGINE_API MergeVectors< helper::vector<bool> >;
template class SOFA_ENGINE_API MergeVectors< helper::vector<std::string> >;

#ifndef SOFA_FLOAT
template class SOFA_ENGINE_API MergeVectors< helper::vector<double> >;
template class SOFA_ENGINE_API MergeVectors< helper::vector<defaulttype::Vec2d> >;
template class SOFA_ENGINE_API MergeVectors< helper::vector<defaulttype::Vec3d> >;
template class SOFA_ENGINE_API MergeVectors< defaulttype::Rigid2dTypes::VecCoord >;
template class SOFA_ENGINE_API MergeVectors< defaulttype::Rigid2dTypes::VecDeriv >;
template class SOFA_ENGINE_API MergeVectors< defaulttype::Rigid3dTypes::VecCoord >;
template class SOFA_ENGINE_API MergeVectors< defaulttype::Rigid3dTypes::VecDeriv >;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
template class SOFA_ENGINE_API MergeVectors< helper::vector<float> >;
template class SOFA_ENGINE_API MergeVectors< helper::vector<defaulttype::Vec2f> >;
template class SOFA_ENGINE_API MergeVectors< helper::vector<defaulttype::Vec3f> >;
template class SOFA_ENGINE_API MergeVectors< defaulttype::Rigid2fTypes::VecCoord >;
template class SOFA_ENGINE_API MergeVectors< defaulttype::Rigid2fTypes::VecDeriv >;
template class SOFA_ENGINE_API MergeVectors< defaulttype::Rigid3fTypes::VecCoord >;
template class SOFA_ENGINE_API MergeVectors< defaulttype::Rigid3fTypes::VecDeriv >;
#endif //SOFA_DOUBLE


} // namespace constraint

} // namespace component

} // namespace sofa

