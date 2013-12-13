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
#define SOFA_COMPONENT_ENGINE_PAIRBOXROI_CPP
#include <sofa/component/engine/PairBoxROI.inl>
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/Vec3Types.h>

namespace sofa
{

namespace component
{

namespace engine
{

SOFA_DECL_CLASS(PairBoxROI)

int PairBoxROIClass = core::RegisterObject("Find the primitives (vertex/edge/triangle/tetrahedron) inside a given box")
#ifndef SOFA_FLOAT
        .add< PairBoxROI<Vec3dTypes> >()
        .add< PairBoxROI<Rigid3dTypes> >()
        .add< PairBoxROI<Vec6dTypes> >() //Phuoc
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
        .add< PairBoxROI<Vec3fTypes> >()
        .add< PairBoxROI<Rigid3fTypes> >()
        .add< PairBoxROI<Vec6fTypes> >() //Phuoc
#endif //SOFA_DOUBLE
        ;

#ifndef SOFA_FLOAT
template class SOFA_ENGINE_API PairBoxROI<Vec3dTypes>;
template class SOFA_ENGINE_API PairBoxROI<Rigid3dTypes>;
template class SOFA_ENGINE_API PairBoxROI<Vec6dTypes>; //Phuoc
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
template class SOFA_ENGINE_API PairBoxROI<Vec3fTypes>;
template class SOFA_ENGINE_API PairBoxROI<Rigid3fTypes>;
template class SOFA_ENGINE_API PairBoxROI<Vec6fTypes>; //Phuoc
#endif //SOFA_DOUBLE


} // namespace constraint

} // namespace component

} // namespace sofa

