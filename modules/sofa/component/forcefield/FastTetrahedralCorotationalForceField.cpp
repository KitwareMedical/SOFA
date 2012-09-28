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
#define SOFA_COMPONENT_INTERACTIONFORCEFIELD_FASTTETRAHEDRALCOROTATIONALFORCEFIELD_CPP
#include <sofa/component/forcefield/FastTetrahedralCorotationalForceField.inl>
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/Vec3Types.h>


namespace sofa
{

namespace component
{

namespace forcefield
{

using namespace sofa::defaulttype;


SOFA_DECL_CLASS(FastTetrahedralCorotationalForceField)

using namespace sofa::defaulttype;


// Register in the Factory
int FastTetrahedralCorotationalForceFieldClass = core::RegisterObject("Fast Corotational Tetrahedral Mesh")
#ifndef SOFA_FLOAT
        .add< FastTetrahedralCorotationalForceField<Vec3dTypes> >()
#endif
#ifndef SOFA_DOUBLE
        .add< FastTetrahedralCorotationalForceField<Vec3fTypes> >()
#endif
        ;

#ifndef SOFA_FLOAT
template class SOFA_MISC_FEM_API FastTetrahedralCorotationalForceField<Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
template class SOFA_MISC_FEM_API FastTetrahedralCorotationalForceField<Vec3fTypes>;
#endif

} // namespace forcefield

} // namespace component

} // namespace sofa

