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
#define CLOSESTPOINTREGISTRATIONFORCEFIELD_CPP

#include "ClosestPointRegistrationForceField.inl"
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace forcefield
{

    using namespace sofa::defaulttype;

      SOFA_DECL_CLASS(ClosestPointRegistrationForceField)

      // Register in the Factory
      int ClosestPointRegistrationForceFieldClass = core::RegisterObject("Compute forces based on closest points from/to a target surface/point set")
    #ifndef SOFA_FLOAT
        .add< ClosestPointRegistrationForceField<Vec3dTypes> >()
    #endif
    #ifndef SOFA_DOUBLE
        .add< ClosestPointRegistrationForceField<Vec3fTypes> >()
    #endif
    ;

    #ifndef SOFA_FLOAT
      template class SOFA_REGISTRATION_API ClosestPointRegistrationForceField<Vec3dTypes>;
    #endif
    #ifndef SOFA_DOUBLE
      template class SOFA_REGISTRATION_API ClosestPointRegistrationForceField<Vec3fTypes>;
    #endif


            
}
}

} // namespace sofa


