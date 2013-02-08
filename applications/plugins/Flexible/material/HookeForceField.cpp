/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 4      *
*                (c) 2006-2009 MGH, INRIA, USTL, UJF, CNRS                    *
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
#define SOFA_HookeFORCEFIELD_CPP

#include "../initFlexible.h"
#include "HookeForceField.h"
#include "../types/StrainTypes.h"
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/State.inl>

namespace sofa
{
namespace component
{
namespace forcefield
{

SOFA_DECL_CLASS(HookeForceField);

using namespace defaulttype;

// Register in the Factory
int HookeForceFieldClass = core::RegisterObject("Hooke's Law for isotropic homogeneous materials")

        .add< HookeForceField< E331Types > >(true)
        .add< HookeForceField< E321Types > >()
        .add< HookeForceField< E311Types > >()
        .add< HookeForceField< E332Types > >()
        .add< HookeForceField< E333Types > >()

        .add< HookeForceField< D331Types > >()
        .add< HookeForceField< D321Types > >()
        .add< HookeForceField< D332Types > >()

        .add< HookeForceField< U331Types > >()
        .add< HookeForceField< U321Types > >()
        ;

template class SOFA_Flexible_API HookeForceField< E331Types >;
template class SOFA_Flexible_API HookeForceField< E321Types >;
template class SOFA_Flexible_API HookeForceField< E311Types >;
template class SOFA_Flexible_API HookeForceField< E332Types >;
template class SOFA_Flexible_API HookeForceField< E333Types >;

template class SOFA_Flexible_API HookeForceField< D331Types >;
template class SOFA_Flexible_API HookeForceField< D321Types >;
template class SOFA_Flexible_API HookeForceField< D332Types >;

template class SOFA_Flexible_API HookeForceField< U331Types >;
template class SOFA_Flexible_API HookeForceField< U321Types >;

}
}
}

