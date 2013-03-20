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
#include "initFlexible.h"

namespace sofa
{

namespace component
{

//Here are just several convenient functions to help user to know what contains the plugin

extern "C" {
    SOFA_Flexible_API void initExternalModule();
    SOFA_Flexible_API const char* getModuleName();
    SOFA_Flexible_API const char* getModuleVersion();
    SOFA_Flexible_API const char* getModuleLicense();
    SOFA_Flexible_API const char* getModuleDescription();
    SOFA_Flexible_API const char* getModuleComponentList();
}

void initExternalModule()
{
    static bool first = true;
    if (first)
    {
        first = false;
    }
}

const char* getModuleName()
{
    return "Flexible";
}

const char* getModuleVersion()
{
    return "0.2";
}

const char* getModuleLicense()
{
    return "LGPL";
}


const char* getModuleDescription()
{
    return "TODO: replace this with the description of your plugin";
}

const char* getModuleComponentList()
{
    /// string containing the names of the classes provided by the plugin
    return  "TopologyGaussPointSampler, ImageGaussPointSampler, ShepardShapeFunction, BarycentricShapeFunction, VoronoiShapeFunction, DefGradientMechanicalObject, LinearMapping, StrainMechanicalObject, GreenStrainMapping, CorotationalStrainMapping, DiagonalStrainMapping, PrincipalStretchesMapping, HookeForceField, AffineMechanicalObject, QuadraticMechanicalObject, DistanceMapping, ExtensionMapping";
}
}
}

/// Use the SOFA_LINK_CLASS macro for each class, to enable linking on all platforms

SOFA_LINK_CLASS(TopologyGaussPointSampler)

#ifdef SOFA_HAVE_IMAGE
SOFA_LINK_CLASS(VoronoiShapeFunction)
SOFA_LINK_CLASS(ImageGaussPointSampler)
SOFA_LINK_CLASS(ImageDensityMass)
#endif

SOFA_LINK_CLASS(ShepardShapeFunction)
SOFA_LINK_CLASS(BarycentricShapeFunction)
SOFA_LINK_CLASS(DefGradientMechanicalObject)
SOFA_LINK_CLASS(LinearMapping_point)
SOFA_LINK_CLASS(LinearMapping_affine)
SOFA_LINK_CLASS(LinearMapping_rigid)
SOFA_LINK_CLASS(LinearMapping_quadratic)
SOFA_LINK_CLASS(StrainMechanicalObject)
SOFA_LINK_CLASS(GreenStrainMapping)
SOFA_LINK_CLASS(CorotationalStrainMapping)
SOFA_LINK_CLASS(PrincipalStretchesMapping)
SOFA_LINK_CLASS(InvariantMapping)
SOFA_LINK_CLASS(HookeForceField)
SOFA_LINK_CLASS(MooneyRivlinForceField)
SOFA_LINK_CLASS(OgdenForceField)
SOFA_LINK_CLASS(VolumePreservationForceField)
SOFA_LINK_CLASS(AffineMechanicalObject)
SOFA_LINK_CLASS(QuadraticMechanicalObject)
SOFA_LINK_CLASS(DistanceMapping)
SOFA_LINK_CLASS(ExtensionMapping)
SOFA_LINK_CLASS(AffineUniformMass)


