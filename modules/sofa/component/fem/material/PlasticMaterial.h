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
#ifndef SOFA_COMPONENT_FEM_MATERIAL_PLASTICMATERIAL_H
#define SOFA_COMPONENT_FEM_MATERIAL_PLASTICMATERIAL_H

#if !defined(__GNUC__) || (__GNUC__ > 3 || (_GNUC__ == 3 && __GNUC_MINOR__ > 3))
#pragma once
#endif

//#define SOFA_FLOAT
#ifndef SOFA_DOUBLE
#define SOFA_DOUBLE
#endif

#include <sofa/component/fem/material/BaseMaterial.h>
#include <newmat/newmat.h>
#include <newmat/newmatap.h>
#include <sofa/defaulttype/Vec.h>

namespace sofa
{

namespace component
{

namespace fem
{

namespace material
{

/**
 * Plastic material (proof of principle)
 */
class PlasticMaterial : public component::fem::BaseMaterial
{

public:
    SOFA_CLASS(PlasticMaterial, component::fem::BaseMaterial);

    typedef sofa::defaulttype::Vector3 Vector3;
    typedef sofa::helper::vector<double> VecDouble;
    typedef sofa::helper::vector<Vector3> VecStress;

	// Material properties
    Data<double> _poissonRatio;
    Data<double> _youngModulus;

    // Stress-strain curve description
    VecDouble _E;
    VecDouble _epsilon;
    VecStress _sigma;

    // Strain of the previous iteration
    VecDouble _previousVonMisesStrain;

    PlasticMaterial();
    void computeStress (Vector3& stress, Vector3 strain, int elementIndex);
    void computeDStress (Vector3& dstress, Vector3 dstrain);

    double computeVonMisesStrain(Vector3 &strain);
    void computeStressOnSection(Vector3& Stress, Vector3 Strain, int section);	// computes the stress on a given section of the piecewise function

    virtual void computeStress (unsigned int /*iElement*/){}

};


} // namespace material

} // namespace fem

} // namespace component

} // namespace sofa
#endif
