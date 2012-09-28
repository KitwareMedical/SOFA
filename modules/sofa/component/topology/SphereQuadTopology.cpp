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
#include <sofa/component/topology/SphereQuadTopology.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/rmath.h>

namespace sofa
{

namespace component
{

namespace topology
{

using namespace sofa::defaulttype;



SOFA_DECL_CLASS(SphereQuadTopology)

int SphereQuadTopologyClass = core::RegisterObject("Sphere topology constructed with deformed quads")
        .addAlias("SphereQuad")
        .add< SphereQuadTopology >()
        ;

SphereQuadTopology::SphereQuadTopology(int nx, int ny, int nz)
    : CubeTopology(nx, ny, nz),
      center(initData(&center,Vector3(0.0f,0.0f,0.0f),"center", "Center of the sphere")),
      radius(initData(&radius,(SReal)1.0,"radius", "Radius of the sphere"))
{
}

SphereQuadTopology::SphereQuadTopology()
    : center(initData(&center,Vector3(0.0f,0.0f,0.0f),"center", "Center of the sphere")),
      radius(initData(&radius,(SReal)1.0,"radius", "Radius of the sphere"))
{
}

Vector3 SphereQuadTopology::getPoint(int x, int y, int z) const
{
    Vector3 p((2*x)/(SReal)(nx.getValue()-1) - 1, (2*y)/(SReal)(ny.getValue()-1) - 1, (2*z)/(SReal)(nz.getValue()-1) - 1);
    p.normalize();
    return center.getValue()+p*radius.getValue();
}

} // namespace topology

} // namespace component

} // namespace sofa

