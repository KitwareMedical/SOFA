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
#include <sofa/component/topology/CylinderGridTopology.h>
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



SOFA_DECL_CLASS(CylinderGridTopology)

int CylinderGridTopologyClass = core::RegisterObject("Cylinder grid in 3D")
        .addAlias("CylinderGrid")
        .add< CylinderGridTopology >()
        ;

CylinderGridTopology::CylinderGridTopology(int nx, int ny, int nz)
    : GridTopology(nx, ny, nz),
      center(initData(&center,Vector3(0.0f,0.0f,0.0f),"center", "Center of the cylinder")),
      axis(initData(&axis,Vector3(0.0f,0.0f,1.0f),"axis", "Main direction of the cylinder")),
      radius(initData(&radius,(SReal)1.0,"radius", "Radius of the cylinder")),
      length(initData(&length,(SReal)1.0,"length", "Length of the cylinder along its axis"))
{
}

CylinderGridTopology::CylinderGridTopology()
    : center(initData(&center,Vector3(0.0f,0.0f,0.0f),"center", "Center of the cylinder")),
      axis(initData(&axis,Vector3(0.0f,0.0f,1.0f),"axis", "Main direction of the cylinder")),
      radius(initData(&radius,(SReal)1.0,"radius", "Radius of the cylinder")),
      length(initData(&length,(SReal)1.0,"length", "Length of the cylinder along its axis"))
{
}

unsigned CylinderGridTopology::getIndex( int i, int j, int k ) const
{
    return n.getValue()[0]* ( n.getValue()[1]*k + j ) + i;
}

Vector3 CylinderGridTopology::getPoint(int i) const
{
    int x = i%n.getValue()[0]; i/=n.getValue()[0];
    int y = i%n.getValue()[1]; i/=n.getValue()[1];
    int z = i;
    return getPoint(x,y,z);
}

Vector3 CylinderGridTopology::getPoint(int x, int y, int z) const
{
    //return p0+dx*x+dy*y+dz*z;
    SReal r = radius.getValue();
    SReal l = length.getValue();
    Vector3 axisZ = axis.getValue();
    axisZ.normalize();
    Vector3 axisX = ((axisZ-Vector3(1,0,0)).norm() < 0.000001 ? Vector3(0,1,0) : Vector3(1,0,0));
    Vector3 axisY = cross(axisZ,axisX);
    axisX = cross(axisY,axisZ);
    axisX.normalize();
    axisY.normalize();
    axisZ.normalize();
    int nx = getNx();
    int ny = getNy();
    int nz = getNz();
    // coordonate on a square
    Vector3 p(x*2*r/(nx-1) - r, y*2*r/(ny-1) - r, 0);
    // scale it to be on a circle
    if (p.norm() > 0.0000001)
        p *= helper::rmax(helper::rabs(p[0]),helper::rabs(p[1]))/p.norm();
    if (nz>1)
        p[2] = z*l/(nz-1);
    return center.getValue()+axisX*p[0] + axisY*p[1] + axisZ * p[2];
}

} // namespace topology

} // namespace component

} // namespace sofa

