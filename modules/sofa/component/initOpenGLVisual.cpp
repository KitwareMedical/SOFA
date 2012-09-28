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
#include <sofa/helper/system/config.h>
#include <sofa/component/initOpenGLVisual.h>


namespace sofa
{

namespace component
{


void initOpenGLVisual()
{
    static bool first = true;
    if (first)
    {
        first = false;
    }
}

SOFA_LINK_CLASS(OglModel)
SOFA_LINK_CLASS(OglViewport)
SOFA_LINK_CLASS(Light)
SOFA_LINK_CLASS(LightManager)
SOFA_LINK_CLASS(PointSplatModel)
SOFA_LINK_CLASS(OglRenderingSRGB)
SOFA_LINK_CLASS(ClipPlane)
#ifdef SOFA_HAVE_GLEW
SOFA_LINK_CLASS(OglShader)
SOFA_LINK_CLASS(OglShaderVisualModel)
SOFA_LINK_CLASS(OglShadowShader)
SOFA_LINK_CLASS(OglTetrahedralModel)
SOFA_LINK_CLASS(OglTexture)
#endif


} // namespace component

} // namespace sofa
