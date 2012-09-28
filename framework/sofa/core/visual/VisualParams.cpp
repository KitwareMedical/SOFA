/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *
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
*                              SOFA :: Framework                              *
*                                                                             *
* Authors: The SOFA Team (see Authors.txt)                                    *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/system/thread/thread_specific_ptr.h>
#include <sofa/helper/BackTrace.h>
#include <cassert>
#include <iostream>

namespace sofa
{

namespace core
{

namespace visual
{

VisualParams::VisualParams()
    : m_viewport(sofa::helper::make_array(0,0,0,0))
    , m_zNear(0)
    , m_zFar(0)
    , m_cameraType(PERSPECTIVE_TYPE)
    , m_pass(Std)
    , m_drawTool(NULL)
    , m_boundFrameBuffer(NULL)
    , m_x (ConstVecCoordId::position())
    , m_v (ConstVecDerivId::velocity())
    , m_supportedAPIs(0)
{
    m_displayFlags.setShowVisualModels(true); // BUGFIX: visual models are visible by default
}

/// Get the default VisualParams, to be used to provide a default values for method parameters
VisualParams* VisualParams::defaultInstance()
{
    SOFA_THREAD_SPECIFIC_PTR(VisualParams, threadParams);
    VisualParams* ptr = threadParams;
    if (!ptr)
    {
        ptr = new VisualParams;
        threadParams = ptr;
    }
    return ptr;
}

} // namespace visual

} // namespace core

} // namespace sofa
