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
#ifndef SOFA_HELPER_SYSTEM_GL_H
#define SOFA_HELPER_SYSTEM_GL_H

#ifndef SOFA_NO_OPENGL

#include <sofa/helper/system/config.h>

#if defined (SOFA_HAVE_GLEW) && !defined(PS3)
#include <GL/glew.h>
#elif defined(PS3)
#include <sofa/helper/gl/ps3gl_compat.h>
#elif defined (__APPLE__)
#include <OpenGL/gl.h>
#else
#define GL_GLEXT_PROTOTYPES // for glext.h : necessary to use glBindBuffer without glew and make GLSLShader file
#include <GL/gl.h>
#include <GL/glext.h> // necessary when you havn't glew
#endif

extern const char* GetGlExtensionsList();

extern bool CanUseGlExtension(char* ext);

#endif /* SOFA_NO_OPENGL */

#endif
