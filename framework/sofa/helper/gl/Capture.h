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
#ifndef SOFA_HELPER_GL_CAPTURE_H
#define SOFA_HELPER_GL_CAPTURE_H

#ifndef SOFA_NO_OPENGL

#include <sofa/helper/system/gl.h>

#include <sofa/helper/io/Image.h>

#include <sofa/helper/helper.h>

namespace sofa
{

namespace helper
{

namespace gl
{

//using namespace sofa::defaulttype;

class SOFA_HELPER_API Capture
{
protected:
    std::string prefix;
    int counter;

public:

    Capture();

    const std::string& getPrefix() const { return prefix; }
    int getCounter() const { return counter; }

    void setPrefix(const std::string v) { prefix=v; }
    void setCounter(int v=-1) { counter = v; }

    std::string findFilename();
    bool saveScreen(const std::string& filename, int compression_level = -1);

    bool saveScreen(int compression_level = -1);
};

} // namespace gl

} // namespace helper

} // namespace sofa

#endif /* SOFA_NO_OPENGL */

#endif
