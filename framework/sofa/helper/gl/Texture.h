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
#ifndef SOFA_HELPER_GL_TEXTURE_H
#define SOFA_HELPER_GL_TEXTURE_H

#ifndef SOFA_NO_OPENGL

#include <sofa/helper/system/gl.h>
#include <sofa/helper/helper.h>
#include <sofa/helper/io/Image.h>

namespace sofa
{

namespace helper
{

namespace gl
{

//using namespace sofa::defaulttype;

class SOFA_HELPER_API Texture
{
private:
    io::Image *image;
    GLuint id, target;
    bool repeat, linearInterpolation, generateMipmaps, srgbColorspace;
    float minLod, maxLod;
public:

    Texture()
        :image(NULL),id(0),repeat(true), linearInterpolation(true), generateMipmaps(false),
         srgbColorspace(false), minLod(-1000), maxLod(1000)
    {
    }

    Texture (io::Image *img, bool repeat = true, bool linearInterpolation = true, bool generateMipmaps = false,
            bool srgbColorspace = false, float minLod = -1000, float maxLod = 1000)
        :image(img),id(0),repeat(repeat), linearInterpolation(linearInterpolation), generateMipmaps(generateMipmaps),
         srgbColorspace(srgbColorspace), minLod(minLod), maxLod(maxLod)
    {}

    io::Image* getImage(void);
    GLuint getTarget() const { return target; }
    void   bind(void);
    void   unbind(void);
    void   init ();
    void   update ();	// to use for dynamic change of the texture image (no memory overhead due to multiple texture creation)
    ~Texture();

    GLuint getId() const { return id; }

private:
    Texture(const Texture& ) {};
    Texture operator=(const Texture& ) { return Texture(); }
};

} // namespace gl

} // namespace helper

} // namespace sofa

#endif /* SOFA_NO_OPENGL */

#endif
