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
#ifndef SOFA_HELPER_GL_CYLINDER_H
#define SOFA_HELPER_GL_CYLINDER_H

#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/Quat.h>

#include <sofa/helper/system/gl.h>
#include <sofa/helper/system/glu.h>

#include <sofa/helper/helper.h>

#include <map>

namespace sofa
{

namespace helper
{

namespace gl
{

using namespace sofa::defaulttype;

class SOFA_HELPER_API Cylinder
{
public:

    Cylinder(SReal len=(SReal)1);
    Cylinder(const Vector3& len);
    Cylinder(const Vector3& center, const Quaternion &orient, const Vector3& length);
    Cylinder(const Vector3& center, const double orient[4][4], const Vector3& length);
    Cylinder(const double *mat, const Vector3& length);
    Cylinder(const Vector3& center, const Quaternion &orient, SReal length=(SReal)1);
    Cylinder(const Vector3& center, const double orient[4][4], SReal length=(SReal)1);
    Cylinder(const double *mat, SReal length=(SReal)1.0);

    ~Cylinder();

    void update(const Vector3& center, const Quaternion& orient = Quaternion());
    void update(const Vector3& center, const double orient[4][4]);
    void update(const double *mat);

    void draw();

    static void draw(const Vector3& center, const Quaternion& orient, const Vector3& length);
    static void draw(const Vector3& center, const double orient[4][4], const Vector3& length);
    static void draw(const double *mat, const Vector3& length);
    static void draw(const Vector3& center, const Quaternion& orient, SReal length=(SReal)1);
    static void draw(const Vector3& center, const double orient[4][4], SReal length=(SReal)1);
    static void draw(const double *mat, SReal length=(SReal)1.0);

private:

    Vector3 length;
    double matTransOpenGL[16];

    GLUquadricObj *quadratic;
    GLuint displayList;

    void initDraw();

    static std::map < std::pair<std::pair<float,float>,float>, Cylinder* > CylinderMap;
    static Cylinder* get(const Vector3& len);

};

} // namespace gl

} // namespace helper

} // namespace sofa

#endif
