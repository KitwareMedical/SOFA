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
#include <sofa/helper/io/bvh/BVHMotion.h>
#include <iostream>

namespace sofa
{

namespace helper
{

namespace io
{

namespace bvh
{

void BVHMotion::init(double _fTime, unsigned int _fCount, unsigned int _fSize)
{
    frameTime = _fTime;
    frameCount = _fCount;

    frames.resize(frameCount);

    for (int i=0; i<frameCount; i++)
        frames[i].resize(_fSize);
}

void BVHMotion::debug(void)
{
    for (unsigned int i=0; i<frames.size(); i++)
    {
        for (unsigned int j=0; j<frames[i].size(); j++)
            std::cout << frames[i][j] << " ";
        std::cout << "\n";
    }
}

} // namespace bvh

} // namespace io

} // namespace helper

} // namespace sofa
