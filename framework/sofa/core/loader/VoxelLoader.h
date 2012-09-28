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
#ifndef SOFA_CORE_VOXELLOADER_H
#define SOFA_CORE_VOXELLOADER_H

#include <sofa/core/loader/BaseLoader.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/core/objectmodel/DataFileName.h>
#include <sofa/helper/fixed_array.h>

namespace sofa
{
namespace helper { namespace io { class Image; }}

namespace core
{

namespace loader
{


class SOFA_CORE_API VoxelLoader : public sofa::core::loader::BaseLoader
{
public:
    SOFA_ABSTRACT_CLASS(VoxelLoader,BaseLoader);

    typedef defaulttype::Vec<3, int> Vec3i;
    typedef defaulttype::Vec<6, int> Vec6i;
    typedef helper::fixed_array<unsigned int,8> Hexahedron;
protected:
    VoxelLoader();
    virtual ~VoxelLoader();
public:

    Data< helper::vector<sofa::defaulttype::Vec<3,SReal> > > positions;
    Data< helper::vector<Hexahedron > > hexahedra;


    void addHexahedron(helper::vector< Hexahedron >* pHexahedra, const helper::fixed_array<unsigned int,8> &p);
    void addHexahedron(helper::vector< Hexahedron >* pHexahedra,
            unsigned int p0, unsigned int p1, unsigned int p2, unsigned int p3,
            unsigned int p4, unsigned int p5, unsigned int p6, unsigned int p7);

    virtual defaulttype::Vector3 getVoxelSize () const = 0;

    virtual helper::vector<unsigned int> getHexaIndicesInGrid() const=0;

    virtual int getDataSize() const = 0;

    virtual Vec6i getROI() const = 0;

    virtual unsigned char * getData() = 0;

    virtual unsigned char * getSegmentID() = 0;

    // fill the texture by 'image' only where there is the 'segmentation' of 'activeValue' and give the 3D texture sizes
    virtual void createSegmentation3DTexture( unsigned char **textureData, int& width, int& height, int& depth) = 0;
};

}

} // namespace component

} // namespace sofa

#endif
