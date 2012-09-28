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
#include <sofa/helper/system/config.h>
#include <sofa/helper/io/Image.h>
#include <sofa/helper/Factory.inl>
#include <stdio.h>
#include <string.h>

namespace sofa
{

namespace helper
{

template class SOFA_HELPER_API Factory<std::string, sofa::helper::io::Image, std::string>;

namespace io
{

SOFA_LINK_CLASS(ImageBMP)
SOFA_LINK_CLASS(ImagePNG)

const char *Image::strFromDataType[COUNT_OF_DATA_TYPES+1] =
{
    "UNORM8",
    "UNORM16",
    "UINT32",
    "HALF",
    "FLOAT",

    "UCOMPRESSED",

    "COUNT_OF_DATA_TYPES"
};

const char *Image::strFromChannelFormat[COUNT_OF_CHANNEL_FORMATS+1] =
{
    "L",
    "LA",
    "R",
    "RG",
    "RGB",
    "RGBA",
    "BGR",
    "BGRA",

    "COUNT_OF_CHANNEL_FORMATS"
};

const char *Image::strFromTextureType[TEXTURE_INVALID+1] =
{
    "TEXTURE_2D",
    "TEXTURE_3D",
    "TEXTURE_CUBE",

    "TEXTURE_INVALID"
};

static unsigned tableBytes[Image::COUNT_OF_DATA_TYPES][Image::COUNT_OF_CHANNEL_FORMATS] =
{
    // Bytes per pixel
    // UNORM8
    {
        1,  // L
        2,  // LA
        1,  // R
        2,  // RG
        3,  // RGB
        4,  // RGBA
        3,  // BGR
        4   // BGRA
    },
    // UNORM16
    {
        2,  // L
        4,  // LA
        2,  // R
        4,  // RG
        6,  // RGB
        8,  // RGBA
        6,  // BGR
        8   // BGRA
    },
    // UINT32
    {
        4,  // L
        8,  // LA
        4,  // R
        8,  // RG
        12, // RGB
        16, // RGBA
        12, // BGR
        16  // BGRA
    },
    // HALF
    {
        2,  // L
        4,  // LA
        2,  // R
        4,  // RG
        6,  // RGB
        8,  // RGBA
        6,  // BGR
        8   // BGRA
    },
    // FLOAT
    {
        4,  // L
        8,  // LA
        4,  // R
        8,  // RG
        12, // RGB
        16, // RGBA
        12, // BGR
        16  // BGRA
    },
    // Bytes per block
    // UCOMPRESSED
    {
        8,  // L
        16, // LA
        8,  // R
        16, // RG
        8,  // RGB
        16, // RGBA
        0,  // BGR
        0   // BGRA
    }
};

Image::Image()
    : data(NULL)
{
}

Image::~Image()
{
    clear();
}

Image::Image(const Image& rhs)
    :data(NULL)
{
    init(rhs.width,rhs.height,rhs.depth,rhs.mipmaps,rhs.dataType,rhs.channelFormat);
    memcpy(data,rhs.data,getImageSize());
}

Image& Image::operator=(const Image& rhs)
{
    clear();
    init(rhs.width,rhs.height,rhs.depth,rhs.mipmaps,rhs.dataType,rhs.channelFormat);
    memcpy(data,rhs.data,getImageSize());
    return *this;
}

unsigned Image::getWidth(unsigned mipmap) const
{
    unsigned result = width >> mipmap;
    return result ? result : 1;
}

unsigned Image::getHeight(unsigned mipmap) const
{
    unsigned result = height >> mipmap;
    return result ? result : 1;
}

unsigned Image::getDepth(unsigned mipmap) const
{
    unsigned result = depth >> mipmap;
    return result ? result : 1;
}

unsigned Image::getBytesPerPixel() const
{
    return dataType != UCOMPRESSED ? tableBytes[dataType][channelFormat] : 0;
}

unsigned Image::getBytesPerBlock() const
{
    return dataType == UCOMPRESSED ? tableBytes[UCOMPRESSED][channelFormat] : 0;
}

unsigned Image::getBytesPerChannel() const
{
    return getBytesPerPixel() / getChannelCount();
}

unsigned Image::getChannelCount() const
{
    static unsigned table[COUNT_OF_CHANNEL_FORMATS] =
    {
        1,  // L
        2,  // LA
        1,  // R
        2,  // RG
        3,  // RGB
        4,  // RGBA
        3,  // BGR
        4   // BGRA
    };

    return table[channelFormat];
}

unsigned Image::getMipmapCount() const
{
    return mipmaps;
}

unsigned Image::getPixelCount() const
{
    return getImageSize() / getBytesPerPixel();
}

unsigned Image::getLineSize(unsigned mipmap) const
{
    return getWidth(mipmap) * getBytesPerPixel();
}

unsigned Image::getMipmapSize(unsigned mipmap) const
{
    // Return the size of one mipmap in bytes. For cubemaps, the size of all six faces is returned.
    unsigned width = getWidth(mipmap);
    unsigned height = getHeight(mipmap);
    unsigned depth = (getTextureType() == TEXTURE_CUBE)? 6 : getDepth(mipmap);

    if (dataType == UCOMPRESSED)
        return ((width + 3) >> 2) * ((height + 3) >> 2) * depth * getBytesPerBlock();

    return width * height * depth * getBytesPerPixel();
}

unsigned Image::getMipmapRangeSize(unsigned firstMipmap, unsigned mipmaps) const
{
    // Return the size of mipmap range (multiple consecutive mipmaps) in bytes.
    unsigned lastMipmap = firstMipmap + mipmaps;
    if (lastMipmap > mipmaps) lastMipmap = mipmaps;
    unsigned size = 0;
    for (unsigned i = firstMipmap; i < lastMipmap; i++)
        size += getMipmapSize(i);
    return size;
}

unsigned Image::getImageSize() const
{
    return getMipmapRangeSize(0, mipmaps);
}

Image::DataType Image::getDataType() const
{
    return dataType;
}

Image::ChannelFormat Image::getChannelFormat() const
{
    return channelFormat;
}

Image::TextureType Image::getTextureType() const
{
    if (dataType == UCOMPRESSED && channelFormat >= BGR)
        return TEXTURE_INVALID;

    if (depth == 0 && width == height)
        return TEXTURE_CUBE;
    else if(depth > 1)
        return TEXTURE_3D;
    else if (depth == 1)
        return TEXTURE_2D;
    else
        return TEXTURE_INVALID;
}

unsigned char *Image::getPixels()
{
    return data;
}

unsigned char *Image::getMipmapPixels(unsigned mipmap)
{
    if (getTextureType() == TEXTURE_CUBE)
        return 0;
    return data + getMipmapRangeSize(0, mipmap);
}

unsigned char *Image::getCubeMipmapPixels(unsigned cubeside, unsigned mipmap)
{
    if (getTextureType() != TEXTURE_CUBE)
        return 0;
    return data + (cubeside * getImageSize() + getMipmapRangeSize(0, mipmap)) / 6;
}

unsigned char *Image::get3DSliceMipmapPixels(unsigned slice, unsigned mipmap)
{
    if (getTextureType() != TEXTURE_3D)
        return 0;
    return getMipmapPixels(mipmap) + getWidth(mipmap) * getHeight(mipmap) * slice;
}

void Image::clear()
{
    if (data) free(data);
    data = NULL;
}

void Image::init(unsigned width, unsigned height, unsigned depth, unsigned mipmaps,
        DataType dataType, ChannelFormat channelFormat)
{
    clear();
    this->width = width;
    this->height = height;
    this->depth = depth;
    this->mipmaps = mipmaps;
    this->dataType = dataType;
    this->channelFormat = channelFormat;
#if 0
    printf("init: w=%i, h=%i, d=%i, mipmaps=%i, type=%s, channels=%s, textype=%s\n",
            width, height, depth, mipmaps, strFromDataType[dataType],
            strFromChannelFormat[channelFormat], strFromTextureType[getTextureType()]);
#endif
    data = (unsigned char*)malloc(getImageSize());
}

void Image::init(unsigned width, unsigned height, unsigned bpp)
{
    ChannelFormat channels;
    DataType type;

    // Guess the real format.
    switch (bpp)
    {
    case 8:
        type = UNORM8;
        channels = L;
        break;
    case 16:
        type = UNORM8;
        channels = LA;
        break;
    case 24:
        type = UNORM8;
        channels = RGB;
        break;
    case 32:
        type = UNORM8;
        channels = RGBA;
        break;
    case 48:
        type = UNORM16;
        channels = RGB;
        break;
    case 64:
        type = UNORM16;
        channels = RGBA;
        break;
    case 96:
        type = UINT32;
        channels = RGB;
        break;
    case 128:
        type = UINT32;
        channels = RGBA;
        break;
    default:
        std::cerr << "Image::init: Unsupported bpp: " << bpp << std::endl;
        return;
    }

    init(width, height, 1, 1, type, channels);
}

Image* Image::Create(std::string filename)
{
    std::string loader="default";
    std::string::size_type p = filename.rfind('.');
    if (p!=std::string::npos)
        loader = std::string(filename, p+1);
    return FactoryImage::CreateObject(loader, filename);
}

} // namespace io

} // namespace helper

} // namespace sofa

