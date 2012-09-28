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
#include <sofa/helper/io/ImageRAW.h>
#include <sofa/helper/system/FileRepository.h>
#include <iostream>

#include <cstdio>		// fopen

namespace sofa
{

namespace helper
{

namespace io
{

ImageRAW::ImageRAW ()
    : headerSize(0)
{}

void ImageRAW::initHeader(unsigned hsize)
{
    headerSize = hsize;
    header = (unsigned char*) malloc(headerSize);
}

bool ImageRAW::load(std::string filename)
{
    if (!sofa::helper::system::DataRepository.findFile(filename))
    {
        std::cerr << "File " << filename << " not found " << std::endl;
        return false;
    }
    FILE *file;
    /* make sure the file is there and open it read-only (binary) */
    if ((file = fopen(filename.c_str(), "rb")) == NULL)
    {
        std::cerr << "File not found : " << filename << std::endl;
        return false;
    }

    // read header and ignore it as we don't know how to interpret it
    for ( unsigned i=0; i<headerSize; ++i )
    {
        int c = getc ( file );

        if ( c == EOF )
        {
            fclose ( file );
            return false;
        }
        else
            header[i] = ( unsigned char ) c;
    }

    const unsigned int numVoxels = getImageSize();

    // get the voxels from the file
    unsigned char *data = getPixels();
    for ( unsigned int i=0; i<numVoxels; ++i )
    {
        int c = getc ( file );

        if ( c == EOF )
        {
            fclose ( file );
            return false;
        }
        else
            data[i] = ( unsigned char ) c;
    }

    fclose(file);

    return true;
}

bool ImageRAW::save(std::string filename, int)
{
    FILE *file;
#ifndef NDEBUG
    std::cout << "Writing RAW file " << filename << std::endl;
#endif
    /* make sure the file is there and open it read-only (binary) */
    if ((file = fopen(filename.c_str(), "wb")) == NULL)
    {
        std::cerr << "File write access failed : " << filename << std::endl;
        return false;
    }

    bool isWriteOk = true;
    if (headerSize > 0)
    {
        isWriteOk = isWriteOk && fwrite(header, 1, headerSize, file) == headerSize;
    }
    isWriteOk = isWriteOk && fwrite(getPixels(), 1, getImageSize(), file) == getImageSize();
    fclose(file);
    return isWriteOk;
}

} // namespace io

} // namespace helper

} // namespace sofa

