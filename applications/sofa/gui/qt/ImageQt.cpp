/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU General Public License as published by the Free  *
* Software Foundation; either version 2 of the License, or (at your option)   *
* any later version.                                                          *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    *
* more details.                                                               *
*                                                                             *
* You should have received a copy of the GNU General Public License along     *
* with this program; if not, write to the Free Software Foundation, Inc., 51  *
* Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.                   *
*******************************************************************************
*                            SOFA :: Applications                             *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#include "ImageQt.h"

#include <sofa/helper/system/FileRepository.h>
#include <sofa/helper/system/SetDirectory.h>
#include <sofa/helper/Factory.inl>

#ifdef SOFA_QT4
#include <QImage>
#include <QList>
#include <QImageReader>
#include <QApplication>
#else
#include <qimage.h>
#endif

#include <string.h>

namespace sofa
{

namespace gui
{

namespace qt
{

SOFA_DECL_CLASS(ImageQt)

class ImageQtCreators
{
public:
    std::vector<sofa::helper::io::Image::FactoryImage::Creator*> creators;
    ImageQtCreators()
    {
#ifdef SOFA_QT4
#ifdef WIN32
        std::string plugdir = sofa::helper::system::SetDirectory::GetRelativeFromProcess("../tools/qt4win/plugins");
        QCoreApplication::addLibraryPath ( QString(plugdir.c_str()) );
#endif
        QList<QByteArray> formats = QImageReader::supportedImageFormats();
        for (QList<QByteArray>::const_iterator it = formats.begin(); it != formats.end(); ++it)
        {
            const char* format = it->data();
#else
        QStrList formats = QImageIO::inputFormats();
        for (const char* format = formats.first(); format; format = formats.next())
        {
#endif
            // ignore if format already supported
            if (sofa::helper::io::Image::FactoryImage::HasKey(format)) continue;
            //std::cout << "ImageQt: supporting format "<<format<<std::endl;
            creators.push_back(new sofa::helper::Creator<sofa::helper::io::Image::FactoryImage, ImageQt>(format));
        }
    }
};


bool ImageQt::Init()
{
    static ImageQtCreators imageQtCreators;
    return true;
}

bool ImageQt::load(std::string filename)
{
    QImage qin(filename.c_str());
    if (qin.isNull())
    {
        return false;
    }
    if (qin.depth() <= 8 && !qin.isGrayscale()) // paletted -> not supported by sofa
    {
        QImage orig = qin;
        qin = orig.convertDepth(32);
    }
    if (qin.depth() >= 24)
    {
        QImage orig = qin;
        qin = orig.swapRGB();
    }
    Image::init(qin.width(), qin.height(), qin.depth());
    int lsize = ((int)qin.bytesPerLine() < (int)getLineSize() ? (int)qin.bytesPerLine() : (int)getLineSize());
    int height = getHeight();
    for (int y=0; y<height; ++y)
    {
        const uchar* in = qin.scanLine(y);
        unsigned char* out = getPixels() + y*getLineSize();
        memcpy(out, in, lsize);
    }
    return true;
}

} //namespace qt

} //namespace gui

} //namespace sofa
