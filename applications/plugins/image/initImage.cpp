/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 4      *
*                (c) 2006-2009 MGH, INRIA, USTL, UJF, CNRS                    *
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
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#include <sofa/helper/system/config.h>
#include "initImage.h"

namespace sofa
{

namespace component
{

//Here are just several convenient functions to help user to know what contains the plugin

extern "C" {
    SOFA_IMAGE_API void initExternalModule();
    SOFA_IMAGE_API const char* getModuleName();
    SOFA_IMAGE_API const char* getModuleVersion();
    SOFA_IMAGE_API const char* getModuleLicense();
    SOFA_IMAGE_API const char* getModuleDescription();
    SOFA_IMAGE_API const char* getModuleComponentList();
}

void initExternalModule()
{
    static bool first = true;
    if (first)
    {
        first = false;
    }
}

const char* getModuleName()
{
    return "Image Plugin";
}

const char* getModuleVersion()
{
    return "0.1";
}

const char* getModuleLicense()
{
    return "LGPL";
}


const char* getModuleDescription()
{
    return "Image support in SOFA";
}

const char* getModuleComponentList()
{
    return "ImageContainer,ImageExporter,ImageViewer,ImageFilter,ImageToMeshEngine";
}

} // namespace image

} // namespace sofa

////////// BEGIN CLASS LIST //////////
SOFA_LINK_CLASS(ImageContainer)
SOFA_LINK_CLASS(ImageExporter)
SOFA_LINK_CLASS(ImageViewer)
SOFA_LINK_CLASS(ImageFilter)
SOFA_LINK_CLASS(MergeImages)
SOFA_LINK_CLASS(DepthMapToMeshEngine)
SOFA_LINK_CLASS(MeshToImageEngine)
SOFA_LINK_CLASS(ImageAccumulator)
SOFA_LINK_CLASS(ImageSampler)
#ifdef SOFA_HAVE_LIBFREENECT
SOFA_LINK_CLASS(Kinect)
#endif

