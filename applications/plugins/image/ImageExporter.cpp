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
#define SOFA_IMAGE_IMAGEEXPORTER_CPP


#include "ImageExporter.h"
#include <sofa/core/ObjectFactory.h>


namespace sofa
{

namespace component
{

namespace misc
{

using namespace defaulttype;

SOFA_DECL_CLASS(ImageExporter)

int ImageExporterClass = core::RegisterObject("Save an image")
        .add<ImageExporter<ImageUC> >(true)
        .add<ImageExporter<ImageD> >()
#ifdef BUILD_ALL_IMAGE_TYPES
        .add<ImageExporter<ImageC> >()
        .add<ImageExporter<ImageI> >()
        .add<ImageExporter<ImageUI> >()
        .add<ImageExporter<ImageS> >()
        .add<ImageExporter<ImageUS> >()
        .add<ImageExporter<ImageL> >()
        .add<ImageExporter<ImageUL> >()
        .add<ImageExporter<ImageF> >()
        .add<ImageExporter<ImageB> >()
#endif
        .add<ImageExporter<BranchingImageUC> >()
        .add<ImageExporter<BranchingImageD> >()
#ifdef BUILD_ALL_IMAGE_TYPES
        .add<ImageExporter<BranchingImageC> >()
        .add<ImageExporter<BranchingImageI> >()
        .add<ImageExporter<BranchingImageUI> >()
        .add<ImageExporter<BranchingImageS> >()
        .add<ImageExporter<BranchingImageUS> >()
        .add<ImageExporter<BranchingImageL> >()
        .add<ImageExporter<BranchingImageUL> >()
        .add<ImageExporter<BranchingImageF> >()
        .add<ImageExporter<BranchingImageB> >()
#endif
        ;

template class SOFA_IMAGE_API ImageExporter<ImageUC>;
template class SOFA_IMAGE_API ImageExporter<ImageD>;
#ifdef BUILD_ALL_IMAGE_TYPES
template class SOFA_IMAGE_API ImageExporter<ImageC>;
template class SOFA_IMAGE_API ImageExporter<ImageI>;
template class SOFA_IMAGE_API ImageExporter<ImageUI>;
template class SOFA_IMAGE_API ImageExporter<ImageS>;
template class SOFA_IMAGE_API ImageExporter<ImageUS>;
template class SOFA_IMAGE_API ImageExporter<ImageL>;
template class SOFA_IMAGE_API ImageExporter<ImageUL>;
template class SOFA_IMAGE_API ImageExporter<ImageF>;
template class SOFA_IMAGE_API ImageExporter<ImageB>;
#endif
template class SOFA_IMAGE_API ImageExporter<BranchingImageUC>;
template class SOFA_IMAGE_API ImageExporter<BranchingImageD>;
#ifdef BUILD_ALL_IMAGE_TYPES
template class SOFA_IMAGE_API ImageExporter<BranchingImageC>;
template class SOFA_IMAGE_API ImageExporter<BranchingImageI>;
template class SOFA_IMAGE_API ImageExporter<BranchingImageUI>;
template class SOFA_IMAGE_API ImageExporter<BranchingImageS>;
template class SOFA_IMAGE_API ImageExporter<BranchingImageUS>;
template class SOFA_IMAGE_API ImageExporter<BranchingImageL>;
template class SOFA_IMAGE_API ImageExporter<BranchingImageUL>;
template class SOFA_IMAGE_API ImageExporter<BranchingImageF>;
template class SOFA_IMAGE_API ImageExporter<BranchingImageB>;
#endif


} // namespace misc

} // namespace component

} // namespace sofa
