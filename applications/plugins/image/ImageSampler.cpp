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
#define SOFA_IMAGE_IMAGESAMPLER_CPP

#include "ImageSampler.h"
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace engine
{

using namespace defaulttype;

SOFA_DECL_CLASS(ImageSampler)

int ImageSamplerClass = core::RegisterObject("Samples an object represented by an image")
        .add<ImageSampler<ImageB> >(true)
        .add<ImageSampler<ImageUC> >()
        .add<ImageSampler<ImageD> >()
#ifdef BUILD_ALL_IMAGE_TYPES
//.add<ImageSampler<ImageC> >()
//.add<ImageSampler<ImageI> >()
//.add<ImageSampler<ImageUI> >()
//.add<ImageSampler<ImageS> >()
//.add<ImageSampler<ImageUS> >()
//.add<ImageSampler<ImageL> >()
//.add<ImageSampler<ImageUL> >()
//.add<ImageSampler<ImageF> >()
#endif
        .add<ImageSampler<BranchingImageB> >()
        .add<ImageSampler<BranchingImageUC> >()
        .add<ImageSampler<BranchingImageD> >()
        ;

template class SOFA_IMAGE_API ImageSampler<ImageB>;
template class SOFA_IMAGE_API ImageSampler<ImageUC>;
template class SOFA_IMAGE_API ImageSampler<ImageD>;
#ifdef BUILD_ALL_IMAGE_TYPES
//template class SOFA_IMAGE_API ImageSampler<ImageC>;
//template class SOFA_IMAGE_API ImageSampler<ImageI>;
//template class SOFA_IMAGE_API ImageSampler<ImageUI>;
//template class SOFA_IMAGE_API ImageSampler<ImageS>;
//template class SOFA_IMAGE_API ImageSampler<ImageUS>;
//template class SOFA_IMAGE_API ImageSampler<ImageL>;
//template class SOFA_IMAGE_API ImageSampler<ImageUL>;
//template class SOFA_IMAGE_API ImageSampler<ImageF>;
#endif
template class SOFA_IMAGE_API ImageSampler<BranchingImageB>;
template class SOFA_IMAGE_API ImageSampler<BranchingImageUC>;
template class SOFA_IMAGE_API ImageSampler<BranchingImageD>;



} //
} // namespace component

} // namespace sofa

