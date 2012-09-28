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
#define SOFA_IMAGE_IMAGECONTAINER_CPP

#include "ImageContainer.h"
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace container
{

using namespace defaulttype;


SOFA_DECL_CLASS (ImageContainer);
// Register in the Factory

int ImageContainerClass = core::RegisterObject ( "Image Container" )
        .add<ImageContainer<ImageUC> >(true)
        .add<ImageContainer<ImageD> >()
#ifdef BUILD_ALL_IMAGE_TYPES
        .add<ImageContainer<ImageC> >()
        .add<ImageContainer<ImageI> >()
        .add<ImageContainer<ImageUI> >()
        .add<ImageContainer<ImageS> >()
        .add<ImageContainer<ImageUS> >()
        .add<ImageContainer<ImageL> >()
        .add<ImageContainer<ImageUL> >()
        .add<ImageContainer<ImageF> >()
        .add<ImageContainer<ImageB> >()
#endif
        ;

template class SOFA_IMAGE_API ImageContainer<ImageUC>;
template class SOFA_IMAGE_API ImageContainer<ImageD>;
#ifdef BUILD_ALL_IMAGE_TYPES
template class SOFA_IMAGE_API ImageContainer<ImageC>;
template class SOFA_IMAGE_API ImageContainer<ImageI>;
template class SOFA_IMAGE_API ImageContainer<ImageUI>;
template class SOFA_IMAGE_API ImageContainer<ImageS>;
template class SOFA_IMAGE_API ImageContainer<ImageUS>;
template class SOFA_IMAGE_API ImageContainer<ImageL>;
template class SOFA_IMAGE_API ImageContainer<ImageUL>;
template class SOFA_IMAGE_API ImageContainer<ImageF>;
template class SOFA_IMAGE_API ImageContainer<ImageB>;
#endif

} // namespace container

} // namespace component

} // namespace sofa
