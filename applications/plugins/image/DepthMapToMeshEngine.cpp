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
#define SOFA_IMAGE_DEPTHMAPTOMESHENGINE_CPP

#include "DepthMapToMeshEngine.h"
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace engine
{

using namespace defaulttype;

SOFA_DECL_CLASS(DepthMapToMeshEngine)

int DepthMapToMeshEngineClass = core::RegisterObject("Compute a mesh from a depth map image ")
        .add<DepthMapToMeshEngine<ImageUC> >(true)
        .add<DepthMapToMeshEngine<ImageD> >()
#ifdef BUILD_ALL_IMAGE_TYPES
        .add<DepthMapToMeshEngine<ImageC> >()
        .add<DepthMapToMeshEngine<ImageI> >()
        .add<DepthMapToMeshEngine<ImageUI> >()
        .add<DepthMapToMeshEngine<ImageS> >()
        .add<DepthMapToMeshEngine<ImageUS> >()
        .add<DepthMapToMeshEngine<ImageL> >()
        .add<DepthMapToMeshEngine<ImageUL> >()
        .add<DepthMapToMeshEngine<ImageF> >()
        .add<DepthMapToMeshEngine<ImageB> >()
#endif
        ;

template class SOFA_IMAGE_API DepthMapToMeshEngine<ImageUC>;
template class SOFA_IMAGE_API DepthMapToMeshEngine<ImageD>;
#ifdef BUILD_ALL_IMAGE_TYPES
template class SOFA_IMAGE_API DepthMapToMeshEngine<ImageC>;
template class SOFA_IMAGE_API DepthMapToMeshEngine<ImageI>;
template class SOFA_IMAGE_API DepthMapToMeshEngine<ImageUI>;
template class SOFA_IMAGE_API DepthMapToMeshEngine<ImageS>;
template class SOFA_IMAGE_API DepthMapToMeshEngine<ImageUS>;
template class SOFA_IMAGE_API DepthMapToMeshEngine<ImageL>;
template class SOFA_IMAGE_API DepthMapToMeshEngine<ImageUL>;
template class SOFA_IMAGE_API DepthMapToMeshEngine<ImageF>;
template class SOFA_IMAGE_API DepthMapToMeshEngine<ImageB>;
#endif

} //
} // namespace component

} // namespace sofa

