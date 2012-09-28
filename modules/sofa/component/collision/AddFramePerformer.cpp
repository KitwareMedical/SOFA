/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 MGH, INRIA, USTL, UJF, CNRS                    *
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
#ifndef SOFA_COMPONENT_COLLISION_ADDFRAMEPERFORMER_CPP
#define SOFA_COMPONENT_COLLISION_ADDFRAMEPERFORMER_CPP
#include <sofa/component/component.h>
#include <sofa/component/collision/AddFramePerformer.inl>
#include <sofa/defaulttype/Vec3Types.h>
#include <sofa/helper/Factory.inl>

namespace sofa
{

namespace component
{

namespace collision
{

#ifndef SOFA_DOUBLE
template class SOFA_USER_INTERACTION_API  AddFramePerformer<defaulttype::Vec3fTypes>;
#endif
#ifndef SOFA_FLOAT
template class SOFA_USER_INTERACTION_API  AddFramePerformer<defaulttype::Vec3dTypes>;
#endif


#ifndef SOFA_DOUBLE
helper::Creator<InteractionPerformer::InteractionPerformerFactory, AddFramePerformer<defaulttype::Vec3fTypes> >  AddFramePerformerVec3fClass("AddFrame",true);
#endif
#ifndef SOFA_FLOAT
helper::Creator<InteractionPerformer::InteractionPerformerFactory, AddFramePerformer<defaulttype::Vec3dTypes> >  AddFramePerformerVec3dClass("AddFrame",true);
#endif
}
}
}
#endif
