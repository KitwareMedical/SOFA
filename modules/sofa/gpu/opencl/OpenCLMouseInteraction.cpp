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
#include "OpenCLTypes.h"

#include <sofa/gpu/opencl/OpenCLMechanicalObject.h>
#include <sofa/gpu/opencl/OpenCLIdentityMapping.h>
#include <sofa/gpu/opencl/OpenCLFixedConstraint.h>
#include <sofa/gpu/opencl/OpenCLSpringForceField.h>

#include <sofa/component/collision/MouseInteractor.inl>
#include <sofa/component/collision/ComponentMouseInteraction.inl>
#include <sofa/component/collision/AttachBodyPerformer.inl>
#include <sofa/component/collision/FixParticlePerformer.inl>
#include <sofa/helper/Factory.inl>
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace collision
{

using namespace sofa::gpu::opencl;

template class MouseInteractor<OpenCLVec3fTypes>;
template class TComponentMouseInteraction< OpenCLVec3fTypes >;
template class AttachBodyPerformer< OpenCLVec3fTypes >;
template class FixParticlePerformer< OpenCLVec3fTypes >;

#ifdef SOFA_GPU_OPENCL_DOUBLE
template class MouseInteractor<OpenCLVec3dTypes>;
template class TComponentMouseInteraction< OpenCLVec3dTypes >;
template class AttachBodyPerformer< OpenCLVec3dTypes >;
template class FixParticlePerformer< OpenCLVec3dTypes >;
#endif

helper::Creator<ComponentMouseInteraction::ComponentMouseInteractionFactory, TComponentMouseInteraction<OpenCLVec3fTypes> > ComponentMouseInteractionOpenCLVec3fClass ("MouseSpringOpenCLVec3f",true);
helper::Creator<InteractionPerformer::InteractionPerformerFactory, AttachBodyPerformer <OpenCLVec3fTypes> >  AttachBodyPerformerOpenCLVec3fClass("AttachBody",true);
helper::Creator<InteractionPerformer::InteractionPerformerFactory, FixParticlePerformer<OpenCLVec3fTypes> >  FixParticlePerformerOpenCLVec3fClass("FixParticle",true);

#ifdef SOFA_GPU_OPENCL_DOUBLE
helper::Creator<ComponentMouseInteraction::ComponentMouseInteractionFactory, TComponentMouseInteraction<OpenCLVec3dTypes> > ComponentMouseInteractionOpenCLVec3dClass ("MouseSpringOpenCLVec3d",true);
helper::Creator<InteractionPerformer::InteractionPerformerFactory, AttachBodyPerformer <OpenCLVec3dTypes> >  AttachBodyPerformerOpenCLVec3dClass("AttachBody",true);
helper::Creator<InteractionPerformer::InteractionPerformerFactory, FixParticlePerformer<OpenCLVec3dTypes> >  FixParticlePerformerOpenCLVec3dClass("FixParticle",true);
#endif

} //namespace collision

} //namespace component


namespace gpu
{

namespace opencl
{

SOFA_DECL_CLASS(OpenCLMouseInteractor)

int MouseInteractorOpenCLClass = core::RegisterObject("Supports Mouse Interaction using OPENCL")
        .add< component::collision::MouseInteractor<OpenCLVec3fTypes> >()
#ifdef SOFA_GPU_OPENCL_DOUBLE
        .add< component::collision::MouseInteractor<OpenCLVec3dTypes> >()
#endif
        ;

}

}

}
