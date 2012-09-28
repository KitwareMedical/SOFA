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
#include "CudaTypes.h"
#include "CudaSpringForceField.inl"
#include "CudaMechanicalObject.inl"
#include "CudaIdentityMapping.inl"
#include "CudaContactMapper.h"
#include "CudaPenalityContactForceField.h"
#include "CudaSpringForceField.h"
#include "CudaSphereModel.h"
#include "CudaPointModel.h"

#include <sofa/component/collision/MouseInteractor.inl>
#include <sofa/component/collision/NewProximityIntersection.inl>
#include <sofa/component/collision/MeshNewProximityIntersection.inl>
#include <sofa/component/collision/RayDiscreteIntersection.inl>
#include <sofa/component/collision/DiscreteIntersection.inl>
#include <sofa/component/collision/ComponentMouseInteraction.inl>
#include <sofa/component/collision/AttachBodyPerformer.inl>
#include <sofa/component/collision/FixParticlePerformer.inl>
#include <sofa/component/collision/RayContact.h>
#include <sofa/component/collision/BarycentricPenalityContact.inl>
#include <sofa/component/collision/BarycentricContactMapper.inl>
#include <sofa/component/interactionforcefield/PenalityContactForceField.h>
#include <sofa/component/interactionforcefield/VectorSpringForceField.h>
#include <sofa/helper/system/gl.h>
#include <sofa/helper/Factory.inl>
#include <fstream>

namespace sofa
{

namespace component
{

namespace collision
{

using namespace sofa::gpu::cuda;


template class MouseInteractor<CudaVec3fTypes>;
template class TComponentMouseInteraction< CudaVec3fTypes >;
template class AttachBodyPerformer< CudaVec3fTypes >;
template class FixParticlePerformer< CudaVec3fTypes >;

#ifdef SOFA_GPU_CUDA_DOUBLE
template class MouseInteractor<CudaVec3dTypes>;
template class TComponentMouseInteraction< CudaVec3dTypes >;
template class AttachBodyPerformer< CudaVec3dTypes >;
template class FixParticlePerformer< CudaVec3dTypes >;
#endif

ContactMapperCreator< ContactMapper<CudaSphereModel> > CudaSphereContactMapperClass("default",true);

helper::Creator<ComponentMouseInteraction::ComponentMouseInteractionFactory, TComponentMouseInteraction<CudaVec3fTypes> > ComponentMouseInteractionCudaVec3fClass ("MouseSpringCudaVec3f",true);
helper::Creator<InteractionPerformer::InteractionPerformerFactory, AttachBodyPerformer <CudaVec3fTypes> >  AttachBodyPerformerCudaVec3fClass("AttachBody",true);
helper::Creator<InteractionPerformer::InteractionPerformerFactory, FixParticlePerformer<CudaVec3fTypes> >  FixParticlePerformerCudaVec3fClass("FixParticle",true);

#ifdef SOFA_GPU_CUDA_DOUBLE
helper::Creator<ComponentMouseInteraction::ComponentMouseInteractionFactory, TComponentMouseInteraction<CudaVec3dTypes> > ComponentMouseInteractionCudaVec3dClass ("MouseSpringCudaVec3d",true);
helper::Creator<InteractionPerformer::InteractionPerformerFactory, AttachBodyPerformer <CudaVec3dTypes> >  AttachBodyPerformerCudaVec3dClass("AttachBody",true);
helper::Creator<InteractionPerformer::InteractionPerformerFactory, FixParticlePerformer<CudaVec3dTypes> >  FixParticlePerformerCudaVec3dClass("FixParticle",true);
#endif

} //namespace collision


} //namespace component


namespace gpu
{

namespace cuda
{


SOFA_DECL_CLASS(CudaMouseInteractor)

int MouseInteractorCudaClass = core::RegisterObject("Supports Mouse Interaction using CUDA")
        .add< component::collision::MouseInteractor<CudaVec3fTypes> >()
#ifdef SOFA_GPU_CUDA_DOUBLE
        .add< component::collision::MouseInteractor<CudaVec3dTypes> >()
#endif
        ;


SOFA_DECL_CLASS(CudaCollision)

using namespace sofa::component::collision;

class CudaProximityIntersection : public sofa::component::collision::NewProximityIntersection
{
public:
    SOFA_CLASS(CudaProximityIntersection,sofa::component::collision::NewProximityIntersection);

    virtual void init()
    {
        sofa::component::collision::NewProximityIntersection::init();
        intersectors.add<CudaSphereModel, CudaSphereModel,   DiscreteIntersection>(this);
        RayDiscreteIntersection* rayIntersector = new RayDiscreteIntersection(this, false);
        intersectors.add<RayModel,        CudaSphereModel,   RayDiscreteIntersection>(rayIntersector);
        MeshNewProximityIntersection* meshIntersector = new MeshNewProximityIntersection(this, false);
        intersectors.add<TriangleModel,   CudaSphereModel,   MeshNewProximityIntersection>(meshIntersector);
    }

};


int CudaProximityIntersectionClass = core::RegisterObject("GPGPU Proximity Intersection based on CUDA")
        .add< CudaProximityIntersection >()
        ;

sofa::helper::Creator<core::collision::Contact::Factory, component::collision::RayContact<CudaSphereModel> > RayCudaSphereContactClass("ray",true);

} // namespace cuda

} // namespace gpu

} // namespace sofa
