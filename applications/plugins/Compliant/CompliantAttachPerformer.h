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
#ifndef SOFA_COMPONENT_COLLISION_CompliantAttachPerformer_H
#define SOFA_COMPONENT_COLLISION_CompliantAttachPerformer_H

#include "initCompliant.h"
#include <sofa/component/collision/InteractionPerformer.h>
#include <sofa/component/collision/BaseContactMapper.h>
#include <../applications/plugins/Flexible/deformationMapping/DistanceMapping.h>
#include <sofa/component/container/MechanicalObject.h>
#include <sofa/defaulttype/Vec3Types.h>
#include <sofa/simulation/common/Node.h>
#include <sofa/gui/MouseOperations.h>


namespace sofa
{
using defaulttype::Vec;

namespace gui
{
class SOFA_Compliant_API CompliantAttachOperation : public Operation
{
public:
    virtual void start() ;
    virtual void execution() ;
    virtual void end() ;
    virtual void endOperation() ;
    static std::string getDescription() {return "CompliantAttach";}
};
}

namespace component
{

namespace collision
{
struct BodyPicked;

/** Mouse interaction using a compliance forcefield, for an object animated using a compliance solver.

  @author Francois Faure, 2012
  */
template <class DataTypes>
class SOFA_Compliant_API CompliantAttachPerformer: public TInteractionPerformer<DataTypes>
{
    typedef typename DataTypes::Real                                  Real;
    typedef defaulttype::StdVectorTypes< Vec<1,Real>, Vec<1,Real>  >  DataTypes1;
    typedef mapping::DistanceMapping< DataTypes,DataTypes1 >          DistanceMapping31;
    typedef sofa::component::container::MechanicalObject< DataTypes > Point3dState;

    simulation::Node::SPtr pickedNode;       ///< Node containing the picked MechanicalState
    int pickedParticleIndex;                 ///< Index of the picked particle in the picked state
    simulation::Node::SPtr interactionNode;  ///< Node used to create the interaction components to constrain the picked point
    core::BaseMapping::SPtr mouseMapping;   ///< Mapping from the mouse position to the 3D point on the ray
    sofa::component::collision::BaseContactMapper< DataTypes >  *mapper;
    Point3dState* mouseState;                  ///< Mouse state container  (position, velocity)
    typename DistanceMapping31::SPtr distanceMapping; ///< computes the distance from the picked point to its target

    void clear();                             ///< release the current interaction

public:
    CompliantAttachPerformer(BaseMouseInteractor *i);
    ~CompliantAttachPerformer();

    void start();
    void execute();

};



#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_COLLISION_CompliantAttachPerformer_CPP)
#ifndef SOFA_DOUBLE
extern template class SOFA_Compliant_API  CompliantAttachPerformer<defaulttype::Vec3fTypes>;
#endif
#ifndef SOFA_FLOAT
extern template class SOFA_Compliant_API  CompliantAttachPerformer<defaulttype::Vec3dTypes>;
#endif
#endif


}
}
}

#endif
