//
// C++ Interface: LDIPenalityContactForceField
//
// Description: Interaction ForceField used to create a repulsion, and separate two colliding objects
//
//
// Author: Francois Faure, Sebastien Barbier, Jeremie Allard, Florent Falipou
//
// Licence: QPL, See LICENCE.txt file that comes with this distribution
// Copyright: See COPYING file that comes with this distribution
//
//
// #include <sofa/component/forcefield/LDIPenalityContactForceField.inl>
#include "LDIPenalityContactForceField.inl"
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/component/container/MechanicalObject.h>
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace forcefield
{

using namespace sofa::defaulttype;

SOFA_DECL_CLASS(LDIPenalityContactForceField)

// Register in the Factory
int LDIPenalityContactForceFieldClass = core::RegisterObject("LDI Contact")
#ifndef SOFA_FLOAT
.add< LDIPenalityContactForceField<Vec3dTypes,Vec3dTypes,Vec3dTypes> >()
//.add< LDIPenalityContactForceField<Rigid3dTypes,Vec3dTypes,Vec3dTypes> >()
//.add< LDIPenalityContactForceField<Vec3dTypes,Rigid3dTypes,Vec3dTypes> >()
//.add< LDIPenalityContactForceField<Rigid3dTypes,Rigid3dTypes,Vec3dTypes> >()
#endif
#ifndef SOFA_DOUBLE
.add< LDIPenalityContactForceField<Vec3fTypes,Vec3fTypes,Vec3fTypes> >()
//.add< LDIPenalityContactForceField<Rigid3fTypes,Vec3fTypes,Vec3fTypes> >()
//.add< LDIPenalityContactForceField<Vec3fTypes,Rigid3fTypes,Vec3fTypes> >()
//.add< LDIPenalityContactForceField<Rigid3fTypes,Rigid3fTypes,Vec3fTypes> >()
#endif

    .addLicense("QPL")
    .addAuthor("Francois Faure, Sebastien Barbier, Jeremie Allard, Florent Falipou");
;

#ifndef SOFA_FLOAT
template class LDIPenalityContactForceField<Vec3dTypes,Vec3dTypes,Vec3dTypes>;
//template class LDIPenalityContactForceField<Rigid3dTypes,Vec3dTypes,Vec3dTypes>;
//template class LDIPenalityContactForceField<Vec3dTypes,Rigid3dTypes,Vec3dTypes>;
//template class LDIPenalityContactForceField<Rigid3dTypes,Rigid3dTypes,Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
template class LDIPenalityContactForceField<Vec3fTypes,Vec3fTypes,Vec3fTypes>;
//template class LDIPenalityContactForceField<Rigid3fTypes,Vec3fTypes,Vec3fTypes>;
//template class LDIPenalityContactForceField<Vec3fTypes,Rigid3fTypes,Vec3fTypes>;
//template class LDIPenalityContactForceField<Rigid3fTypes,Rigid3fTypes,Vec3fTypes>;
#endif

/*
#ifndef SOFA_FLOAT
#ifndef SOFA_DOUBLE
template class LDIPenalityContactForceField<Vec3dTypes,Vec3dTypes,Vec3fTypes>;
template class LDIPenalityContactForceField<Rigid3dTypes,Vec3dTypes,Vec3fTypes>;
template class LDIPenalityContactForceField<Vec3dTypes,Rigid3dTypes,Vec3fTypes>;
template class LDIPenalityContactForceField<Rigid3dTypes,Rigid3dTypes,Vec3fTypes>;
template class LDIPenalityContactForceField<Vec3fTypes,Vec3fTypes,Vec3dTypes>;
template class LDIPenalityContactForceField<Rigid3fTypes,Vec3fTypes,Vec3dTypes>;
template class LDIPenalityContactForceField<Vec3fTypes,Rigid3fTypes,Vec3dTypes>;
template class LDIPenalityContactForceField<Rigid3fTypes,Rigid3fTypes,Vec3dTypes>;
#endif
#endif*/


} // namespace forcefield

} // namespace component

} // namespace sofa

