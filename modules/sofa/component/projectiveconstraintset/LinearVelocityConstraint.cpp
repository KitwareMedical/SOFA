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
#define SOFA_COMPONENT_PROJECTIVECONSTRAINTSET_LINEARVELOCITYCONSTRAINT_CPP
#include <sofa/component/projectiveconstraintset/LinearVelocityConstraint.inl>
#include <sofa/core/behavior/ProjectiveConstraintSet.inl>
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/Vec3Types.h>
#include <sofa/defaulttype/RigidTypes.h>

#include <sofa/simulation/common/Node.h>

namespace sofa
{

namespace component
{

namespace projectiveconstraintset
{

using namespace sofa::defaulttype;
using namespace sofa::helper;


//display specialisation for rigid types
#ifndef SOFA_FLOAT
template <>
void LinearVelocityConstraint<Rigid3dTypes>::draw(const core::visual::VisualParams* vparams)
{
    const SetIndexArray & indices = m_indices.getValue();
    if (!vparams->displayFlags().getShowBehaviorModels()) return;
    const VecCoord& x = *mstate->getX();
    glDisable (GL_LIGHTING);
    glPointSize(10);
    glColor4f (0,0.5,0.5,0);
    glBegin (GL_POINTS);
    for (SetIndex::const_iterator it = indices.begin(); it != indices.end(); ++it)
        gl::glVertexT(x[*it].getCenter());
    glEnd();
}
template <>
void LinearVelocityConstraint<Vec6dTypes>::draw(const core::visual::VisualParams* vparams)
{
    const SetIndexArray & indices = m_indices.getValue();
    if (!vparams->displayFlags().getShowBehaviorModels()) return;
    const VecCoord& x = *mstate->getX();
    glDisable (GL_LIGHTING);
    glPointSize(10);
    glColor4f (0,0.5,0.5,0);
    glBegin (GL_POINTS);
    for (SetIndex::const_iterator it = indices.begin(); it != indices.end(); ++it)
        gl::glVertexT(Vec<3,double>(x[*it][0], x[*it][1], x[*it][2]));
    glEnd();
}
#endif
#ifndef SOFA_DOUBLE
template <>
void LinearVelocityConstraint<Rigid3fTypes>::draw(const core::visual::VisualParams* vparams)
{
    const SetIndexArray & indices = m_indices.getValue();
    if (!vparams->displayFlags().getShowBehaviorModels()) return;
    const VecCoord& x = *mstate->getX();
    glDisable (GL_LIGHTING);
    glPointSize(10);
    glColor4f (0,0.5,0.5,0);
    glBegin (GL_POINTS);
    for (SetIndex::const_iterator it = indices.begin(); it != indices.end(); ++it)
        gl::glVertexT(x[*it].getCenter());
    glEnd();
}
template <>
void LinearVelocityConstraint<Vec6fTypes>::draw(const core::visual::VisualParams* vparams)
{
    const SetIndexArray & indices = m_indices.getValue();
    if (!vparams->displayFlags().getShowBehaviorModels()) return;
    const VecCoord& x = *mstate->getX();
    glDisable (GL_LIGHTING);
    glPointSize(10);
    glColor4f (0,0.5,0.5,0);
    glBegin (GL_POINTS);
    for (SetIndex::const_iterator it = indices.begin(); it != indices.end(); ++it)
        gl::glVertexT(Vec<3,float>(x[*it][0], x[*it][1], x[*it][2]));
    glEnd();
}
#endif

//declaration of the class, for the factory
SOFA_DECL_CLASS(LinearVelocityConstraint)


int LinearVelocityConstraintClass = core::RegisterObject("apply velocity to given particles")
#ifndef SOFA_FLOAT
        .add< LinearVelocityConstraint<Vec3dTypes> >()
        .add< LinearVelocityConstraint<Vec2dTypes> >()
        .add< LinearVelocityConstraint<Vec1dTypes> >()
        .add< LinearVelocityConstraint<Vec6dTypes> >()
        .add< LinearVelocityConstraint<Rigid3dTypes> >()
#endif
#ifndef SOFA_DOUBLE
        .add< LinearVelocityConstraint<Vec3fTypes> >()
        .add< LinearVelocityConstraint<Vec2fTypes> >()
        .add< LinearVelocityConstraint<Vec1fTypes> >()
        .add< LinearVelocityConstraint<Vec6fTypes> >()
        .add< LinearVelocityConstraint<Rigid3fTypes> >()
#endif
        ;

#ifndef SOFA_FLOAT
template class SOFA_BOUNDARY_CONDITION_API LinearVelocityConstraint<Vec3dTypes>;
template class SOFA_BOUNDARY_CONDITION_API LinearVelocityConstraint<Vec2dTypes>;
template class SOFA_BOUNDARY_CONDITION_API LinearVelocityConstraint<Vec1dTypes>;
template class SOFA_BOUNDARY_CONDITION_API LinearVelocityConstraint<Vec6dTypes>;
template class SOFA_BOUNDARY_CONDITION_API LinearVelocityConstraint<Rigid3dTypes>;
#endif
#ifndef SOFA_DOUBLE
template class SOFA_BOUNDARY_CONDITION_API LinearVelocityConstraint<Vec3fTypes>;
template class SOFA_BOUNDARY_CONDITION_API LinearVelocityConstraint<Vec2fTypes>;
template class SOFA_BOUNDARY_CONDITION_API LinearVelocityConstraint<Vec1fTypes>;
template class SOFA_BOUNDARY_CONDITION_API LinearVelocityConstraint<Vec6fTypes>;
template class SOFA_BOUNDARY_CONDITION_API LinearVelocityConstraint<Rigid3fTypes>;
#endif

} // namespace projectiveconstraintset

} // namespace component

} // namespace sofa

