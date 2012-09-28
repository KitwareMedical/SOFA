/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *
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
*                              SOFA :: Framework                              *
*                                                                             *
* Authors: The SOFA Team (see Authors.txt)                                    *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_CORE_BEHAVIOR_BASEPROJECTIVECONSTRAINTSET_H
#define SOFA_CORE_BEHAVIOR_BASEPROJECTIVECONSTRAINTSET_H

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/behavior/MultiMatrixAccessor.h>

namespace sofa
{

namespace core
{

namespace behavior
{

/**
*  \brief Component computing projective constraints within a simulated body.
*
*  This class define the abstract API common to all projective constraints.
*  A BaseConstraint computes constraints applied to one or more simulated body
*  given its current position and velocity.
*
*  Constraints can be internal to a given body (attached to one MechanicalState,
*  see the Constraint class), or link several bodies together (such as contacts,
*  see the InteractionConstraint class).
*
*/
class SOFA_CORE_API BaseProjectiveConstraintSet : public virtual objectmodel::BaseObject
{
public:
    SOFA_ABSTRACT_CLASS(BaseProjectiveConstraintSet, objectmodel::BaseObject);
protected:
    BaseProjectiveConstraintSet()
        : group(initData(&group, 0, "group", "ID of the group containing this constraint. This ID is used to specify which constraints are solved by which solver, by specifying in each solver which groups of constraints it should handle."))
    {
    }

    virtual ~BaseProjectiveConstraintSet() {}
public:
    /// Get the ID of the group containing this constraint.
    /// This ID is used to specify which constraints are solved by which solver, by specifying in each solver which groups of constraints it should handle.
    int getGroup() const { return group.getValue(); }

    /// Set the ID of the group containing this constraint.
    /// This ID is used to specify which constraints are solved by which solver, by specifying in each solver which groups of constraints it should handle.
    void setGroup(int g) { group.setValue(g); }

    /// @name Vector operations
    /// @{

    /// Project dx to constrained space (dx models an acceleration).
    /// \param dxId output vector
    virtual void projectResponse(const MechanicalParams* mparams /* PARAMS FIRST */, MultiVecDerivId dxId) = 0;

    /// Project the L matrix of the Lagrange Multiplier equation system.
    /// \param cId output vector
    virtual void projectJacobianMatrix(const MechanicalParams* mparams /* PARAMS FIRST */, MultiMatrixDerivId cId) = 0;

    /// Project v to constrained space (v models a velocity).
    /// \param vId output vector
    virtual void projectVelocity(const MechanicalParams* mparams /* PARAMS FIRST */, MultiVecDerivId vId) = 0;

    /// Project x to constrained space (x models a position).
    /// \param xId output vector
    virtual void projectPosition(const MechanicalParams* mparams /* PARAMS FIRST */, MultiVecCoordId xId) = 0;

    /// @}


    /// @name Matrix operations
    /// @{

    /// Project the compliance Matrix to constrained space.
    virtual void projectResponse(const MechanicalParams* /*mparams*/ /* PARAMS FIRST */, double **) {};

    /// Project the global Mechanical Matrix to constrained space using offset parameter
    virtual void applyConstraint(const MechanicalParams* /*mparams*/ /* PARAMS FIRST */, const behavior::MultiMatrixAccessor* /*matrix*/) {};

    /// Project the global Mechanical Vector to constrained space using offset parameter
    virtual void applyConstraint(const MechanicalParams* /*mparams*/ /* PARAMS FIRST */, defaulttype::BaseVector* /*vector*/, const behavior::MultiMatrixAccessor* /*matrix*/) {};

    /** Project the the given matrix (Experimental API).
      Replace M with PMP, where P is the projection matrix corresponding to the projectResponse method. Contrary to applyConstraint(), the diagonal blocks of the result are not reset to the identity.
      Typically, M is the (generalized) mass matrix of the whole system, offset is the starting index of the local state in this global matrix, and P is the identity matrix with a block on the diagonal replaced by the projection matrix.
      If M is the matrix of the local state, then offset should be 0.
      */
    virtual void projectMatrix( sofa::defaulttype::BaseMatrix* /*M*/, unsigned /*offset*/ ) { serr<<"projectMatrix not implemented, projection will not be handled appropriately"<<sendl; }

    /// @}


    /// If the constraint is applied only on a subset of particles.
    /// That way, we can optimize the time spent traversing the mappings
    /// Deactivated by default. The constraints using only a subset of particles should activate the mask,
    /// and during projectResponse(), insert the indices of the particles modified
    virtual bool useMask() const {return false;}

protected:
    Data<int> group;
};

} // namespace behavior

} // namespace core

} // namespace sofa

#endif
