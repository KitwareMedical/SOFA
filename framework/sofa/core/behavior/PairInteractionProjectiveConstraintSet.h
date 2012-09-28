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
#ifndef SOFA_CORE_BEHAVIOR_PAIRINTERACTIONPROJECTIVECONSTRAINTSET_H
#define SOFA_CORE_BEHAVIOR_PAIRINTERACTIONPROJECTIVECONSTRAINTSET_H

#include <sofa/core/core.h>
#include <sofa/core/behavior/BaseInteractionProjectiveConstraintSet.h>
#include <sofa/core/behavior/MechanicalState.h>

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>

namespace sofa
{

namespace core
{

namespace behavior
{

/**
 *  \brief Component computing constraints between a pair of simulated body.
 *
 *  This class define the abstract API common to interaction constraints
 *  between a pair of bodies using a given type of DOFs.
 */
template<class TDataTypes>
class PairInteractionProjectiveConstraintSet : public BaseInteractionProjectiveConstraintSet
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(PairInteractionProjectiveConstraintSet,TDataTypes), BaseInteractionProjectiveConstraintSet);

    typedef TDataTypes DataTypes;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::MatrixDeriv MatrixDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef objectmodel::Data<VecCoord> DataVecCoord;
    typedef objectmodel::Data<VecDeriv> DataVecDeriv;
    typedef typename DataTypes::Real Real;
    typedef helper::ParticleMask ParticleMask;
protected:
    PairInteractionProjectiveConstraintSet(MechanicalState<DataTypes> *mm1 = NULL, MechanicalState<DataTypes> *mm2 = NULL);

    virtual ~PairInteractionProjectiveConstraintSet();
public:
    Data<double> endTime;  ///< Time when the constraint becomes inactive (-1 for infinitely active)
    virtual bool isActive() const; ///< if false, the constraint does nothing

    virtual void init();

    /// Retrieve the associated MechanicalState
    MechanicalState<DataTypes>* getMState1() { return mstate1.get(); }
    BaseMechanicalState* getMechModel1() { return mstate1.get(); }
    /// Retrieve the associated MechanicalState
    MechanicalState<DataTypes>* getMState2() { return mstate2.get(); }
    BaseMechanicalState* getMechModel2() { return mstate2.get(); }

    /// @name Vector operations
    /// @{

    /// Project dx to constrained space (dx models an acceleration).
    ///
    /// This method retrieves the dx vector from the MechanicalState and call
    /// the internal projectResponse(VecDeriv&,VecDeriv&) method implemented by
    /// the component.
    virtual void projectResponse(const MechanicalParams* mparams /* PARAMS FIRST */, MultiVecDerivId dxId);

    /// Project the L matrix of the Lagrange Multiplier equation system.
    ///
    /// This method retrieves the lines of the Jacobian Matrix from the MechanicalState and call
    /// the internal projectResponse(MatrixDeriv&) method implemented by
    /// the component.
    virtual void projectJacobianMatrix(const MechanicalParams* mparams /* PARAMS FIRST */, MultiMatrixDerivId cId);

    /// Project v to constrained space (v models a velocity).
    ///
    /// This method retrieves the v vector from the MechanicalState and call
    /// the internal projectVelocity(VecDeriv&,VecDeriv&) method implemented by
    /// the component.
    virtual void projectVelocity(const MechanicalParams* mparams /* PARAMS FIRST */, MultiVecDerivId vId);

    /// Project x to constrained space (x models a position).
    ///
    /// This method retrieves the x vector from the MechanicalState and call
    /// the internal projectPosition(VecCoord&,VecCoord&) method implemented by
    /// the component.
    virtual void projectPosition(const MechanicalParams* mparams /* PARAMS FIRST */, MultiVecCoordId xId);

    /// Project dx to constrained space (dx models an acceleration).
    virtual void projectResponse(const MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecDeriv& dx1, DataVecDeriv& dx2)
#ifdef SOFA_DEPRECATE_OLD_API
        = 0;
#else
    {
        projectResponse(*dx1.beginEdit(), *dx2.beginEdit());
        dx1.endEdit();
        dx2.endEdit();
    }

    /// Project dx to constrained space (dx models an acceleration).
    ///
    /// This method must be implemented by the component, and is usually called
    /// by the generic Constraint::projectResponse() method.
    virtual void projectResponse(VecDeriv& /*dx1*/, VecDeriv& /*dx2*/) {};
#endif

    /// Project v to constrained space (v models a velocity).
    virtual void projectVelocity(const MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecDeriv& v1, DataVecDeriv& v2)
#ifdef SOFA_DEPRECATE_OLD_API
        = 0;
#else
    {
        projectVelocity(*v1.beginEdit(), *v2.beginEdit());
        v1.endEdit();
        v2.endEdit();
    }

    /// Project v to constrained space (v models a velocity).
    ///
    /// This method must be implemented by the component, and is usually called
    /// by the generic Constraint::projectVelocity() method.
    virtual void projectVelocity(VecDeriv& /*v1*/, VecDeriv& /*v2*/) {};
#endif

    /// Project x to constrained space (x models a position).
    virtual void projectPosition(const MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecCoord& x1, DataVecCoord& x2)
#ifdef SOFA_DEPRECATE_OLD_API
        = 0;
#else
    {
        projectPosition(*x1.beginEdit(), *x2.beginEdit());
        x1.endEdit();
        x2.endEdit();
    }

    /// Project x to constrained space (x models a position).
    ///
    /// This method must be implemented by the component, and is usually called
    /// by the generic Constraint::projectPosition() method.
    virtual void projectPosition(VecCoord& /*x1*/, VecCoord& /*x2*/) {};
#endif

    /// @}

    /// Project the global Mechanical Matrix to constrained space using offset parameter
    virtual void applyConstraint(const MechanicalParams* /*mparams*/ /* PARAMS FIRST */, const sofa::core::behavior::MultiMatrixAccessor* /*matrix*/)
    {

    }

    /// Project the global Mechanical Vector to constrained space using offset parameter
    virtual void applyConstraint(const MechanicalParams* /*mparams*/ /* PARAMS FIRST */, defaulttype::BaseVector* /*vector*/, const sofa::core::behavior::MultiMatrixAccessor* /*matrix*/)
    {

    }

    /// Pre-construction check method called by ObjectFactory.
    /// Check that DataTypes matches the MechanicalState.
    template<class T>
    static bool canCreate(T* obj, objectmodel::BaseContext* context, objectmodel::BaseObjectDescription* arg)
    {
        MechanicalState<DataTypes>* mstate1 = NULL;
        MechanicalState<DataTypes>* mstate2 = NULL;
        std::string object1 = arg->getAttribute("object1","@./");
        std::string object2 = arg->getAttribute("object2","@./");
        if (object1.empty()) object1 = "@./";
        if (object2.empty()) object2 = "@./";
        if (object1[0] != '@')
            object1 = BaseLink::ConvertOldPath(object1, "object1", "object1", context, false);
        if (object2[0] != '@')
            object2 = BaseLink::ConvertOldPath(object2, "object2", "object2", context, false);
        context->findLinkDest(mstate1, object1, NULL);
        context->findLinkDest(mstate2, object2, NULL);

        if (!mstate1 || !mstate2)
            return false;

        return BaseInteractionProjectiveConstraintSet::canCreate(obj, context, arg);
    }

    /// Construction method called by ObjectFactory.
    template<class T>
    static typename T::SPtr create(T* /*p0*/, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg)
    {
        typename T::SPtr obj = sofa::core::objectmodel::New<T>();

        if (context)
            context->addObject(obj);

        if (arg)
        {
            std::string object1 = arg->getAttribute("object1","");
            std::string object2 = arg->getAttribute("object2","");
            if (!object1.empty() && object1[0] != '@')
            {
                object1 = BaseLink::ConvertOldPath(object1, "object1", "object1", context, false);
                arg->setAttribute("object1", object1.c_str());
            }
            if (!object2.empty() && object2[0] != '@')
            {
                object2 = BaseLink::ConvertOldPath(object2, "object2", "object2", context, false);
                arg->setAttribute("object2", object2.c_str());
            }

            obj->parse(arg);
        }

        return obj;
    }

    virtual std::string getTemplateName() const
    {
        return templateName(this);
    }

    static std::string templateName(const PairInteractionProjectiveConstraintSet<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }

protected:
    SingleLink<PairInteractionProjectiveConstraintSet<DataTypes>, MechanicalState<DataTypes>, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> mstate1;
    SingleLink<PairInteractionProjectiveConstraintSet<DataTypes>, MechanicalState<DataTypes>, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> mstate2;
    ParticleMask *mask1;
    ParticleMask *mask2;
};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_BUILD_CORE)
extern template class SOFA_CORE_API PairInteractionProjectiveConstraintSet<defaulttype::Vec3dTypes>;
extern template class SOFA_CORE_API PairInteractionProjectiveConstraintSet<defaulttype::Vec2dTypes>;
extern template class SOFA_CORE_API PairInteractionProjectiveConstraintSet<defaulttype::Vec1dTypes>;
extern template class SOFA_CORE_API PairInteractionProjectiveConstraintSet<defaulttype::Rigid3dTypes>;
extern template class SOFA_CORE_API PairInteractionProjectiveConstraintSet<defaulttype::Rigid2dTypes>;

extern template class SOFA_CORE_API PairInteractionProjectiveConstraintSet<defaulttype::Vec3fTypes>;
extern template class SOFA_CORE_API PairInteractionProjectiveConstraintSet<defaulttype::Vec2fTypes>;
extern template class SOFA_CORE_API PairInteractionProjectiveConstraintSet<defaulttype::Vec1fTypes>;
extern template class SOFA_CORE_API PairInteractionProjectiveConstraintSet<defaulttype::Rigid3fTypes>;
extern template class SOFA_CORE_API PairInteractionProjectiveConstraintSet<defaulttype::Rigid2fTypes>;
#endif

} // namespace behavior

} // namespace core

} // namespace sofa

#endif
