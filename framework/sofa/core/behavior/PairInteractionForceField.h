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
#ifndef SOFA_CORE_BEHAVIOR_PAIRINTERACTIONFORCEFIELD_H
#define SOFA_CORE_BEHAVIOR_PAIRINTERACTIONFORCEFIELD_H

#include <sofa/core/core.h>
#include <sofa/core/MultiVecId.h>
#include <sofa/core/MechanicalParams.h>
#include <sofa/core/behavior/BaseInteractionForceField.h>
#include <sofa/core/behavior/MechanicalState.h>

#include <sofa/core/objectmodel/Data.h>
#include <sofa/defaulttype/Vec.h>

#include <sofa/defaulttype/Vec3Types.h>
#include <sofa/defaulttype/RigidTypes.h>

namespace sofa
{

namespace core
{

namespace behavior
{

/**
 *  \brief Component computing forces between a pair of simulated body.
 *
 *  This class define the abstract API common to interaction force fields
 *  between a pair of bodies using a given type of DOFs.
 */
template<class TDataTypes>
class PairInteractionForceField : public BaseInteractionForceField
{
public:
    SOFA_ABSTRACT_CLASS(SOFA_TEMPLATE(PairInteractionForceField,TDataTypes), BaseInteractionForceField);

    typedef TDataTypes DataTypes;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef helper::ParticleMask ParticleMask;

    typedef core::objectmodel::Data<VecCoord>    DataVecCoord;
    typedef core::objectmodel::Data<VecDeriv>    DataVecDeriv;

protected:
    PairInteractionForceField(MechanicalState<DataTypes> *mm1 = NULL, MechanicalState<DataTypes> *mm2 = NULL);

    virtual ~PairInteractionForceField();
public:
    virtual void init();

    /// Retrieve the associated MechanicalState
    MechanicalState<DataTypes>* getMState1() { return mstate1.get(); }
    BaseMechanicalState* getMechModel1() { return mstate1.get(); }
    /// Retrieve the associated MechanicalState
    MechanicalState<DataTypes>* getMState2() { return mstate2.get(); }
    BaseMechanicalState* getMechModel2() { return mstate2.get(); }

    /// Set the Object1 path
    void setPathObject1(const std::string & path) { mstate1.setPath(path); }
    /// Set the Object2 path
    void setPathObject2(const std::string & path) { mstate2.setPath(path); }
    /// Retrieve the Object1 path
    std::string getPathObject1() const { return mstate1.getPath(); }
    /// Retrieve the Object2 path
    std::string getPathObject2() const { return mstate2.getPath(); }


    /// @name Vector operations
    /// @{

    /// Given the current position and velocity states, update the current force
    /// vector by computing and adding the forces associated with this
    /// ForceField.
    ///
    /// If the ForceField can be represented as a matrix, this method computes
    /// $ f += B v + K x $
    ///
    /// This method retrieves the force, x and v vector from the two MechanicalState
    /// and call the internal addForce(VecDeriv&,VecDeriv&,const VecCoord&,const VecCoord&,const VecDeriv&,const VecDeriv&)
    /// method implemented by the component.
    virtual void addForce(const MechanicalParams* mparams /* PARAMS FIRST */, MultiVecDerivId fId );

    /// Compute the force derivative given a small displacement from the
    /// position and velocity used in the previous call to addForce().
    ///
    /// The derivative should be directly derived from the computations
    /// done by addForce. Any forces neglected in addDForce will be integrated
    /// explicitly (i.e. using its value at the beginning of the timestep).
    ///
    /// If the ForceField can be represented as a matrix, this method computes
    /// $ df += kFactor K dx + bFactor B dx $
    ///
    /// This method retrieves the force and dx vector from the two MechanicalState
    /// and call the internal addDForce(VecDeriv&,VecDeriv&,const VecDeriv&,const VecDeriv&,double,double)
    /// method implemented by the component.
    virtual void addDForce(const MechanicalParams* mparams /* PARAMS FIRST */, MultiVecDerivId dfId );


    /// Get the potential energy associated to this ForceField.
    ///
    /// Used to extimate the total energy of the system by some
    /// post-stabilization techniques.
    ///
    /// This method retrieves the x vector from the MechanicalState and call
    /// the internal getPotentialEnergy(const VecCoord&,const VecCoord&) method implemented by
    /// the component.
    virtual double getPotentialEnergy(const MechanicalParams* mparams) const;

    /// Given the current position and velocity states, update the current force
    /// vector by computing and adding the forces associated with this
    /// ForceField.
    ///
    /// If the ForceField can be represented as a matrix, this method computes
    /// $ f += B v + K x $
    ///
    /// This method must be implemented by the component, and is usually called
    /// by the generic ForceField::addForce() method.

    virtual void addForce(const MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& f1, DataVecDeriv& f2, const DataVecCoord& x1, const DataVecCoord& x2, const DataVecDeriv& v1, const DataVecDeriv& v2 )=0;

    /// @deprecated
    //virtual void addForce(VecDeriv& f1, VecDeriv& f2, const VecCoord& x1, const VecCoord& x2, const VecDeriv& v1, const VecDeriv& v2);


    /// Compute the force derivative given a small displacement from the
    /// position and velocity used in the previous call to addForce().
    ///
    /// The derivative should be directly derived from the computations
    /// done by addForce. Any forces neglected in addDForce will be integrated
    /// explicitly (i.e. using its value at the beginning of the timestep).
    ///
    /// If the ForceField can be represented as a matrix, this method computes
    /// $ df += kFactor K dx + bFactor B dx $
    ///
    /// This method must be implemented by the component, and is usually called
    /// by the generic PairInteractionForceField::addDForce() method.
    ///
    /// To support old components that implement the deprecated addForce method
    /// without scalar coefficients, it defaults to using a temporaty vector to
    /// compute $ K dx $ and then manually scaling all values by kFactor.

    virtual void addDForce(const MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& df1, DataVecDeriv& df2, const DataVecDeriv& dx1, const DataVecDeriv& dx2)=0;

    /// @deprecated
    //virtual void addDForce(VecDeriv& df1, VecDeriv& df2, const VecDeriv& dx1, const VecDeriv& dx2, double kFactor, double /*bFactor*/);

    /// Compute the force derivative given a small displacement from the
    /// position and velocity used in the previous call to addForce().
    ///
    /// The derivative should be directly derived from the computations
    /// done by addForce. Any forces neglected in addDForce will be integrated
    /// explicitly (i.e. using its value at the beginning of the timestep).
    ///
    /// If the ForceField can be represented as a matrix, this method computes
    /// $ df += K dx $
    ///
    /// This method must be implemented by the component, and is usually called
    /// by the generic PairInteractionForceField::addDForce() method.
    ///
    /// @deprecated to more efficiently accumulate contributions from all terms
    ///   of the system equation, a new addDForce method allowing to pass two
    ///   coefficients for the stiffness and damping terms should now be used.
    /// @deprecated
    //virtual void addDForce(VecDeriv& df1, VecDeriv& df2, const VecDeriv& dx1, const VecDeriv& dx2);

    /// Get the potential energy associated to this ForceField.
    ///
    /// Used to extimate the total energy of the system by some
    /// post-stabilization techniques.
    ///
    /// This method must be implemented by the component, and is usually called
    /// by the generic ForceField::getPotentialEnergy() method.

    virtual double getPotentialEnergy(const MechanicalParams* mparams /* PARAMS FIRST */, const DataVecCoord& x1, const DataVecCoord& x2) const=0;

    /// @deprecated
    //virtual double getPotentialEnergy(const VecCoord& x1, const VecCoord& x2) const;

    /// @}

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

        return BaseInteractionForceField::canCreate(obj, context, arg);
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

    static std::string templateName(const PairInteractionForceField<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }

    template<class T>
    static std::string shortName(const T* ptr = NULL, objectmodel::BaseObjectDescription* arg = NULL)
    {
        std::string name = Inherit1::shortName(ptr, arg);
        sofa::helper::replaceAll(name, "InteractionForceField", "IFF");
        sofa::helper::replaceAll(name, "ForceField", "FF");
        return name;
    }

protected:
    SingleLink<PairInteractionForceField<DataTypes>, MechanicalState<DataTypes>, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> mstate1;
    SingleLink<PairInteractionForceField<DataTypes>, MechanicalState<DataTypes>, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> mstate2;

    ParticleMask *mask1;
    ParticleMask *mask2;
};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_BUILD_CORE)
#ifndef SOFA_FLOAT
extern template class SOFA_CORE_API PairInteractionForceField<defaulttype::Vec3dTypes>;
extern template class SOFA_CORE_API PairInteractionForceField<defaulttype::Vec2dTypes>;
extern template class SOFA_CORE_API PairInteractionForceField<defaulttype::Vec1dTypes>;
extern template class SOFA_CORE_API PairInteractionForceField<defaulttype::Rigid3dTypes>;
extern template class SOFA_CORE_API PairInteractionForceField<defaulttype::Rigid2dTypes>;
#endif

#ifndef SOFA_DOUBLE
extern template class SOFA_CORE_API PairInteractionForceField<defaulttype::Vec3fTypes>;
extern template class SOFA_CORE_API PairInteractionForceField<defaulttype::Vec2fTypes>;
extern template class SOFA_CORE_API PairInteractionForceField<defaulttype::Vec1fTypes>;
extern template class SOFA_CORE_API PairInteractionForceField<defaulttype::Rigid3fTypes>;
extern template class SOFA_CORE_API PairInteractionForceField<defaulttype::Rigid2fTypes>;
#endif
#endif

} // namespace behavior

} // namespace core

} // namespace sofa

#endif
