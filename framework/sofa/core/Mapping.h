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
#ifndef SOFA_CORE_MAPPING_H
#define SOFA_CORE_MAPPING_H

#include <sofa/core/BaseMapping.h>
#include <sofa/core/State.h>

#include <sofa/defaulttype/Vec3Types.h>
#include <sofa/defaulttype/RigidTypes.h>

namespace sofa
{

namespace core
{

/**
*  \brief Specialized interface to convert a model state of type TIn to a model state of type TOut.
* This is basically a sofa::core::BaseMapping with given input and output types.
*
*
*/

template <class TIn, class TOut>
class Mapping : public BaseMapping
{
public:
    SOFA_ABSTRACT_CLASS(SOFA_TEMPLATE2(Mapping,TIn,TOut), BaseMapping);

    /// Input Data Type
    typedef TIn In;
    /// Output Data Type
    typedef TOut Out;

    typedef typename In::VecCoord InVecCoord;
    typedef typename In::VecDeriv InVecDeriv;
    typedef typename In::MatrixDeriv InMatrixDeriv;
    typedef Data<InVecCoord> InDataVecCoord;
    typedef Data<InVecDeriv> InDataVecDeriv;
    typedef Data<InMatrixDeriv> InDataMatrixDeriv;

    typedef typename Out::VecCoord OutVecCoord;
    typedef typename Out::VecDeriv OutVecDeriv;
    typedef typename Out::MatrixDeriv OutMatrixDeriv;
    typedef Data<OutVecCoord> OutDataVecCoord;
    typedef Data<OutVecDeriv> OutDataVecDeriv;
    typedef Data<OutMatrixDeriv> OutDataMatrixDeriv;

protected:
    /// Input Model, also called parent
    SingleLink<Mapping<In,Out>, State< In >, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> fromModel;
    /// Output Model, also called child
    SingleLink<Mapping<In,Out>, State< Out >, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> toModel;
public:

    Data<bool> f_applyRestPosition; ///< @todo document this
    Data<bool> f_checkJacobian;     ///< @todo document this
protected:
    /// Constructor, taking input and output models as parameters.
    ///
    /// Note that if you do not specify these models here, you must call
    /// setModels with non-NULL value before the intialization (i.e. before
    /// init() is called).
    Mapping(State< In >* from=NULL, State< Out >* to=NULL);
    /// Destructor
    virtual ~Mapping();
public:
    /// Specify the input and output models.
    virtual void setModels(State< In > * from, State< Out >* to);
    /// If the type is compatible set the output model and return true, otherwise do nothing and return false.
    virtual bool setTo( BaseState* to );

    /// Set the path to the objects mapped in the scene graph
    void setPathInputObject(const std::string &o) {fromModel.setPath(o);}
    void setPathOutputObject(const std::string &o) {toModel.setPath(o);}

    /// Return the pointer to the input model.
    State< In >* getFromModel();
    /// Return the pointer to the output model.
    State< Out >* getToModel();

    /// Return the pointer to the input model.
    helper::vector<BaseState*> getFrom();
    /// Return the pointer to the output model.
    helper::vector<BaseState*> getTo();

    /// Apply ///
    /// Apply the mapping to position vectors.
    ///
    /// If the Mapping can be represented as a matrix J, this method computes
    /// $ out = J in $
    virtual void apply (const MechanicalParams* mparams /* PARAMS FIRST  = MechanicalParams::defaultInstance()*/, MultiVecCoordId outPos, ConstMultiVecCoordId inPos ) ;

    /// This method must be reimplemented by all mappings.
    virtual void apply( const MechanicalParams* mparams /* PARAMS FIRST */, OutDataVecCoord& out, const InDataVecCoord& in)
#ifdef SOFA_DEPRECATE_OLD_API
        = 0;
#else
    {
        this->apply(*out.beginEdit(mparams), in.getValue(mparams));
        out.endEdit(mparams);
    }
    /// Compat Method
    /// @deprecated
    virtual void apply( OutVecCoord& /* out */, const InVecCoord& /* in */) { };
#endif //SOFA_DEPRECATE_OLD_API

    /// ApplyJ ///
    /// Apply the mapping to derived (velocity, displacement) vectors.
    /// $ out = J in $
    /// where J is the tangent operator (the linear approximation) of the mapping
    virtual void applyJ(const MechanicalParams* mparams /* PARAMS FIRST  = MechanicalParams::defaultInstance()*/, MultiVecDerivId outVel, ConstMultiVecDerivId inVel );

    /// This method must be reimplemented by all mappings.
    virtual void applyJ( const MechanicalParams* mparams /* PARAMS FIRST */, OutDataVecDeriv& out, const InDataVecDeriv& in)
#ifdef SOFA_DEPRECATE_OLD_API
        = 0;
#else
    {
        this->applyJ(*out.beginEdit(mparams), in.getValue(mparams));
        out.endEdit(mparams);
    }
    /// Compat Method
    /// @deprecated
    virtual void applyJ( OutVecDeriv& /* out */, const InVecDeriv& /* in */) { }
#endif //SOFA_DEPRECATE_OLD_API

    /// ApplyJT (Force)///
    /// Apply the reverse mapping to force vectors.
    /// $ out += J^t in $
    /// where J is the tangent operator (the linear approximation) of the mapping
    virtual void applyJT(const MechanicalParams* mparams /* PARAMS FIRST  = MechanicalParams::defaultInstance()*/, MultiVecDerivId inForce, ConstMultiVecDerivId outForce );

    /// This method must be reimplemented by all mappings.
    virtual void applyJT( const MechanicalParams* mparams /* PARAMS FIRST */, InDataVecDeriv& out, const OutDataVecDeriv& in)
#ifdef SOFA_DEPRECATE_OLD_API
        = 0;
#else
    {
        this->applyJT(*out.beginEdit(mparams), in.getValue(mparams));
        out.endEdit(mparams);
    }
    /// Compat Method
    /// @deprecated
    virtual void applyJT( InVecDeriv& /* out */, const OutVecDeriv& /* in */) { }
#endif //SOFA_DEPRECATE_OLD_API

    /// ApplyDJT (Force)///
    /// Apply the change of force due to the nonlinearity of the mapping and the last propagated displacement. Also called geometric stiffness.
    /// The default implementation does nothing, assuming a linear mapping.
    ///
    /// This method computes
    /// \f$ f_p += dJ^t f_c \f$, where \f$ f_p \f$ is the parent force and  \f$ f_c \f$ is the child force.
    /// where J is the tangent operator (the linear approximation) of the mapping
    /// The child force is accessed in the child state using mparams->readF() .  This requires that the child force vector is used by the solver to compute the force \f$ f(x,v)\f$ corresponding to the current positions and velocities, and not to store auxiliary values.
    /// The displacement is accessed in the parent state using mparams->readDx() .
    /// This method generally corresponds to a symmetric stiffness matrix, but with rotations (which are not a commutative group) it is not the case.
    /// Since some solvers (including the Conjugate Gradient) require symmetric matrices, a flag is set in the MechanicalParams to say if symmetric matrices are required. If so, non-symmetric geometric stiffness should not be applied.
    virtual void applyDJT(const MechanicalParams* /*mparams = MechanicalParams::defaultInstance()*/ , MultiVecDerivId /*parentForce*/, ConstMultiVecDerivId  /*childForce*/ );

    /// ApplyJT (Constraint)///
    virtual void applyJT(const ConstraintParams* cparams /* PARAMS FIRST  = ConstraintParams::defaultInstance()*/, MultiMatrixDerivId inConst, ConstMultiMatrixDerivId outConst );

    /// This method must be reimplemented by all mappings if they need to support constraints.
    virtual void applyJT( const ConstraintParams* mparams /* PARAMS FIRST */, InDataMatrixDeriv& out, const OutDataMatrixDeriv& in)
#ifdef SOFA_DEPRECATE_OLD_API
    {
        serr << "This mapping does not support constraints because Mapping::applyJT( const ConstraintParams* , InDataMatrixDeriv&, const OutDataMatrixDeriv&) is not overloaded." << sendl;
    }
#else
    {
        this->applyJT(*out.beginEdit(mparams), in.getValue(mparams));
        out.endEdit(mparams);
    }
    /// Compat Method
    /// @deprecated
    virtual void applyJT( InMatrixDeriv& /*out*/, const OutMatrixDeriv& /*in*/ )
    {
        serr << "This mapping does not support constraints because Mapping::applyJT( InMatrixDeriv&, const OutMatrixDeriv& ) is not overloaded. " << sendl;
    }
#endif //SOFA_DEPRECATE_OLD_API

    /// computeAccFromMapping
    /// Compute the acceleration of the child, based on the acceleration and the velocity of the parent.
    /// Let \f$ v_c = J v_p \f$ be the velocity of the child given the velocity of the parent, then the acceleration is \f$ a_c = J a_p + dJ v_p \f$.
    /// The second term is null in linear mappings, otherwise it encodes the acceleration due to the change of mapping at constant parent velocity.
    /// For instance, in a rigid mapping with angular velocity\f$ w \f$,  the second term is $ w^(w^rel_pos) $
    virtual void computeAccFromMapping(const MechanicalParams* mparams /* PARAMS FIRST  = MechanicalParams::defaultInstance()*/, MultiVecDerivId outAcc, ConstMultiVecDerivId inVel, ConstMultiVecDerivId inAcc );

    /// This method must be reimplemented by all mappings if they need to support composite accelerations
    virtual void computeAccFromMapping(const MechanicalParams* mparams /* PARAMS FIRST */, OutDataVecDeriv& accOut, const InDataVecDeriv& vIn, const InDataVecDeriv& accIn)
#ifdef SOFA_DEPRECATE_OLD_API
    {
    }
#else
    {
        this->computeAccFromMapping(*accOut.beginEdit(mparams), vIn.getValue(mparams), accIn.getValue(mparams));
        accOut.endEdit(mparams);
    }
    /// Compat Method
    /// @deprecated
    virtual void computeAccFromMapping( OutVecDeriv& /*acc_out*/, const InVecDeriv& /*v_in*/, const InVecDeriv& /*acc_in*/)
    {}
#endif //SOFA_DEPRECATE_OLD_API

    /// Propagate positions and velocities to the output
    virtual void init();

    ///<TO REMOVE>  FF:why would we remove this, is there any alternative function ?
    // Useful ?
    /// Get the source (upper) model.
    virtual helper::vector<behavior::BaseMechanicalState*> getMechFrom();

    /// Get the destination (lower, mapped) model.
    virtual helper::vector<behavior::BaseMechanicalState*> getMechTo();

    //Create a matrix for mapped mechanical objects
    //If the two mechanical objects is identical, create a new stiffness matrix for this mapped objects
    //If the two mechanical objects is different, create a new interaction matrix
    virtual sofa::defaulttype::BaseMatrix* createMappedMatrix(const behavior::BaseMechanicalState* state1, const behavior::BaseMechanicalState* state2, func_createMappedMatrix);

    ///<TO REMOVE>
    /// Apply the mapping to position and velocity vectors.
    ///
    /// This method call the internal apply(Out::VecCoord&,const In::VecCoord&)
    /// and applyJ(Out::VecDeriv&,const In::VecDeriv&) methods.
    //virtual void updateMapping();

    /// Disable the mapping to get the original coordinates of the mapped model.
    ///
    /// It is for instance used in RigidMapping to get the local coordinates of the object.
    virtual void disable();

    /// Pre-construction check method called by ObjectFactory.
    ///
    /// This implementation read the input and output attributes and check
    /// if they are compatible with the input and output model types of this
    /// mapping.
    template<class T>
    static bool canCreate(T*& obj, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg)
    {
        State<In>* stin = NULL;
        State<Out>* stout = NULL;

        std::string inPath, outPath;

        if (arg->getAttribute("input"))
            inPath = arg->getAttribute("input");
#ifndef SOFA_DEPRECATE_OLD_API
        else if (arg->getAttribute("object1"))
            inPath = BaseLink::ConvertOldPath(arg->getAttribute("object1"), "object1", "input", context, false);
#endif
        else
            inPath = "@../";

        context->findLinkDest(stin, inPath, NULL);

        if (arg->getAttribute("output"))
            outPath = arg->getAttribute("output");
#ifndef SOFA_DEPRECATE_OLD_API
        else if (arg->getAttribute("object2"))
            outPath = BaseLink::ConvertOldPath(arg->getAttribute("object2"), "object2", "output", context, false);
#endif
        else
            outPath = "@./";

        context->findLinkDest(stout, outPath, NULL);

        if (stin == NULL)
        {
//            This warning seems irrelevant, as it is raised multiple times while the creation works fine (Francois Faure, Feb. 2012)
//            context->serr << "Cannot create "<<className(obj)<<" as input model "<< inPath << " is missing or invalid." << context->sendl;
            return false;
        }

        if (stout == NULL)
        {
//            This warning seems irrelevant, as it is raised multiple times while the creation works fine (Francois Faure, Feb. 2012)
//            context->serr << "Cannot create "<<className(obj)<<" as output model "<< outPath << " is missing or invalid." << context->sendl;
            return false;
        }

        if (static_cast<BaseObject*>(stin) == static_cast<BaseObject*>(stout))
        {
            // we should refuse to create mappings with the same input and output model, which may happen if a State object is missing in the child node
            context->serr << "Creation of " << className(obj) << " mapping failed because the same object \"" << stin->getName() << "\" is linked as both input and output." << context->sendl;
            context->serr << "  Maybe a MechanicalObject should be added before this mapping." << context->sendl;
            return false;
        }

        return BaseMapping::canCreate(obj, context, arg);
    }

    /// Construction method called by ObjectFactory.
    ///
    /// This implementation read the input and output attributes to
    /// find the input and output models of this mapping.
    template<class T>
    static typename T::SPtr create(T*, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg)
    {
        typename T::SPtr obj = sofa::core::objectmodel::New<T>();

        if (context)
            context->addObject(obj);

        if (arg)
        {
            std::string inPath, outPath;
            if (arg->getAttribute("input"))
                inPath = arg->getAttribute("input");
#ifndef SOFA_DEPRECATE_OLD_API
            else if (arg->getAttribute("object1"))
                inPath = BaseLink::ConvertOldPath(arg->getAttribute("object1"), "object1", "input", obj.get());
#endif
            else
                inPath = "@../";

            if (arg->getAttribute("output"))
                outPath = arg->getAttribute("output");
#ifndef SOFA_DEPRECATE_OLD_API
            else if (arg->getAttribute("object2"))
                outPath = BaseLink::ConvertOldPath(arg->getAttribute("object2"), "object2", "output", obj.get());
#endif
            else
                outPath = "@./";

            obj->fromModel.setPath( inPath );
            obj->toModel.setPath( outPath );

            obj->parse(arg);
        }

        return obj;
    }

    virtual std::string getTemplateName() const
    {
        return templateName(this);
    }

    static std::string templateName(const Mapping<TIn, TOut>* = NULL);


    template<class T>
    static std::string shortName(const T* ptr = NULL, objectmodel::BaseObjectDescription* arg = NULL)
    {
        std::string name = Inherit1::shortName(ptr, arg);
        sofa::helper::replaceAll(name, "Mapping", "Map");
        return name;
    }

protected:

    void matrixApplyJ( OutVecDeriv& /* out */, const InVecDeriv& /* in */, const sofa::defaulttype::BaseMatrix* /* J */);
    void matrixApplyJT( InVecDeriv& /* out */, const OutVecDeriv& /* in */, const sofa::defaulttype::BaseMatrix* /* J */);
    void matrixApplyJT( InMatrixDeriv& /* out */, const OutMatrixDeriv& /* in */, const sofa::defaulttype::BaseMatrix* /* J */);
    bool checkApplyJ( OutVecDeriv& /* out */, const InVecDeriv& /* in */, const sofa::defaulttype::BaseMatrix* /* J */);
    bool checkApplyJT( InVecDeriv& /* out */, const OutVecDeriv& /* in */, const sofa::defaulttype::BaseMatrix* /* J */);
    bool checkApplyJT( InMatrixDeriv& /* out */, const OutMatrixDeriv& /* in */, const sofa::defaulttype::BaseMatrix* /* J */);

};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_BUILD_CORE)

using namespace sofa::defaulttype;

#ifndef SOFA_FLOAT
extern template class SOFA_CORE_API Mapping< Vec3dTypes, Vec3dTypes >;
extern template class SOFA_CORE_API Mapping< Rigid3dTypes, Vec3dTypes >;
extern template class SOFA_CORE_API Mapping< Vec3dTypes, ExtVec3dTypes >;
extern template class SOFA_CORE_API Mapping< Vec3dTypes, Vec1dTypes >;
extern template class SOFA_CORE_API Mapping< Rigid2dTypes, Vec2dTypes >;
#endif

#ifndef SOFA_DOUBLE
extern template class SOFA_CORE_API Mapping< Vec3fTypes, Vec3fTypes >;
extern template class SOFA_CORE_API Mapping< Rigid3fTypes, Vec3fTypes >;
extern template class SOFA_CORE_API Mapping< Vec3fTypes, ExtVec3fTypes >;
extern template class SOFA_CORE_API Mapping< Vec3fTypes, Vec1fTypes >;
extern template class SOFA_CORE_API Mapping< Rigid2fTypes, Vec2fTypes >;
#endif

#ifndef SOFA_FLOAT
#ifndef SOFA_DOUBLE
extern template class SOFA_CORE_API Mapping< Vec3dTypes, Vec3fTypes >;
extern template class SOFA_CORE_API Mapping< Rigid3dTypes, Vec3fTypes >;
extern template class SOFA_CORE_API Mapping< Vec3dTypes, ExtVec3fTypes >;
extern template class SOFA_CORE_API Mapping< Vec3dTypes, Vec1fTypes >;
extern template class SOFA_CORE_API Mapping< Rigid2dTypes, Vec2fTypes >;

extern template class SOFA_CORE_API Mapping< Vec3fTypes, Vec3dTypes >;
extern template class SOFA_CORE_API Mapping< Rigid3fTypes, Vec3dTypes >;
extern template class SOFA_CORE_API Mapping< Vec3fTypes, ExtVec3dTypes >;
extern template class SOFA_CORE_API Mapping< Vec3fTypes, Vec1dTypes >;
extern template class SOFA_CORE_API Mapping< Rigid2fTypes, Vec2dTypes >;
#endif
#endif
#endif

} // namespace core

} // namespace sofa

#endif // SOFA_CORE_MAPPING_H
