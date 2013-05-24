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
#ifndef SOFA_COMPONENT_COLLISION_SPHEREMODEL_H
#define SOFA_COMPONENT_COLLISION_SPHEREMODEL_H

#include <sofa/core/CollisionModel.h>
#include <sofa/component/container/MechanicalObject.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/objectmodel/DataFileName.h>
#include <sofa/component/component.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/helper/accessor.h>
//#include <sofa/component/collision/RigidContactMapper.h>

namespace sofa
{

namespace component
{

namespace collision
{

using namespace sofa::defaulttype;

template<class T>
struct _SphereDataTypes{
    typedef T DataTypes;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;

    static const Coord & position(const typename DataTypes::Coord & c){return c;}
    static const Deriv & lvelocity(const typename DataTypes::Deriv & d){return d;}
};

template <class MyReal>
struct _SphereDataTypes<StdRigidTypes<3,MyReal> >{
    typedef defaulttype::StdRigidTypes<3,MyReal> DataTypes;
    typedef typename StdRigidTypes<3,MyReal>::Coord::Pos Coord;

    static const Coord & position(const typename DataTypes::Coord & c){return c.getCenter();}
    static const Coord & lvelocity(const typename DataTypes::Deriv & d){return d.getLinear();}
};

//template <>
//struct _SphereDataTypes<defaulttype::Vec3dTypes>{
//    typedef defaulttype::Vec3dTypes DataTypes;
//    typedef DataTypes::Coord Coord;
//    typedef DataTypes::Deriv Deriv;

//    static const Coord & position(const DataTypes::Coord & c){return c;}
//    static const Deriv & lvelocity(const DataTypes::Deriv & d){return d;}
//};

//template <>
//struct _SphereDataTypes<defaulttype::Vec3fTypes>{
//    typedef defaulttype::Vec3fTypes DataTypes;
//    typedef DataTypes::Coord Coord;
//    typedef DataTypes::Deriv Deriv;

//    static const Coord & position(const DataTypes::Coord & c){return c;}
//    static const Deriv & lvelocity(const DataTypes::Deriv & d){return d;}
//};


template<class DataTypes>
class TSphereModel;

template <class TDataTypes>
class TSphere;

template<class TDataTypes>
class TSphere : public core::TCollisionElementIterator< TSphereModel<TDataTypes> >
{
public:
    typedef TDataTypes DataTypes;
    typedef typename DataTypes::Real   Real;
    typedef typename _SphereDataTypes<TDataTypes>::Coord Coord;

    typedef TSphereModel<DataTypes> ParentModel;

    TSphere(ParentModel* model, int index);

    explicit TSphere(core::CollisionElementIterator& i);

    const Coord& center() const;
    const typename TDataTypes::Coord & rigidCenter() const;
    const Coord& p() const;
    const Coord& pFree() const;
    const Coord& v() const;

    /// Return true if the element stores a free position vector
    bool hasFreePosition() const;

    Real r() const;
};

template< class TDataTypes>
class TSphereModel : public core::CollisionModel
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(TSphereModel, TDataTypes), core::CollisionModel);

    typedef TDataTypes DataTypes;
    typedef DataTypes InDataTypes;

    typedef typename _SphereDataTypes<DataTypes>::Coord Coord;
    //typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::VecReal VecReal;
    typedef TSphere<DataTypes> Element;
    friend class TSphere<DataTypes>;
protected:
    TSphereModel();

    TSphereModel(core::behavior::MechanicalState<TDataTypes>* _mstate );
public:
    virtual void init();

    // -- CollisionModel interface

    virtual void resize(int size);

    virtual void computeBoundingTree(int maxDepth=0);

    virtual void computeContinuousBoundingTree(double dt, int maxDepth=0);

    void draw(const core::visual::VisualParams*,int index);

    void draw(const core::visual::VisualParams* vparams);


    core::behavior::MechanicalState<DataTypes>* getMechanicalState() { return mstate; }

    const VecReal& getR() const { return this->radius.getValue(); }

    Real getRadius(const int i) const;

    const Coord & velocity(int index)const;

    /// Pre-construction check method called by ObjectFactory.
    /// Check that DataTypes matches the MechanicalState.
    template<class T>
    static bool canCreate(T*& obj, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg)
    {
        if (dynamic_cast<core::behavior::MechanicalState<TDataTypes>*>(context->getMechanicalState()) == NULL && context->getMechanicalState() != NULL)
            return false;

        return BaseObject::canCreate(obj, context, arg);
    }

    template<class T>
    static typename T::SPtr create(T*, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg)
    {
        typename T::SPtr obj;
        core::behavior::MechanicalState<TDataTypes>* _mstate = NULL;

        if( context)
        {
            _mstate = dynamic_cast<core::behavior::MechanicalState<TDataTypes>*>(context->getMechanicalState());
            if (_mstate)
                obj = sofa::core::objectmodel::New<T>(_mstate);
            else
                obj = sofa::core::objectmodel::New<T>();

            context->addObject(obj);
        }

        if (arg) obj->parse(arg);

        return obj;
    }


    virtual std::string getTemplateName() const
    {
        return templateName(this);
    }

    static std::string templateName(const TSphereModel<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }

    Data< VecReal > radius;
    Data< SReal > defaultRadius;

protected:
    core::behavior::MechanicalState<DataTypes>* mstate;
};

////////////////////////////////////////////////////////////////////////////
//Specialisation for rigid spheres//////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

//template< class TReal>
//class TSphere<StdRigidTypes<3,TReal> > : public core::TCollisionElementIterator< TSphereModel<StdRigidTypes<3,TReal> > >
//{
//public:
//    typedef StdRigidTypes<3,TReal>  DataTypes;
//    typedef typename DataTypes::Real   Real;
//    typedef typename DataTypes::Coord::Pos Coord;

//    typedef TSphereModel<DataTypes> ParentModel;

//    TSphere(ParentModel* model, int index);

//    explicit TSphere(core::CollisionElementIterator& i);

//    const Coord& center() const;
//    const Coord& p() const;
//    const Coord& pFree() const;
//    const Coord& v() const;

//    /// Return true if the element stores a free position vector
//    bool hasFreePosition() const;

//    Real r() const;
//};


//template< class TReal>
//class TSphereModel<StdRigidTypes<3,TReal> > : public core::CollisionModel
//{
//public:
//    SOFA_CLASS(SOFA_TEMPLATE(TSphereModel, SOFA_TEMPLATE2(StdRigidTypes,3,TReal) ), core::CollisionModel);

//    typedef StdRigidTypes<3,TReal> DataTypes;
//    typedef DataTypes InDataTypes;
//    typedef typename DataTypes::Coord::Pos Coord;
//    typedef TReal Real;
//    typedef typename DataTypes::VecReal VecReal;
//    typedef TSphere<DataTypes> Element;
//    friend class TSphere<DataTypes>;
//protected:
//    TSphereModel();

//    TSphereModel(core::behavior::MechanicalState<DataTypes>* _mstate );
//public:
//    virtual void init();

//    // -- CollisionModel interface

//    virtual void resize(int size);

//    virtual void computeBoundingTree(int maxDepth=0);

//    virtual void computeContinuousBoundingTree(double dt, int maxDepth=0);

//    void draw(const core::visual::VisualParams*,int index);

//    void draw(const core::visual::VisualParams* vparams);


//    core::behavior::MechanicalState<DataTypes>* getMechanicalState() { return mstate; }

//    const VecReal& getR() const { return this->radius.getValue(); }

//    Real getRadius(const int i) const;

//    const Coord & velocity(int index)const;

//    /// Pre-construction check method called by ObjectFactory.
//    /// Check that DataTypes matches the MechanicalState.
//    template<class T>
//    static bool canCreate(T*& obj, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg)
//    {
//        if (dynamic_cast<core::behavior::MechanicalState<DataTypes>*>(context->getMechanicalState()) == NULL && context->getMechanicalState() != NULL)
//            return false;

//        return BaseObject::canCreate(obj, context, arg);
//    }

//    template<class T>
//    static typename T::SPtr create(T*, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg)
//    {
//        typename T::SPtr obj;
//        core::behavior::MechanicalState<DataTypes>* _mstate = NULL;

//        if( context)
//        {
//            _mstate = dynamic_cast<core::behavior::MechanicalState<DataTypes>*>(context->getMechanicalState());
//            if (_mstate)
//                obj = sofa::core::objectmodel::New<T>(_mstate);
//            else
//                obj = sofa::core::objectmodel::New<T>();

//            context->addObject(obj);
//        }

//        if (arg) obj->parse(arg);

//        return obj;
//    }


//    virtual std::string getTemplateName() const
//    {
//        return templateName(this);
//    }

//    static std::string templateName(const TSphereModel<DataTypes>* = NULL)
//    {
//        return DataTypes::Name();
//    }

//    Data< VecReal > radius;
//    Data< SReal > defaultRadius;

//protected:
//    core::behavior::MechanicalState<DataTypes>* mstate;
//};


template<class DataTypes>
inline TSphere<DataTypes>::TSphere(ParentModel* model, int index)
    : core::TCollisionElementIterator<ParentModel>(model, index)
{}

template<class DataTypes>
inline TSphere<DataTypes>::TSphere(core::CollisionElementIterator& i)
    : core::TCollisionElementIterator<ParentModel>(static_cast<ParentModel*>(i.getCollisionModel()), i.getIndex())
{
}

template<class DataTypes>
inline const typename TSphere<DataTypes>::Coord& TSphere<DataTypes>::center() const { return _SphereDataTypes<DataTypes>::position((*this->model->mstate->getX())[this->index]); }

template<class DataTypes>
inline const typename DataTypes::Coord & TSphere<DataTypes>::rigidCenter() const { return (*this->model->mstate->getX())[this->index];}

template<class DataTypes>
inline const typename TSphere<DataTypes>::Coord& TSphere<DataTypes>::p() const { return _SphereDataTypes<DataTypes>::position((*this->model->mstate->getX())[this->index]); }

template<class DataTypes>
inline const typename TSphere<DataTypes>::Coord& TSphere<DataTypes>::pFree() const { return (*this->model->mstate->read(core::ConstVecCoordId::freePosition())).getValue()[this->index]; }

template<class DataTypes>
inline const typename TSphereModel<DataTypes>::Coord& TSphereModel<DataTypes>::velocity(int index) const { return _SphereDataTypes<DataTypes>::lvelocity((*mstate->getV())[index]);}

template<class DataTypes>
inline const typename TSphere<DataTypes>::Coord& TSphere<DataTypes>::v() const { return _SphereDataTypes<DataTypes>::lvelocity((*this->model->mstate->getV())[this->index]); }

template<class DataTypes>
inline typename DataTypes::Real TSphere<DataTypes>::r() const { return (Real) this->model->getRadius((unsigned)this->index); }

template<class DataTypes>
inline bool TSphere<DataTypes>::hasFreePosition() const { return this->model->mstate->read(core::ConstVecCoordId::freePosition())->isSet(); }


typedef TSphereModel<Vec3Types> SphereModel;
typedef TSphere<Vec3Types> Sphere;

typedef TSphereModel<Rigid3Types> RigidSphereModel;
typedef TSphere<Rigid3Types> RigidSphere;

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_BUILD_BASE_COLLISION)
#ifndef SOFA_FLOAT
extern template class SOFA_BASE_COLLISION_API TSphereModel<defaulttype::Vec3dTypes>;
extern template class SOFA_BASE_COLLISION_API TSphereModel<defaulttype::Rigid3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_BASE_COLLISION_API TSphereModel<defaulttype::Vec3fTypes>;
extern template class SOFA_BASE_COLLISION_API TSphereModel<defaulttype::Rigid3fTypes>;
#endif
#endif

} // namespace collision

} // namespace component

} // namespace sofa

#endif
