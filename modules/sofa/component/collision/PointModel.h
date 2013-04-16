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
#ifndef SOFA_COMPONENT_COLLISION_POINTMODEL_H
#define SOFA_COMPONENT_COLLISION_POINTMODEL_H

#include <sofa/core/CollisionModel.h>
#include <sofa/component/collision/LocalMinDistanceFilter.h>
#include <sofa/component/container/MechanicalObject.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/defaulttype/Vec3Types.h>
#include <vector>

namespace sofa
{

namespace component
{

namespace collision
{

using namespace sofa::defaulttype;

template<class DataTypes>
class TPointModel;

class PointLocalMinDistanceFilter;

template<class TDataTypes>
class TPoint : public core::TCollisionElementIterator<TPointModel<TDataTypes> >
{
public:
    typedef TDataTypes DataTypes;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef TPointModel<DataTypes> ParentModel;

    TPoint(ParentModel* model, int index);
    TPoint() {}

    explicit TPoint(core::CollisionElementIterator& i);

    const Coord& p() const;
    const Coord& pFree() const;
    const Deriv& v() const;
    Deriv n() const;

    /// Return true if the element stores a free position vector
    bool hasFreePosition() const;

    bool testLMD(const Vector3 &, double &, double &);

    bool activated(core::CollisionModel *cm = 0) const;
};

class PointActiver
{
public:
    PointActiver() {}
    virtual ~PointActiver() {}
    virtual bool activePoint(int /*index*/, core::CollisionModel * /*cm*/ = 0) {return true;}
};

template<class TDataTypes>
class SOFA_MESH_COLLISION_API TPointModel : public core::CollisionModel
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(TPointModel, TDataTypes), core::CollisionModel);

//    typedef Vec3Types InDataTypes;
//    typedef Vec3Types DataTypes;
    typedef TDataTypes DataTypes;
    typedef DataTypes InDataTypes;
    typedef TPointModel<DataTypes> ParentModel;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef TPoint<DataTypes> Element;
    typedef helper::vector<unsigned int> VecIndex;

    friend class TPoint<DataTypes>;
protected:
    TPointModel();
public:
    virtual void init();

    // -- CollisionModel interface

    virtual void resize(int size);

    virtual void computeBoundingTree(int maxDepth=0);

    virtual void computeContinuousBoundingTree(double dt, int maxDepth=0);

    void draw(const core::visual::VisualParams*,int index);

    void draw(const core::visual::VisualParams* vparams);

    virtual bool canCollideWithElement(int index, CollisionModel* model2, int index2);

    core::behavior::MechanicalState<DataTypes>* getMechanicalState() { return mstate; }

    //virtual const char* getTypeName() const { return "Point"; }

    PointLocalMinDistanceFilter *getFilter() const;

    //template< class TFilter >
    //TFilter *getFilter() const
    //{
    //	if (m_lmdFilter != 0)
    //		return m_lmdFilter;
    //	else
    //		return &m_emptyFilter;
    //}

    void setFilter(PointLocalMinDistanceFilter * /*lmdFilter*/);

    const Deriv& velocity(int index) const;

    Data<bool> bothSide; // to activate collision on both side of the point model (when surface normals are defined on these points)

    /// Pre-construction check method called by ObjectFactory.
    /// Check that DataTypes matches the MechanicalState.
    template<class T>
    static bool canCreate(T*& obj, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg)
    {
        if (dynamic_cast<core::behavior::MechanicalState<DataTypes>*>(context->getMechanicalState()) == NULL)
            return false;
        return BaseObject::canCreate(obj, context, arg);
    }

    virtual std::string getTemplateName() const
    {
        return templateName(this);
    }

    static std::string templateName(const TPointModel<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }

protected:

    core::behavior::MechanicalState<DataTypes>* mstate;

    Data<bool> computeNormals;

    Data<std::string> PointActiverPath;

    VecDeriv normals;

    PointLocalMinDistanceFilter *m_lmdFilter;
    EmptyFilter m_emptyFilter;

    Data<bool> m_displayFreePosition;

    void updateNormals();

    PointActiver *myActiver;
};

template<class DataTypes>
inline TPoint<DataTypes>::TPoint(ParentModel* model, int index)
    : core::TCollisionElementIterator<ParentModel>(model, index)
{

}

template<class DataTypes>
inline TPoint<DataTypes>::TPoint(core::CollisionElementIterator& i)
    : core::TCollisionElementIterator<ParentModel>(static_cast<ParentModel*>(i.getCollisionModel()), i.getIndex())
{

}

template<class DataTypes>
inline const typename DataTypes::Coord& TPoint<DataTypes>::p() const { return (*this->model->mstate->getX())[this->index]; }

template<class DataTypes>
inline const typename DataTypes::Coord& TPoint<DataTypes>::pFree() const
{
    if (hasFreePosition())
        return this->model->mstate->read(core::ConstVecCoordId::freePosition())->getValue()[this->index];
    else
        return p();
}

template<class DataTypes>
inline const typename DataTypes::Deriv& TPoint<DataTypes>::v() const { return (*this->model->mstate->getV())[this->index]; }

template<class DataTypes>
inline const typename DataTypes::Deriv& TPointModel<DataTypes>::velocity(int index) const { return (*mstate->getV())[index]; }

template<class DataTypes>
inline typename DataTypes::Deriv TPoint<DataTypes>::n() const { return ((unsigned)this->index<this->model->normals.size()) ? this->model->normals[this->index] : Deriv(); }

template<class DataTypes>
inline bool TPoint<DataTypes>::hasFreePosition() const { return this->model->mstate->read(core::ConstVecCoordId::freePosition())->isSet(); }

template<class DataTypes>
inline bool TPoint<DataTypes>::activated(core::CollisionModel *cm) const
{
    return this->model->myActiver->activePoint(this->index, cm);
}

typedef TPointModel<Vec3Types> PointModel;
typedef TPoint<Vec3Types> Point;

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_BUILD_MESH_COLLISION)
#ifndef SOFA_FLOAT
extern template class SOFA_MESH_COLLISION_API TPointModel<defaulttype::Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_MESH_COLLISION_API TPointModel<defaulttype::Vec3fTypes>;
#endif
#endif

//bool Point::testLMD(const Vector3 &PQ, double &coneFactor, double &coneExtension);

} // namespace collision

} // namespace component

} // namespace sofa

#endif
