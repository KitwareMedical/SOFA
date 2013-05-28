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
#ifndef SOFA_COMPONENT_COLLISION_LINEMODEL_H
#define SOFA_COMPONENT_COLLISION_LINEMODEL_H

#include <sofa/core/CollisionModel.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/component/collision/LocalMinDistanceFilter.h>
#include <sofa/component/container/MechanicalObject.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/defaulttype/Vec3Types.h>
#include <sofa/component/collision/PointModel.h>

namespace sofa
{

namespace component
{

namespace collision
{

using namespace sofa::defaulttype;

template<class DataTypes>
class TLineModel;

class LineLocalMinDistanceFilter;

template<class TDataTypes>
class TLine : public core::TCollisionElementIterator<TLineModel<TDataTypes> >
{
public:
    typedef TDataTypes DataTypes;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef TLineModel<DataTypes> ParentModel;

    TLine(ParentModel* model, int index);
    TLine() {}

    explicit TLine(const core::CollisionElementIterator& i);

    unsigned i1() const;
    unsigned i2() const;

    const Coord& p1() const;
    const Coord& p2() const;

    const Coord& p1Free() const;
    const Coord& p2Free() const;

    const Deriv& v1() const;
    const Deriv& v2() const;

    /// Return true if the element stores a free position vector
    bool hasFreePosition() const;

    bool activated(core::CollisionModel *cm = 0) const;

    // Return respectively the Vertex composing the neighbor Rigt and Left Triangle
//	const Vector3* tRight() const;
//	const Vector3* tLeft() const;
};

class LineActiver
{
public:
    LineActiver() {}
    virtual ~LineActiver() {}
    virtual bool activeLine(int /*index*/, core::CollisionModel * /*cm*/ = 0) {return true;}
};

template<class TDataTypes>
class TLineModel : public core::CollisionModel
{
public :
    SOFA_CLASS(SOFA_TEMPLATE(TLineModel, TDataTypes), core::CollisionModel);

protected:
    struct LineData
    {
        int i1,i2;
        // Triangles neighborhood
//		int tRight, tLeft;
    };

    sofa::helper::vector<LineData> elems;
    bool needsUpdate;
    virtual void updateFromTopology();

    TLineModel();

public:
//    typedef Vec3Types InDataTypes;
//    typedef Vec3Types DataTypes;
//    typedef DataTypes DataTypes;
    typedef TDataTypes DataTypes;
    typedef DataTypes InDataTypes;
    typedef TLineModel<DataTypes> ParentModel;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef TLine<DataTypes> Element;
    friend class TLine<DataTypes>;

    virtual void init();

    // -- CollisionModel interface

    virtual void resize(int size);

    virtual void computeBoundingTree(int maxDepth=0);

    virtual void computeContinuousBoundingTree(double dt, int maxDepth=0);

    void draw(const core::visual::VisualParams*,int index);

    void draw(const core::visual::VisualParams* vparams);

    virtual void handleTopologyChange();

    bool canCollideWithElement(int index, CollisionModel* model2, int index2);

    core::behavior::MechanicalState<DataTypes>* getMechanicalState() { return mstate; }

    //virtual const char* getTypeName() const { return "Line"; }

    Deriv velocity(int index)const;

    LineLocalMinDistanceFilter *getFilter() const;

    virtual int getElemEdgeIndex(int index) const { return index; }


    //template< class TFilter >
    //TFilter *getFilter() const
    //{
    //	if (m_lmdFilter != 0)
    //		return m_lmdFilter;
    //	else
    //		return &m_emptyFilter;
    //}

    void setFilter(LineLocalMinDistanceFilter * /*lmdFilter*/);

    Data<bool> bothSide; // to activate collision on both-side of the both side of the line model (when surface normals are defined on these lines)

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

    static std::string templateName(const TLineModel<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }


protected:

    core::behavior::MechanicalState<DataTypes>* mstate;
    Topology* topology;
    PointModel* mpoints;
    int meshRevision;
    LineLocalMinDistanceFilter *m_lmdFilter;

    Data< std::string  > LineActiverPath;
    Data<bool> m_displayFreePosition;

    LineActiver *myActiver;


};

template<class DataTypes>
inline TLine<DataTypes>::TLine(ParentModel* model, int index)
    : core::TCollisionElementIterator<ParentModel>(model, index)
{
//	activated = model->myActiver->activeLine(index);
}

template<class DataTypes>
inline TLine<DataTypes>::TLine(const core::CollisionElementIterator& i)
    : core::TCollisionElementIterator<ParentModel>(static_cast<ParentModel*>(i.getCollisionModel()), i.getIndex())
{
//	LineModel* CM = static_cast<LineModel*>(i.getCollisionModel());
//	activated = CM->myActiver->activeLine(i.getIndex());
}

template<class DataTypes>
inline unsigned TLine<DataTypes>::i1() const { return this->model->elems[this->index].i1; }

template<class DataTypes>
inline unsigned TLine<DataTypes>::i2() const { return this->model->elems[this->index].i2; }

template<class DataTypes>
inline const typename DataTypes::Coord& TLine<DataTypes>::p1() const { return (*this->model->mstate->getX())[this->model->elems[this->index].i1]; }

template<class DataTypes>
inline const typename DataTypes::Coord& TLine<DataTypes>::p2() const { return (*this->model->mstate->getX())[this->model->elems[this->index].i2]; }

template<class DataTypes>
inline const typename DataTypes::Coord& TLine<DataTypes>::p1Free() const
{
    if (hasFreePosition())
        return this->model->mstate->read(core::ConstVecCoordId::freePosition())->getValue()[this->model->elems[this->index].i1];
    else
        return p1();
}

template<class DataTypes>
inline const typename DataTypes::Coord& TLine<DataTypes>::p2Free() const
{
    if (hasFreePosition())
        return this->model->mstate->read(core::ConstVecCoordId::freePosition())->getValue()[this->model->elems[this->index].i2];
    else
        return p2();
}

template<class DataTypes>
inline const typename DataTypes::Deriv& TLine<DataTypes>::v1() const { return (*this->model->mstate->getV())[this->model->elems[this->index].i1]; }

template<class DataTypes>
inline const typename DataTypes::Deriv& TLine<DataTypes>::v2() const { return (*this->model->mstate->getV())[this->model->elems[this->index].i2]; }


template<class DataTypes>
inline typename TLineModel<DataTypes>::Deriv TLineModel<DataTypes>::velocity(int index) const { return ((*mstate->getV())[elems[index].i1] + (*mstate->getV())[elems[index].i2])/((Real)(2.0)); }

template<class DataTypes>
inline bool TLine<DataTypes>::hasFreePosition() const { return this->model->mstate->read(core::ConstVecCoordId::freePosition())->isSet(); }

template<class DataTypes>
inline bool TLine<DataTypes>::activated(core::CollisionModel *cm) const
{
    return this->model->myActiver->activeLine(this->index, cm);
}

//inline const Vector3* Line::tRight() const {
//	if (model->elems[index].tRight != -1)
//		return &(*model->mstate->getX())[model->elems[index].tRight];
//	else
//		return NULL;
//}

//inline const Vector3* Line::tLeft() const {
//	if (model->elems[index].tLeft != -1)
//		return &(*model->mstate->getX())[model->elems[index].tLeft];
//	else
//		return NULL;
//}

//class LineMeshModel : public LineModel
//{
//protected:
//    int meshRevision;
//    void updateFromTopology();
//
//public:
//    typedef topology::MeshTopology Topology;
//
//    LineMeshModel();
//
//    virtual void init();
//
//    Topology* getMeshTopology() { return mesh; }
//
//protected:
//    Topology* mesh;
//};

//class LineSetModel : public LineModel
//{
//public:
//
//    LineSetModel();
//
//    ///\Todo
//    virtual void init();
//
//    Topology* getMeshTopology() { return mesh; }
//
//protected:
//
//    Topology* mesh;
//
//    void updateFromTopology();
//};

typedef TLineModel<Vec3Types> LineModel;
typedef TLine<Vec3Types> Line;

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_BUILD_MESH_COLLISION)
#ifndef SOFA_FLOAT
extern template class SOFA_MESH_COLLISION_API TLineModel<defaulttype::Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_MESH_COLLISION_API TLineModel<defaulttype::Vec3fTypes>;
#endif
#endif

} // namespace collision

} // namespace component

} // namespace sofa

#endif
