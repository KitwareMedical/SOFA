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
#ifndef SOFA_GPU_CUDA_CUDAPOINTMODEL_H
#define SOFA_GPU_CUDA_CUDAPOINTMODEL_H

#include "CudaTypes.h"

#include <sofa/core/CollisionModel.h>
#include <sofa/component/container/MechanicalObject.h>
#include <sofa/defaulttype/Vec3Types.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/helper/io/Mesh.h>
#include <sofa/core/topology/BaseMeshTopology.h>

namespace sofa
{

namespace gpu
{

namespace cuda
{

using namespace sofa::defaulttype;
//using namespace sofa::component::collision;

class CudaPointModel;

class CudaPoint : public core::TCollisionElementIterator<CudaPointModel>
{
public:
    CudaPoint(CudaPointModel* model, int index);

    int i0();
    int getSize();

    explicit CudaPoint(const core::CollisionElementIterator& i);
};

class CudaPointModel : public core::CollisionModel
{
public:
    SOFA_CLASS(CudaPointModel,core::CollisionModel);
    typedef CudaVec3fTypes InDataTypes;
    typedef CudaVec3fTypes DataTypes;
    typedef DataTypes::VecCoord VecCoord;
    typedef DataTypes::VecDeriv VecDeriv;
    typedef DataTypes::Coord Coord;
    typedef DataTypes::Deriv Deriv;
    typedef CudaPoint Element;
    friend class CudaPoint;

    Data<int> groupSize;

    CudaPointModel();

    virtual void init();

    // -- CollisionModel interface

    virtual void resize(int size);

    virtual void computeBoundingTree(int maxDepth=0);

    //virtual void computeContinuousBoundingTree(double dt, int maxDepth=0);

    void draw(const core::visual::VisualParams*,int index);

    void draw(const core::visual::VisualParams*);

    core::behavior::MechanicalState<InDataTypes>* getMechanicalState() { return mstate; }

protected:

    core::behavior::MechanicalState<InDataTypes>* mstate;
};

inline CudaPoint::CudaPoint(CudaPointModel* model, int index)
    : core::TCollisionElementIterator<CudaPointModel>(model, index)
{}

inline CudaPoint::CudaPoint(const core::CollisionElementIterator& i)
    : core::TCollisionElementIterator<CudaPointModel>(static_cast<CudaPointModel*>(i.getCollisionModel()), i.getIndex())
{
}

inline int CudaPoint::i0()
{
    return model->groupSize.getValue()*index;
}

inline int CudaPoint::getSize()
{
    if (index == model->getSize()-1)
        return model->getMechanicalState()->getX()->size();
    else
        return model->groupSize.getValue();
}

} // namespace cuda

} // namespace gpu

} // namespace sofa

#endif
