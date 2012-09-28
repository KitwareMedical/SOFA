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
#ifndef SOFA_COMPONENT_TOPOLOGY_DYNAMICSPARSEGRIDGEOMETRYALGORITHMS_H
#define SOFA_COMPONENT_TOPOLOGY_DYNAMICSPARSEGRIDGEOMETRYALGORITHMS_H

#include <sofa/component/topology/HexahedronSetGeometryAlgorithms.h>
#include <sofa/component/topology/DynamicSparseGridTopologyContainer.h>
#include <sofa/component/container/MechanicalObject.h>

namespace sofa
{
namespace component
{
namespace topology
{
/**
* A class that provides geometry information on an HexahedronSet.
*/
template < class DataTypes >
class DynamicSparseGridGeometryAlgorithms : public HexahedronSetGeometryAlgorithms<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(DynamicSparseGridGeometryAlgorithms,DataTypes),SOFA_TEMPLATE(HexahedronSetGeometryAlgorithms,DataTypes));

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
    typedef sofa::component::container::MechanicalObject<DataTypes> MObject;
protected:
    DynamicSparseGridGeometryAlgorithms()
        : HexahedronSetGeometryAlgorithms<DataTypes>()
    {}

    virtual ~DynamicSparseGridGeometryAlgorithms() {}
public:
    /// finds a hexahedron, in its rest position, which is nearest to a given point. Computes barycentric coordinates and a distance measure.
    virtual int findNearestElementInRestPos(const Coord& pos, defaulttype::Vector3& baryC, Real& distance) const;

    virtual void init();

    HexaID getTopoIndexFromRegularGridIndex ( unsigned int index, bool& existing );
    unsigned int getRegularGridIndexFromTopoIndex ( HexaID index );

protected:
    DynamicSparseGridTopologyContainer* topoContainer;
    MObject* dof;
};


#ifndef SOFA_FLOAT
template <>
int SOFA_NON_UNIFORM_FEM_API DynamicSparseGridGeometryAlgorithms<defaulttype::Vec2dTypes>::findNearestElementInRestPos(const Coord& pos, sofa::defaulttype::Vector3& baryC, Real& distance) const;

template <>
int SOFA_NON_UNIFORM_FEM_API DynamicSparseGridGeometryAlgorithms<defaulttype::Vec1dTypes>::findNearestElementInRestPos(const Coord& pos, sofa::defaulttype::Vector3& baryC, Real& distance) const;
#endif

#ifndef SOFA_DOUBLE
template <>
int SOFA_NON_UNIFORM_FEM_API DynamicSparseGridGeometryAlgorithms<defaulttype::Vec2fTypes>::findNearestElementInRestPos(const Coord& pos, sofa::defaulttype::Vector3& baryC, Real& distance) const;

template <>
int SOFA_NON_UNIFORM_FEM_API DynamicSparseGridGeometryAlgorithms<defaulttype::Vec1fTypes>::findNearestElementInRestPos(const Coord& pos, sofa::defaulttype::Vector3& baryC, Real& distance) const;
#endif

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_BUILD_NON_UNIFORM_FEM)
#ifndef SOFA_FLOAT
extern template class SOFA_NON_UNIFORM_FEM_API DynamicSparseGridGeometryAlgorithms<defaulttype::Vec3dTypes>;
extern template class SOFA_NON_UNIFORM_FEM_API DynamicSparseGridGeometryAlgorithms<defaulttype::Vec2dTypes>;
extern template class SOFA_NON_UNIFORM_FEM_API DynamicSparseGridGeometryAlgorithms<defaulttype::Vec1dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_NON_UNIFORM_FEM_API DynamicSparseGridGeometryAlgorithms<defaulttype::Vec3fTypes>;
extern template class SOFA_NON_UNIFORM_FEM_API DynamicSparseGridGeometryAlgorithms<defaulttype::Vec2fTypes>;
extern template class SOFA_NON_UNIFORM_FEM_API DynamicSparseGridGeometryAlgorithms<defaulttype::Vec1fTypes>;
#endif
#endif

} // namespace topology

} // namespace component

} // namespace sofa

#endif
