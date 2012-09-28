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
#ifndef SOFA_COMPONENT_COLLISION_SPATIALGRIDPOINTMODEL_H
#define SOFA_COMPONENT_COLLISION_SPATIALGRIDPOINTMODEL_H

#include <sofa/component/collision/PointModel.h>
#include <sofa/component/container/SpatialGridContainer.h>
#include <sofa/defaulttype/Vec3Types.h>
#include <vector>

namespace sofa
{

namespace component
{

namespace collision
{

using namespace sofa::defaulttype;

class SOFA_MISC_COLLISION_API SpatialGridPointModel : public PointModel
{
public:
    SOFA_CLASS(SpatialGridPointModel, PointModel);

    typedef container::SpatialGridContainer<Vec3Types> GridContainer;
    typedef GridContainer::Grid Grid;

    Data<int> d_leafScale;
protected:
    SpatialGridPointModel();
public:
    virtual void init();

    // -- CollisionModel interface

    virtual void computeBoundingTree(int maxDepth=0);

    Grid* getGrid() { return grid->getGrid(); }
protected:

    GridContainer* grid;
    class OctreeCell
    {
    public:
        Grid::Key k;
        int pfirst;
        int plast;
        OctreeCell( Grid::Key k = Grid::Key(), int pfirst = 0, int plast = -1)
            : k(k), pfirst(pfirst), plast(plast)
        {
        }
    };
    class OctreeSorter
    {
    public:
        int root_shift;
        OctreeSorter(int root_shift=8) : root_shift(root_shift) {}
        bool operator()(const Grid::Key& k1, const Grid::Key &k2);
        bool operator()(const OctreeCell& c1, const OctreeCell &c2)
        {
            return (*this)(c1.k,c2.k);
        }
    };
};

} // namespace collision

} // namespace component

} // namespace sofa

#endif
