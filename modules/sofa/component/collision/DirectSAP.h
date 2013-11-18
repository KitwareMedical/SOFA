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

#ifndef DIRECTSAP_H
#define DIRECTSAP_H

#include <sofa/core/collision/BroadPhaseDetection.h>
#include <sofa/core/collision/NarrowPhaseDetection.h>
#include <sofa/core/CollisionElement.h>
#include <sofa/core/CollisionModel.h>
#include <sofa/component/component.h>
#include <sofa/component/collision/CubeModel.h>
#include <sofa/component/collision/EndPoint.h>
#include <sofa/defaulttype/Vec.h>
#include <set>
#include <map>
#include <deque>
#include <sofa/helper/AdvancedTimer.h>

namespace sofa
{

namespace component
{

namespace collision
{

class EndPoint;

/**
  *SAPBox is a simple bounding box. It contains a Cube which contains only one final
  *CollisionElement and pointers to min and max EndPoints. min and max end points
  *are respectively min and max coordinates of the cube on a coordinate axis.
  *min and max are updated with the method update(int i), so min and max have
  *min/max values on the i-th axis after the method update(int i).
  */
class SOFA_MESH_COLLISION_API DSAPBox{
public:
    DSAPBox(Cube c,EndPoint * mi = 0x0,EndPoint * ma = 0x0) : cube(c),min(mi),max(ma){}

    void update(int axis,double alarmDist);

    bool overlaps(const DSAPBox & other,int axis,double alarmDist)const;

    bool overlaps(const DSAPBox &other,double alarmDist)const;

    bool sqOverlaps(const DSAPBox &other,double squaredAlarmDist)const;

    double squaredDistance(const DSAPBox & other)const;

    double squaredDistance(const DSAPBox & other,int axis)const;


    inline void show()const{
        std::cout<<"MIN "<<cube.minVect()<<std::endl;
        std::cout<<"MAX "<<cube.maxVect()<<std::endl;
    }

    Cube cube;
    EndPoint * min;
    EndPoint * max;
};

using namespace sofa::defaulttype;

/**
  *This class is an implementation of sweep and prune in its "direct" version, i.e. at each step
  *it sorts all the primitives along an axis (not checking the moving ones) and computes overlaping pairs without
  *saving it. But the memory used to save these primitives is created just once, the first time we add CollisionModels.
  */
template <template<class T,class Allocator> class List,template <class T> class Allocator = std::allocator>
class TDirectSAP :
    public core::collision::BroadPhaseDetection,
    public core::collision::NarrowPhaseDetection
{
public:
    SOFA_CLASS2(SOFA_TEMPLATE2(TDirectSAP,List,Allocator), core::collision::BroadPhaseDetection, core::collision::NarrowPhaseDetection);

    typedef List<EndPoint*,Allocator<EndPoint*> > EndPointList;

    typedef DSAPBox SAPBox;

    //void collidingCubes(std::vector<std::pair<Cube,Cube> > & col_cubes)const;
private:
    /**
      *Returns the axis number which have the greatest variance for the primitive end points.
      *This axis is used when updating and sorting end points. The greatest variance means
      *that this axis have the most chance to eliminate a maximum of not overlaping SAPBox pairs
      *because along this axis, SAPBoxes are the sparsest.
      */
    int greatestVarianceAxis()const;

    bool added(core::CollisionModel * cm)const;

    void add(core::CollisionModel * cm);

    /**
      *Updates values of end points. These values are coordinates of AABB on axis that maximazes the variance for the AABBs.
      */
    void update();

    Data<bool> bDraw;

    Data< helper::fixed_array<Vector3,2> > box;

    CubeModel::SPtr boxModel;

    std::vector<DSAPBox> _boxes;//boxes
    EndPointList _end_points;//end points of _boxes
    int _cur_axis;//the current greatest variance axis

    std::set<core::CollisionModel*> collisionModels;//used to check if a collision model is added
    std::vector<core::CollisionModel*> _new_cm;//eventual new collision models to  add at a step

    double _alarmDist;
    double _alarmDist_d2;
    double _sq_alarmDist;
protected:
    TDirectSAP();

    ~TDirectSAP();

    std::vector<EndPoint*> _to_del;//EndPoint arrays to delete when deleting DirectSAP
public:
    void setDraw(bool val) { bDraw.setValue(val); }

    void init();
    void reinit();

    void addCollisionModel (core::CollisionModel *cm);    

    /**
      *Unuseful methods because all is done in addCollisionModel
      */
    void addCollisionPair (const std::pair<core::CollisionModel*, core::CollisionModel*>& ){}
    void addCollisionPairs (std::vector<std::pair<core::CollisionModel*, core::CollisionModel*> >&){}

    virtual void endBroadPhase();
    virtual void beginNarrowPhase();


    /* for debugging */
    inline void draw(const core::visual::VisualParams*){}

    inline virtual bool needsDeepBoundingTree()const{return false;}
};

typedef TDirectSAP<std::vector,std::allocator> DirectSAP;

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_BUILD_MESH_COLLISION)
extern template class SOFA_MESH_COLLISION_API TDirectSAP<helper::vector,helper::CPUMemoryManager>;
extern template class SOFA_MESH_COLLISION_API TDirectSAP<std::vector,std::allocator>;
#endif

} // namespace collision

} // namespace component

} // namespace sofa

#endif // BRUTESAP_H
