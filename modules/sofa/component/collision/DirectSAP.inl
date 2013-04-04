#ifndef DIRECTSAP_INL
#define DIRECTSAP_INL


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
#include <sofa/component/collision/DirectSAP.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/component/collision/CapsuleModel.h>
#include <sofa/component/collision/Sphere.h>
#include <sofa/component/collision/Triangle.h>
#include <sofa/component/collision/Line.h>
#include <sofa/component/collision/Point.h>
#include <sofa/helper/FnDispatcher.h>
#include <sofa/core/ObjectFactory.h>
#include <map>
#include <queue>
#include <stack>

#include <sofa/helper/system/gl.h>
#include <sofa/helper/system/glut.h>

namespace sofa
{

namespace component
{

namespace collision
{

inline void SAPBox::update(int axis){
    min->value = (cube.minVect())[axis];
    max->value = (cube.maxVect())[axis];
}

inline bool SAPBox::overlaps(const SAPBox &other, int axis) const{
    const Vector3 & min0 = this->cube.minVect();
    const Vector3 & max0 = this->cube.maxVect();
    const Vector3 & min1 = other.cube.minVect();
    const Vector3 & max1 = other.cube.maxVect();

    if(min0[axis] >= max1[axis] || min1[axis] >= max0[axis])
        return false;

    return true;
}

using namespace sofa::defaulttype;
using namespace sofa::helper;
using namespace collision;

using namespace core::objectmodel;

template <template<class T,class Allocator> class List,template <class T> class Allocator>
TDirectSAP<List,Allocator>::TDirectSAP()
    : bDraw(initData(&bDraw, false, "draw", "enable/disable display of results"))
    , box(initData(&box, "box", "if not empty, objects that do not intersect this bounding-box will be ignored"))
{
}

template <template<class T,class Allocator> class List,template <class T> class Allocator>
TDirectSAP<List,Allocator>::~TDirectSAP()
{
//    for(typename EndPointList::iterator it = _end_points.begin() ; it != _end_points.end() ; ++it){
//        delete (*it);
//    }

    for(unsigned int i = 0 ; i < _to_del.size() ; ++i)
        delete _to_del[i];
}

template <template<class T,class Allocator> class List,template <class T> class Allocator>
void TDirectSAP<List,Allocator>::init()
{
    reinit();
}

template <template<class T,class Allocator> class List,template <class T> class Allocator>
void TDirectSAP<List,Allocator>::reinit()
{
    if (box.getValue()[0][0] >= box.getValue()[1][0])
    {
        boxModel.reset();
    }
    else
    {
        if (!boxModel) boxModel = sofa::core::objectmodel::New<CubeModel>();
        boxModel->resize(1);
        boxModel->setParentOf(0, box.getValue()[0], box.getValue()[1]);
    }
}

template <template<class T,class Allocator> class List,template <class T> class Allocator>
inline bool TDirectSAP<List,Allocator>::added(core::CollisionModel *cm) const
{
    return collisionModels.count(cm->getLast()) >= 1;
}

//template <template<class T,class Allocator> class List,template <class T> class Allocator>
//inline void TDirectSAP<List,Allocator>::add(core::CollisionModel *cm)
//{
//    collisionModels.insert(cm->getLast());
//}

template <template<class T,class Allocator> class List,template <class T> class Allocator>
inline void TDirectSAP<List,Allocator>::add(core::CollisionModel *cm)
{
    collisionModels.insert(cm->getLast());
    _new_cm.push_back(cm->getLast());
}


template <template<class T,class Allocator> class List,template <class T> class Allocator>
void TDirectSAP<List,Allocator>::endBroadPhase()
{
    BroadPhaseDetection::endBroadPhase();

    if(_new_cm.size() == 0)
        return;

    //to gain time, we create at the same time all SAPboxes so as to allocate
    //memory the less times
    std::vector<CubeModel*> cube_models;
    cube_models.reserve(_new_cm.size());

    int n = 0;
    for(unsigned int i = 0 ; i < _new_cm.size() ; ++i){
        n += _new_cm[i]->getSize();
        cube_models.push_back(dynamic_cast<CubeModel*>(_new_cm[i]->getPrevious()));
    }

    _boxes.reserve(_boxes.size() + n);
    EndPoint * end_pts = new EndPoint[2*n];
    _to_del.push_back(end_pts);

    int cur_EndPtID = 0;
    int cur_boxID = _boxes.size();
    for(unsigned int i = 0 ; i < cube_models.size() ; ++i){
        CubeModel * cm = cube_models[i];
        for(int j = 0 ; j < cm->getSize() ; ++j){
            EndPoint * min = &end_pts[cur_EndPtID];
            ++cur_EndPtID;
            EndPoint * max = &end_pts[cur_EndPtID];
            ++cur_EndPtID;

            min->setBoxID(cur_boxID);
            max->setBoxID(cur_boxID);
            max->setMax();

            _end_points.push_back(min);
            _end_points.push_back(max);

            _boxes.push_back(SAPBox(Cube(cm,j),min,max));
            ++cur_boxID;
        }
    }

    _new_cm.clear();
}


template <template<class T,class Allocator> class List,template <class T> class Allocator>
void TDirectSAP<List,Allocator>::addCollisionModel(core::CollisionModel *cm){
    if(!added(cm))
        add(cm);
}

//template <template<class T,class Allocator> class List,template <class T> class Allocator>
//void TDirectSAP<List,Allocator>::addCollisionModel(core::CollisionModel *cm)
//{
//    if(!added(cm)){
//        add(cm);
//        CubeModel * cube_model = dynamic_cast<CubeModel*>(cm->getLast()->getPrevious());

//        int old_size = _boxes.size();
//        int cm_size = cube_model->getSize();
//        int new_size = old_size + cm_size;

//        _boxes.reserve(new_size);
//        EndPoint * end_pts = new EndPoint[2*cm_size];
//        _to_del.push_back(end_pts);

//        for(int i = old_size ; i < new_size ; ++i){
//            //EndPoint * min = new EndPoint;
//            //EndPoint * max = new EndPoint;
//            end_pts[2*(i - old_size)].setBoxID(i);
//            end_pts[2*(i - old_size) + 1].setBoxID(i);
//            end_pts[2*(i - old_size) + 1].setMax();
//            _end_points.push_back(&end_pts[2*(i - old_size)]);
//            _end_points.push_back(&end_pts[2*(i - old_size) + 1]);
//            _boxes.push_back(SAPBox(Cube(cube_model,i - old_size),&end_pts[2*(i - old_size)],&end_pts[2*(i - old_size) + 1]));
//        }
//    }
//}

//template <template<class T,class Allocator> class List,template <class T> class Allocator>
//void TDirectSAP<List,Allocator>::addCollisionModel(core::CollisionModel *cm)
//{
//    if(!added(cm)){
//        add(cm);
//        CubeModel * cube_model = dynamic_cast<CubeModel*>(cm->getLast()->getPrevious());

//        int old_size = _boxes.size();
//        int cm_size = cube_model->getSize();
//        int new_size = old_size + cm_size;

//        _boxes.reserve(new_size);
//        EndPoint * end_pts = new EndPoint[2*cm_size];
//        _to_del.push_back(end_pts);

//        for(int i = old_size ; i < new_size ; ++i){
//            EndPoint * min = new EndPoint;
//            EndPoint * max = new EndPoint;
//            min->setBoxID(i);
//            max->setBoxID(i);
//            max->setMax();
//            _end_points.push_back(min);
//            _end_points.push_back(max);
//            _boxes.push_back(SAPBox(Cube(cube_model,i - old_size),min,max));
//        }
//    }
//}

template <template<class T,class Allocator> class List,template <class T> class Allocator>
int TDirectSAP<List,Allocator>::greatestVarianceAxis()const{
    double diff;
    double v[3];//variances for each axis
    double m[3];//means for each axis
    for(int i = 0 ; i < 3 ; ++i)
        v[i] = m[i] = 0;

    //computing the mean value of end points on each axis
    for(unsigned int i = 0 ; i < _boxes.size() ; ++i){
        const Vector3 & min = _boxes[i].cube.minVect();
        const Vector3 & max = _boxes[i].cube.maxVect();
        m[0] += min[0] + max[0];
        m[1] += min[1] + max[1];
        m[2] += min[2] + max[2];
    }

    m[0] /= 2*_boxes.size();
    m[1] /= 2*_boxes.size();
    m[2] /= 2*_boxes.size();

    //computing the variance of end points on each axis
    for(unsigned int i = 0 ; i < _boxes.size() ; ++i){
        const Vector3 & min = _boxes[i].cube.minVect();
        const Vector3 & max = _boxes[i].cube.maxVect();

        diff = min[0] - m[0];
        v[0] += diff*diff;
        diff = max[0] - m[0];
        v[0] += diff*diff;

        diff = min[1] - m[1];
        v[1] += diff*diff;
        diff = max[1] - m[1];
        v[1] += diff*diff;

        diff = min[2] - m[2];
        v[2] += diff*diff;
        diff = max[2] - m[2];
        v[2] += diff*diff;
    }

    if(v[0] >= v[1] && v[0] >= v[2])
        return 0;
    else if(v[1] >= v[2])
        return 1;
    else
        return 2;
}


template <template<class T,class Allocator> class List,template <class T> class Allocator>
void TDirectSAP<List,Allocator>::update(){
    _cur_axis = greatestVarianceAxis();
    for(unsigned int i = 0 ; i < _boxes.size() ; ++i){
        _boxes[i].update(_cur_axis);
    }
}

struct CompPEndPoint{
    bool operator()(const EndPoint * ep1,const EndPoint * ep2)const{
        if(ep1->value != ep2->value)
            return ep1->value < ep2->value;
        else if(ep1->boxID() == ep2->boxID())
            return ep2->max();
        else
            return ep1->boxID() < ep2->boxID();
    }
};

template <template<class T,class Allocator> class List,template <class T> class Allocator>
void TDirectSAP<List,Allocator>::beginNarrowPhase()
{
    core::collision::NarrowPhaseDetection::beginNarrowPhase();
    update();

    CompPEndPoint comp;
    std::sort(_end_points.begin(),_end_points.end(),comp);

    int axis1 = (1  << _cur_axis) & 3;
    int axis2 = (1  << axis1) & 3;

//    std::cout<<"sorted"<<std::endl;
//    for(int i = 0 ; i < _end_points.size() ; ++i){
//        //std::cout<<"boxID "<<(_end_points[i]->boxID())<<std::endl;
//        std::cout<<"nice index "<<(_boxes[_end_points[i]->boxID()].cube.getIndex())<<std::endl;
//        std::cout<<"collision model "<<(_boxes[_end_points[i]->boxID()].cube.getCollisionModel())<<std::endl;
//        if(_end_points[i]->max())
//            std::cout<<"max"<<std::endl;
//        else
//            std::cout<<"min"<<std::endl;

//        std::cout<<"value "<<(_end_points[i]->value)<<std::endl;
//        std::cout<<"==="<<std::endl;
//    }

    std::deque<int> active_boxes;//active boxes are the one that we encoutered only their min (end point), so if there are two boxes b0 and b1,
                                 //if we encounter b1_min as b0_min < b1_min, on the current axis, the two boxes intersect :  b0_min--------------------b0_max
                                 //                                                                                                      b1_min---------------------b1_max
                                 //once we encouter b0_max, b0 will not intersect with nothing (trivial), so we delete it from active_boxes.
                                 //so the rule is : -every time we encounter a box min end point, we check if it is overlapping with other active_boxes and add the owner (a box) of this end point to
                                 //                  the active boxes.
                                 //                 -every time we encounter a max end point of a box, we are sure that we encountered min end point of a box because _end_points is sorted,
                                 //                  so, we delete the owner box, of this max end point from the active boxes
    for(typename EndPointList::iterator it = _end_points.begin() ; it != _end_points.end() ; ++it){
        if((**it).max()){//erase it from the active_boxes
            assert(std::find(active_boxes.begin(),active_boxes.end(),(**it).boxID()) != active_boxes.end());
            active_boxes.erase(std::find(active_boxes.begin(),active_boxes.end(),(**it).boxID()));
        }
        else{//we encounter a min possible intersection between it and active_boxes
            int new_box = (**it).boxID();            

            SAPBox & box0 = _boxes[new_box];
            for(unsigned int i = 0 ; i < active_boxes.size() ; ++i){
                SAPBox & box1 = _boxes[active_boxes[i]];

                core::CollisionModel *finalcm1 = box0.cube.getCollisionModel()->getLast();//get the finnest CollisionModel which is not a CubeModel
                core::CollisionModel *finalcm2 = box1.cube.getCollisionModel()->getLast();
                if((finalcm1->isSimulated() || finalcm2->isSimulated()) &&
                        (((finalcm1->getContext() != finalcm2->getContext()) || finalcm1->canCollideWith(finalcm2)) && box0.overlaps(box1,axis1) && box0.overlaps(box1,axis2))){//intersection on all axes
                    //sout << "Final phase "<<gettypename(typeid(*finalcm1))<<" - "<<gettypename(typeid(*finalcm2))<<sendl;
//                    std::cout<<"finalcm1 finalcm2 "<<finalcm1<<" "<<finalcm2<<std::endl;
//                    std::cout<<"intersectionMethod "<<intersectionMethod->getClass()->className<<std::endl;
//                    std::cout<<"Final phase "<<finalcm1->getClass()->className<<" - "<<finalcm2->getClass()->className<<std::endl;

                    bool swapModels = false;
                    core::collision::ElementIntersector* finalintersector = intersectionMethod->findIntersector(finalcm1, finalcm2, swapModels);//find the method for the finnest CollisionModels

                    assert(box0.cube.getExternalChildren().first.getIndex() == box0.cube.getIndex());
                    assert(box1.cube.getExternalChildren().first.getIndex() == box1.cube.getIndex());

                    if((!swapModels) && finalcm1->getClass() == finalcm2->getClass() && finalcm1 > finalcm2)//we do that to have only pair (p1,p2) without having (p2,p1)
                        swapModels = true;

//                    std::cout<<"COLLISION"<<std::endl;
//                    std::cout<<"\t"<<finalcm1<<" "<<finalcm2<<std::endl;
//                    std::cout<<"\t"<<box0.cube.getIndex()<<" "<<box1.cube.getIndex()<<std::endl;
//                    std::cout<<"\t nice indices "<<new_box<<" "<<active_boxes[i]<<std::endl;

                    if(finalintersector != 0x0){
//                        std::cout<<"Final phase "<<finalcm1->getClass()->className<<" - "<<finalcm2->getClass()->className<<std::endl;
                        if(swapModels){
                            sofa::core::collision::DetectionOutputVector*& outputs = this->getDetectionOutputs(finalcm2, finalcm1);
                            finalintersector->beginIntersect(finalcm2, finalcm1, outputs);//creates outputs if null

                            if(finalintersector->intersect(box1.cube.getExternalChildren().first,box0.cube.getExternalChildren().first,outputs)){
                                //std::cout<<"\tREAL contact"<<std::endl;
                            }
                            else{
                                //std::cout<<"\tFALSE contact"<<std::endl;
                            }

                        }
                        else{
                            sofa::core::collision::DetectionOutputVector*& outputs = this->getDetectionOutputs(finalcm1, finalcm2);

                            finalintersector->beginIntersect(finalcm1, finalcm2, outputs);//creates outputs if null

                            if(finalintersector->intersect(box0.cube.getExternalChildren().first,box1.cube.getExternalChildren().first,outputs)){
                                //std::cout<<"\tREAL contact"<<std::endl;
                            }
                            else{
                                //std::cout<<"\tFALSE contact"<<std::endl;
                            }
                        }
                    }
                    else{
//                        std::cout<<"Final phase "<<finalcm1->getClass()->className<<" - "<<finalcm2->getClass()->className<<std::endl;
//                        std::cout<<"not found with intersectionMethod : "<<intersectionMethod->getClass()->className<<std::endl;
                    }
                }
            }
            active_boxes.push_back(new_box);
        }
    }
}


} // namespace collision

} // namespace component

} // namespace sofa







#endif // DIRECTSAP_INL
