#ifndef OBBMODEL_INL
#define OBBMODEL_INL

#include <sofa/component/collision/OBBModel.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/collision/Intersection.inl>

#include <sofa/helper/system/config.h>
#include <sofa/helper/proximity.h>
#include <sofa/defaulttype/Mat.h>
#include <sofa/defaulttype/Vec.h>
#include <iostream>
#include <algorithm>

#include <sofa/helper/system/FileRepository.h>
#include <sofa/component/collision/CubeModel.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/system/gl.h>
#include <sofa/helper/system/glut.h>

#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/simulation/common/Simulation.h>

namespace sofa
{

namespace component
{

namespace collision
{

using namespace sofa::defaulttype;
using namespace sofa::core::collision;
using namespace helper;


template<class DataTypes>
TOBBModel<DataTypes>::TOBBModel():
    ext(initData(&ext,"extents","Extents in x,y and z directions")),
    default_ext(initData(&default_ext,(Real)(1.0), "defaultExtent","Default extent")),
    _mstate(NULL)
{
    enum_type = OBB_TYPE;
}

template<class DataTypes>
TOBBModel<DataTypes>::TOBBModel(core::behavior::MechanicalState<DataTypes>* mstate):
    ext(initData(&ext, "extents","Extents in x,y and z directions")),
    default_ext(initData(&default_ext,(Real)(1.0), "defaultExtent","Default extent")),
    _mstate(mstate)
{
    enum_type = OBB_TYPE;
}


template<class DataTypes>
void TOBBModel<DataTypes>::init()
{
    this->CollisionModel::init();
    _mstate = dynamic_cast< core::behavior::MechanicalState<DataTypes>* > (getContext()->getMechanicalState());
    if (_mstate==NULL)
    {
        serr<<"TOBBModel requires a Rigid Mechanical Model" << sendl;
        return;
    }

    const int npoints = _mstate->getX()->size();
    resize(npoints);
}


template<class DataTypes>
void TOBBModel<DataTypes>::resize(int size){
    this->core::CollisionModel::resize(size);

    VecCoord & vext = *(ext.beginEdit());

    if ((int)vext.size() < size)
    {
        while((int)vext.size() < size)
            vext.push_back(Coord(default_ext.getValue(),default_ext.getValue(),default_ext.getValue()));
    }
    else
    {
        vext.reserve(size);
    }

    ext.endEdit();
}


template<class DataTypes>
void TOBBModel<DataTypes>::computeBoundingTree(int maxDepth){
    CubeModel* cubeModel = createPrevious<CubeModel>();
    const int npoints = _mstate->getX()->size();
    bool updated = false;
    if (npoints != size)
    {
        resize(npoints);
        updated = true;
        cubeModel->resize(0);
    }

    if (!isMoving() && !cubeModel->empty() && !updated)
        return; // No need to recompute BBox if immobile

    cubeModel->resize(size);
    if (!empty())
    {
        const typename TOBBModel<DataTypes>::Real distance = (typename TOBBModel<DataTypes>::Real)this->proximity.getValue();

        std::vector<Coord> vs;
        vs.reserve(8);
        for (int i=0; i<size; i++)
        {
            vs.clear();
            vertices(i,vs);

            Coord minElem = vs[0];
            Coord maxElem = vs[0];

            for(int j = 1 ; j < 8 ; ++j){
                for(int jj = 0 ; jj < 3 ; ++jj){
                    if(minElem[jj] > vs[j][jj])
                        minElem[jj] = vs[j][jj];
                    else if(maxElem[jj] < vs[j][jj])
                        maxElem[jj] = vs[j][jj];;
                }
            }

            for(int jj = 0 ; jj < 3 ; ++jj){
                minElem[jj] -= distance;
                maxElem[jj] += distance;
            }

            cubeModel->setParentOf(i, minElem, maxElem);
        }
        cubeModel->computeBoundingTree(maxDepth);
    }
}


template<class DataTypes>
void TOBBModel<DataTypes>::draw(const core::visual::VisualParams* vparams,int index){
    std::vector<Coord> p;
    vertices(index,p);

    Vec4f col4f(getColor4f());

    std::vector<Vector3> n;
    n.push_back(axis(index,1));
    //n.push_back(n.back());
    n.push_back(axis(index,0));
    //n.push_back(n.back());
    n.push_back(-n[0]);
//    n.push_back(n.back());
    n.push_back(-n[2]);
//    n.push_back(n.back());
    n.push_back(n.front());

    std::vector<Vector3> points;
    points.push_back(p[3]);
    points.push_back(p[0]);
    points.push_back(p[2]);
    points.push_back(p[1]);
    points.push_back(p[6]);
    points.push_back(p[5]);
    points.push_back(p[7]);
    points.push_back(p[4]);
    points.push_back(p[3]);
    points.push_back(p[0]);

    vparams->drawTool()->drawTriangleStrip(points,n,col4f);

    n.clear();
    points.clear();

    points.push_back(p[0]);
    points.push_back(p[1]);
    points.push_back(p[4]);
    points.push_back(p[5]);

    n.push_back(-axis(index,2));
    n.push_back(n.back());

    vparams->drawTool()->drawTriangleStrip(points,n,col4f);

    n.clear();
    points.clear();

    points.push_back(p[6]);
    points.push_back(p[7]);
    points.push_back(p[2]);
    points.push_back(p[3]);

    n.push_back(axis(index,2));
    n.push_back(n.back());

    vparams->drawTool()->drawTriangleStrip(points,n,col4f);
}

template<class DataTypes>
void TOBBModel<DataTypes>::draw(const core::visual::VisualParams* vparams){
    if (vparams->displayFlags().getShowCollisionModels())
    {
        vparams->drawTool()->setPolygonMode(0,vparams->displayFlags().getShowWireFrame());

        const int npoints = _mstate->getX()->size();
        vparams->drawTool()->setLightingEnabled(true); //Enable lightning
        for(int i = 0 ; i < npoints ; ++i )
            draw(vparams,i);
        vparams->drawTool()->setLightingEnabled(false); //Disable lightning
    }

    if (getPrevious()!=NULL && vparams->displayFlags().getShowBoundingCollisionModels())
        getPrevious()->draw(vparams);

    vparams->drawTool()->setPolygonMode(0,false);
}


template <class DataTypes>
inline typename TOBB<DataTypes>::Coord TOBB<DataTypes>::generalCoordinates(const Coord &c)const{
    return this->model->generalCoordinates(c,this->index);
}

template <class DataTypes>
inline bool TOBB<DataTypes>::onSurface(const Coord &c)const{
    Coord loc = this->localCoordinates(c);
    std::cout<<"loc "<<loc<<std::endl;
    for(int i = 0 ; i < 3 ; ++i){
        std::cout<<"val on surface "<<fabs(this->extent(i) - fabs(loc[i]))<<std::endl;
        if(loc[i] > extent(i) + 1e-6 || loc[i] < - extent(i) - 1e-6)
            return false;
        else if(fabs(this->extent(i) - fabs(loc[i])) < 1e-6)
            return true;
    }

    return false;
}

template <class DataTypes>
inline typename TOBBModel<DataTypes>::Coord TOBBModel<DataTypes>::generalCoordinates(const Coord &c,int index)const{
    return orientation(index).rotate(c) + center(index);
}


template <class DataTypes>
inline typename TOBBModel<DataTypes>::Coord TOBBModel<DataTypes>::localCoordinates(const Coord &c,int index)const{
    return orientation(index).inverseRotate(c - center(index));
}

template <class DataTypes>
inline typename TOBB<DataTypes>::Coord TOBB<DataTypes>::localCoordinates(const Coord & c)const{
    return this->model->localCoordinates(c,this->index);
}

template <class DataTypes>
inline const typename TOBBModel<DataTypes>::Coord & TOBBModel<DataTypes>::lvelocity(int index)const{
    return (*(_mstate->getV()))[index].getLinear();
}

template <class DataTypes>
inline const typename TOBB<DataTypes>::Coord & TOBB<DataTypes>::v()const{
    return this->model->lvelocity(this->index);
}

template<class DataTypes>
inline typename TOBBModel<DataTypes>::Coord TOBBModel<DataTypes>::axis(int index,int dim)const{
    Coord unit;
    if(dim == 0){
        unit[0] = 1;
        unit[1] = 0;
        unit[2] = 0;
    }
    else if(dim == 1){
        unit[0] = 0;
        unit[1] = 1;
        unit[2] = 0;
    }
    else{
        unit[0] = 0;
        unit[1] = 0;
        unit[2] = 1;
    }

    return orientation(index).rotate(unit);
}

template<class DataTypes>
inline typename TOBBModel<DataTypes>::Coord TOBBModel<DataTypes>::vertex(int index,int num)const{
    Real s0 = extent(index,0);
    Real s1 = extent(index,1);
    Real s2 = extent(index,2);

    if(num == 0){
        s0*=-1;
        s1*=1;
        s2*=-1;
    }
    else if(num == 1){
        s0*=1;
        s1*=1;
        s2*=-1;
    }
    else if(num == 2){
        s0*=1;
        s1*=1;
        s2*=1;
    }
    else if(num == 3){
        s0*=-1;
        s1*=1;
        s2*=1;
    }
    else if(num == 4){
        s0*=-1;
        s1*=-1;
        s2*=-1;
    }
    else if(num == 5){
        s0*=1;
        s1*=-1;
        s2*=-1;
    }
    else if(num == 6){
        s0*=1;
        s1*=-1;
        s2*=1;
    }
    else{
        s0*=-1;
        s1*=-1;
        s2*=1;
    }

    return center(index) + s0 * axis(index,0) + s1 * axis(index,1) + s2 * axis(index,2);
}

template<class DataTypes>
inline void TOBBModel<DataTypes>::axes(int index,Coord * v_axes)const{
    v_axes[0] = axis(index,0);
    v_axes[1] = axis(index,1);
    v_axes[2] = axis(index,2);
}

template<class DataTypes>
inline void TOBB<DataTypes>::axes(Coord * v_axes)const{
    this->model->axes(this->index,v_axes);
}

template<class DataTypes>
inline void TOBBModel<DataTypes>::vertices(int index,std::vector<Coord> & vs)const{
    Coord a0(axis(index,0) * extent(index,0));
    Coord a1(axis(index,1) * extent(index,1));
    Coord a2(axis(index,2) * extent(index,2));

    const Coord & c = center(index);

    vs.push_back(c + a0 - a1 - a2);
    vs.push_back(c + a0 + a1 - a2);
    vs.push_back(c + a0 + a1 + a2);
    vs.push_back(c + a0 - a1 + a2);
    vs.push_back(c - a0 - a1 - a2);
    vs.push_back(c - a0 + a1 - a2);
    vs.push_back(c - a0 + a1 + a2);
    vs.push_back(c - a0 - a1 + a2);
}

template<class DataTypes>
inline const typename TOBBModel<DataTypes>::Coord & TOBBModel<DataTypes>::center(int index)const{
    return (*_mstate->getX())[index].getCenter();
}

template<class DataTypes>
inline const typename TOBBModel<DataTypes>::Quaternion & TOBBModel<DataTypes>::orientation(int index)const{
    return (*_mstate->getX())[index].getOrientation();
}

template<class DataTypes>
inline typename TOBBModel<DataTypes>::Real TOBBModel<DataTypes>::extent(int index,int dim)const{
    return ((ext.getValue())[index])[dim];
}

template<class DataTypes>
inline const typename TOBBModel<DataTypes>::Coord & TOBBModel<DataTypes>::extents(int index)const{
    return (ext.getValue())[index];
}

template<class DataTypes>
inline typename TOBB<DataTypes>::Coord TOBB<DataTypes>::axis(int dim)const{
    return this->model->axis(this->index,dim);
}

template<class DataTypes>
inline typename TOBB<DataTypes>::Real TOBB<DataTypes>::extent(int dim)const{
    return this->model->extent(this->index,dim);
}

template<class DataTypes>
inline const typename TOBB<DataTypes>::Coord & TOBB<DataTypes>::extents()const{
    return this->model->extents(this->index);
}

template<class DataTypes>
inline const typename TOBB<DataTypes>::Coord & TOBB<DataTypes>::center()const{
    return this->model->center(this->index);
}

template<class DataTypes>
inline const typename TOBB<DataTypes>::Quaternion & TOBB<DataTypes>::orientation()const{
    return this->model->orientation(this->index);
}

template <class DataTypes>
inline Data<typename TOBBModel<DataTypes>::VecCoord> & TOBBModel<DataTypes>::writeExtents(){return ext;}

template <class DataTypes>
inline void TOBB<DataTypes>::vertices(std::vector<Coord> & vs)const{return this->model->vertices(this->index,vs);}

//template <class DataTypes>
//bool TOBBModel<DataTypes>::bidon()const{return false;}

}
}
}
#endif // OBBMODEL_INL
