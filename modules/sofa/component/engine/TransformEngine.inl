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
#ifndef SOFA_COMPONENT_ENGINE_TRANSFORMENGINE_INL
#define SOFA_COMPONENT_ENGINE_TRANSFORMENGINE_INL

#include <sofa/component/engine/TransformEngine.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/rmath.h> //M_PI

namespace sofa
{

namespace component
{

namespace engine
{

template <class DataTypes>
TransformEngine<DataTypes>::TransformEngine()
    : f_inputX ( initData (&f_inputX, "input_position", "input array of 3d points") )
    , f_outputX( initData (&f_outputX, "output_position", "output array of 3d points") )
    , translation(initData(&translation, defaulttype::Vector3(0,0,0),"translation", "translation vector ") )
    , rotation(initData(&rotation, defaulttype::Vector3(0,0,0), "rotation", "rotation vector ") )
    , scale(initData(&scale, defaulttype::Vector3(1,1,1),"scale", "scale factor") )
    , inverse(initData(&inverse, false, "inverse", "true to apply inverse transformation"))
{
}


template <class DataTypes>
void TransformEngine<DataTypes>::init()
{
    addInput(&f_inputX);
    addInput(&translation);
    addInput(&rotation);
    addInput(&scale);
    addInput(&inverse);
    addOutput(&f_outputX);
    setDirtyValue();
}

template <class DataTypes>
void TransformEngine<DataTypes>::reinit()
{
    update();
}

//Declare a TransformOperation class able to do an operation on a Coord
template <class DataTypes>
struct TransformOperation
{
    virtual ~TransformOperation() {}
    virtual void execute(typename DataTypes::Coord &v) const =0;
};

//*****************************************************************
//Scale Operation
template <class DataTypes>
struct Scale : public TransformOperation<DataTypes>
{
    typedef typename DataTypes::Real Real;
    Scale():sx(0),sy(0),sz(0) {}

    void execute(typename DataTypes::Coord &p) const
    {
        Real x,y,z;
        DataTypes::get(x,y,z,p);
        DataTypes::set(p,x*sx,y*sy,z*sz);
    }

    void configure(const defaulttype::Vector3 &s, bool inverse)
    {
        if (inverse)
        {
            sx=(Real)(1.0/s[0]); sy=(Real)(1.0/s[1]); sz=(Real)(1.0/s[2]);
        }
        else
        {
            sx=(Real)s[0]; sy=(Real)s[1]; sz=(Real)s[2];
        }
    }
private:
    Real sx,sy,sz;
};


//*****************************************************************
//Rotation Operation
template <class DataTypes>
struct Rotation : public TransformOperation<DataTypes>
{
    typedef typename DataTypes::Real Real;

    void execute(typename DataTypes::Coord &p) const
    {
        defaulttype::Vector3 pos;
        DataTypes::get(pos[0],pos[1],pos[2],p);
        pos=q.rotate(pos);
        DataTypes::set(p,pos[0],pos[1],pos[2]);
    }

    void configure(const defaulttype::Vector3 &r, bool inverse)
    {
        q=helper::Quater<Real>::createQuaterFromEuler( r*(M_PI/180.0));
        if (inverse)
            q = q.inverse();
    }
private:
    defaulttype::Quaternion q;
};

#ifndef SOFA_DOUBLE
template<>
void Rotation<defaulttype::Rigid3fTypes>::execute(defaulttype::Rigid3fTypes::Coord &p) const
{
    p.getCenter() = q.rotate(p.getCenter());
    p.getOrientation() = q*p.getOrientation();
}

#endif
#ifndef SOFA_FLOAT
template<>
void Rotation<defaulttype::Rigid3dTypes>::execute(defaulttype::Rigid3dTypes::Coord &p) const
{
    p.getCenter() = q.rotate(p.getCenter());
    p.getOrientation() = q*p.getOrientation();
}
#endif


//*****************************************************************
//Translation Operation
template <class DataTypes>
struct Translation : public TransformOperation<DataTypes>
{
    typedef typename DataTypes::Real Real;
    Translation():tx(0),ty(0),tz(0) {}
    void execute(typename DataTypes::Coord &p) const
    {
        Real x,y,z;
        DataTypes::get(x,y,z,p);
        DataTypes::set(p,x+tx,y+ty,z+tz);
    }
    void configure(const defaulttype::Vector3 &t, bool inverse)
    {
        if (inverse)
        {
            tx=(Real)-t[0]; ty=(Real)-t[1]; tz=(Real)-t[2];
        }
        else
        {
            tx=(Real)t[0]; ty=(Real)t[1]; tz=(Real)t[2];
        }
    }
private:
    Real tx,ty,tz;
};


//*****************************************************************
//Functor to apply the operations wanted
template <class DataTypes>
struct Transform
{
    typedef TransformOperation<DataTypes> Op;

    template <class  Operation>
    Operation* add(Operation *op, bool inverse)
    {
//     Operation *op=new Operation();
        if (inverse)
            list.push_front(op);
        else
            list.push_back(op);
        return op;
    }

    std::list< Op* > &getOperations() {return list;}

    void operator()(typename DataTypes::Coord &v) const
    {
        for (typename std::list< Op* >::const_iterator it=list.begin(); it != list.end() ; ++it)
        {
            (*it)->execute(v);
        }
    }
private:
    std::list< Op* > list;
};


template <class DataTypes>
void TransformEngine<DataTypes>::update()
{
    cleanDirty();

    const defaulttype::Vector3 &s=scale.getValue();
    const defaulttype::Vector3 &r=rotation.getValue();
    const defaulttype::Vector3 &t=translation.getValue();

    //Create the object responsible for the transformations
    Transform<DataTypes> transformation;
    const bool inv = inverse.getValue();
    if (s != defaulttype::Vector3(1,1,1))  transformation.add(new Scale<DataTypes>, inv)->configure(s, inv);
    if (r != defaulttype::Vector3(0,0,0))  transformation.add(new Rotation<DataTypes>, inv)->configure(r, inv);
    if (t != defaulttype::Vector3(0,0,0))  transformation.add(new Translation<DataTypes>, inv)->configure(t, inv);

    //Get input
    const VecCoord& in = f_inputX.getValue();
    VecCoord& out = *(f_outputX.beginEdit());

    //Set Output
    out.resize(in.size());
    //Set the output to the input
    std::copy(in.begin(),in.end(), out.begin());
    //Apply the transformation of the output
    std::for_each(out.begin(), out.end(), transformation);

    //Deleting operations
    std::list< TransformOperation<DataTypes>* > operations=transformation.getOperations();
    while (!operations.empty())
    {
        delete operations.back();
        operations.pop_back();
    }

    f_outputX.endEdit();
}



} // namespace engine

} // namespace component

} // namespace sofa

#endif
