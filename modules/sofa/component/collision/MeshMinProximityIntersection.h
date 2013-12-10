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
#ifndef SOFA_COMPONENT_COLLISION_MESHMINPROXIMITYINTERSECTION_H
#define SOFA_COMPONENT_COLLISION_MESHMINPROXIMITYINTERSECTION_H

#include <sofa/component/collision/MinProximityIntersection.h>
#include <sofa/helper/FnDispatcher.h>
#include <sofa/component/collision/MeshIntTool.h>
#include <sofa/component/collision/CapsuleModel.h>
#include <sofa/component/collision/SphereModel.h>
#include <sofa/component/collision/TriangleModel.h>
#include <sofa/component/collision/LineModel.h>
#include <sofa/component/collision/PointModel.h>
#include <sofa/component/collision/CubeModel.h>
#include <sofa/component/collision/IntrUtility3.h>

namespace sofa
{

namespace component
{

namespace collision
{

class SOFA_MESH_COLLISION_API MeshMinProximityIntersection : public core::collision::BaseIntersector
{
    typedef MinProximityIntersection::OutputVector OutputVector;

public:
    MeshMinProximityIntersection(MinProximityIntersection* object, bool addSelf=true);

    bool testIntersection(Point&, Point&);
    template<class T> bool testIntersection(TSphere<T>&, Point&);
    bool testIntersection(Line&, Point&);
    template<class T> bool testIntersection(Line&, TSphere<T>&);
    bool testIntersection(Line&, Line&);
    bool testIntersection(Triangle&, Point&);
    template<class T> bool testIntersection(Triangle&, TSphere<T>&);
    bool testIntersection(Capsule&,Triangle&);
    bool testIntersection(Capsule&,Line&);

    int computeIntersection(Point&, Point&, OutputVector*);        
    template<class T> int computeIntersection(TSphere<T>&, Point&, OutputVector*);
    int computeIntersection(Line&, Point&, OutputVector*);
    template<class T> int computeIntersection(Line&, TSphere<T>&, OutputVector*);
    int computeIntersection(Line&, Line&, OutputVector*);
    int computeIntersection(Triangle&, Point&, OutputVector*);
    template<class T> int computeIntersection(Triangle&, TSphere<T>&, OutputVector*);
    int computeIntersection(Capsule & cap,Triangle & tri,OutputVector* contacts);
    int computeIntersection(Capsule & cap,Line & lin,OutputVector* contacts);

protected:

    MinProximityIntersection* intersection;
};




template <class T>
bool MeshMinProximityIntersection::testIntersection(Triangle& e2, TSphere<T>& e1)
{
    const SReal alarmDist = intersection->getAlarmDistance() + e1.r() + e1.getProximity() + e2.getProximity();

    const Vector3 x13 = e2.p1()-e2.p2();
    const Vector3 x23 = e2.p1()-e2.p3();
    const Vector3 x03 = e2.p1()-e1.center();
    Matrix2 A;
    Vector2 b;
    A[0][0] = x13*x13;
    A[1][1] = x23*x23;
    A[0][1] = A[1][0] = x13*x23;
    b[0] = x13*x03;
    b[1] = x23*x03;
    const SReal det = determinant(A);

    SReal alpha = 0.5;
    SReal beta = 0.5;

    //if (det < -0.000001 || det > 0.000001)
    {
        alpha = (b[0]*A[1][1] - b[1]*A[0][1])/det;
        beta  = (b[1]*A[0][0] - b[0]*A[1][0])/det;
        if (alpha < 0.000001 ||
            beta  < 0.000001 ||
            alpha + beta  > 0.999999)
            return false;
    }

    Vector3 P,Q,PQ;
    P = e1.center();
    Q = e2.p1() - x13 * alpha - x23 * beta;
    PQ = Q-P;

    if (PQ.norm2() < alarmDist*alarmDist)
    {
        return true;
    }
    else
        return false;
}

template <class T>
int MeshMinProximityIntersection::computeIntersection(Triangle& e2, TSphere<T>& e1, OutputVector* contacts)
{
    const SReal alarmDist = intersection->getAlarmDistance() + e1.r() + e1.getProximity() + e2.getProximity();

    const Vector3 x13 = e2.p1()-e2.p2();
    const Vector3 x23 = e2.p1()-e2.p3();
    const Vector3 x03 = e2.p1()-e1.center();
    Matrix2 A;
    Vector2 b;
    A[0][0] = x13*x13;
    A[1][1] = x23*x23;
    A[0][1] = A[1][0] = x13*x23;
    b[0] = x13*x03;
    b[1] = x23*x03;
    const SReal det = determinant(A);

    SReal alpha = 0.5;
    SReal beta = 0.5;

    //if (det < -0.000001 || det > 0.000001)
    {
        alpha = (b[0]*A[1][1] - b[1]*A[0][1])/det;
        beta  = (b[1]*A[0][0] - b[0]*A[1][0])/det;
        if (alpha < 0.000001 ||
            beta  < 0.000001 ||
            alpha + beta  > 0.999999)
            return 0;
    }

    Vector3 P = e1.center();
    Vector3 Q = e2.p1() - x13 * alpha - x23 * beta;
    Vector3 QP = P-Q;
    //Vector3 PQ = Q-P;

    if (QP.norm2() >= alarmDist*alarmDist)
        return 0;

    const SReal contactDist = intersection->getContactDistance() + e1.r() + e1.getProximity() + e2.getProximity();

    contacts->resize(contacts->size()+1);
    sofa::core::collision::DetectionOutput *detection = &*(contacts->end()-1);
    detection->elem = std::pair<core::CollisionElementIterator, core::CollisionElementIterator>(e2, e1);
    detection->id = e1.getIndex();
    detection->normal=QP;
    detection->value = detection->normal.norm();
    if(detection->value>1e-15)
    {
        detection->normal /= detection->value;
    }
    else
    {
        intersection->serr<<"WARNING: null distance between contact detected"<<intersection->sendl;
        detection->normal= Vector3(1,0,0);
    }
    detection->value -= contactDist;
    detection->point[0]=Q;
    detection->point[1]=e1.getContactPointByNormal( detection->normal );
    return 1;
}

template <class T>
bool MeshMinProximityIntersection::testIntersection(Line& e2, TSphere<T>& e1)
{
    const SReal alarmDist = intersection->getAlarmDistance() + e1.r() + e1.getProximity() + e2.getProximity();

    const Vector3 x32 = e2.p1()-e2.p2();
    const Vector3 x31 = e1.center()-e2.p2();
    SReal A;
    SReal b;
    A = x32*x32;
    b = x32*x31;

    SReal alpha = 0.5;

    //if (A < -0.000001 || A > 0.000001)
    {
        alpha = b/A;
        if (alpha < 0.000001 || alpha > 0.999999)
            return false;
    }

    Vector3 P,Q,PQ;
    P = e1.center();
    Q = e2.p1() - x32 * alpha;
    PQ = Q-P;

    if (PQ.norm2() < alarmDist*alarmDist)
    {
        return true;
    }
    else
        return false;
}

template <class T>
int MeshMinProximityIntersection::computeIntersection(Line& e2, TSphere<T>& e1, OutputVector* contacts)
{
    const SReal alarmDist = intersection->getAlarmDistance() + e1.r() + e1.getProximity() + e2.getProximity();

    const Vector3 x32 = e2.p1()-e2.p2();
    const Vector3 x31 = e1.center()-e2.p2();
    SReal A;
    SReal b;
    A = x32*x32;
    b = x32*x31;

    SReal alpha = 0.5;

    //if (A < -0.000001 || A > 0.000001)
    {
        alpha = b/A;
        if (alpha < 0.000001 || alpha > 0.999999)
            return 0;
    }

    Vector3 P = e1.center();
    Vector3 Q = e2.p1() - x32 * alpha;
    Vector3 QP = P-Q;
    //Vector3 PQ = Q-P;

    if (QP.norm2() >= alarmDist*alarmDist)
        return 0;

    const SReal contactDist = intersection->getContactDistance() + e1.r() + e1.getProximity() + e2.getProximity();

    contacts->resize(contacts->size()+1);
    sofa::core::collision::DetectionOutput *detection = &*(contacts->end()-1);
    detection->elem = std::pair<core::CollisionElementIterator, core::CollisionElementIterator>(e2, e1);
    detection->id = e1.getIndex();
    detection->normal=QP;
    detection->value = detection->normal.norm();
    if(detection->value>1e-15)
    {
        detection->normal /= detection->value;
    }
    else
    {
        intersection->serr<<"WARNING: null distance between contact detected"<<intersection->sendl;
        detection->normal= Vector3(1,0,0);
    }
    detection->point[0]=Q;
    detection->point[1]=e1.getContactPointByNormal( detection->normal );
    detection->value -= contactDist;
    return 1;
}

template <class T>
bool MeshMinProximityIntersection::testIntersection(TSphere<T>& e1, Point& e2)
{
    const SReal alarmDist = intersection->getAlarmDistance() + e1.r() + e1.getProximity() + e2.getProximity();

    Vector3 P,Q,PQ;
    P = e1.center();
    Q = e2.p();
    PQ = Q-P;

    if (PQ.norm2() < alarmDist*alarmDist)
    {
        return true;
    }
    else
        return false;
}

template <class T>
int MeshMinProximityIntersection::computeIntersection(TSphere<T>& e1, Point& e2, OutputVector* contacts)
{
    const SReal alarmDist = intersection->getAlarmDistance() + e1.r() + e1.getProximity() + e2.getProximity();

    Vector3 P,Q,PQ;
    P = e1.center();
    Q = e2.p();
    PQ = Q-P;
    if (PQ.norm2() >= alarmDist*alarmDist)
        return 0;

    const SReal contactDist = intersection->getContactDistance() + e1.r() + e1.getProximity() + e2.getProximity();

    contacts->resize(contacts->size()+1);
    sofa::core::collision::DetectionOutput *detection = &*(contacts->end()-1);
    detection->elem = std::pair<core::CollisionElementIterator, core::CollisionElementIterator>(e1, e2);
    //detection->id = (e1.getCollisionModel()->getSize() > e2.getCollisionModel()->getSize()) ? e1.getIndex() : e2.getIndex();
    detection->id = e1.getIndex();
    detection->normal=PQ;
    detection->value = detection->normal.norm();
    if(detection->value>1e-15)
    {
        detection->normal /= detection->value;
    }
    else
    {
        intersection->serr<<"WARNING: null distance between contact detected"<<intersection->sendl;
        detection->normal= Vector3(1,0,0);
    }
    detection->value -= contactDist;
    detection->point[0]=e1.getContactPointByNormal( -detection->normal );
    detection->point[1]=Q;
    return 1;
}



} // namespace collision

} // namespace component

} // namespace sofa

#endif
