#include <sofa/component/collision/CapsuleIntTool.h>

namespace sofa
{
namespace component
{
namespace collision
{
using namespace sofa::defaulttype;
using namespace sofa::core::collision;

int CapsuleIntTool::computeIntersection(Capsule & e1,Capsule & e2,SReal alarmDist,SReal contactDist,OutputVector * contacts){
    SReal contact_exists = e1.radius() + e2.radius() + alarmDist;

    Vector3 A = e1.point1();
    Vector3 B = e1.point2();
    Vector3 C = e2.point1();
    Vector3 D = e2.point2();
    const Vector3 AB = B-A;//segment of the capsule e1
    const Vector3 CD = D-C;//segment of the capsule e2
    const Vector3 AC = C-A;
    Matrix2 Amat;//matrix helping us to find the two nearest points lying on the segments of the two capsules
    Vector2 b;

    Amat[0][0] = AB*AB;
    Amat[1][1] = CD*CD;
    Amat[0][1] = Amat[1][0] = -CD*AB;
    b[0] = AB*AC;
    b[1] = -CD*AC;
    const SReal det = determinant(Amat);

    SReal AB_norm2 = AB.norm2();
    SReal CD_norm2 = CD.norm2();
    SReal alpha = 0.5;
    SReal beta = 0.5;
    Vector3 P,Q,PQ;
    //Check that the determinant is not null which would mean that the capsule segments are lying on a same plane.
    //in this case we can solve the little system which gives us
    //the two coefficients alpha and beta. We obtain the two nearest points P and Q lying on the segments of the two capsules.
    //P = A + AB * alpha;
    //Q = C + CD * beta;
    if (det < -1.0e-15 || det > 1.0e-15)
    {
        alpha = (b[0]*Amat[1][1] - b[1]*Amat[0][1])/det;
        beta  = (b[1]*Amat[0][0] - b[0]*Amat[1][0])/det;

        if (alpha < 0)
            alpha = 0;
        else if(alpha > 1)
            alpha = 1;

        if (beta < 0)
            beta = 0;
        else if(beta > 1)
            beta = 1;
    }
    else{//Capsule segments on a same plane. Here the idea to find the nearest points
        //is to project segment apexes on the other segment.
        //Visual example with semgents AB and CD :
        //            A----------------B
        //                     C----------------D
        //After projection :
        //            A--------c-------B
        //                     C-------b--------D
        //So the nearest points are p and q which are respecively in the middle of cB and Cb:
        //            A--------c---p---B
        //                     C---q---b--------D

        Vector3 AD = D - A;
        Vector3 CB = B - C;

        SReal c_proj= b[0]/AB_norm2;//alpha = (AB * AC)/AB_norm2
        SReal d_proj = (AB * AD)/AB_norm2;
        SReal a_proj = b[1]/CD_norm2;//beta = (-CD*AC)/CD_norm2
        SReal b_proj= (CD*CB)/CD_norm2;

        if(c_proj >= 0 && c_proj <= 1){//projection of C on AB is lying on AB
            if(d_proj > 1){//case :
                           //             A----------------B
                           //                      C---------------D
                alpha = (1.0 + c_proj)/2.0;
                beta = b_proj/2.0;
            }
            else if(d_proj < 0){//case :
                                //             A----------------B
                                //     D----------------C
                alpha = c_proj/2.0;
                beta = (1 + a_proj)/2.0;
            }
            else{//case :
                //             A----------------B
                //                 C------D
                alpha = (c_proj + d_proj)/2.0;
                beta  = 0.5;
            }
        }
        else if(d_proj >= 0 && d_proj <= 1){
            if(c_proj < 0){//case :
                           //             A----------------B
                           //     C----------------D
                alpha = d_proj /2.0;
                beta = (1 + a_proj)/2.0;
            }
            else{//case :
                 //          A---------------B
                 //                 D-------------C
                alpha = (1 + d_proj)/2.0;
                beta = b_proj/2.0;
            }
        }
        else{
            if(c_proj * d_proj < 0){//case :
                                    //           A--------B
                                    //       D-----------------C
                alpha = 0.5;
                beta = (a_proj + b_proj)/2.0;
            }
            else{
                if(c_proj < 0){//case :
                               //                    A---------------B
                               // C-------------D
                    alpha = 0;
                }
                else{
                    alpha = 1;
                }

                if(a_proj < 0){//case :
                               // A---------------B
                               //                     C-------------D
                    beta = 0;
                }
                else{//case :
                     //                     A---------------B
                     //   C-------------D
                    beta = 1;
                }
            }
        }
    }

    if(alpha < 0){
        alpha = 0;
        beta = (CD * (A - C))/CD_norm2;
    }
    else if(alpha > 1){
        alpha = 1;
        beta = (CD * (B - C))/CD_norm2;
    }

    if(beta < 0){
        beta = 0;
        alpha = (AB * (C - A))/AB_norm2;
    }
    else if(beta > 1){
        beta = 1;
        alpha = (AB * (D - A))/AB_norm2;
    }

    if(alpha < 0)
        alpha = 0;
    else if (alpha > 1)
        alpha = 1;

    P = A + AB * alpha;
    Q = C + CD * beta;
    PQ = Q-P;

    SReal norm2 = PQ.norm2();

    if (norm2 > contact_exists * contact_exists)
        return 0;

    contacts->resize(contacts->size()+1);
    DetectionOutput *detection = &*(contacts->end()-1);

    const SReal theory_contactDist = e1.radius() + e2.radius() + contactDist;

    detection->elem = std::pair<core::CollisionElementIterator, core::CollisionElementIterator>(e1, e2);
    detection->id = (e1.getCollisionModel()->getSize() > e2.getCollisionModel()->getSize()) ? e1.getIndex() : e2.getIndex();
    detection->value = helper::rsqrt( norm2 );
    detection->normal = PQ / detection->value;
    detection->point[0] = P + e1.radius() * detection->normal;
    detection->point[1] = Q - e2.radius() * detection->normal;

    detection->value -= theory_contactDist;
    return 1;
}


int CapsuleIntTool::computeIntersection(Capsule& cap, OBB& obb,SReal alarmDist,SReal contactDist,OutputVector* contacts){
    IntrCapsuleOBB intr(cap,obb);
    if(intr.Find(alarmDist)){
        OBB::Real dist2 = (intr.pointOnFirst() - intr.pointOnSecond()).norm2();
        if((!intr.colliding()) && dist2 > alarmDist * alarmDist)
            return 0;

        contacts->resize(contacts->size()+1);
        DetectionOutput *detection = &*(contacts->end()-1);

        detection->normal = intr.separatingAxis();
        detection->point[0] = intr.pointOnFirst();
        detection->point[1] = intr.pointOnSecond();

        if(intr.colliding())
            detection->value = -helper::rsqrt(dist2) - contactDist;
        else
            detection->value = helper::rsqrt(dist2) - contactDist;

        detection->elem.first = cap;
        detection->elem.second = obb;
        //detection->id = (cap.getCollisionModel()->getSize() > obb.getCollisionModel()->getSize()) ? cap.getIndex() : obb.getIndex();
        detection->id = cap.getIndex();

        return 1;
    }

    return 0;
}

}
}
}
