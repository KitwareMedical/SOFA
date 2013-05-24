#ifndef INTRSPHEREOBB_H
#define INTRSPHEREOBB_H
#include <sofa/component/collision/OBBModel.h>
#include <sofa/component/collision/SphereModel.h>
#include <sofa/component/collision/IntrUtility3.h>
#include <sofa/component/collision/Intersector.h>

namespace sofa{
namespace component{
namespace collision{

/**
  *TDataTypes is the sphere type and TDataTypes2 the OBB type.
  */
template <typename TDataTypes,typename TDataTypes2>
class TIntrSphereOBB : public Intersector<typename TDataTypes::Real>
{
public:
    typedef TSphere<TDataTypes> IntrSph;
    typedef typename IntrSph::Real Real;
    typedef typename IntrSph::Coord Coord;
    typedef TOBB<TDataTypes2> Box;
    typedef Vec<3,Real> Vec3;

    TIntrSphereOBB (const IntrSph& sphere, const Box & box);

    /**
      *The idea of finding contact points is simple : project
      *the sphere center on the OBB and find the intersection point
      *on the OBB. Once we have this point we project it on the sphere.
      */
    bool Find ();

    Real distance()const;
private:
    using Intersector<Real>::_is_colliding;
    using Intersector<Real>::_pt_on_first;
    using Intersector<Real>::_pt_on_second;
    using Intersector<Real>::mContactTime;
    using Intersector<Real>::_sep_axis;

    // The objects to intersect.
    const IntrSph* _sph;
    const Box * mBox;
};

typedef TIntrSphereOBB<Vec3Types,Rigid3Types> IntrSphereOBB;

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_BUILD_BASE_COLLISION)
#ifndef SOFA_FLOAT
extern template class SOFA_BASE_COLLISION_API TIntrSphereOBB<Vec3dTypes,Rigid3dTypes>;
extern template class SOFA_BASE_COLLISION_API TIntrSphereOBB<Rigid3dTypes,Rigid3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_BASE_COLLISION_API TIntrSphereOBB<Vec3fTypes,Rigid3fTypes>;
extern template class SOFA_BASE_COLLISION_API TIntrSphereOBB<Rigid3fTypes,Rigid3fTypes>;
#endif
#endif

}
}
}
#endif // INTRSPHEREOBB_H
