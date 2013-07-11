#define SOFA_FLEXIBLE_ImageDensityMass_CPP

#include "ImageDensityMass.inl"
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/Mat.h>
#include <sofa/defaulttype/Vec3Types.h>

namespace sofa
{

namespace component
{

namespace mass
{

using namespace sofa::defaulttype;



//#ifndef SOFA_FLOAT
//template <> SOFA_BASE_MECHANICS_API
//double ImageDensityMass<Rigid3dTypes,core::behavior::ShapeFunction3d,Rigid3dMass>::getPotentialEnergy( const core::MechanicalParams* /* PARAMS FIRST */, const DataVecCoord& vx ) const
//{
//    const VecCoord& _x = vx.getValue();

//    VecCoord Mx = m_massMatrix * _x;

//    double e = 0;
//    // gravity
//    Vec3d g ( this->getContext()->getGravity() );
//    for( unsigned int i=0 ; i<_x.size() ; i++ )
//    {
//        e -= g*Mx[i].getCenter();
//    }
//    return e;
//}


//#endif
//#ifndef SOFA_DOUBLE
//template <> SOFA_BASE_MECHANICS_API
//double ImageDensityMass<Rigid3fTypes,core::behavior::ShapeFunction3f,Rigid3fMass>::getPotentialEnergy( const core::MechanicalParams* /* PARAMS FIRST */, const DataVecCoord& vx ) const
//{
//    const VecCoord& _x = vx.getValue();

//    VecCoord Mx = m_massMatrix * _x;

//    float e = 0;
//    // gravity
//    Vec3f g ( this->getContext()->getGravity() );
//    for( unsigned int i=0 ; i<_x.size() ; i++ )
//    {
//        e -= g*Mx[i].getCenter();
//    }
//    return e;
//}


//#endif




SOFA_DECL_CLASS(ImageDensityMass)


// Register in the Factory
int ImageDensityMassClass = core::RegisterObject("Define a global mass matrix including non diagonal terms")
#ifndef SOFA_FLOAT
        .add< ImageDensityMass<Vec3dTypes,core::behavior::ShapeFunction3d,Mat3x3d> >( true )
        .add< ImageDensityMass<Vec3dTypes,core::behavior::ShapeFunction2d,Mat3x3d> >()
//        .add< ImageDensityMass<Rigid3dTypes,core::behavior::ShapeFunction3d,Rigid3dMass> >()
#endif
#ifndef SOFA_DOUBLE
        .add< ImageDensityMass<Vec3fTypes,core::behavior::ShapeFunction3f,Mat3x3f> >()
        .add< ImageDensityMass<Vec3fTypes,core::behavior::ShapeFunction2f,Mat3x3f> >()
//        .add< ImageDensityMass<Rigid3fTypes,core::behavior::ShapeFunction3f,Rigid3fMass> >()
#endif
        ;
#ifndef SOFA_FLOAT
template class SOFA_Flexible_API ImageDensityMass<Vec3dTypes,core::behavior::ShapeFunction3d,Mat3x3d>;
template class SOFA_Flexible_API ImageDensityMass<Vec3dTypes,core::behavior::ShapeFunction2d,Mat3x3d>;
//template class SOFA_Flexible_API ImageDensityMass<Rigid3dTypes,core::behavior::ShapeFunction3d,Rigid3dMass>;
#endif
#ifndef SOFA_DOUBLE
template class SOFA_Flexible_API ImageDensityMass<Vec3fTypes,core::behavior::ShapeFunction3f,Mat3x3f>;
template class SOFA_Flexible_API ImageDensityMass<Vec3fTypes,core::behavior::ShapeFunction2f,Mat3x3f>;
//template class SOFA_Flexible_API ImageDensityMass<Rigid3fTypes,core::behavior::ShapeFunction3f,Rigid3fMass>;
#endif




} // namespace mass

} // namespace component

} // namespace sofa

