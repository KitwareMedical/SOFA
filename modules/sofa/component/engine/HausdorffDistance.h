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
#ifndef SOFA_COMPONENT_ENGINE_HAUSDORFFDISTANCE_H
#define SOFA_COMPONENT_ENGINE_HAUSDORFFDISTANCE_H

#include <sofa/core/DataEngine.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/defaulttype/Vec3Types.h>
#include <sofa/defaulttype/RigidTypes.h>

#include <sofa/component/component.h>

namespace sofa
{

namespace component
{

namespace engine
{

/**
 * This class compute the Hausdorff distance of two point clouds
 * \todo: mean and mean square error
 */
template <class DataTypes>
class HausdorffDistance : public core::DataEngine
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(HausdorffDistance,DataTypes),core::DataEngine);
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::VecCoord VecCoord;

protected:

    HausdorffDistance();

    virtual ~HausdorffDistance() {}
public:
    void init();

    void reinit();

    void update();

    //Input
    Data<VecCoord> f_points_1;
    Data<VecCoord> f_points_2;

    //Output
    Data<Real> d12;
    Data<Real> d21;
    Data<Real> max;

    Data<bool> f_update;

    virtual std::string getTemplateName() const
    {
        return templateName(this);
    }


    static std::string templateName(const HausdorffDistance<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }

protected:

    void computeDistances();

    Real distance(Coord p, VecCoord S);

};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_ENGINE_HAUSDORFFDISTANCE_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_ENGINE_API HausdorffDistance<defaulttype::Vec2dTypes>;
extern template class SOFA_ENGINE_API HausdorffDistance<defaulttype::Vec3dTypes>;
extern template class SOFA_ENGINE_API HausdorffDistance<defaulttype::Rigid2dTypes>;
extern template class SOFA_ENGINE_API HausdorffDistance<defaulttype::Rigid3dTypes>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
extern template class SOFA_ENGINE_API HausdorffDistance<defaulttype::Vec2fTypes>;
extern template class SOFA_ENGINE_API HausdorffDistance<defaulttype::Vec3fTypes>;
extern template class SOFA_ENGINE_API HausdorffDistance<defaulttype::Rigid2fTypes>;
extern template class SOFA_ENGINE_API HausdorffDistance<defaulttype::Rigid3fTypes>;
#endif //SOFA_DOUBLE
#endif

} // namespace engine

} // namespace component

} // namespace sofa

#endif //SOFA_COMPONENT_ENGINE_HAUSDORFFDISTANCE_H
