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
#define SOFA_COMPONENT_PROJECTIVECONSTRAINTSET_FIXEDPLANECONSTRAINT_CPP
#include <sofa/component/projectiveconstraintset/FixedPlaneConstraint.inl>
#include <sofa/core/behavior/ProjectiveConstraintSet.inl>
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/Vec3Types.h>
#include <sofa/defaulttype/RigidTypes.h>

namespace sofa
{

namespace component
{

namespace projectiveconstraintset
{

using namespace sofa::defaulttype;
using namespace sofa::helper;

SOFA_DECL_CLASS(FixedPlaneConstraint)

int FixedPlaneConstraintClass = core::RegisterObject("Project particles on a given plane")
#ifndef SOFA_FLOAT
        .add< FixedPlaneConstraint<Vec3dTypes> >()
        .add< FixedPlaneConstraint<Vec6dTypes> >()
        .add< FixedPlaneConstraint<Rigid3dTypes> >()
#endif
#ifndef SOFA_DOUBLE
        .add< FixedPlaneConstraint<Vec3fTypes> >()
        .add< FixedPlaneConstraint<Vec6fTypes> >()
        .add< FixedPlaneConstraint<Rigid3fTypes> >()
#endif
        ;


#ifndef SOFA_FLOAT
template class SOFA_BOUNDARY_CONDITION_API FixedPlaneConstraint<Rigid3dTypes>;
template class SOFA_BOUNDARY_CONDITION_API FixedPlaneConstraint<Vec3dTypes>;
template class SOFA_BOUNDARY_CONDITION_API FixedPlaneConstraint<Vec6dTypes>;
#endif
#ifndef SOFA_DOUBLE
template class SOFA_BOUNDARY_CONDITION_API FixedPlaneConstraint<Rigid3fTypes>;
template class SOFA_BOUNDARY_CONDITION_API FixedPlaneConstraint<Vec3fTypes>;
template class SOFA_BOUNDARY_CONDITION_API FixedPlaneConstraint<Vec6fTypes>;
#endif

#ifndef SOFA_FLOAT
template <> template <class DataDeriv>
void FixedPlaneConstraint<Rigid3dTypes>::projectResponseT(const core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataDeriv& res)
{
    const int N=Coord::spatial_dimensions;
    Vec<N,Real> dir=direction.getValue().getCenter();

    for (helper::vector< unsigned int > ::const_iterator it = this->indices.getValue().begin(); it != this->indices.getValue().end(); ++it)
    {
        getVCenter(res[*it]) -= dir*(dir*(getVCenter(res[*it])));
    }
}

template <>
bool FixedPlaneConstraint<Rigid3dTypes>::isPointInPlane(Rigid3dTypes::Coord p)
{
    const int N=Coord::spatial_dimensions;
    Vec<N,Real> pos = p.getCenter();
    Real d=pos*direction.getValue().getCenter();
    if ((d>dmin.getValue())&& (d<dmax.getValue()))
        return true;
    else
        return false;
}
#endif

#ifndef SOFA_DOUBLE
template <> template <class DataDeriv>
void FixedPlaneConstraint<Rigid3fTypes>::projectResponseT(const core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataDeriv& res)
{
    const int N=Coord::spatial_dimensions;
    Vec<N,Real> dir=direction.getValue().getCenter();

    for (helper::vector< unsigned int > ::const_iterator it = this->indices.getValue().begin(); it != this->indices.getValue().end(); ++it)
    {
        getVCenter(res[*it]) -= dir*(dir*(getVCenter(res[*it])));
    }
}

template <>
bool FixedPlaneConstraint<Rigid3fTypes>::isPointInPlane(Coord p)
{
    const int N=Coord::spatial_dimensions;
    Vec<N,Real> pos = p.getCenter();
    Real d=pos*direction.getValue().getCenter();
    if ((d>dmin.getValue())&& (d<dmax.getValue()))
        return true;
    else
        return false;
}
#endif

} // namespace projectiveconstraintset

} // namespace component

} // namespace sofa

