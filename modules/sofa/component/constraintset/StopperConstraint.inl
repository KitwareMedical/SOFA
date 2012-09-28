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
#ifndef SOFA_COMPONENT_CONSTRAINTSET_STOPPERCONSTRAINT_INL
#define SOFA_COMPONENT_CONSTRAINTSET_STOPPERCONSTRAINT_INL

#include <sofa/component/constraintset/StopperConstraint.h>
#include <sofa/core/visual/VisualParams.h>

#include <sofa/defaulttype/Vec.h>
#include <sofa/helper/gl/template.h>
namespace sofa
{

namespace component
{

namespace constraintset
{

template<class DataTypes>
void StopperConstraint<DataTypes>::init()
{
    this->mstate = dynamic_cast<MechanicalState*>(this->getContext()->getMechanicalState());
    assert(this->mstate);

    helper::WriteAccessor<Data<VecCoord> > xData = *this->mstate->write(core::VecCoordId::position());
    VecCoord& x = xData.wref();
    if (x[index.getValue()].x() < min.getValue())
        x[index.getValue()].x() = (Real) min.getValue();
    if (x[index.getValue()].x() > max.getValue())
        x[index.getValue()].x() = (Real) max.getValue();
}

template<class DataTypes>
void StopperConstraint<DataTypes>::buildConstraintMatrix(const core::ConstraintParams* /*cParams*/ /* PARAMS FIRST */, DataMatrixDeriv &c_d, unsigned int &cIndex, const DataVecCoord &/*x*/)
{
    cid = cIndex;
//	assert(mstate);

    MatrixDeriv& c = *c_d.beginEdit();

    MatrixDerivRowIterator c_it = c.writeLine(cid);
    c_it.setCol(index.getValue(), Coord(1));

    cIndex++;
    c_d.endEdit();
}

template<class DataTypes>
void StopperConstraint<DataTypes>::getConstraintViolation(const core::ConstraintParams* /*cParams*/ /* PARAMS FIRST */, defaulttype::BaseVector *resV, const DataVecCoord &x, const DataVecDeriv &/*v*/)
{
    resV->set(cid, x.getValue()[index.getValue()][0]);
}

template<class DataTypes>
void StopperConstraint<DataTypes>::getConstraintResolution(std::vector<core::behavior::ConstraintResolution*>& resTab, unsigned int& offset)
{
    for(int i=0; i<1; i++)
        resTab[offset++] = new StopperConstraintResolution1Dof(min.getValue(), max.getValue());
}

} // namespace constraintset

} // namespace component

} // namespace sofa

#endif
