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
#ifndef SOFA_COMPONENT_ENGINE_POINTSFROMINDICES_INL
#define SOFA_COMPONENT_ENGINE_POINTSFROMINDICES_INL

#if !defined(__GNUC__) || (__GNUC__ > 3 || (_GNUC__ == 3 && __GNUC_MINOR__ > 3))
#pragma once
#endif

#include <sofa/component/engine/PointsFromIndices.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/gl/template.h>

namespace sofa
{

namespace component
{

namespace engine
{

using namespace sofa::helper;
using namespace sofa::defaulttype;
using namespace core::objectmodel;

template <class DataTypes>
PointsFromIndices<DataTypes>::PointsFromIndices()
    : f_X( initData (&f_X, "position", "Position coordinates of the degrees of freedom") )
    , f_indices( initData(&f_indices,"indices","Indices of the points") )
    , f_indices_position( initData (&f_indices_position, "indices_position", "Coordinates of the points contained in indices"))
{
}

template <class DataTypes>
void PointsFromIndices<DataTypes>::init()
{
    if (f_X.getValue().empty())
    {
        MechanicalState<DataTypes>* mstate;
        this->getContext()->get(mstate);
        if (mstate)
        {
            BaseData* parent = mstate->findField("position");
            if (parent)
            {
                f_X.setParent(parent);
                f_X.setReadOnly(true);
            }
        }
    }
    addInput(&f_X);
    addInput(&f_indices);
    addOutput(&f_indices_position);
    setDirtyValue();
}

template <class DataTypes>
void PointsFromIndices<DataTypes>::reinit()
{
    update();
}

template <class DataTypes>
bool PointsFromIndices<DataTypes>::contains(VecCoord& v, Coord c)
{
    for( unsigned i=0; i<v.size(); ++i )
    {
        if (c == v[i]) return true;
    }
    return false;
}

template <class DataTypes>
void PointsFromIndices<DataTypes>::update()
{
    cleanDirty();
    VecCoord& indices_position = *(f_indices_position.beginEdit());
    const SetIndex& indices = f_indices.getValue();
    const VecCoord& x = f_X.getValue();

    if(!x.empty())
    {
        indices_position.clear();
        //indices_position.resize(indices.size());

        for( unsigned i=0; i<indices.size(); ++i )
        {
            if (!contains(indices_position, x[indices[i]]))
                indices_position.push_back(x[indices[i]]);
        }

    }

    f_indices_position.endEdit();
}

} // namespace engine

} // namespace component

} // namespace sofa

#endif
