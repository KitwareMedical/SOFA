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
#ifndef SOFA_COMPONENT_ENGINE_INDICESFROMVALUES_INL
#define SOFA_COMPONENT_ENGINE_INDICESFROMVALUES_INL

#include <sofa/component/engine/IndicesFromValues.h>
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

template <class T>
IndicesFromValues<T>::IndicesFromValues()
    : f_values( initData (&f_values, "values", "input values") )
    , f_global( initData (&f_global, "global", "Global values, in which the input values are searched") )
    , f_indices( initData(&f_indices, "indices","Output indices of the given values, searched in global") )
{
}

template <class T>
IndicesFromValues<T>::~IndicesFromValues()
{
}

template <class T>
void IndicesFromValues<T>::init()
{
    addInput(&f_values);
    addInput(&f_global);
    addOutput(&f_indices);
    setDirtyValue();
}

template <class T>
void IndicesFromValues<T>::reinit()
{
    update();
}

template <class T>
void IndicesFromValues<T>::update()
{
    cleanDirty();
    helper::ReadAccessor<Data<VecValue> > global = f_global;
    helper::ReadAccessor<Data<VecValue> > values = f_values;
    helper::WriteAccessor<Data<VecIndex> > indices = f_indices;

    indices.clear();
    indices.reserve(values.size());
    for (unsigned int i=0; i<values.size(); ++i)
    {
        const Value v = values[i];
        int index=-1;
        for (unsigned int j=0; j<global.size(); ++j)
        {
            //if (global[j] == v)
            /// @TODO: add operator== to helper::fixed_array and defaulttype::RididCoord/Deriv
            if (!(global[j] < v) && !(v < global[j]))
            {
                index = j;
                break;
            }
        }
        if (index >= 0)
            indices.push_back(index);
        else
            serr << "Input value " << i <<" not found : " << v << sendl;
    }
}

} // namespace engine

} // namespace component

} // namespace sofa

#endif
