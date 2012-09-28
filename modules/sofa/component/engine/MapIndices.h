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
#ifndef SOFA_COMPONENT_ENGINE_MAPINDICES_H
#define SOFA_COMPONENT_ENGINE_MAPINDICES_H

#include <sofa/core/DataEngine.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/topology/BaseMeshTopology.h>
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
 * This class apply a permutation to a set of indices
 */
template <class T>
class MapIndices : public core::DataEngine
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(MapIndices,T),core::DataEngine);
    typedef T Value;
    typedef sofa::helper::vector<T> VecValue;
    typedef unsigned int Index;
    typedef sofa::helper::vector<Index> VecIndex;
    typedef std::map<Index, Index> MapIndex;
protected:
    MapIndices();

    virtual ~MapIndices();
public:
    void init();

    void reinit();

    void update();

    core::objectmodel::Data<VecValue> f_in;
    core::objectmodel::Data<VecIndex> f_indices;
    core::objectmodel::Data<VecValue> f_out;
    core::objectmodel::Data<std::string> f_outStr;
    core::objectmodel::Data<bool> f_transpose;

    template<class V>
    void applyIndex(V& v, const MapIndex& m)
    {
        typename MapIndex::const_iterator it = m.find(v);
        if (it != m.end())
            v = it->second;
    }

    void apply(Value& v, const MapIndex& m);
};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_ENGINE_MAPINDICES_CPP)
extern template class SOFA_ENGINE_API MapIndices<int>;
extern template class SOFA_ENGINE_API MapIndices<unsigned int>;
extern template class SOFA_ENGINE_API MapIndices< helper::fixed_array<unsigned int, 2> >;
extern template class SOFA_ENGINE_API MapIndices< helper::fixed_array<unsigned int, 3> >;
extern template class SOFA_ENGINE_API MapIndices< helper::fixed_array<unsigned int, 4> >;
extern template class SOFA_ENGINE_API MapIndices< helper::fixed_array<unsigned int, 8> >;
#endif

} // namespace engine

} // namespace component

} // namespace sofa

#endif
