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
#ifndef SOFA_COMPONENT_ENGINE_RIGIDTOQUATENGINE_INL
#define SOFA_COMPONENT_ENGINE_RIGIDTOQUATENGINE_INL

#include <sofa/component/engine/RigidToQuatEngine.h>
#include <sofa/core/visual/VisualParams.h>

namespace sofa
{

namespace component
{

namespace engine
{
template <class DataTypes>
RigidToQuatEngine<DataTypes>::RigidToQuatEngine()
    : f_positions( initData (&f_positions, "positions", "Positions (Vector of 3)") )
    , f_orientations( initData (&f_orientations, "orientations", "Orientations (Quaternion)") )
    , f_rigids( initData (&f_rigids, "rigids", "Rigid (Position + Orientation)") )
{
    //
    addAlias(&f_positions,"position");
    addAlias(&f_orientations,"orientation");
    addAlias(&f_rigids,"rigid");
}

template <class DataTypes>
RigidToQuatEngine<DataTypes>::~RigidToQuatEngine()
{

}

template <class DataTypes>
void RigidToQuatEngine<DataTypes>::init()
{
    addInput(&f_rigids);

    addOutput(&f_positions);
    addOutput(&f_orientations);

    setDirtyValue();
}

template <class DataTypes>
void RigidToQuatEngine<DataTypes>::reinit()
{
    update();
}

template <class DataTypes>
void RigidToQuatEngine<DataTypes>::update()
{
    cleanDirty();

    helper::ReadAccessor< Data< helper::vector<RigidVec3> > > rigids = f_rigids;
    helper::WriteAccessor< Data< helper::vector<Vec3> > > positions = f_positions;
    helper::WriteAccessor< Data< helper::vector<Quat> > > orientations = f_orientations;

    unsigned int sizeRigids = rigids.size();
    positions.resize(sizeRigids);
    orientations.resize(sizeRigids);
    for (unsigned int i=0 ; i< sizeRigids ; i++)
    {
        RigidVec3 r = rigids[i];
        positions[i] = r.getCenter();
        orientations[i] = r.getOrientation();
    }
}

} // namespace engine

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_ENGINE_RIGIDTOQUATENGINE_INL
