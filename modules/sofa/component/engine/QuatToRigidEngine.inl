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
#ifndef SOFA_COMPONENT_ENGINE_QUATTORIGIDENGINE_INL
#define SOFA_COMPONENT_ENGINE_QUATTORIGIDENGINE_INL

#include <sofa/component/engine/QuatToRigidEngine.h>
#include <sofa/core/visual/VisualParams.h>

namespace sofa
{

namespace component
{

namespace engine
{
template <class DataTypes>
QuatToRigidEngine<DataTypes>::QuatToRigidEngine()
    : f_positions( initData (&f_positions, "positions", "Positions (Vector of 3)") )
    , f_orientations( initData (&f_orientations, "orientations", "Orientations (Quaternion)") )
    , f_colinearPositions( initData (&f_colinearPositions, "colinearPositions", "Optional positions to restrict output to be colinear in the quaternion Z direction") )
    , f_rigids( initData (&f_rigids, "rigids", "Rigid (Position + Orientation)") )
{
    //
    addAlias(&f_positions,"position");
    addAlias(&f_orientations,"orientation");
    addAlias(&f_colinearPositions,"colinearPosition");
    addAlias(&f_rigids,"rigid");
}

template <class DataTypes>
QuatToRigidEngine<DataTypes>::~QuatToRigidEngine()
{

}

template <class DataTypes>
void QuatToRigidEngine<DataTypes>::init()
{
    addInput(&f_positions);
    addInput(&f_orientations);
    addInput(&f_colinearPositions);

    addOutput(&f_rigids);

    setDirtyValue();
}

template <class DataTypes>
void QuatToRigidEngine<DataTypes>::reinit()
{
    update();
}

template <class DataTypes>
void QuatToRigidEngine<DataTypes>::update()
{
    cleanDirty();

    const helper::vector<Vec3>& positions = f_positions.getValue();
    const helper::vector<Quat>& orientations = f_orientations.getValue();
    const helper::vector<Vec3>& colinearPositions = f_colinearPositions.getValue();

    helper::vector<RigidVec3>& rigids = *(f_rigids.beginEdit());

    unsigned int sizeRigids = positions.size();

    if(positions.size() != orientations.size())
    {
        serr << "Warnings : size of positions and orientations are not equal" << sendl;
        sizeRigids = ( positions.size() > orientations.size() ) ?  orientations.size() :  positions.size() ;
    }

    rigids.clear();
    for (unsigned int i=0 ; i< sizeRigids ; i++)
    {
        Vec3 pos = positions[i];
        Quat q = orientations[i];
        if (i < colinearPositions.size())
        {
            Vec3 colP = colinearPositions[i];
            Vec3 dir = q.rotate(Vec3(0,0,1));
            pos = colP + dir*(dot(dir,(pos-colP))/dir.norm());
        }
        RigidVec3 r(pos, q);
        rigids.push_back(r);
    }

    f_rigids.endEdit();
}

} // namespace engine

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_ENGINE_QUATTORIGIDENGINE_INL
