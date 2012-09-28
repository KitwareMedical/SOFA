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
#ifndef SOFA_COMPONENT_COLLISION_CUTTINGPOINT_H
#define SOFA_COMPONENT_COLLISION_CUTTINGPOINT_H

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/collision/Intersection.h>
#include <sofa/core/collision/NarrowPhaseDetection.h>
#include <sofa/component/collision/TriangleModel.h>
#include <sofa/component/topology/TopologySubsetData.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/simulation/common/AnimateBeginEvent.h>
#include <sofa/simulation/common/AnimateEndEvent.h>

#include <sofa/core/behavior/BaseController.h>
#include <fstream>

#include <sofa/component/topology/TriangleSetTopologyAlgorithms.h>
#include <sofa/component/topology/TriangleSetGeometryAlgorithms.h>

namespace sofa
{

namespace component
{

namespace collision
{

class SOFA_USER_INTERACTION_API CuttingPoint : public virtual core::objectmodel::BaseObject
{
public:
    SOFA_CLASS(CuttingPoint,sofa::core::objectmodel::BaseObject);

    typedef TriangleModel::DataTypes DataTypes;
    typedef DataTypes::VecCoord VecCoord;
    typedef DataTypes::Coord Coord;
    typedef DataTypes::Real Real;
    typedef helper::vector<unsigned int> SetIndexArray;
    typedef sofa::component::topology::PointSubsetData< SetIndexArray > SetIndex;
    Coord lastPos, pos, newPos;
    int toolElemIndex;
public:

    SetIndex f_points;
    Data<int> id;
    Data<int> prevID;
    Data<int> nextID;
    Data<bool> cutInProgress;
    core::behavior::MechanicalState<DataTypes>* mstate;
    sofa::core::topology::BaseMeshTopology* topology;
protected:
    CuttingPoint(int myid = -1)
        : toolElemIndex(-1)
        , f_points( initData(&f_points, "points", "cut-related points") )
        , id( initData(&id, myid, "id", "ID of this cutting point") )
        , prevID( initData(&prevID, -1, "prevID", "ID of the previous cutting point (if this point is in the middle or at the end of a cut line)") )
        , nextID( initData(&nextID, -1, "nextID", "ID of the next cutting point (if this point is in the middle or at the start of a cut line)") )
        , cutInProgress( initData(&cutInProgress, false, "cutInProgress", "True if this point is currently being cut") )
        , mstate(NULL)
    {
    }
    virtual ~CuttingPoint()
    {
    }
public:
    virtual void init();

    const SetIndexArray& getPointSubSet() const {return f_points.getValue();}

    void setPointSubSet(SetIndexArray pSub);
    void setPoint(unsigned int p);

    bool isValid() const { return !f_points.getValue().empty(); }
    bool canBeCut() const { return isValid() && (prevID.getValue() == -1 || nextID.getValue() == -1); }
    bool keepAfterCut() const { return isValid() && (prevID.getValue() == -1 && nextID.getValue() == -1); }
    bool canBeFractured() const { return canBeCut() && !cutInProgress.getValue(); }

    virtual void draw(const core::visual::VisualParams* vparams);
};

} //collision
} //component
} // sofa

#endif // SOFA_COMPONENT_COLLISION_CUTTINGPOINT_H
