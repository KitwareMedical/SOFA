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
#ifndef SOFA_COMPONENT_INTERACTIONFORCEFIELD_QUADBENDINGSPRINGS_H
#define SOFA_COMPONENT_INTERACTIONFORCEFIELD_QUADBENDINGSPRINGS_H

#include <sofa/component/interactionforcefield/StiffSpringForceField.h>
#include <map>
#include <set>

namespace sofa
{

namespace component
{

namespace interactionforcefield
{

/**
Bending springs added between vertices of quads sharing a common edge.
The springs connect the vertices not belonging to the common edge. It compresses when the surface bends along the common edge.
*/
template<class DataTypes>
class QuadBendingSprings : public sofa::component::interactionforcefield::StiffSpringForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(QuadBendingSprings, DataTypes), SOFA_TEMPLATE(StiffSpringForceField, DataTypes));

    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef core::behavior::MechanicalState<DataTypes> MechanicalState;
    /// optional range of local DOF indices. Any computation involving only indices outside of this range are discarded (useful for parallelization using mesh partitionning)
    Data< defaulttype::Vec<2,int> > localRange;
protected:
    QuadBendingSprings();

    ~QuadBendingSprings();
public:
    /// Searches triangle topology and creates the bending springs
    virtual void init();

    virtual void draw()
    {
    }
    void setObject1(MechanicalState* object1) {this->mstate1=object1;}
    void setObject2(MechanicalState* object2) {this->mstate2=object2;}

protected:
    typedef std::pair<unsigned,unsigned> IndexPair;
    void addSpring( unsigned, unsigned, std::set<IndexPair>& );
    void registerEdge( IndexPair, IndexPair, std::map<IndexPair, IndexPair>&, std::set<IndexPair>&);

};

} // namespace interactionforcefield

} // namespace component

} // namespace sofa

#endif  /* SOFA_COMPONENT_INTERACTIONFORCEFIELD_QUADBENDINGSPRINGS_H */
