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
//
// C++ Implementation: LineBendingSprings
//
// Description:
//
//
// Author: The SOFA team </www.sofa-framework.org>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SOFA_COMPONENT_INTERACTIONFORCEFIELD_LINEBENDINGSPRINGS_INL
#define SOFA_COMPONENT_INTERACTIONFORCEFIELD_LINEBENDINGSPRINGS_INL

#include <sofa/component/interactionforcefield/LineBendingSprings.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <iostream>

namespace sofa
{

namespace component
{

namespace interactionforcefield
{

using namespace core::behavior;

template<class DataTypes>
LineBendingSprings<DataTypes>::LineBendingSprings()
{
    //serr<<"LineBendingSprings<DataTypes>::LineBendingSprings"<<sendl;
}


template<class DataTypes>
LineBendingSprings<DataTypes>::~LineBendingSprings()
{}

template<class DataTypes>
void LineBendingSprings<DataTypes>::addSpring( unsigned a, unsigned b )
{
    const VecCoord& x = *this->mstate1->getX();
    Real s = (Real)this->ks.getValue();
    Real d = (Real)this->kd.getValue();
    Real l = (x[a]-x[b]).norm();
    this->SpringForceField<DataTypes>::addSpring(a,b, s, d, l );
    //sout<<"=================================LineBendingSprings<DataTypes>::addSpring "<<a<<", "<<b<<sendl;
}

template<class DataTypes>
void LineBendingSprings<DataTypes>::registerLine( unsigned a, unsigned b, std::map<Index, unsigned>& ptMap)
{
    //sout<<"=================================LineBendingSprings<DataTypes>::registerLine "<<a<<", "<<b<<sendl;
    using namespace std;
    {
        if( ptMap.find( a ) != ptMap.end() )
        {
            // create a spring between the opposite
            this->addSpring(b,ptMap[a]);
        }
        ptMap[a] = b;
        ptMap[b] = a;
    }

}



template<class DataTypes>
void LineBendingSprings<DataTypes>::init()
{
    this->mstate1 = this->mstate2 = dynamic_cast<core::behavior::MechanicalState<DataTypes>*>( this->getContext()->getMechanicalState() );
    StiffSpringForceField<DataTypes>::clear();

    // Set the bending springs

    std::map< Index, unsigned > ptMap;
    sofa::core::topology::BaseMeshTopology* topology = this->getContext()->getMeshTopology();
    assert( topology );

    const sofa::core::topology::BaseMeshTopology::SeqEdges& Edges = topology->getEdges();
    //sout<<"==================================LineBendingSprings<DataTypes>::init(), Lines size = "<<Lines.size()<<sendl;
    for( unsigned i= 0; i<Edges.size(); ++i )
    {
        const sofa::core::topology::BaseMeshTopology::Edge& edge = Edges[i];
        {
            registerLine( edge[0], edge[1] , ptMap );
        }

    }

    // init the parent class
    StiffSpringForceField<DataTypes>::init();

}


} // namespace interactionforcefield

} // namespace component

} // namespace sofa

#endif  /* SOFA_COMPONENT_INTERACTIONFORCEFIELD_LINEBENDINGSPRINGS_INL */
