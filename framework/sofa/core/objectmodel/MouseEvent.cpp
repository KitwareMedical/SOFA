/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *
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
*                              SOFA :: Framework                              *
*                                                                             *
* Authors: The SOFA Team (see Authors.txt)                                    *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/

// Author: Pierre-Jean Bensoussan, Digtal Trainers, (C) 2008

#include <sofa/core/objectmodel/MouseEvent.h>

namespace sofa
{

namespace core
{

namespace objectmodel
{


MouseEvent::MouseEvent(State state, int wheelDelta)
    : sofa::core::objectmodel::Event()
    , m_state(state)
    , m_wheelDelta(wheelDelta)
{
    m_posX = 0;
    m_posY = 0;
}



MouseEvent::MouseEvent(State state, int posX, int posY)
    : sofa::core::objectmodel::Event()
    , m_state(state)
    , m_posX(posX)
    , m_posY(posY)
{

}



MouseEvent::~MouseEvent()
{

}

} // namespace tree

} // namespace simulation

} // namespace sofa
