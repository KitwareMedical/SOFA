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
// C++ Implementation: GetAssembledSizeVisitor
//
// Description:
//
//
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <sofa/simulation/common/GetAssembledSizeVisitor.h>
#include <sofa/defaulttype/Vec.h>

namespace sofa
{

namespace simulation
{


GetAssembledSizeVisitor::GetAssembledSizeVisitor( const sofa::core::ExecParams* params )
    : Visitor(params)
    , xsize(0)
    , vsize(0)
{}

GetAssembledSizeVisitor::~GetAssembledSizeVisitor()
{}

Visitor::Result GetAssembledSizeVisitor::processNodeTopDown( simulation::Node* gnode )
{
    if (gnode->mechanicalState != NULL) // independent DOFs
    {
        xsize += gnode->mechanicalState->getSize() * gnode->mechanicalState->getCoordDimension();
        vsize += gnode->mechanicalState->getMatrixSize();
    }
    return Visitor::RESULT_CONTINUE;
}

} // namespace simulation

} // namespace sofa

