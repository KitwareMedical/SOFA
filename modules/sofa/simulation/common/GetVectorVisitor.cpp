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
// C++ Implementation: GetVectorVisitor
//
// Description:
//
//
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <sofa/simulation/common/GetVectorVisitor.h>
#include <sofa/defaulttype/Vec.h>
#include <iostream>
using std::cerr;
using std::endl;

namespace sofa
{

namespace simulation
{


GetVectorVisitor::GetVectorVisitor( const sofa::core::ExecParams* params, defaulttype::BaseVector* vec, core::ConstVecId src )
    : Visitor(params), vec(vec), src(src), offset(0)
{}

GetVectorVisitor::~GetVectorVisitor()
{}

Visitor::Result GetVectorVisitor::processNodeTopDown( simulation::Node* gnode )
{
//    cerr << "GetVectorVisitor::processNodeTopDown, node "<< gnode->getName() << endl;
    if (gnode->mechanicalState != NULL) // independent DOFs
    {
//        cerr << "GetVectorVisitor::processNodeTopDown, node has mechanical state "<< endl;
        gnode->mechanicalState->copyToBaseVector(vec,src,offset);
    }
    return Visitor::RESULT_CONTINUE;
}

} // namespace simulation

} // namespace sofa

