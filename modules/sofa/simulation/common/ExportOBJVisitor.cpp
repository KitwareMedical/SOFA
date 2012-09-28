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
#include <sofa/simulation/common/ExportOBJVisitor.h>
#include <sofa/helper/system/config.h>
#include <sofa/helper/Factory.h>
#include <sofa/simulation/common/Node.h>
#include <sofa/core/objectmodel/BaseContext.h>
namespace sofa
{

namespace simulation
{


ExportOBJVisitor::ExportOBJVisitor(const core::ExecParams* params /* PARAMS FIRST */, std::ostream* out)
    : Visitor(params) , out(out), mtl(NULL), ID(0), vindex(0), nindex(0), tindex(0), count(0)
{
}

ExportOBJVisitor::ExportOBJVisitor(const core::ExecParams* params /* PARAMS FIRST */, std::ostream* out,std::ostream* mtl)
    : Visitor(params) , out(out), mtl(mtl), ID(0), vindex(0), nindex(0), tindex(0), count(0)
{
}

ExportOBJVisitor::~ExportOBJVisitor()
{
}

void ExportOBJVisitor::processVisualModel(Node* /*node*/, core::visual::VisualModel* vm)
{
    std::ostringstream oname;
    oname << ++ID << "_" << vm->getName();

    vm->exportOBJ(oname.str(),out,mtl,vindex,nindex,tindex, ++count);
}

simulation::Visitor::Result ExportOBJVisitor::processNodeTopDown(Node* node)
{
    //simulation::Node* node = static_cast<simulation::Node*>(n);
    for_each(this, node, node->visualModel,              &ExportOBJVisitor::processVisualModel);
    count = 0;
    return RESULT_CONTINUE;
}

void ExportOBJVisitor::processNodeBottomUp(Node* /*node*/)
{
}

} // namespace simulation

} // namespace sofa

