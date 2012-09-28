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
#ifndef SOFA_SIMULATION_COMMON_VISITOREXECUTE_H
#define SOFA_SIMULATION_COMMON_VISITOREXECUTE_H


#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/simulation/common/Visitor.h>
#include <sofa/simulation/common/MechanicalVisitor.h>

namespace sofa
{
namespace simulation
{
namespace common
{

struct VisitorExecuteFunc
{
protected:
    core::objectmodel::BaseContext& ctx;
public:
    VisitorExecuteFunc(core::objectmodel::BaseContext& ctx):ctx(ctx) {};

    template< class Visitor >
    void operator()(Visitor* pv)
    {
        prepareVisitor(pv);
        pv->execute(&ctx);
    }
    template< class Visitor >
    void operator()(Visitor v)
    {
        prepareVisitor(&v);
        v.execute(&ctx);
    }
protected:
    void prepareVisitor( sofa::simulation::Visitor* v)
    {
        v->setTags( ctx.getTags() );
    }
    void prepareVisitor( sofa::simulation::BaseMechanicalVisitor* mv)
    {
        prepareVisitor( (sofa::simulation::Visitor*)mv );
    }
};
}
}
}

#endif // SOFA_SIMULATION_COMMON_VISITOREXECUTE_H
