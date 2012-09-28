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
#ifndef SOFA_COMPONENT_COLLISION_SUBSETCONTACTMAPPER_INL
#define SOFA_COMPONENT_COLLISION_SUBSETCONTACTMAPPER_INL

#include <sofa/component/collision/SubsetContactMapper.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/simulation/common/Node.h>
#include <sofa/simulation/common/Simulation.h>
#include <sofa/simulation/common/DeleteVisitor.h>
#include <iostream>

namespace sofa
{

namespace component
{

namespace collision
{

template < class TCollisionModel, class DataTypes >
void SubsetContactMapper<TCollisionModel,DataTypes>::cleanup()
{
    if (child!=NULL)
    {
        child->detachFromGraph();
        child->execute<simulation::DeleteVisitor>(sofa::core::ExecParams::defaultInstance());
        child.reset();
    }
}

template < class TCollisionModel, class DataTypes >
typename SubsetContactMapper<TCollisionModel,DataTypes>::MMechanicalState* SubsetContactMapper<TCollisionModel,DataTypes>::createMapping(const char* name)
{
    if (model==NULL) return NULL;
    InMechanicalState* instate = model->getMechanicalState();
    if (instate!=NULL)
    {
        simulation::Node* parent = dynamic_cast<simulation::Node*>(instate->getContext());
        if (parent==NULL)
        {
            std::cerr << "ERROR: SubsetContactMapper only works for scenegraph scenes.\n";
            return NULL;
        }
        child = parent->createChild(name);
        outmodel = sofa::core::objectmodel::New<MMechanicalObject>(); child->addObject(outmodel);
        outmodel->useMask.setValue(true);
        mapping = sofa::core::objectmodel::New<MMapping>();
        child->addObject(mapping);
        mapping->setModels(instate, outmodel.get());
    }
    else
    {
        simulation::Node* parent = dynamic_cast<simulation::Node*>(model->getContext());
        if (parent==NULL)
        {
            std::cerr << "ERROR: SubsetContactMapper only works for scenegraph scenes.\n";
            return NULL;
        }
        child = parent->createChild(name);
        outmodel = sofa::core::objectmodel::New<MMechanicalObject>(); child->addObject(outmodel);
        outmodel->useMask.setValue(true);
        mapping = NULL;
    }
    return outmodel.get();
}
} // namespace collision

} // namespace component

} // namespace sofa

#endif /* SOFA_COMPONENT_COLLISION_SUBSETCONTACTMAPPER_INL */
