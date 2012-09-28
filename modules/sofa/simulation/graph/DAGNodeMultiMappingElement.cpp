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
#include <sofa/simulation/graph/DAGNodeMultiMappingElement.h>
#include <sofa/core/BaseMapping.h>
#include <sofa/simulation/common/Node.h>
#include <sofa/helper/Factory.h>

namespace sofa
{

namespace simulation
{

namespace graph
{
DAGNodeMultiMappingElement::DAGNodeMultiMappingElement(const std::string &name,
        const std::string &type, BaseElement *parent /*= 0*/)
    :BaseMultiMappingElement(name,type,parent)
{

}

void DAGNodeMultiMappingElement::updateSceneGraph(
    sofa::core::BaseMapping* multiMapping,
    const helper::vector<simulation::Node*>& /*ancestorInputs*/,
    helper::vector<simulation::Node*>& otherInputs,
    helper::vector<simulation::Node*>& /*outputs*/)
{

    helper::vector<simulation::Node*>::const_iterator it;
    for( it = otherInputs.begin(); it != otherInputs.end(); ++it)
    {
        multiMapping->serr << "Node: " << (*it)->getName() << " does not belong to "
                << multiMapping->getContext()->getName() << "ancestors" << multiMapping->sendl;
    }
}


SOFA_DECL_CLASS(DAGNodeMultiMappingElement)

helper::Creator<sofa::simulation::xml::BaseElement::NodeFactory, DAGNodeMultiMappingElement> DAGNodeMultiMappingClass("DAGNodeMultiMapping");

const char* DAGNodeMultiMappingElement::getClass() const
{
    return DAGNodeMultiMappingClass.c_str();
}


} // namespace graph

} // namespace simulation

} // namespace sofa
