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
#include <sofa/core/topology/TopologyElementHandler.h>

namespace sofa
{

namespace core
{

namespace topology
{

/// Apply swap between indices elements.
template<class TopologyElementType>
void TopologyElementHandler<TopologyElementType>::ApplyTopologyChange(const EIndicesSwap* event)
{
    this->swap(event->index[0], event->index[1]);
}
/// Apply adding elements.
template<class TopologyElementType>
void TopologyElementHandler<TopologyElementType>::ApplyTopologyChange(const EAdded* event)
{
    //this->add(event->getNbAddedElements(), event->getElementArray(),
    //    event->ancestorsList, event->coefs);
    this->add(event->getIndexArray(), event->getElementArray(),
        event->ancestorsList, event->coefs, event->ancestorElems);
}

/// Apply removing elements.
template<class TopologyElementType>
void TopologyElementHandler<TopologyElementType>::ApplyTopologyChange(const ERemoved* event)
{
    this->remove(event->getArray());
}

/// Apply renumbering on elements.
template<class TopologyElementType>
void TopologyElementHandler<TopologyElementType>::ApplyTopologyChange(const ERenumbering* event)
{
    this->renumber(event->getIndexArray());
}

/// Apply moving elements.
template<class TopologyElementType>
void TopologyElementHandler<TopologyElementType>::ApplyTopologyChange(const EMoved* /*event*/)
{
    std::cerr << "ERROR: MOVED topology event not handled on " << ElementInfo::name()
        << " (it should not even exist!)" << std::endl;
}

/// Apply moving elements on points
template<>
void TopologyElementHandler<Point>::ApplyTopologyChange(const EMoved* event)
{
    this->move(event->getIndexArray(), event->ancestorsList, event->baryCoefsList);
}

/// Apply adding function on moved elements.
template<class TopologyElementType>
void TopologyElementHandler<TopologyElementType>::ApplyTopologyChange(const EMoved_Adding* event)
{
    this->addOnMovedPosition(event->getIndexArray(), event->getElementArray());
}

/// Apply adding function on moved point.
template<>
void TopologyElementHandler<Point>::ApplyTopologyChange(const EMoved_Adding* /* event */)
{
    std::cerr << "ERROR: MOVED_ADDING topology event not handled on " << ElementInfo::name()
        << " (it should not even exist!)" << std::endl;
}

/// Apply removing function on moved elements.
template<class TopologyElementType>
void TopologyElementHandler<TopologyElementType>::ApplyTopologyChange(const EMoved_Removing* event)
{
    this->removeOnMovedPosition(event->getIndexArray());
}

/// Apply removing function on moved elements.
template<>
void TopologyElementHandler<Point>::ApplyTopologyChange(const EMoved_Removing* /* event */)
{
    std::cerr << "ERROR: MOVED_REMOVING topology event not handled on " << ElementInfo::name()
        << " (it should not even exist!)" << std::endl;
}

template class SOFA_CORE_API TopologyElementHandler<Point>;
template class SOFA_CORE_API TopologyElementHandler<Edge>;
template class SOFA_CORE_API TopologyElementHandler<Triangle>;
template class SOFA_CORE_API TopologyElementHandler<Quad>;
template class SOFA_CORE_API TopologyElementHandler<Tetrahedron>;
template class SOFA_CORE_API TopologyElementHandler<Hexahedron>;

} // namespace topology

} // namespace core

} // namespace sofa
