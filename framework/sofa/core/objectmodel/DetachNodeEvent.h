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
#ifndef SOFA_CORE_OBJECTMODEL_DETACHNODEEVENT_H
#define SOFA_CORE_OBJECTMODEL_DETACHNODEEVENT_H

#include <sofa/core/objectmodel/Event.h>
#include <sofa/core/objectmodel/BaseNode.h>
#include <sofa/core/objectmodel/BaseObject.h>

namespace sofa
{

namespace core
{

namespace objectmodel
{

/**
 *  Event indicating that a child node is being detached from the scene.
 *  Any reference to ony of its descendant (such as active contacts) should be removed.
*/
class SOFA_CORE_API DetachNodeEvent : public Event
{
public:
    DetachNodeEvent( BaseNode* n );

    ~DetachNodeEvent();

    BaseNode* getNode() const;

    bool contains(BaseNode* n) const;

    bool contains(BaseObject* o) const;

    virtual const char* getClassName() const { return "DetachNodeEvent"; }
protected:
    BaseNode* node;
};

} // namespace objectmodel

} // namespace core

} // namespace sofa

#endif
