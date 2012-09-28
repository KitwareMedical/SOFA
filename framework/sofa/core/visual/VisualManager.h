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
#ifndef SOFA_CORE_VISUAL_VISUALMANAGER_H
#define SOFA_CORE_VISUAL_VISUALMANAGER_H

#include <sofa/core/visual/VisualModel.h>
#include <sofa/core/visual/VisualParams.h>

namespace sofa
{

namespace core
{

namespace visual
{

class VisualManager : public virtual VisualModel
{
public:
    SOFA_CLASS(VisualManager, VisualModel);
protected:
    /// Destructor
    virtual ~VisualManager() { }
public:
    /**
     *  \brief Called before rendering the scene
     */
    virtual void preDrawScene(VisualParams* /*vparams*/) {}

    /**
     *  \brief Called after rendering the scene
     */
    virtual void postDrawScene(VisualParams* /*vparams*/) {}

    /**
     *  \brief Called instead of rendering the scene
     *
     *  Return true if this object actually did the rendering, or false if it wasn't done.
     */
    virtual bool drawScene(VisualParams* /*vparams*/) { return false; }
};

} // namespace visual

} // namespace core

} // namespace sofa

#endif // SOFA_CORE_VISUAL_VISUALMANAGER_H
