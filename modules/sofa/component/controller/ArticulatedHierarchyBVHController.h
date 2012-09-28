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
// C++ Implementation: ArticulatedHierarchyBVHController
//
// Description:
//
//
// Author: The SOFA team </www.sofa-framework.org>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SOFA_COMPONENT_CONTROLLER_ARTICULATEDHIERARCHYBVHCONTROLLER_H
#define SOFA_COMPONENT_CONTROLLER_ARTICULATEDHIERARCHYBVHCONTROLLER_H

#include <sofa/component/controller/ArticulatedHierarchyController.h>
#include <sofa/helper/io/bvh/BVHLoader.h>

namespace sofa
{

namespace component
{

namespace controller
{

/**
 * @brief ArticulatedHierarchyController Class.
 *
 * Implements a handler that controls the values of the
 * articulations of an articulated hierarchy container.
 * .bvh files are controlling the value.
 */
class SOFA_USER_INTERACTION_API ArticulatedHierarchyBVHController : public ArticulatedHierarchyController
{
public:
    SOFA_CLASS(ArticulatedHierarchyBVHController,ArticulatedHierarchyController);
protected:
    /**
    * @brief Default Constructor.
     */
    ArticulatedHierarchyBVHController()
        : useExternalTime( initData(&useExternalTime, false, "useExternalTime", "use the external time line"))
        , externalTime( initData(&externalTime, 0.0, "externalTime", " value of the External Time") )
    {
        this->f_listening.setValue(true);
    };

    /**
     * @brief Default Destructor.
     */
    virtual ~ArticulatedHierarchyBVHController() {};
public:
    /**
     * @brief Init method called during the scene graph initialization.
     */
    virtual void init();

    /**
     * @brief Reset to initial state
     */
    virtual void reset();

    /**
     * @brief Apply the controller current modifications to its controled component.
     */
    virtual void applyController(void);

protected:
    Data< bool > useExternalTime;
    Data< double > externalTime;
    ArtCenterVec m_artCenterVec; ///< List of ArticulationCenters controlled by the controller.
    ArticulatedHierarchyContainer* ahc;
    int frame;
    int n;
};

} // namespace controller

} // namespace component

} // namespace sofa

#endif //SOFA_COMPONENT_CONTROLLER_ARTICULATEDHIERARCHYBVHCONTROLLER_H
