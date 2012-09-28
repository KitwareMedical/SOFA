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
*                               SOFA :: Plugins                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_CORE_OBJECTMODEL_ARTRACKEVENT_H
#define SOFA_CORE_OBJECTMODEL_ARTRACKEVENT_H

#include <sofa/core/objectmodel/Event.h>
#include <sofa/defaulttype/Vec3Types.h>
#include <sofa/defaulttype/Quat.h>
#include <sofa/helper/system/config.h>

namespace sofa
{

namespace core
{

namespace objectmodel
{

using namespace sofa::defaulttype;

/**
 * @brief This event notifies about ARTrack device interaction.
 */
class ARTrackEvent : public sofa::core::objectmodel::Event
{
public:

    /**
     * @brief Constructor.
     */
    ARTrackEvent(const Vector3& position, const Quat& orientation, const sofa::helper::fixed_array<double,3>& angles, const sofa::helper::fixed_array<Vector3,3>& fingersPosition);

    /**
     * @brief Destructor.
     */
    virtual ~ARTrackEvent() {}

    const Vector3 getPosition() const;
    const Quat getOrientation() const;
    const sofa::helper::fixed_array<double,3> getAngles() const;
    const Vector3 getFingerposition(const unsigned int i) const;

private:
    Vector3 m_position; ///< ARTrack coordinates in a Vec3d type.
    Quat m_orientation; ///< ARTrack orientation.
    sofa::helper::fixed_array<double,3> m_angles; ///< ARTrack finger angles.
    sofa::helper::fixed_array<Vector3,3> m_fingersPosition; ///< ARTrack fingers position.
};

} // namespace objectmodel

} // namespace core

} // namespace sofa

#endif
