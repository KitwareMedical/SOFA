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
#ifndef SOFA_CORE_BEHAVIOR_BASEANIMATIONLOOP_H
#define SOFA_CORE_BEHAVIOR_BASEANIMATIONLOOP_H

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/ExecParams.h>

namespace sofa
{

namespace core
{

namespace behavior
{

/**
 *  \brief Component responsible for main animation algorithms, managing how
 *  and when mechanical computation happens in one animation step.
 *
 *  This class can optionally replace the default computation scheme of computing
 *  collisions then doing an integration step.
 *
 *  Note that it is in a preliminary stage, hence its fonctionnalities and API will
 *  certainly change soon.
 *
 */

class SOFA_CORE_API BaseAnimationLoop : public virtual objectmodel::BaseObject
{

public:
    SOFA_ABSTRACT_CLASS(BaseAnimationLoop, objectmodel::BaseObject);

protected:
    BaseAnimationLoop();

    virtual ~BaseAnimationLoop();

    /// Stores starting time of the simulation
    double m_resetTime;

    /// Save the initial state for later uses in reset()
    virtual void storeResetState();

public:
    /// Main computation method.
    ///
    /// Specify and execute all computations for computing a timestep, such
    /// as one or more collisions and integrations stages.
    virtual void step(const core::ExecParams* params /* PARAMS FIRST =ExecParams::defaultInstance()*/, double dt) = 0;

    /// Returns starting time of the simulation
    double getResetTime() const;
};

} // namespace behavior

} // namespace core

} // namespace sofa

#endif /* SOFA_CORE_BEHAVIOR_BASEANIMATIONLOOP_H */
