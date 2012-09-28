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
#ifndef SOFA_CORE_BEHAVIORMODEL_H
#define SOFA_CORE_BEHAVIORMODEL_H

#include <sofa/core/objectmodel/BaseObject.h>

namespace sofa
{

namespace core
{

/**
 *  \brief Abstract Interface of components defining the behavior of a simulated object.
 *
 *  This Interface is used by "black-box" objects (such as some fluid algorithms)
 *  that are present in a SOFA simulation but which do not use the internal
 *  behavior components (MechanicalState, ForceField, etc) defined in the
 *  sofa::core::behavior namespace.
 *
 *  A BehaviorModel simply has to implement the updatePosition method
 *  to compute a new simulation step.
 *
 */
class BehaviorModel : public virtual sofa::core::objectmodel::BaseObject
{
public:
    SOFA_ABSTRACT_CLASS(BehaviorModel, objectmodel::BaseObject);
protected:
    /// Destructor
    virtual ~BehaviorModel() {}
public:
    /// Computation of a new simulation step.
    virtual void updatePosition(double dt) = 0;

    virtual bool addBBox(double* /*minBBox*/, double* /*maxBBox*/)
    {
        return false;
    }
};

} // namespace core

} // namespace sofa

#endif
