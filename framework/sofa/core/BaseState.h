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
#ifndef SOFA_CORE_BASESTATE_H
#define SOFA_CORE_BASESTATE_H

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/VecId.h>

namespace sofa
{

namespace core
{

/**
 *  \brief Component storing position and velocity vectors.
 *
 *  This class define the interface of components used as source and
 *  destination of regular (non mechanical) mapping. It is then specialized as
 *  MechanicalState (storing other mechanical data) or MappedModel (if no
 *  mechanical data is used, such as for VisualModel).
 */
class SOFA_CORE_API BaseState : public virtual objectmodel::BaseObject
{
public:
    SOFA_ABSTRACT_CLASS(BaseState, objectmodel::BaseObject);
protected:
    virtual ~BaseState() { }
public:
    /// Current size of all stored vectors
    virtual int getSize() const = 0;

    /// Resize all stored vector
    virtual void resize(int vsize) = 0;

    /// @name BaseData vectors access API based on VecId
    /// @{

    virtual objectmodel::BaseData* baseWrite(VecId v) = 0;
    virtual const objectmodel::BaseData* baseRead(ConstVecId v) const = 0;

    /// @}

};

} // namespace core

} // namespace sofa

#endif
