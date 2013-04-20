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
// C++ Interface: GetAssembledSizeVisitor
//
// Description:
//
//
// Author: Francois Faure, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SOFA_SIMULATION_GetAssembledSizeVisitor_H
#define SOFA_SIMULATION_GetAssembledSizeVisitor_H

#if !defined(__GNUC__) || (__GNUC__ > 3 || (_GNUC__ == 3 && __GNUC_MINOR__ > 3))
#pragma once
#endif

#include <sofa/simulation/common/Visitor.h>
#include <sofa/core/MultiVecId.h>
#include <sofa/defaulttype/BaseVector.h>


namespace sofa
{

namespace simulation
{

/** Compute the size of the assembled position vector and velocity vector.
  Only the independent DOFs are considered.
  The two values may be different, such as for rigid objects.
    Francois Faure, 2013
*/
class SOFA_SIMULATION_COMMON_API GetAssembledSizeVisitor: public Visitor
{
public:
    GetAssembledSizeVisitor( const sofa::core::ExecParams* params /* PARAMS FIRST */=core::MechanicalParams::defaultInstance() );
    virtual ~GetAssembledSizeVisitor();

    virtual Result processNodeTopDown( simulation::Node*  );
    virtual const char* getClassName() const { return "GetAssembledSizeVisitor"; }

    unsigned positionSize() const { return xsize; }
    unsigned velocitySize() const { return vsize; }

protected:
    unsigned xsize;
    unsigned vsize;
};

} // namespace simulation
} // namespace sofa

#endif
