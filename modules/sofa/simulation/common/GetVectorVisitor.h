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
// C++ Interface: GetVectorVisitor
//
// Description:
//
//
// Author: Francois Faure, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SOFA_SIMULATION_GetVectorVisitor_H
#define SOFA_SIMULATION_GetVectorVisitor_H

#if !defined(__GNUC__) || (__GNUC__ > 3 || (_GNUC__ == 3 && __GNUC_MINOR__ > 3))
#pragma once
#endif

#include <sofa/simulation/common/Visitor.h>
#include <sofa/core/MultiVecId.h>
#include <sofa/defaulttype/BaseVector.h>
#include <Eigen/Dense>


namespace sofa
{

namespace simulation
{

/** Copy a given MultiVector (generally spread across the MechanicalStates) to a BaseVector
    Only the independent DOFs are used.
    Francois Faure, 2013
*/
class SOFA_SIMULATION_COMMON_API GetVectorVisitor: public Visitor
{
public:
//    typedef Eigen::Matrix<SReal, Eigen::Dynamic, 1> Vector;
    typedef defaulttype::BaseVector Vector;
    GetVectorVisitor( const sofa::core::ExecParams* params /* PARAMS FIRST */, Vector* vec, core::ConstVecId src );
    virtual ~GetVectorVisitor();

    virtual Result processNodeTopDown( simulation::Node*  );
    virtual const char* getClassName() const { return "GetVectorVisitor"; }

protected:
    Vector* vec;
    core::ConstVecId src;
    unsigned offset;
};

} // namespace simulation
} // namespace sofa

#endif
