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
#ifndef SOFA_CORE_BEHAVIOR_SingleMatrixAccessor_H
#define SOFA_CORE_BEHAVIOR_SingleMatrixAccessor_H

#include <sofa/component/component.h>
#include <sofa/core/behavior/MultiMatrixAccessor.h>
#include <sofa/defaulttype/BaseMatrix.h>

namespace sofa
{

namespace component
{

namespace linearsolver
{

/** Special case to access a single square matrix.
*/
class SOFA_BASE_LINEAR_SOLVER_API SingleMatrixAccessor : public core::behavior::MultiMatrixAccessor
{
public:
    typedef defaulttype::BaseMatrix BaseMatrix;

    SingleMatrixAccessor( BaseMatrix* m=0 ) { setMatrix(m); }
    virtual ~SingleMatrixAccessor();

    void setMatrix( BaseMatrix* m );
    BaseMatrix* getMatrix() { return matrix; }
    const BaseMatrix* getMatrix() const { return matrix; }


    virtual int getGlobalDimension() const { return matrix->rowSize(); }
    virtual int getGlobalOffset(const core::behavior::BaseMechanicalState*) const { return 0; }
    virtual MatrixRef getMatrix(const core::behavior::BaseMechanicalState*) const;


    virtual InteractionMatrixRef getMatrix(const core::behavior::BaseMechanicalState* mstate1, const core::behavior::BaseMechanicalState* mstate2) const;

protected:
    BaseMatrix* matrix;   ///< The single matrix
    MatrixRef matRef; ///< The accessor to the single matrix

};

} // namespace behavior

} // namespace core

} // namespace sofa

#endif
