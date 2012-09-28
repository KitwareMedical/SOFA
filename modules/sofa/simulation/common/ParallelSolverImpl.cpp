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
#include "ParallelSolverImpl.h"
#include "ParallelMechanicalVisitor.h"
#include <sofa/simulation/common/MechanicalMatrixVisitor.h>
#include <sofa/simulation/common/MechanicalVPrintVisitor.h>
#include <sofa/simulation/common/VelocityThresholdVisitor.h>
#include <sofa/core/behavior/LinearSolver.h>


#include <stdlib.h>
#include <math.h>

namespace sofa
{

namespace simulation
{
namespace common
{


ParallelSolverImpl::ParallelSolverImpl()
{}

ParallelSolverImpl::~ParallelSolverImpl()
{}



void ParallelSolverImpl::v_op(VecId v, VecId a, VecId b, Shared<double>  &f) ///< v=a+b*f
{
    ParallelMechanicalVOpVisitor(v,a,b,1.0,&f).execute( getContext() );
}
void ParallelSolverImpl::v_peq(VecId v, VecId a, Shared<double> &fSh,double f) ///< v+=f*a
{
    ParallelMechanicalVOpVisitor(v,v,a,f,&fSh).execute( getContext() );
}
void ParallelSolverImpl::v_peq(VecId v, VecId a, double f) ///< v+=f*a
{
    ParallelMechanicalVOpVisitor(v,v,a,f).execute( getContext() );
}
void ParallelSolverImpl::v_meq(VecId v, VecId a, Shared<double> &fSh) ///< v+=f*a
{
    ParallelMechanicalVOpMecVisitor(v,a,&fSh).execute( getContext() );
}

void ParallelSolverImpl::v_dot(Shared<double> &result,VecId a, VecId b) ///< a dot b ( get result using finish )
{

    ParallelMechanicalVDotVisitor(&result,a,b).execute( getContext() );
}



} // namespace simulation
} // namespace simulation

} // namespace sofa
