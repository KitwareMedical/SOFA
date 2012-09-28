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
#ifndef SOFA_COMPONENT_LINEARSOLVER_GRAPHSCATTEREDTYPES_H
#define SOFA_COMPONENT_LINEARSOLVER_GRAPHSCATTEREDTYPES_H

#include <sofa/simulation/common/MechanicalOperations.h>
#include <sofa/core/behavior/MultiVec.h>
#include <sofa/core/behavior/LinearSolver.h>
#include <sofa/component/component.h>
#include <sofa/component/linearsolver/SparseMatrix.h>
#include <sofa/component/linearsolver/FullMatrix.h>
#ifdef SOFA_SMP
#include <sofa/core/behavior/ParallelMultiVec.h>
#endif

namespace sofa
{

namespace component
{

namespace linearsolver
{

class GraphScatteredMatrix;
class GraphScatteredVector;
#ifdef SOFA_SMP
class ParallelGraphScatteredVector;
#endif
template <class T1, class T2>
class MultExpr
{
public:
    T1& a;
    T2& b;
    MultExpr(T1& a, T2& b) : a(a), b(b) {}
};

class SOFA_BASE_LINEAR_SOLVER_API GraphScatteredMatrix
{
public:
    typedef SReal Real;
    //simulation::SolverImpl* parent;
    //double mFact, bFact, kFact;
    core::MechanicalParams mparams;
    simulation::common::MechanicalOperations* parent;
public:
    GraphScatteredMatrix()
        : parent(NULL) //, mFact(0.0), bFact(0.0), kFact(0.0)
    {
    }
    void setMBKFacts(const core::MechanicalParams* mparams)
    {
        this->mparams = *mparams;
    }
    MultExpr<GraphScatteredMatrix,GraphScatteredVector> operator*(GraphScatteredVector& v)
    {
        return MultExpr<GraphScatteredMatrix,GraphScatteredVector>(*this, v);
    }
#ifdef SOFA_SMP
    MultExpr<GraphScatteredMatrix,ParallelGraphScatteredVector> operator*(ParallelGraphScatteredVector& v)
    {
        return MultExpr<GraphScatteredMatrix,ParallelGraphScatteredVector>(*this, v);
    }
#endif
    void apply(GraphScatteredVector& res, GraphScatteredVector& x);
#ifdef SOFA_SMP
    void apply(ParallelGraphScatteredVector& res, ParallelGraphScatteredVector& x);
#endif

    // compatibility with baseMatrix
    unsigned int rowSize()
    {
        unsigned int nbRow=0, nbCol=0;
        this->parent->getMatrixDimension(&nbRow, &nbCol);
        return nbRow;

    }
    int colSize()
    {
        unsigned int nbRow=0, nbCol=0;
        this->parent->getMatrixDimension(&nbRow, &nbCol);
        return nbCol;
    }

    //void papply(GraphScatteredVector& res, GraphScatteredVector& x);

    static const char* Name() { return "GraphScattered"; }
};

class SOFA_BASE_LINEAR_SOLVER_API GraphScatteredVector : public sofa::core::behavior::MultiVecDeriv
{
public:
    typedef sofa::core::behavior::MultiVecDeriv Inherit;
    typedef SReal Real;
    GraphScatteredVector(core::behavior::BaseVectorOperations* p, core::VecDerivId id)
        : Inherit(p, id)
    {
    }
    GraphScatteredVector(core::behavior::BaseVectorOperations* p)
        : Inherit(p)
    {
    }
    void set(core::MultiVecDerivId id)
    {
        this->v = id;
    }
    void reset()
    {
        this->v = core::VecDerivId::null();
    }


    /// TO IMPLEMENT
    void add(int /*row*/, SReal /*v*/)
    {
        std::cerr<<"WARNING : add an element is not supported in MultiVector"<<std::endl;
    }

    /// TO IMPLEMENT
    void set(int /*row*/, SReal /*v*/)
    {
        std::cerr<<"WARNING : set an element is not supported in MultiVector"<<std::endl;
    }

    SReal element(int /*i*/)
    {
        std::cerr<<"WARNING : get a single element is not supported in MultiVector"<<std::endl;
        return 0;
    }

    friend class GraphScatteredMatrix;

    void operator=(const MultExpr<GraphScatteredMatrix,GraphScatteredVector>& expr)
    {
        expr.a.apply(*this,expr.b);
    }

    //void operator+=(const MultExpr<GraphScatteredMatrix,GraphScatteredVector>& expr)
    //{
    //    expr.a.papply(*this,expr.b);
    //}

    static const char* Name() { return "GraphScattered"; }
};

#ifdef SOFA_SMP
class SOFA_BASE_LINEAR_SOLVER_API ParallelGraphScatteredVector : public sofa::core::behavior::ParallelMultiVecDeriv
{
public:
    typedef sofa::core::behavior::ParallelMultiVecDeriv Inherit;
    ParallelGraphScatteredVector(core::behavior::BaseVectorOperations* p, core::VecDerivId id)
        : Inherit(p, id)
    {
    }
    ParallelGraphScatteredVector(core::behavior::BaseVectorOperations* p)
        : Inherit(p)
    {
    }
    void set(core::MultiVecDerivId id)
    {
        this->v = id;
    }
    void reset()
    {
        this->v = core::VecDerivId::null();
    }

    /// TO IMPLEMENT
    void add(int /*row*/, SReal /*v*/)
    {
        std::cerr<<"WARNING : add an element is not supported in ParallelMultiVector"<<std::endl;
    }

    SReal element(int /*i*/)
    {
        std::cerr<<"WARNING : get a single element is not supported in ParallelMultiVector"<<std::endl;
        return 0;
    }

    friend class GraphScatteredMatrix;

    void operator=(const MultExpr<GraphScatteredMatrix,ParallelGraphScatteredVector>& expr)
    {
        expr.a.apply(*this,expr.b);
    }

    static const char* Name() { return "ParallelGraphScattered"; }
};
#endif /* SOFA_SMP */

} // namespace linearsolver

} // namespace component

} // namespace sofa

#endif
