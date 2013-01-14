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
#ifndef SOFA_CORE_COLLISION_LINEARSOLVERCONTACTCORRECTION_H
#define SOFA_CORE_COLLISION_LINEARSOLVERCONTACTCORRECTION_H

#include <sofa/core/behavior/ConstraintCorrection.h>

#include <sofa/core/behavior/OdeSolver.h>
#include <sofa/core/behavior/LinearSolver.h>

#include <sofa/defaulttype/Mat.h>
#include <sofa/defaulttype/Vec.h>

#include <sofa/component/linearsolver/SparseMatrix.h>
#include <sofa/component/linearsolver/FullMatrix.h>

namespace sofa
{

namespace component
{

namespace constraintset
{

using namespace sofa::core;
using namespace sofa::defaulttype;


/// to avoid compilation problem under gcc3.3
extern inline behavior::OdeSolver* getOdeSolver(objectmodel::BaseContext* context)
{
    return context->get<behavior::OdeSolver>();
}

/**
 *  \brief Component computing contact forces within a simulated body using the compliance method.
 */
template<class TDataTypes>
class LinearSolverConstraintCorrection : public behavior::ConstraintCorrection< TDataTypes >
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(LinearSolverConstraintCorrection, TDataTypes), SOFA_TEMPLATE(behavior::ConstraintCorrection, TDataTypes));

    typedef TDataTypes DataTypes;
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::MatrixDeriv MatrixDeriv;
    typedef typename DataTypes::MatrixDeriv::RowConstIterator MatrixDerivRowConstIterator;
    typedef typename DataTypes::MatrixDeriv::ColConstIterator MatrixDerivColConstIterator;
    typedef typename DataTypes::MatrixDeriv::RowIterator MatrixDerivRowIterator;
    typedef typename DataTypes::MatrixDeriv::ColIterator MatrixDerivColIterator;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;

    typedef std::list<int> ListIndex;
    typedef behavior::ConstraintCorrection< TDataTypes > Inherit;
protected:
    LinearSolverConstraintCorrection(behavior::MechanicalState<DataTypes> *mm = NULL);

    virtual ~LinearSolverConstraintCorrection();
public:
    virtual void init();

    virtual void addComplianceInConstraintSpace(const ConstraintParams *cparams, defaulttype::BaseMatrix* W);

    virtual void getComplianceMatrix(defaulttype::BaseMatrix* ) const;

    virtual void computeAndApplyMotionCorrection(const ConstraintParams *cparams, MultiVecCoordId x, MultiVecDerivId v, MultiVecDerivId f, const BaseVector * lambda);

    virtual void computeAndApplyPositionCorrection(const ConstraintParams *cparams, MultiVecCoordId x, MultiVecDerivId f, const BaseVector *lambda);

    virtual void computeAndApplyVelocityCorrection(const ConstraintParams *cparams, MultiVecDerivId v, MultiVecDerivId f, const BaseVector *lambda);

    virtual void computeAndApplyMotionCorrection(const ConstraintParams * /*cparams*/, Data< VecCoord > &/*x*/, Data< VecDeriv > &/*v*/, Data< VecDeriv > &/*f*/, const BaseVector * /*lambda*/) {};

    virtual void computeAndApplyPositionCorrection(const ConstraintParams * /*cparams*/, Data< VecCoord > &/*x*/, Data< VecDeriv > &/*f*/, const BaseVector * /*lambda*/) {};

    virtual void computeAndApplyVelocityCorrection(const ConstraintParams * /*cparams*/, Data< VecDeriv > &/*v*/, Data< VecDeriv > &/*f*/, const BaseVector * /*lambda*/) {};

    virtual void applyPredictiveConstraintForce(const ConstraintParams *cparams, Data< VecDeriv > &f, const defaulttype::BaseVector *lambda);

    virtual void rebuildSystem(double massFactor, double forceFactor);

    /// @name Deprecated API
    /// @{

    virtual void applyContactForce(const defaulttype::BaseVector *f);

    virtual void resetContactForce();

    /// @}

    /// @name Unbuilt constraint system during resolution
    /// @{

    Data< bool > wire_optimization;
    Data< helper::vector< std::string > >  solverName;

    void verify_constraints();

    virtual bool hasConstraintNumber(int index) ;  // virtual ???

    virtual void resetForUnbuiltResolution(double * f, std::list<unsigned int>& renumbering);

    virtual void addConstraintDisplacement(double *d, int begin,int end) ;

    virtual void setConstraintDForce(double *df, int begin, int end, bool update) ;

    virtual void getBlockDiagonalCompliance(defaulttype::BaseMatrix* W, int begin, int end) ;

    /// @}


    /// Pre-construction check method called by ObjectFactory.
    template<class T>
    static bool canCreate(T*& obj, objectmodel::BaseContext* context, objectmodel::BaseObjectDescription* arg)
    {
        if (getOdeSolver(context) == NULL)
            return false;

        return Inherit::canCreate(obj, context, arg);
    }

protected:

    behavior::OdeSolver* odesolver;
    std::vector<sofa::core::behavior::LinearSolver*> linearsolvers;

    linearsolver::SparseMatrix<SReal> J; ///< constraint matrix
    linearsolver::FullVector<SReal> F; ///< forces computed from the constraints
#if 0 // refMinv is not use in normal case    
    linearsolver::FullMatrix<SReal> refMinv; ///< reference inverse matrix
#endif

    /**
     * @brief Compute dx correction from motion space force vector.
     */
    void computeDx(MultiVecDerivId f);

private:
    // new :  for non building the constraint system during solving process //
    VecDeriv constraint_disp, constraint_force;
    std::list<int> constraint_dofs;		// list of indices of each point which is involve with constraint // TODO : verify if useful !!
    defaulttype::BaseMatrix* systemMatrix_buf;
    defaulttype::BaseVector* systemRHVector_buf;
    defaulttype::BaseVector* systemLHVector_buf;


    // par un vecteur de listes precaclues pour chaque contrainte
    std::vector< ListIndex > Vec_I_list_dof;   // vecteur donnant la liste des indices des dofs par block de contrainte
    int last_force, last_disp; //last_force indice du dof le plus petit portant la force/le dpt qui a ?t? modifi? pour la derni?re fois (wire optimisation only?)
    bool _new_force; // if true, a "new" force was added in setConstraintDForce which is not yet integrated by a new computation in addConstraintDisplacements

};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_CONSTRAINT_LINEARSOLVERCONSTRAINTCORRECTION_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_CONSTRAINT_API LinearSolverConstraintCorrection<Vec3dTypes>;
extern template class SOFA_CONSTRAINT_API LinearSolverConstraintCorrection<Vec1dTypes>;
extern template class SOFA_CONSTRAINT_API LinearSolverConstraintCorrection<Rigid3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_CONSTRAINT_API LinearSolverConstraintCorrection<Vec3fTypes>;
extern template class SOFA_CONSTRAINT_API LinearSolverConstraintCorrection<Vec1fTypes>;
extern template class SOFA_CONSTRAINT_API LinearSolverConstraintCorrection<Rigid3fTypes>;
#endif
#endif

} // namespace collision

} // namespace component

} // namespace sofa

#endif
