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
#ifndef SOFA_CORE_BEHAVIOR_BASEMECHANICALSTATE_H
#define SOFA_CORE_BEHAVIOR_BASEMECHANICALSTATE_H

#include <sofa/core/BaseState.h>
#include <sofa/core/MultiVecId.h>
#include <sofa/defaulttype/BaseMatrix.h>
#include <sofa/defaulttype/BaseVector.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/Quat.h>
#include <sofa/helper/ParticleMask.h>

#include <iostream>


namespace sofa
{

namespace core
{

namespace behavior
{

/**
 *  \brief Component storing all state vectors of a simulated body (position, velocity, etc).
 *
 *  This class only contains the data of the body and not any of its
 *  <i>active</i> computations, which are handled by the Mass, ForceField, and
 *  Constraint components.
 *
 *  Two types of vectors are used :
 *  \li \code VecCoord \endcode : containing positions.
 *  \li \code VecDeriv \endcode : derivative values, i.e. velocity, forces, displacements.
 *  In most cases they are the same (i.e. 3D/2D point particles), but they can
 *  be different (rigid frames for instance).
 *
 *  Several pre-defined vectors are stored :
 *  \li \code position \endcode
 *  \li \code velocity \endcode
 *  \li \code force \endcode
 *  \li \code dx \endcode (displacement)
 *
 *  Other vectors can be allocated to store other temporary values.
 *  Vectors can be assigned efficiently by just swapping pointers.
 *
 *  In addition to state vectors, the current constraint system matrix is also
 *  stored, containing the coefficient of each constraint defined over the DOFs
 *  in this body.
 *
 */
class SOFA_CORE_API BaseMechanicalState : public virtual BaseState
{
public:
    SOFA_ABSTRACT_CLASS(BaseMechanicalState, BaseState);
protected:
    BaseMechanicalState();

    virtual ~BaseMechanicalState();
public:
    /// @name Methods allowing to have access to the geometry without a template class (generic but not efficient)
    /// @{
    virtual double getPX(int /*i*/) const { return 0.0; }
    virtual double getPY(int /*i*/) const { return 0.0; }
    virtual double getPZ(int /*i*/) const { return 0.0; }

    /// @}

    /// @name Vectors allocation and generic operations (based on VecId)
    /// @{
    /// Increment the index of the given VecCoordId, so that all 'allocated' vectors in this state have a lower index
    virtual void vAvail(const ExecParams* params /* PARAMS FIRST  = ExecParams::defaultInstance()*/, VecCoordId& v) = 0;
    /// Increment the index of the given VecDerivId, so that all 'allocated' vectors in this state have a lower index
    virtual void vAvail(const ExecParams* params /* PARAMS FIRST  = ExecParams::defaultInstance()*/, VecDerivId& v) = 0;
    /// Increment the index of the given MatrixDerivId, so that all 'allocated' vectors in this state have a lower index
    //virtual void vAvail(MatrixDerivId& v) = 0;

    /// Allocate a new temporary vector
    virtual void vAlloc(const ExecParams* params /* PARAMS FIRST  = ExecParams::defaultInstance()*/, VecCoordId v) = 0;
    /// Allocate a new temporary vector
    virtual void vAlloc(const ExecParams* params /* PARAMS FIRST  = ExecParams::defaultInstance()*/, VecDerivId v) = 0;
    /// Allocate a new temporary vector
    //virtual void vAlloc(MatrixDerivId v) = 0;

    /// Reallocate a new temporary vector
    virtual void vRealloc(const ExecParams* params /* PARAMS FIRST  = ExecParams::defaultInstance()*/, VecCoordId v) = 0;
    /// Reallocate a new temporary vector
    virtual void vRealloc(const ExecParams* params /* PARAMS FIRST  = ExecParams::defaultInstance()*/, VecDerivId v) = 0;


    /// Free a temporary vector
    virtual void vFree(const ExecParams* params /* PARAMS FIRST  = ExecParams::defaultInstance()*/, VecCoordId v) = 0;
    /// Free a temporary vector
    virtual void vFree(const ExecParams* params /* PARAMS FIRST  = ExecParams::defaultInstance()*/, VecDerivId v) = 0;
    /// Free a temporary vector
    //virtual void vFree(MatrixDerivId v) = 0;

    /// Initialize an unset vector
    virtual void vInit(const ExecParams* params /* PARAMS FIRST  = ExecParams::defaultInstance()*/, VecCoordId v, ConstVecCoordId vSrc) = 0;
    /// Initialize an unset vector
    virtual void vInit(const ExecParams* params /* PARAMS FIRST  = ExecParams::defaultInstance()*/, VecDerivId v, ConstVecDerivId vSrc) = 0;
    /// Initialize an unset vector
    //virtual void vInit(const ExecParams* params /* PARAMS FIRST  = ExecParams::defaultInstance()*/, MatrixDerivId v, ConstMatrixDerivId vSrc) = 0;

    /// Compute a linear operation on vectors : v = a + b * f.
    ///
    /// This generic operation can be used for many simpler cases :
    /// \li v = 0
    /// \li v = a
    /// \li v = a + b
    /// \li v = b * f
    virtual void vOp(const ExecParams* params /* PARAMS FIRST  = ExecParams::defaultInstance()*/, VecId v, ConstVecId a = ConstVecId::null(), ConstVecId b = ConstVecId::null(), double f = 1.0 ) = 0;
#ifdef SOFA_SMP
    virtual void vOp(const ExecParams* params /* PARAMS FIRST  = ExecParams::defaultInstance()*/, VecId v, ConstVecId a, ConstVecId b, double f, a1::Shared<double> * fSh ) = 0;
    virtual void vOpMEq(const ExecParams* params /* PARAMS FIRST  = ExecParams::defaultInstance()*/, VecId v, ConstVecId a = ConstVecId::null(), a1::Shared<double> * fSh=NULL ) = 0;
    virtual void vDot(const ExecParams* params /* PARAMS FIRST  = ExecParams::defaultInstance()*/, a1::Shared<double> *result,ConstVecId a, ConstVecId b ) = 0;
#endif
    /// Data structure describing a set of linear operation on vectors
    /// \see vMultiOp
    class VMultiOpEntry : public std::pair< MultiVecId, helper::vector< std::pair< ConstMultiVecId, double > > >
    {
    public:
        typedef std::pair< ConstMultiVecId, double > Fact;
        typedef helper::vector< Fact > VecFact;
        typedef std::pair< MultiVecId, VecFact > Inherit;
        VMultiOpEntry() : Inherit(MultiVecId::null(), VecFact()) {}
        VMultiOpEntry(MultiVecId v) : Inherit(v, VecFact()) {}
        VMultiOpEntry(MultiVecId v, ConstMultiVecId a, double af = 1.0) : Inherit(v, VecFact())
        { this->second.push_back(Fact(a, af)); }
        VMultiOpEntry(MultiVecId v, ConstMultiVecId a, ConstMultiVecId b, double bf = 1.0) : Inherit(v, VecFact())
        { this->second.push_back(Fact(a,1.0));  this->second.push_back(Fact(b, bf)); }
        VMultiOpEntry(MultiVecId v, ConstMultiVecId a, double af, ConstMultiVecId b, double bf = 1.0) : Inherit(v, VecFact())
        { this->second.push_back(Fact(a, af));  this->second.push_back(Fact(b, bf)); }
    };

    typedef helper::vector< VMultiOpEntry > VMultiOp;

    /// Perform a sequence of linear vector accumulation operation $r_i = sum_j (v_j*f_{ij})$
    ///
    /// This is used to compute in on steps operations such as $v = v + a*dt, x = x + v*dt$.
    /// Note that if the result vector appears inside the expression, it must be the first operand.
    /// By default this method decompose the computation into multiple vOp calls.
    virtual void vMultiOp(const ExecParams* params /* PARAMS FIRST  = ExecParams::defaultInstance()*/, const VMultiOp& ops);

    /// Compute the scalar products between two vectors.
    virtual double vDot(const ExecParams* params /* PARAMS FIRST  = ExecParams::defaultInstance()*/, ConstVecId a, ConstVecId b) = 0; //{ return 0; }

    /// Get vector size
    virtual size_t vSize( const ExecParams* params, ConstVecId v ) = 0;

    /// Apply a threshold to all entries
    virtual void vThreshold( VecId a, double threshold ) = 0;

    /// @}

    /// @name Mechanical integration related methods
    /// Note: all these methods can now be implemented generically using VecId-based operations
    /// @{

    /// Called at the beginning of each integration step.
    virtual void beginIntegration(double /*dt*/)
    {
        // it is no longer necessary to switch forceId to internalForce here...
    }

    /// Called at the end of each integration step.
    virtual void endIntegration(const ExecParams* params /* PARAMS FIRST */, double /*dt*/)
    {
        vOp(params /* PARAMS FIRST */, VecId::externalForce(), ConstVecId::null(), ConstVecId::null(), 1.0); // externalForce = 0
    }

    /// Set F = 0
    virtual void resetForce( const ExecParams* params /* PARAMS FIRST  = ExecParams::defaultInstance()*/, VecId f = VecId::force())
    { vOp( params /* PARAMS FIRST */, f, ConstVecId::null(), ConstVecId::null(), 1.0 ); }

    /// Set Acc =0
    virtual void resetAcc( const ExecParams* params /* PARAMS FIRST   = ExecParams::defaultInstance()*/, VecId a = VecId::dx() /* VecId::accFromFrame() */ )
    { vOp( params /* PARAMS FIRST */, a, ConstVecId::null(), ConstVecId::null(), 1.0 ); }

    /// Add stored external forces to F
    virtual void accumulateForce( const ExecParams* params /* PARAMS FIRST  = ExecParams::defaultInstance()*/, VecId f = VecId::force() )
    {
        vOp( params /* PARAMS FIRST */, f, f, ConstVecId::externalForce(), 1.0 ); // f += externalForce
    }

    /// @}

    /// @name Constraints related methods
    /// @{

    /// Reset the constraint matrix
    virtual void resetConstraint(const ExecParams* params = ExecParams::defaultInstance()) = 0;

    /// build the jacobian of the constraint in a baseMatrix
    virtual void getConstraintJacobian(const ExecParams* params, sofa::defaulttype::BaseMatrix* J,unsigned int & off) = 0;

    /// Renumber the constraint ids with the given permutation vector
    virtual void renumberConstraintId(const sofa::helper::vector<unsigned>& renumbering) = 0;

    class ConstraintBlock
    {
    public:
        ConstraintBlock( unsigned int c, defaulttype::BaseMatrix *m):column(c),matrix(m) {}

        unsigned int getColumn() const {return column;}
        const defaulttype::BaseMatrix &getMatrix() const {return *matrix;};
        defaulttype::BaseMatrix *getMatrix() {return matrix;};
    protected:
        unsigned int column;
        defaulttype::BaseMatrix *matrix;
    };

    /// Express the matrix L in term of block of matrices, using the indices of the lines in the MatrixDeriv container
    virtual std::list<ConstraintBlock> constraintBlocks( const std::list<unsigned int> &/* indices */) const
    {  return std::list<ConstraintBlock>();  }

    /// Compute the error given a state vector and a line of the Jacobian (line in vector C)
    virtual SReal getConstraintJacobianTimesVecDeriv( unsigned int /*line*/, ConstVecId /*id*/)
    {  this->serr << "NOT IMPLEMENTED YET" << this->sendl; return (SReal)0;  }

    /// @}

    /// @name events
    ///   Methods related to Event processing
    /// @{

    /// Handle state Changes
    /// @deprecated topological changes now rely on TopologyEngine
    virtual void handleStateChange() {}

    /// Handle state Changes from a given Topology
    /// @deprecated topological changes now rely on TopologyEngine
    virtual void handleStateChange(core::topology::Topology* t);

    ///@}

    /// @name Misc properties and actions
    /// @{

    /// Write current state to the given output stream
    virtual void writeState( std::ostream& out );

    virtual unsigned int getCoordDimension() const { return 0; }
    virtual unsigned int getDerivDimension() const { return 0; }

    /// Translate the current state
    virtual void applyTranslation(const double dx, const double dy, const double dz)=0;

    /// Rotate the current state
    /// This method is optional, it is used when the user want to interactively change the position of an object using Euler angles
    virtual void applyRotation (const double /*rx*/, const double /*ry*/, const double /*rz*/) {};

    /// Rotate the current state
    virtual void applyRotation(const defaulttype::Quat q)=0;

    /// Scale the current state
    virtual void applyScale(const double /*sx*/,const double /*sy*/,const double /*sz*/)=0;

    virtual defaulttype::Vector3 getScale() const { return defaulttype::Vector3(1.0,1.0,1.0); }

    virtual bool addBBox(double* /*minBBox*/, double* /*maxBBox*/)
    {
        return false;
    }

    /// Find mechanical particles hit by the given ray.
    /// A mechanical particle is defined as a 2D or 3D, position or rigid DOF
    /// Returns false if this object does not support picking
    virtual bool pickParticles(const ExecParams* /* params */ /* PARAMS FIRST */, double /*rayOx*/, double /*rayOy*/, double /*rayOz*/,
            double /*rayDx*/, double /*rayDy*/, double /*rayDz*/,
            double /*radius0*/, double /*dRadius*/,
            std::multimap< double, std::pair<sofa::core::behavior::BaseMechanicalState*, int> >& /*particles*/)
    {
        return false;
    }

    /// @}

    /// @name Mask-based optimized computations (by only updating a subset of the DOFs)
    /// @{

    Data<bool> useMask;
    /// Mask to filter the particles. Used inside MechanicalMappings inside applyJ and applyJT methods.
    helper::ParticleMask forceMask;

    /// @}

    /// @name Interface with BaseMatrix / BaseVector
    /// @{

    /// Get the number of scalars per Deriv value, as necessary to build mechanical matrices and vectors.
    /// If not all Derivs have the same number of scalars, then return 1 here and overload the getMatrixSize() method.
    virtual unsigned int getMatrixBlockSize() const { return getDerivDimension(); }

    /// Get the number of rows necessary to build mechanical matrices and vectors.
    /// In most cases this is equivalent to getSize() * getMatrixBlockSize().
    virtual unsigned int getMatrixSize() const { return getSize() * getMatrixBlockSize(); }

    /// Copy data to a global BaseVector from the state stored in a local vector
    /// @param offset the offset in the BaseVector where the scalar values will be used. It will be updated to the first scalar value after the ones used by this operation when this method returns
    virtual void copyToBaseVector(defaulttype::BaseVector* dest, ConstVecId src, unsigned int &offset) = 0;

    /// Copy data to a local vector from the state stored in a global BaseVector
    /// @param offset the offset in the BaseVector where the scalar values will be used. It will be updated to the first scalar value after the ones used by this operation when this method returns
    virtual void copyFromBaseVector(VecId dest, const defaulttype::BaseVector* src, unsigned int &offset) = 0;

    /// Add data to a global BaseVector from the state stored in a local vector
    /// @param offset the offset in the BaseVector where the scalar values will be used. It will be updated to the first scalar value after the ones used by this operation when this method returns
    virtual void addToBaseVector(defaulttype::BaseVector* dest, ConstVecId src, unsigned int &offset) = 0;

    /// Performs dest[i][j] += src[offset + i][j] 0<= i < src_entries 0<= j < 3 (for 3D objects) 0 <= j < 2 (for 2D objects)
    /// @param offset the offset in the BaseVector where the scalar values will be used. It will be updated to the first scalar value after the ones used by this operation when this method returns
    virtual void addFromBaseVectorSameSize(VecId dest, const defaulttype::BaseVector* src, unsigned int &offset) = 0;


    /// Performs dest[ offset + i ][j] += src[i][j]  0<= i < src_entries  0<= j < 3 (for 3D objects) 0 <= j < 2 (for 2D objects)
    /// @param offset the offset in the MechanicalObject local vector specified by VecId dest. It will be updated to the first scalar value after the ones used by this operation when this method returns.
    virtual void addFromBaseVectorDifferentSize(VecId dest, const defaulttype::BaseVector* src, unsigned int &offset ) = 0;
    /// @}

    /// @name Data output
    /// @{

    virtual void printDOF( ConstVecId v, std::ostream& out = std::cerr, int firstIndex = 0, int range = -1 ) const = 0;
    virtual unsigned printDOFWithElapsedTime(ConstVecId /*v*/, unsigned /*count*/ = 0, unsigned /*time*/ = 0, std::ostream& /*out*/ = std::cerr ) { return 0; }
    virtual void initGnuplot(const std::string /*filepath*/) {}
    virtual void exportGnuplot(double /*time*/) {}

    virtual void writeVec(ConstVecId v, std::ostream &out) = 0;
    virtual void readVec(VecId v, std::istream &in) = 0;
    virtual double compareVec(ConstVecId v, std::istream &in) = 0;

    /// @}getPotent
};

} // namespace behavior

} // namespace core

} // namespace sofa

#endif
