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
#ifndef SOFA_CORE_BEHAVIOR_MECHANICALSTATE_H
#define SOFA_CORE_BEHAVIOR_MECHANICALSTATE_H

#include <sofa/core/behavior/BaseMechanicalState.h>
#include <sofa/core/VecId.h>
#include <sofa/core/State.h>
#include <sofa/defaulttype/DataTypeInfo.h>

namespace sofa
{

namespace core
{

namespace behavior
{

/**
 *  \brief Component storing all state vectors of a simulated body (position,
 *  velocity, etc), using the datatype specified in the templace.
 *
 *  The given DataTypes class should define the following internal types:
 *  \li \code Real \endcode : scalar values (float or double).
 *  \li \code Coord \endcode : position values.
 *  \li \code Deriv \endcode : derivative values (velocity, forces, displacements).
 *  \li \code VecReal \endcode : container of scalar values with the same API as sofa::helper::vector.
 *  \li \code VecCoord \endcode : container of Coord values with the same API as sofa::helper::vector.
 *  \li \code VecDeriv \endcode : container of Deriv values with the same API as sofa::helper::vector.
 *  \li \code MatrixDeriv \endcode : vector of constraints.
 *
 *  Other vectors can be allocated to store other temporary values.
 *  Vectors can be assigned efficiently by just swapping pointers.
 *
 *  In addition to state vectors, the current constraint system matrix is also
 *  stored, containing the coefficient of each constraint defined over the DOFs
 *  in this body.
 *
 */
template<class TDataTypes>
class MechanicalState : public BaseMechanicalState, public State<TDataTypes>
{
public:
    SOFA_CLASS2(SOFA_TEMPLATE(MechanicalState,TDataTypes), BaseMechanicalState, SOFA_TEMPLATE(State,TDataTypes));

    typedef TDataTypes DataTypes;
    /// Scalar values (float or double).
    typedef typename DataTypes::Real Real;
    /// Position values.
    typedef typename DataTypes::Coord Coord;
    /// Derivative values (velocity, forces, displacements).
    typedef typename DataTypes::Deriv Deriv;
    /// Container of scalar values with the same API as sofa::helper::vector.
    typedef typename DataTypes::VecReal VecReal;
    /// Container of Coord values with the same API as sofa::helper::vector.
    typedef typename DataTypes::VecCoord VecCoord;
    /// Container of Deriv values with the same API as sofa::helper::vector.
    typedef typename DataTypes::VecDeriv VecDeriv;
    /// Sparse matrix containing derivative values (constraints)
    typedef typename DataTypes::MatrixDeriv MatrixDeriv;
protected:
    virtual ~MechanicalState() {}
public:
    /// Return the current free-motion position vector.
    /// @deprecated use readVecCoord(ConstVecCoordId::freePosition()) instead.
    virtual const VecCoord* getXfree()  const
    {
        const Data<VecCoord>* v = this->read(ConstVecCoordId::freePosition());
        return (v == NULL) ? NULL : &(v->getValue());
    }

    /// Return the current free-motion velocity vector.
    /// @deprecated use readVecDeriv(ConstVecDerivId::freeVelocity()) instead.
    virtual const VecDeriv* getVfree()  const
    {
        const Data<VecDeriv>* v = this->read(ConstVecDerivId::freeVelocity());
        return (v == NULL) ? NULL : &(v->getValue());
    }

    /// Return the current reset position vector.
    /// (return NULL if the state does not store rest position).
    /// @deprecated use readVecCoord(ConstVecCoordId::resetPosition()) instead.
    virtual const VecCoord* getXReset() const
    {
        const Data<VecCoord>* v = this->read(ConstVecCoordId::resetPosition());
        return (v == NULL) ? NULL : &(v->getValue());
    }

    /// Return the force vector.
    /// @deprecated use readVecDeriv(ConstVecDerivId::force()) instead.
    virtual const VecDeriv* getF() const
    {
        const Data<VecDeriv>* v = this->read(ConstVecDerivId::force());
        return (v == NULL) ? NULL : &(v->getValue());
    }

    /// Return the external forces vector.
    /// @deprecated use readVecDeriv(ConstVecDerivId::externalForce()) instead.
    virtual const VecDeriv* getExternalForces() const
    {
        const Data<VecDeriv>* v = this->read(ConstVecDerivId::externalForce());
        return (v == NULL) ? NULL : &(v->getValue());
    }

    /// Return the dx vector.
    /// @deprecated use readVecDeriv(ConstVecDerivId::dx()) instead.
    virtual const VecDeriv* getDx() const
    {
        const Data<VecDeriv>* v = this->read(ConstVecDerivId::dx());
        return (v == NULL) ? NULL : &(v->getValue());
    }

    /// Return the dx vector.
    /// @deprecated use readVecDeriv(ConstVecDerivId::restVelocity()) instead.
    virtual const VecDeriv* getVReset() const
    {
        const Data<VecDeriv>* v = this->read(ConstVecDerivId::resetVelocity());
        return (v == NULL) ? NULL : &(v->getValue());
    }

    /// Return the constraint sparse matrix.
    /// @deprecated use readVecDeriv(ConstMatrixDerivId::holonomicC()) instead.
    virtual const MatrixDeriv* getC() const
    {
        const Data<MatrixDeriv>* v = this->read(ConstMatrixDerivId::holonomicC());
        return (v == NULL) ? NULL : &(v->getValue());
    }

    virtual unsigned int getCoordDimension() const { return defaulttype::DataTypeInfo<Coord>::size(); }
    virtual unsigned int getDerivDimension() const { return defaulttype::DataTypeInfo<Deriv>::size(); }

    /// Get the indices of the particles located in the given bounding box
    virtual void getIndicesInSpace(sofa::helper::vector<unsigned>& /*indices*/, Real /*xmin*/, Real /*xmax*/,Real /*ymin*/, Real /*ymax*/, Real /*zmin*/, Real /*zmax*/) const=0;

    virtual std::string getTemplateName() const
    {
        return templateName(this);
    }

    static std::string templateName(const MechanicalState<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }

    template<class T>
    static std::string shortName(const T* ptr = NULL, objectmodel::BaseObjectDescription* arg = NULL)
    {
        std::string name = Inherit1::shortName(ptr, arg);
        sofa::helper::replaceAll(name, "Mechanical", "M");
        sofa::helper::replaceAll(name, "mechanical", "m");
        return name;
    }

	virtual void copyToBuffer(SReal* dst, ConstVecId src, unsigned n) const {
		const unsigned size = this->getSize();
		
		switch(src.type) {
		case V_COORD: {
			helper::ReadAccessor< Data<VecCoord> > vec = this->read(ConstVecCoordId(src));
			const unsigned dim = defaulttype::DataTypeInfo<Coord>::size();
			assert( n == dim * size );
			
			for(unsigned i = 0; i < size; ++i) {
				for(unsigned j = 0; j < dim; ++j) {
					defaulttype::DataTypeInfo<Coord>::getValue(vec[i], j, *(dst++));
				}
			}
			
		}; break;
		case V_DERIV: {
			helper::ReadAccessor< Data<VecDeriv> > vec = this->read(ConstVecDerivId(src));
			const unsigned dim = defaulttype::DataTypeInfo<Deriv>::size();
			assert( n == dim * size );
			
			for(unsigned i = 0; i < size; ++i) {
				for(unsigned j = 0; j < dim; ++j) {
					defaulttype::	DataTypeInfo<Deriv>::getValue(vec[i], j, *(dst++));
				}
			}
			
		}; break;
		default: 
			assert( false );
		}
		
		// get rid of unused parameter warnings in release build
		(void) n;
	}

	virtual void copyFromBuffer(VecId dst, SReal* src, unsigned n) const {
		const unsigned size = this->getSize();
		
		switch(dst.type) {
		case V_COORD: {
			helper::WriteAccessor< Data<VecCoord> > vec = this->read(VecCoordId(dst));
			const unsigned dim = defaulttype::DataTypeInfo<Coord>::size();
			assert( n == dim * size );
			
			for(unsigned i = 0; i < size; ++i) {
				for(unsigned j = 0; j < dim; ++j) {
					defaulttype::DataTypeInfo<Coord>::setValue(vec[i], j, *(src++));
				}
			}
			
		}; break;
		case V_DERIV: {
			helper::WriteAccessor< Data<VecDeriv> > vec = this->read(VecDerivId(dst));
			const unsigned dim = defaulttype::DataTypeInfo<Deriv>::size();
			assert( n == dim * size );
			
			for(unsigned i = 0; i < size; ++i) {
				for(unsigned j = 0; j < dim; ++j) {
					defaulttype::	DataTypeInfo<Deriv>::setValue(vec[i], j, *(src++));
				}
			}
			
		}; break;
		default: 
			assert( false );
		}
		
		// get rid of unused parameter warnings in release build
		(void) n;
	}

};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_BUILD_CORE)
#ifndef SOFA_FLOAT
extern template class SOFA_CORE_API MechanicalState<defaulttype::Vec3dTypes>;
extern template class SOFA_CORE_API MechanicalState<defaulttype::Vec2dTypes>;
extern template class SOFA_CORE_API MechanicalState<defaulttype::Vec1dTypes>;
extern template class SOFA_CORE_API MechanicalState<defaulttype::Vec6dTypes>;
extern template class SOFA_CORE_API MechanicalState<defaulttype::Rigid3dTypes>;
extern template class SOFA_CORE_API MechanicalState<defaulttype::Rigid2dTypes>;
extern template class SOFA_CORE_API MechanicalState<defaulttype::ExtVec3dTypes>;
#endif

#ifndef SOFA_DOUBLE
extern template class SOFA_CORE_API MechanicalState<defaulttype::Vec3fTypes>;
extern template class SOFA_CORE_API MechanicalState<defaulttype::Vec2fTypes>;
extern template class SOFA_CORE_API MechanicalState<defaulttype::Vec1fTypes>;
extern template class SOFA_CORE_API MechanicalState<defaulttype::Vec6fTypes>;
extern template class SOFA_CORE_API MechanicalState<defaulttype::Rigid3fTypes>;
extern template class SOFA_CORE_API MechanicalState<defaulttype::Rigid2fTypes>;
extern template class SOFA_CORE_API MechanicalState<defaulttype::ExtVec3fTypes>;
#endif

#endif

} // namespace behavior

} // namespace core

} // namespace sofa

#endif
