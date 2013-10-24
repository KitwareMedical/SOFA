#ifndef RIGIDMASS_H
#define RIGIDMASS_H

#include <sofa/core/behavior/Mass.h>
#include <sofa/core/behavior/MechanicalState.h>

#include "utils/se3.h"
#include "utils/map.h"

#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/gl/Axis.h>

namespace sofa
{

namespace component
{

namespace mass
{

/**
   An actual rigid mass matrix. It applies to center-of-mass-centered,
   principal-axis-aligned rigid frames. 

   It is still unclear to me whether the defaulttype::RigidMass
   class/UniformMass/DiagonalMass accounts for frame changes
   correctly, so i rolled my own.

   Since SOFA uses absolute frame for rotational velocities
   (i.e. spatial velocities), the corresponding spatial inertia
   tensors are Is = R.Ib.R^T, where Ib is the body-fixed inertia
   tensor. It seems that Ib is used as a spatial inertia tensor, but I
   might be wrong.

   @author: maxime.tournier@inria.fr
 */

template <class DataTypes>
class RigidMass : public core::behavior::Mass<DataTypes>
{
public:
	SOFA_CLASS(SOFA_TEMPLATE(RigidMass, DataTypes), SOFA_TEMPLATE(core::behavior::Mass,DataTypes));
	
	typedef typename DataTypes::Real real;
	typedef helper::vector<real> mass_type;
	typedef helper::vector< defaulttype::Vec<3, real> > inertia_type;
	
	Data<mass_type> mass;
	Data<inertia_type> inertia;
	Data<bool> inertia_forces;
	
	typedef SE3<real> se3;

	RigidMass() 
		: mass(initData(&mass, "mass", "mass of each rigid body")),
		  inertia(initData(&inertia, "inertia", "inertia of each rigid body")),
		  inertia_forces(initData(&inertia_forces, false, "inertia_forces", "compute (explicit) inertia forces")) {
		
	}
	
protected:
	
	// clamps an index to the largest index in mass/inertia vectors
	unsigned clamp(unsigned i) const {
		return std::min<unsigned>(i, mass.getValue().size() - 1);
	}
public:

	void init() {
		this->core::behavior::Mass<DataTypes>::init();
		
		if( mass.getValue().size() != inertia.getValue().size() ) throw std::logic_error("mass and inertia must have the same size");
		if( !mass.getValue().size() ) throw std::logic_error("empty mass field");
		if( !this->mstate )  throw std::logic_error("no mstate");

		this->reinit();
	} 
	
	typedef typename DataTypes::VecCoord VecCoord;
	typedef typename DataTypes::VecDeriv VecDeriv;
	
	typedef core::objectmodel::Data<VecCoord> DataVecCoord;
	typedef core::objectmodel::Data<VecDeriv> DataVecDeriv;

	void draw(const core::visual::VisualParams* vparams) {

		if ( !vparams->displayFlags().getShowBehaviorModels() )
			return;
		helper::ReadAccessor<VecCoord> x = *this->mstate->getX();


		for(unsigned i = 0, n = x.size(); i < n; ++i) {
			const unsigned index = clamp(i);

			double m00 = inertia.getValue()[index][0];
			double m11 = inertia.getValue()[index][1];
			double m22 = inertia.getValue()[index][2];
			
			defaulttype::Vec3d len;
			len[0] = sqrt(m11+m22-m00);
			len[1] = sqrt(m00+m22-m11);
			len[2] = sqrt(m00+m11-m22);
			
			helper::gl::Axis::draw(x[i].getCenter(), x[i].getOrientation(), len);
		}
		
	}

	void addForce(const core::MechanicalParams*  /* PARAMS FIRST */, 
	              DataVecDeriv& _f, 
	              const DataVecCoord& _x, const DataVecDeriv& _v) {

		helper::WriteAccessor< DataVecDeriv > f(_f);
		helper::ReadAccessor< DataVecCoord >  x(_x);
    helper::ReadAccessor< DataVecDeriv >  v(_v);

		typename se3::vec3 g = SE3<SReal>::map(this->getContext()->getGravity()).template cast<real>();

		for(unsigned i = 0, n = this->mstate->getSize(); i < n; ++i) {
			const unsigned index = clamp(i);
		
			se3::map( f[i].getVCenter() ).template head<3>() += mass.getValue()[index] * g;
			
			// explicit inertia
			if( inertia_forces.getValue() ) {
				typename se3::mat33 R = se3::rotation( x[i] ).toRotationMatrix();

				// spatial velocity
				typename se3::vec3 omega = se3::map( v[i].getVOrientation() );

				// body inertia tensor
				typename se3::vec3 I = se3::map(inertia.getValue()[ index ]);
				
				se3::map( f[i].getVCenter() ).template head<3>() -= omega.cross( R * I.cwiseProduct( R.transpose() * omega) );
				
			}

		}
		
	
	}

	// perdu: il faut aussi la position pour connaître l'énergie
	// cinétique d'un rigide, sauf si on considère les vitesses en
	// coordonnées locales (ce que sofa ne fait pas). du coup on tape
	// dans this->mstate->getX pour l'obtenir mais l'api devrait gérer ca
	double getKineticEnergy( const core::MechanicalParams* /* PARAMS FIRST */, 
	                         const DataVecDeriv& _v  ) const {
		helper::ReadAccessor< DataVecDeriv >  v(_v);
		
		double res = 0;

		for(unsigned i = 0, n = v.size(); i < n; ++i) {
			const unsigned index = clamp(i);
			
			// body-fixed velocity
			typename se3::vec3 omega_body = se3::rotation( (*this->mstate->getX())[i] ).inverse() * se3::map(v[i].getVOrientation());
			
			res += 
				mass.getValue()[index] * v[i].getVCenter().norm2() +
				se3::map(inertia.getValue()[index]).cwiseProduct( omega_body ).dot(omega_body);
		}
		
		res *= 0.5;
		
		return res;
	}

	// TODO maybe sign is wrong 
	double getPotentialEnergy( const core::MechanicalParams* /* PARAMS FIRST */, 
	                           const DataVecCoord& _x  ) const {
		helper::ReadAccessor< DataVecCoord >  x(_x);
				
		defaulttype::Vec3d g ( this->getContext()->getGravity() );

		double res = 0;

		for(unsigned i = 0, n = x.size(); i < n; ++i) {
			const unsigned index = clamp(i);
	
			res += mass.getValue()[index] * (g * x[i].getCenter()); 
		}
		
		return res;
	}


	virtual void addMDx(const core::MechanicalParams*  /* PARAMS FIRST */, 
	                    DataVecDeriv& _f, 
	                    const DataVecDeriv& _dx, 
	                    double factor) {
		helper::WriteAccessor< DataVecDeriv >  f(_f);
		helper::ReadAccessor< DataVecDeriv >  dx(_dx);

		for(unsigned i = 0, n = this->mstate->getSize(); i < n; ++i) {
			const unsigned index = clamp(i);
			
			map(f[i].getLinear()) += (factor * mass.getValue()[ index ]) * map(dx[i].getLinear());

			typename se3::quat q = se3::rotation( (*this->mstate->getX())[i] );
			map(f[i].getAngular()) += factor * 
				( q * map(inertia.getValue()[ index ]).cwiseProduct( q.conjugate() * map(dx[i].getAngular() ) ));
			
		}
		
	}


    virtual void addMToMatrix(const core::MechanicalParams* mparams,
	                          const sofa::core::behavior::MultiMatrixAccessor* matrix) {
		
		sofa::core::behavior::MultiMatrixAccessor::MatrixRef r = matrix->getMatrix( this->mstate );
		const unsigned size = defaulttype::DataTypeInfo<typename DataTypes::Deriv>::size();

        real mFactor = (real)mparams->mFactorIncludingRayleighDamping(this->rayleighMass.getValue());
		
		for(unsigned i = 0, n = this->mstate->getSize(); i < n; ++i) {

			const unsigned index = clamp(i);
			
			// translation
			for(unsigned j = 0; j < 3; ++j) {
				r.matrix->add(r.offset + size * i + j,
				              r.offset + size * i + j,
                              mass.getValue()[ index ] * mFactor );
			}			              
			
			typename se3::mat33 R = se3::rotation( (*this->mstate->getX())[i] ).toRotationMatrix();
			
			typename se3::mat33 chunk = R * se3::map( inertia.getValue()[ index ] ).asDiagonal() * R.transpose();
			
			// rotation
			for(unsigned j = 0; j < 3; ++j) {
				for(unsigned k = 0; k < 3; ++k) {
					
					r.matrix->add(r.offset + size * i + 3 + j,
					              r.offset + size * i + 3 + k,
                                  chunk(j, k) * mFactor );
				}
			}			              
			
		}
		
	}
};

}
}
}

#endif
