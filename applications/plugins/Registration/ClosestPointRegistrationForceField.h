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
*                               SOFA :: Plugins                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef CLOSESTPOINTREGISTRATIONFORCEFIELD_H
#define CLOSESTPOINTREGISTRATIONFORCEFIELD_H

#include <sofa/core/core.h>
#include <sofa/core/behavior/BaseForceField.h>
#include <sofa/core/behavior/ForceField.h>
#include <sofa/component/interactionforcefield/SpringForceField.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/helper/accessor.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/helper/vector.h>
#include <sofa/component/component.h>
#include <sofa/helper/OptionsGroup.h>
#include <sofa/core/objectmodel/DataFileName.h>

#include <sofa/helper/kdTree.inl>

#include <set>
#include "initRegistration.h"

#ifdef WIN32
#include <process.h>
#else
#include <pthread.h>
#endif

namespace sofa
{

namespace component
{

namespace forcefield
{

using helper::vector;
using namespace sofa::defaulttype;

template<class DataTypes>
class ClosestPointRegistrationForceFieldInternalData
{
public:
};


template<class DataTypes>
class ClosestPointRegistrationForceField : public core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(ClosestPointRegistrationForceField,DataTypes),SOFA_TEMPLATE(core::behavior::ForceField, DataTypes));

    typedef core::behavior::ForceField<DataTypes> Inherit;

    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef Data<typename DataTypes::VecCoord> DataVecCoord;
    typedef Data<typename DataTypes::VecDeriv> DataVecDeriv;

    typedef core::behavior::MechanicalState<DataTypes> MechanicalState;
    enum { N=DataTypes::spatial_dimensions };
    typedef defaulttype::Mat<N,N,Real> Mat;

    typedef typename interactionforcefield::LinearSpring<Real> Spring;
    typedef helper::fixed_array <unsigned int,3> tri;
    typedef helper::kdTree<Coord> KDT;
    typedef typename KDT::distanceSet distanceSet;

public:
    ClosestPointRegistrationForceField(core::behavior::MechanicalState<DataTypes> *mm = NULL);
    virtual ~ClosestPointRegistrationForceField();

    core::behavior::MechanicalState<DataTypes>* getObject() { return this->mstate; }

    const sofa::helper::vector< Spring >& getSprings() const {return springs.getValue();}

    // -- ForceField interface
    void reinit();
    void init();
    void addForce(const core::MechanicalParams* /*mparams*/,DataVecDeriv& f , const DataVecCoord& x , const DataVecDeriv& v);
    void addDForce(const core::MechanicalParams* mparams ,DataVecDeriv&   df , const DataVecDeriv&   dx);
    double getPotentialEnergy(const core::MechanicalParams* ,const DataVecCoord&) const { return m_potentialEnergy; }
    void addKToMatrix( const core::MechanicalParams* mparams,const sofa::core::behavior::MultiMatrixAccessor* matrix);

    Real getStiffness() const{ return ks.getValue(); }
    Real getDamping() const{ return kd.getValue(); }
    void setStiffness(Real _ks){ ks.setValue(_ks); }
    void setDamping(Real _kd){ kd.setValue(_kd); }
    Real getArrowSize() const{return showArrowSize.getValue();}
    void setArrowSize(float s){showArrowSize.setValue(s);}
    int getDrawMode() const{return drawMode.getValue();}
    void setDrawMode(int m){drawMode.setValue(m);}

    void draw(const core::visual::VisualParams* vparams);

    // -- Modifiers

    void clearSprings(int reserve=0)
    {
        sofa::helper::vector<Spring>& springs = *this->springs.beginEdit();
        springs.clear();
        if (reserve) springs.reserve(reserve);
        this->springs.endEdit();
    }

    void removeSpring(unsigned int idSpring)
    {
        if (idSpring >= (this->springs.getValue()).size())
            return;

        sofa::helper::vector<Spring>& springs = *this->springs.beginEdit();
        springs.erase(springs.begin() +idSpring );
        this->springs.endEdit();
    }

    void addSpring(int m1, SReal ks, SReal kd )
    {
        springs.beginEdit()->push_back(Spring(m1,-1,ks,kd,0));
        springs.endEdit();
    }

    void addSpring(const Spring & spring)
    {
        springs.beginEdit()->push_back(spring);
        springs.endEdit();
    }

protected :
    void updateClosestPoints();

    vector<Mat>  dfdx;
    VecCoord closestPos;
    vector<unsigned int>  cnt;
    double m_potentialEnergy;

    Real min,max;

    /// Accumulate the spring force and compute and store its stiffness
    virtual void addSpringForce(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v, int i, const Spring& spring);
    /// Apply the stiffness, i.e. accumulate df given dx
    virtual void addSpringDForce(VecDeriv& df,const  VecDeriv& dx, int i, const Spring& spring, double kFactor, double bFactor);

    Data<Real> ks;
    Data<Real> kd;
    Data<unsigned int> cacheSize;
    Data<Real> blendingFactor;
    Data<Real> outlierThreshold;
    Data<Real> normalThreshold;
    Data<bool> projectToPlane;
    Data<bool> rejectBorders;
    Data<sofa::helper::vector<Spring> > springs;

    // source mesh data
    Data< helper::vector< tri > > sourceTriangles;
    Data< VecCoord > sourceNormals;
    vector< distanceSet >  closestSource; // CacheSize-closest target points from source
    vector< Real > cacheDist;	vector< Real > cacheDist2; VecCoord previousX; // storage for cache acceleration
    KDT sourceKdTree;
    vector< bool > sourceBorder;
    vector< bool > sourceIgnored;  // flag ignored vertices
    vector< bool > targetIgnored;  // flag ignored vertices
    void initSource(); // built k-d tree and identify border vertices

    // target mesh data
    Data< VecCoord > targetPositions;
    Data< VecCoord > targetNormals;
    Data< helper::vector< tri > > targetTriangles;
    vector< distanceSet >  closestTarget; // CacheSize-closest source points from target
    KDT targetKdTree;
    vector< bool > targetBorder;
    void initTarget();  // built k-d tree and identify border vertices

    Data<float> showArrowSize;
    Data<int> drawMode; //Draw Mode: 0=Line - 1=Cylinder - 2=Arrow
    Data<bool> drawColorMap;

    void detectBorder(vector<bool> &border,const helper::vector< tri > &triangles);
};


#if defined(SOFA_EXTERN_TEMPLATE) && !defined(CLOSESTPOINTREGISTRATIONFORCEFIELD_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_REGISTRATION_API ClosestPointRegistrationForceField<defaulttype::Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_REGISTRATION_API ClosestPointRegistrationForceField<defaulttype::Vec3fTypes>;
#endif
#endif


} //

} //

} // namespace sofa

#endif
