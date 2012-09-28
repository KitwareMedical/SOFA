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
#ifndef SOFA_COMPONENT_FORCEFIELD_INTERACTION_ELLIPSOIDFORCEFIELD_H
#define SOFA_COMPONENT_FORCEFIELD_INTERACTION_ELLIPSOIDFORCEFIELD_H

#include <sofa/core/behavior/MixedInteractionForceField.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/MechanicalParams.h>



namespace sofa
{

namespace component
{

namespace interactionforcefield
{
using namespace sofa::defaulttype;
using namespace sofa::core;

/// This class can be overridden if needed for additionnal storage within template specializations.
template<class DataTypes1, class DataTypes2>
class InteractionEllipsoidForceFieldInternalData
{
public:
};

template<typename TDataTypes1, typename TDataTypes2>
class InteractionEllipsoidForceField : public core::behavior::MixedInteractionForceField<TDataTypes1, TDataTypes2>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(InteractionEllipsoidForceField, TDataTypes1, TDataTypes2), SOFA_TEMPLATE2(core::behavior::MixedInteractionForceField, TDataTypes1, TDataTypes2));

    typedef core::behavior::MixedInteractionForceField<TDataTypes1, TDataTypes2> Inherit;
    typedef TDataTypes1 DataTypes1;
    typedef typename DataTypes1::VecCoord VecCoord1;
    typedef typename DataTypes1::VecDeriv VecDeriv1;
    typedef typename DataTypes1::Coord    Coord1;
    typedef typename DataTypes1::Deriv    Deriv1;
    typedef typename DataTypes1::Real     Real1;
    typedef TDataTypes2 DataTypes2;
    typedef typename DataTypes2::VecCoord VecCoord2;
    typedef typename DataTypes2::VecDeriv VecDeriv2;
    typedef typename DataTypes2::Coord    Coord2;
    typedef typename DataTypes2::Deriv    Deriv2;
    typedef typename DataTypes2::Real     Real2;
    typedef sofa::helper::ParticleMask ParticleMask;

    typedef core::objectmodel::Data<VecCoord1>    DataVecCoord1;
    typedef core::objectmodel::Data<VecDeriv1>    DataVecDeriv1;
    typedef core::objectmodel::Data<VecCoord2>    DataVecCoord2;
    typedef core::objectmodel::Data<VecDeriv2>    DataVecDeriv2;



    enum { N=DataTypes1::spatial_dimensions };
    typedef defaulttype::Mat<N,N,Real1> Mat;
protected:
    class Contact
    {
    public:
        int index;
        Deriv1 pos,force;
        Vec<3,SReal> bras_levier;
        Mat m;
        Contact( int index=0, const Mat& m=Mat())
            : index(index), m(m)
        {
        }

        inline friend std::istream& operator >> ( std::istream& in, Contact& c )
        {
            in>>c.index>>c.m;
            return in;
        }

        inline friend std::ostream& operator << ( std::ostream& out, const Contact& c )
        {
            out << c.index << " " << c.m ;
            return out;
        }

    };

    Data<sofa::helper::vector<Contact> > contacts;

    InteractionEllipsoidForceFieldInternalData<DataTypes1, DataTypes2> data;

    bool calcF(const Coord1& p1, const Deriv1 &v1, Deriv1& f1,Mat& dfdx);
    void initCalcF();

public:

    Data<VecCoord1> center;
    Data<VecCoord1> vradius;
    Data<Real1> stiffness;
    Data<Real1> damping;
    Data<defaulttype::Vec3f> color;
    Data<bool> bDraw;
    Data<int> object2_dof_index;
    Data<bool> object2_forces;
    Data<bool> object2_invert;

protected:
    InteractionEllipsoidForceField()
        : contacts(initData(&contacts,"contacts", "Contacts"))
        , center(initData(&center, "center", "ellipsoid center"))
        , vradius(initData(&vradius, "vradius", "ellipsoid radius"))
        , stiffness(initData(&stiffness, (Real1)500, "stiffness", "force stiffness (positive to repulse outward, negative inward)"))
        , damping(initData(&damping, (Real1)5, "damping", "force damping"))
        , color(initData(&color, defaulttype::Vec3f(0.0f,0.5f,1.0f), "color", "ellipsoid color"))
        , bDraw(initData(&bDraw, true, "draw", "enable/disable drawing of the ellipsoid"))
        , object2_dof_index(initData(&object2_dof_index, (int)0, "object2_dof_index", "Dof index of object 2 where the forcefield is attached"))
        , object2_forces(initData(&object2_forces, true, "object2_forces", "enable/disable propagation of forces to object 2"))
        , object2_invert(initData(&object2_invert, false, "object2_invert", "inverse transform from object 2 (use when object 1 is in local coordinates within a frame defined by object 2)"))
    {
    }
public:
    void setStiffness(Real1 stiff)
    {
        stiffness.setValue( stiff );
    }

    void setDamping(Real1 damp)
    {
        damping.setValue( damp );
    }

    virtual void addForce(const MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv1& f1, DataVecDeriv2& f2, const DataVecCoord1& x1, const DataVecCoord2& x2, const DataVecDeriv1& v1, const DataVecDeriv2& v2);
    ///SOFA_DEPRECATED_ForceField <<<virtual void addForce(VecDeriv1& f1, VecDeriv2& f2, const VecCoord1& p1, const VecCoord2& p2, const VecDeriv1& v1, const VecDeriv2& v2);

    virtual void addForce2(DataVecDeriv1& f1, DataVecDeriv2& f2, const DataVecCoord1& p1, const DataVecCoord2& p2, const DataVecDeriv1& v1, const DataVecDeriv2& v2);
    ///SOFA_DEPRECATED_ForceField <<<virtual void addForce2(VecDeriv1& f1, VecDeriv2& f2, const VecCoord1& p1, const VecCoord2& p2, const VecDeriv1& v1, const VecDeriv2& v2);

    virtual void addDForce(const MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv1& df1, DataVecDeriv2& df2, const DataVecDeriv1& dx1, const DataVecDeriv2& dx2);
    ///SOFA_DEPRECATED_ForceField <<<virtual void addDForce(VecDeriv1& df1, VecDeriv2& df2, const VecDeriv1& dx1, const VecDeriv2& dx2);

    virtual double getPotentialEnergy(const MechanicalParams* mparams /* PARAMS FIRST */, const DataVecCoord1& x1, const DataVecCoord2& x2)const;
    ///SOFA_DEPRECATED_ForceField <<<virtual double getPotentialEnergy(const VecCoord1& x1, const VecCoord2& x2) const;

    void init();
    void reinit();

    void draw(const core::visual::VisualParams* vparams);

protected:
    struct TempVars
    {
        unsigned int nelems;
        VecCoord1 vcenter; // center in the local frame ;
        VecCoord1 vr;
        Real1 stiffness;
        Real1 stiffabs;
        Real1 damping;
        VecCoord1 vinv_r2;
        Coord2 pos6D;
    } vars;
};

} // namespace interactionforcefield

} // namespace component

} // namespace sofa


#endif
