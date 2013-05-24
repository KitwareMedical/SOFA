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
// Author: François Faure, INRIA-UJF, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
#ifndef SOFA_COMPONENT_INTERACTIONFORCEFIELD_SPRINGFORCEFIELD_H
#define SOFA_COMPONENT_INTERACTIONFORCEFIELD_SPRINGFORCEFIELD_H

#include <sofa/core/behavior/PairInteractionForceField.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/helper/vector.h>
#include <sofa/helper/accessor.h>
#include <sofa/component/component.h>

#include <sofa/core/objectmodel/DataFileName.h>

namespace sofa
{

namespace component
{

namespace interactionforcefield
{

using namespace sofa::defaulttype;

/// This class contains the description of one linear spring
template<class T>
class LinearSpring
{
public:
    typedef T Real;
    int     m1, m2;  ///< the two extremities of the spring: masses m1 and m2
    Real  ks;      ///< spring stiffness
    Real  kd;      ///< damping factor
    Real  initpos; ///< rest length of the spring

    LinearSpring(int m1=0, int m2=0, double ks=0.0, double kd=0.0, double initpos=0.0)
        : m1(m1), m2(m2), ks((Real)ks), kd((Real)kd), initpos((Real)initpos)
    {
    }

    LinearSpring(int m1, int m2, float ks, float kd=0, float initpos=0)
        : m1(m1), m2(m2), ks((Real)ks), kd((Real)kd), initpos((Real)initpos)
    {
    }

    inline friend std::istream& operator >> ( std::istream& in, LinearSpring<Real>& s )
    {
        in>>s.m1>>s.m2>>s.ks>>s.kd>>s.initpos;
        return in;
    }

    inline friend std::ostream& operator << ( std::ostream& out, const LinearSpring<Real>& s )
    {
        out<<s.m1<<" "<<s.m2<<" "<<s.ks<<" "<<s.kd<<" "<<s.initpos<<"\n";
        return out;
    }

};


/// This class can be overridden if needed for additionnal storage within template specializations.
template<class DataTypes>
class SpringForceFieldInternalData
{
public:
};

/// Set of simple springs between particles
template<class DataTypes>
class SpringForceField : public core::behavior::PairInteractionForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(SpringForceField,DataTypes), SOFA_TEMPLATE(core::behavior::PairInteractionForceField,DataTypes));

    typedef typename core::behavior::PairInteractionForceField<DataTypes> Inherit;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename Coord::value_type Real;
    typedef core::objectmodel::Data<VecDeriv>    DataVecDeriv;
    typedef core::objectmodel::Data<VecCoord>    DataVecCoord;

    typedef helper::ReadAccessor< Data< VecCoord > > RDataRefVecCoord;
    typedef helper::WriteAccessor< Data< VecCoord > > WDataRefVecCoord;
    typedef helper::ReadAccessor< Data< VecDeriv > > RDataRefVecDeriv;
    typedef helper::WriteAccessor< Data< VecDeriv > > WDataRefVecDeriv;

    typedef core::behavior::MechanicalState<DataTypes> MechanicalState;

    typedef LinearSpring<Real> Spring;

protected:
    bool maskInUse;
    Real m_potentialEnergy;
    Data<SReal> ks;
    Data<SReal> kd;
    Data<float> showArrowSize;
    Data<int> drawMode; //Draw Mode: 0=Line - 1=Cylinder - 2=Arrow
    Data<sofa::helper::vector<Spring> > springs;
    core::objectmodel::DataFileName fileSprings;
    class Loader;

    SpringForceFieldInternalData<DataTypes> data;
    friend class SpringForceFieldInternalData<DataTypes>;

    virtual void addSpringForce(Real& potentialEnergy, VecDeriv& f1, const VecCoord& p1, const VecDeriv& v1, VecDeriv& f2, const VecCoord& p2, const VecDeriv& v2, int /*i*/, const Spring& spring);
    void updateMaskStatus();


    SpringForceField(MechanicalState* object1, MechanicalState* object2, SReal _ks=100.0, SReal _kd=5.0);
    SpringForceField(SReal _ks=100.0, SReal _kd=5.0);

public:
    bool load(const char *filename);

    core::behavior::MechanicalState<DataTypes>* getObject1() { return this->mstate1; }
    core::behavior::MechanicalState<DataTypes>* getObject2() { return this->mstate2; }

    const sofa::helper::vector< Spring >& getSprings() const {return springs.getValue();}

    virtual void reinit();
    virtual void init();

    virtual void addForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& f1, DataVecDeriv& f2, const DataVecCoord& x1, const DataVecCoord& x2, const DataVecDeriv& v1, const DataVecDeriv& v2);
    virtual void addDForce(const core::MechanicalParams* /* PARAMS FIRST */, DataVecDeriv& df1, DataVecDeriv& df2, const DataVecDeriv& dx1, const DataVecDeriv& dx2 );
    virtual double getPotentialEnergy(const core::MechanicalParams* /* PARAMS FIRST */, const DataVecCoord&, const DataVecCoord& ) const {return m_potentialEnergy; }


    virtual void addKToMatrix(sofa::defaulttype::BaseMatrix * /*mat*/, double /*kFact*/, unsigned int &/*offset*/);

    SReal getStiffness() const { return ks.getValue(); }
    SReal getDamping() const { return kd.getValue(); }
    void setStiffness(SReal _ks) { ks.setValue(_ks); }
    void setDamping(SReal _kd) { kd.setValue(_kd); }
    SReal getArrowSize() const {return showArrowSize.getValue();}
    void setArrowSize(float s) {showArrowSize.setValue(s);}
    int getDrawMode() const {return drawMode.getValue();}
    void setDrawMode(int m) {drawMode.setValue(m);}

    virtual void draw(const core::visual::VisualParams* vparams);

    // -- Modifiers

    void clear(int reserve=0)
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

        updateMaskStatus();
    }

    void addSpring(int m1, int m2, SReal ks, SReal kd, SReal initlen)
    {
        springs.beginEdit()->push_back(Spring(m1,m2,ks,kd,initlen));
        springs.endEdit();
        updateMaskStatus();
    }

    void addSpring(const Spring & spring)
    {
        springs.beginEdit()->push_back(spring);
        springs.endEdit();
        updateMaskStatus();
    }

    virtual void handleTopologyChange(core::topology::Topology *topo);

    virtual bool useMask() const ;

    /// initialization to export kinetic, potential energy  and force intensity to gnuplot files format
    virtual void initGnuplot(const std::string path);

    /// export kinetic and potential energy state at "time" to a gnuplot file
    virtual void exportGnuplot(double time);

    protected:
    /// stream to export Potential Energy to gnuplot files
    std::ofstream* m_gnuplotFileEnergy;
};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_BUILD_DEFORMABLE)
#ifndef SOFA_FLOAT
extern template class SOFA_DEFORMABLE_API SpringForceField<defaulttype::Vec3dTypes>;
extern template class SOFA_DEFORMABLE_API SpringForceField<defaulttype::Vec2dTypes>;
extern template class SOFA_DEFORMABLE_API SpringForceField<defaulttype::Vec1dTypes>;
extern template class SOFA_DEFORMABLE_API SpringForceField<defaulttype::Vec6dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_DEFORMABLE_API SpringForceField<defaulttype::Vec3fTypes>;
extern template class SOFA_DEFORMABLE_API SpringForceField<defaulttype::Vec2fTypes>;
extern template class SOFA_DEFORMABLE_API SpringForceField<defaulttype::Vec1fTypes>;
extern template class SOFA_DEFORMABLE_API SpringForceField<defaulttype::Vec6fTypes>;
#endif
#endif

} // namespace interactionforcefield

} // namespace component

} // namespace sofa

#endif  /* SOFA_COMPONENT_INTERACTIONFORCEFIELD_SPRINGFORCEFIELD_H */
