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
#ifndef SOFA_COMPONENT_INTERACTIONFORCEFIELD_JOINTSPRINGFORCEFIELD_INL
#define SOFA_COMPONENT_INTERACTIONFORCEFIELD_JOINTSPRINGFORCEFIELD_INL

#include <sofa/component/interactionforcefield/JointSpringForceField.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/behavior/PairInteractionForceField.inl>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/helper/io/MassSpringLoader.h>
#include <sofa/helper/gl/template.h>
#include <sofa/helper/gl/Cylinder.h>
#include <sofa/helper/gl/BasicShapes.h>
#include <sofa/helper/gl/Axis.h>
#include <sofa/helper/system/config.h>
#include <cassert>
#include <iostream>



namespace sofa
{

namespace component
{

namespace interactionforcefield
{


template<class DataTypes>
JointSpringForceField<DataTypes>::JointSpringForceField(MechanicalState* object1, MechanicalState* object2)
    : Inherit(object1, object2)
    , infile(NULL)
    , outfile(NULL)
    , springs(initData(&springs,"spring","pairs of indices, stiffness, damping, rest length"))
    , f_outfilename( initData(&f_outfilename, "outfile", "output file name"))
    , f_infilename( initData(&f_infilename, "infile", "input file containing constant joint force"))
    , f_period( initData(&f_period, (Real)0.0, "period", "period between outputs"))
    , f_reinit( initData(&f_reinit, false, "reinit", "flag enabling reinitialization of the output file at each timestep"))
    , lastTime((Real)0.0)
    , showLawfulTorsion(initData(&showLawfulTorsion, false, "showLawfulTorsion", "display the lawful part of the joint rotation"))
    , showExtraTorsion(initData(&showExtraTorsion, false, "showExtraTorsion", "display the illicit part of the joint rotation"))
    , showFactorSize(initData(&showFactorSize, (Real)1.0, "showFactorSize", "modify the size of the debug information of a given factor" ))
{
}

template<class DataTypes>
JointSpringForceField<DataTypes>::JointSpringForceField()
    : infile(NULL)
    , outfile(NULL)
    , springs(initData(&springs,"spring","pairs of indices, stiffness, damping, rest length"))
    , f_outfilename( initData(&f_outfilename, "outfile", "output file name"))
    , f_infilename( initData(&f_infilename, "infile", "input file containing constant joint force"))
    , f_period( initData(&f_period, (Real)0.0, "period", "period between outputs"))
    , f_reinit( initData(&f_reinit, false, "reinit", "flag enabling reinitialization of the output file at each timestep"))
    , lastTime((Real)0.0)
    , showLawfulTorsion(initData(&showLawfulTorsion, false, "showLawfulTorsion", "display the lawful part of the joint rotation"))
    , showExtraTorsion(initData(&showExtraTorsion, false, "showExtraTorsion", "display the illicit part of the joint rotation"))
    , showFactorSize(initData(&showFactorSize, (Real)1.0, "showFactorSize", "modify the size of the debug information of a given factor" ))
{
}

template<class DataTypes>
JointSpringForceField<DataTypes>::~JointSpringForceField()
{
    if (outfile) 	  delete outfile;
    if (infile) 	  delete infile;
}


template <class DataTypes>
void JointSpringForceField<DataTypes>::init()
{
    this->Inherit::init();

    const std::string& outfilename = f_outfilename.getFullPath();
    if (!outfilename.empty())
    {
        outfile = new std::ofstream(outfilename.c_str());
        if( !outfile->is_open() )
        {
            serr << "Error creating file "<<outfilename<<sendl;
            delete outfile;
            outfile = NULL;
        }
    }

    const std::string& infilename = f_infilename.getFullPath();
    if (!infilename.empty())
    {
        infile = new std::ifstream(infilename.c_str());
        if( !infile->is_open() )
        {
            serr << "Error opening file "<<infilename<<sendl;
            delete infile;
            infile = NULL;
        }
    }
}


template <class DataTypes>
void JointSpringForceField<DataTypes>::bwdInit()
{
    //   this->Inherit::bwdInit();

    const VecCoord& x1= *this->mstate1->getX();

    const VecCoord& x2= *this->mstate2->getX();
    sofa::helper::vector<Spring> &springsVector=*(springs.beginEdit());
    for (unsigned int i=0; i<springs.getValue().size(); ++i)
    {
        Spring &s=springsVector[i];
        if (s.needToInitializeTrans)
        {
            s.initTrans = x2[s.m2].getCenter() - x1[s.m1].getCenter();
        }
        if (s.needToInitializeRot)
        {
            s.initRot = x2[s.m2].getOrientation()*x1[s.m1].getOrientation().inverse();
        }
    }
    springs.endEdit();
}



static const double PI=3.14159265358979323846264338327950288;

template<class DataTypes>
void JointSpringForceField<DataTypes>::projectTorsion(Spring& spring)
{
    Real pi2=(Real)2.*(Real)PI;

    //std::cout<<"torsion:=";

    for (unsigned int i=0; i<3; i++)
    {
        if (!spring.freeMovements[3+i]) // hard constraint
        {
            spring.lawfulTorsion[i]=0;
        }
        else if(spring.torsion[i]>spring.limitAngles[i*2] && spring.torsion[i]<spring.limitAngles[i*2+1]) // inside limits
        {
            spring.lawfulTorsion[i]=spring.torsion[i];
        }
        else // outside limits
        {
            Real d1,d2;
            if(spring.torsion[i]>0)
            {
                d1=spring.torsion[i]-spring.limitAngles[i*2+1];
                d2=spring.limitAngles[i*2]+pi2-spring.torsion[i];
                if(d1<d2) spring.lawfulTorsion[i]=spring.limitAngles[i*2+1];
                else spring.lawfulTorsion[i]=spring.limitAngles[i*2];
            }
            else
            {
                d1=spring.torsion[i]-spring.limitAngles[i*2+1]+pi2;
                d2=spring.limitAngles[i*2]-spring.torsion[i];
                if(d1<d2) spring.lawfulTorsion[i]=spring.limitAngles[i*2+1];
                else spring.lawfulTorsion[i]=spring.limitAngles[i*2];
            }
        }
    }
}

template<class DataTypes>
void JointSpringForceField<DataTypes>::addSpringForce( double& /*potentialEnergy*/, VecDeriv& f1, const VecCoord& p1, const VecDeriv& v1, VecDeriv& f2, const VecCoord& p2, const VecDeriv& v2, int , /*const*/ Spring& spring)
{

    Deriv constantForce;
    if(infile)
    {
        if (infile->eof())	{ infile->clear(); infile->seekg(0); }
        std::string line;  getline(*infile, line);
        std::istringstream str(line);
        str >> constantForce;
    }


    int a = spring.m1;
    int b = spring.m2;

    spring.ref = p1[a].getOrientation();

    //compute p2 position and velocity, relative to p1 referential
    Coord Mp1p2 = p2[b] - p1[a];
    Deriv Vp1p2 = v2[b] - v1[a];

    // offsets
    Mp1p2.getCenter() -= spring.initTrans;
    Mp1p2.getOrientation() = Mp1p2.getOrientation() * spring.initRot;
    Mp1p2.getOrientation().normalize();

    // get relative orientation in axis/angle format
    Real phi; Mp1p2.getOrientation().quatToAxis(spring.torsion,phi);
    Real pi2=(Real)2.*(Real)PI; while(phi<-PI) phi+=pi2; while(phi>PI) phi-=pi2; 		// remove modulo(2PI) from torsion
    spring.torsion*=phi;

    //compute forces
    for (unsigned int i=0; i<3; i++) spring.KT[i]=spring.freeMovements[i]==0?spring.hardStiffnessTrans:spring.softStiffnessTrans;
    Vector fT0 = spring.ref.rotate(getVCenter(constantForce) + spring.KT.linearProduct(spring.ref.inverseRotate(Mp1p2.getCenter()))) + getVCenter(Vp1p2)*spring.kd;

    // comput torques
    for (unsigned int i=0; i<3; i++) spring.KR[i]=spring.freeMovements[3+i]==0?spring.hardStiffnessRot:spring.softStiffnessRot;
    Vector fR0;
    if(!spring.freeMovements[3] && !spring.freeMovements[4] && !spring.freeMovements[5]) // encastrement
    {
        for (unsigned int i=0; i<3; i++) fR0[i]=spring.torsion[i]*spring.KR[i];
    }
    else if(spring.freeMovements[3] && !spring.freeMovements[4] && !spring.freeMovements[5] && spring.torsion[0]>spring.limitAngles[0] && spring.torsion[0]<spring.limitAngles[1]) // pivot /x
    {
        Mat M;
        Mp1p2.writeRotationMatrix(M);
        Real crossnorm=sqrt(M[1][0]*M[1][0]+M[2][0]*M[2][0]);
        Real thet=atan2(crossnorm,M[0][0]);
        fR0[0]=spring.torsion[0]*spring.KR[0]; // soft constraint
        fR0[1]=-M[2][0]*thet*spring.KR[1];
        fR0[2]=M[1][0]*thet*spring.KR[2];
    }
    else if(!spring.freeMovements[3] && spring.freeMovements[4] && !spring.freeMovements[5] && spring.torsion[1]>spring.limitAngles[2] && spring.torsion[1]<spring.limitAngles[3]) // pivot /y
    {
        Mat M;
        Mp1p2.writeRotationMatrix(M);
        Real crossnorm=sqrt(M[0][1]*M[0][1]+M[2][1]*M[2][1]);
        Real thet=atan2(crossnorm,M[1][1]);
        fR0[0]=M[2][1]*thet*spring.KR[0];
        fR0[1]=spring.torsion[1]*spring.KR[1]; // soft constraint
        fR0[2]=-M[0][1]*thet*spring.KR[2];
    }
    else if(!spring.freeMovements[3] && !spring.freeMovements[4] && spring.freeMovements[5] && spring.torsion[2]>spring.limitAngles[4] && spring.torsion[2]<spring.limitAngles[5]) // pivot /z
    {
        Mat M;
        Mp1p2.writeRotationMatrix(M);
        Real crossnorm=sqrt(M[1][2]*M[1][2]+M[0][2]*M[0][2]);
        Real thet=atan2(crossnorm,M[2][2]);
        fR0[0]=-M[1][2]*thet*spring.KR[0];
        fR0[1]=M[0][2]*thet*spring.KR[1];
        fR0[2]=spring.torsion[2]*spring.KR[2]; // soft constraint
    }
    else // general case
    {
        // update lawfull torsion
        projectTorsion(spring);
        Vector extraTorsion=spring.torsion-spring.lawfulTorsion;
        Real psi=extraTorsion.norm(); extraTorsion/=psi;		while(psi<-PI) psi+=pi2; while(psi>PI) psi-=pi2; 		extraTorsion*=psi;

        for (unsigned int i=0; i<3; i++)
            if(spring.freeMovements[3+i] && spring.torsion[i]!=spring.lawfulTorsion[i]) // outside limits
            {
                spring.KR[i]=spring.blocStiffnessRot;
                fR0[i]=extraTorsion[i]*spring.KR[i];
            }
            else fR0[i]=spring.torsion[i]*spring.KR[i]; // hard constraint or soft constraint inside limits
    }


    Vector fR = spring.ref.rotate(getVOrientation(constantForce) + fR0) + getVOrientation(Vp1p2)*spring.kd;

    // add force
    const Deriv force(fT0, fR );
    f1[a] += force;
    this->mask1->insertEntry(a);
    f2[b] -= force;
    this->mask2->insertEntry(b);

    // write output file
    if (outfile)
    {
        if(f_reinit.getValue()) outfile->seekp(std::ios::beg);

        double time = this->getContext()->getTime();
        if (time >= (lastTime + f_period.getValue()))
        {
            lastTime += f_period.getValue();
            (*outfile) << "T= "<< time << "\n";

            const Coord xrel(spring.ref.inverseRotate(Mp1p2.getCenter()), Mp1p2.getOrientation());
            (*outfile) << "  Xrel= " << xrel << "\n";

            (*outfile) << "  Vrel= " << Vp1p2 << "\n";

            const Deriv frel(spring.KT.linearProduct(spring.ref.inverseRotate(Mp1p2.getCenter())) , fR0 );
            (*outfile) << "  Frel= " << frel << "\n";

            const Deriv damp(getVCenter(Vp1p2)*spring.kd , getVOrientation(Vp1p2)*spring.kd );
            (*outfile) << "  Fdamp= " << damp << "\n";

            if(infile) (*outfile) << "  Fconstant= " << constantForce << "\n";

            (*outfile) << "  F= " << force << "\n";

            if(f_reinit.getValue()) (*outfile) << "\n\n\n\n\n";
            outfile->flush();
        }
    }

}



template<class DataTypes>
void JointSpringForceField<DataTypes>::addSpringDForce(VecDeriv& f1, const VecDeriv& dx1, VecDeriv& f2, const VecDeriv& dx2, int , /*const*/ Spring& spring, Real kFactor)
{
    const Deriv Mdx1dx2 = dx2[spring.m2]-dx1[spring.m1];

    Vector df = spring.ref.rotate(spring.KT.linearProduct(spring.ref.inverseRotate(getVCenter(Mdx1dx2))));
    Vector dR = spring.ref.rotate(spring.KR.linearProduct(spring.ref.inverseRotate(getVOrientation(Mdx1dx2))));

    const Deriv dforce(df,dR);

    f1[spring.m1] += dforce * kFactor;
    f2[spring.m2] -= dforce * kFactor;
}

template<class DataTypes>
void JointSpringForceField<DataTypes>::addForce(const core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecDeriv& data_f1, DataVecDeriv& data_f2, const DataVecCoord& data_x1, const DataVecCoord& data_x2, const DataVecDeriv& data_v1, const DataVecDeriv& data_v2 )
{

    VecDeriv&       f1 = *data_f1.beginEdit();
    const VecCoord& x1 =  data_x1.getValue();
    const VecDeriv& v1 =  data_v1.getValue();
    VecDeriv&       f2 = *data_f2.beginEdit();
    const VecCoord& x2 =  data_x2.getValue();
    const VecDeriv& v2 =  data_v2.getValue();

    helper::vector<Spring>& springs = *this->springs.beginEdit();

    f1.resize(x1.size());
    f2.resize(x2.size());
    m_potentialEnergy = 0;
    for (unsigned int i=0; i<springs.size(); i++)
    {
        this->addSpringForce(m_potentialEnergy,f1,x1,v1,f2,x2,v2, i, springs[i]);
    }
    this->springs.endEdit();

    data_f1.endEdit();
    data_f2.endEdit();
}

template<class DataTypes>
void JointSpringForceField<DataTypes>::addDForce(const core::MechanicalParams *mparams /* PARAMS FIRST */, DataVecDeriv& data_df1, DataVecDeriv& data_df2, const DataVecDeriv& data_dx1, const DataVecDeriv& data_dx2)
{
    VecDeriv&        df1 = *data_df1.beginEdit();
    VecDeriv&        df2 = *data_df2.beginEdit();
    const VecDeriv&  dx1 =  data_dx1.getValue();
    const VecDeriv&  dx2 =  data_dx2.getValue();


    df1.resize(dx1.size());
    df2.resize(dx2.size());

    Real kFactor = (Real)mparams->kFactor();

    helper::vector<Spring>& springs = *this->springs.beginEdit();
    for (unsigned int i=0; i<springs.size(); i++)
    {
        this->addSpringDForce(df1, dx1, df2, dx2, i, springs[i], kFactor);
    }
    this->springs.endEdit();

    data_df1.endEdit();
    data_df2.endEdit();
}

template<class DataTypes>
void JointSpringForceField<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    if (!((this->mstate1 == this->mstate2)?vparams->displayFlags().getShowForceFields():vparams->displayFlags().getShowInteractionForceFields())) return;
    const VecCoord& p1 = *this->mstate1->getX();
    const VecCoord& p2 = *this->mstate2->getX();

    glDisable(GL_LIGHTING);
    bool external = (this->mstate1!=this->mstate2);
    const helper::vector<Spring>& springs = this->springs.getValue();

    for (unsigned int i=0; i<springs.size(); i++)
    {
        Real d = (p2[springs[i].m2]-p1[springs[i].m1]).getCenter().norm();
        if (external)
        {
            if (d<springs[i].initTrans.norm()*0.9999)
                glColor4f(1,0,0,1);
            else
                glColor4f(0,1,0,1);
        }
        else
        {
            if (d<springs[i].initTrans.norm()*0.9999)
                glColor4f(1,0.5f,0,1);
            else
                glColor4f(0,1,0.5f,1);
        }
        glBegin(GL_LINES);
        helper::gl::glVertexT(p1[springs[i].m1].getCenter());
        helper::gl::glVertexT(p2[springs[i].m2].getCenter());
        glEnd();

        if(springs[i].freeMovements[3] == 1)
        {
            helper::gl::Cylinder::draw(p1[springs[i].m1].getCenter(), p1[springs[i].m1].getOrientation(), Vector((Real)(1.0*showFactorSize.getValue()),0,0));
        }
        if(springs[i].freeMovements[4] == 1)
        {
            helper::gl::Cylinder::draw(p1[springs[i].m1].getCenter(), p1[springs[i].m1].getOrientation(), Vector(0,(Real)(1.0*showFactorSize.getValue()),0));
        }
        if(springs[i].freeMovements[5] == 1)
        {
            helper::gl::Cylinder::draw(p1[springs[i].m1].getCenter(), p1[springs[i].m1].getOrientation(), Vector(0,0,(Real)(1.0*showFactorSize.getValue())) );
        }

        //---debugging
        if (showLawfulTorsion.getValue())
            helper::gl::drawArrow(p1[springs[i].m1].getCenter(), p1[springs[i].m1].pointToParent(springs[i].lawfulTorsion), (float)(0.5*showFactorSize.getValue()));
        if (showExtraTorsion.getValue())
            helper::gl::drawArrow(p1[springs[i].m1].getCenter(), p1[springs[i].m1].pointToParent(springs[i].torsion-springs[i].lawfulTorsion), (float)(0.5*showFactorSize.getValue()));
    }

}

} // namespace interactionforcefield

} // namespace component

} // namespace sofa

#endif  /* SOFA_COMPONENT_INTERACTIONFORCEFIELD_JOINTSPRINGFORCEFIELD_INL */

