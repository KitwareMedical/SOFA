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
#ifndef SOFA_COMPONENT_CONSTRAINTSET_BILATERALINTERACTIONCONSTRAINT_INL
#define SOFA_COMPONENT_CONSTRAINTSET_BILATERALINTERACTIONCONSTRAINT_INL

#include <sofa/component/constraintset/BilateralInteractionConstraint.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/behavior/Constraint.inl>

#include <sofa/defaulttype/Vec.h>
#include <sofa/helper/gl/template.h>

#include <algorithm> // for std::min

namespace sofa
{

namespace component
{

namespace constraintset
{

template<class DataTypes>
void BilateralInteractionConstraint<DataTypes>::init()
{
    Inherit1::init();
    assert(this->mstate1);
    assert(this->mstate2);
    prevForces.clear();
    iteration = 0;
    activated = (activateAtIteration.getValue() >= 0 && activateAtIteration.getValue() <= iteration);
}

template<class DataTypes>
void BilateralInteractionConstraint<DataTypes>::reinit()
{
    prevForces.clear();
    activated = (activateAtIteration.getValue() >= 0 && activateAtIteration.getValue() <= iteration);
}

template<class DataTypes>
void BilateralInteractionConstraint<DataTypes>::buildConstraintMatrix(const core::ConstraintParams*, DataMatrixDeriv &c1_d, DataMatrixDeriv &c2_d, unsigned int &constraintId
        , const DataVecCoord &/*x1*/, const DataVecCoord &/*x2*/)
{
    if (!activated)
        return;
    
    unsigned minp = std::min(m1.getValue().size(), m2.getValue().size());
    if (minp == 0)
        return;

    const helper::vector<int> &m1Indices = m1.getValue();
    const helper::vector<int> &m2Indices = m2.getValue();

    MatrixDeriv &c1 = *c1_d.beginEdit();
    MatrixDeriv &c2 = *c2_d.beginEdit();

    cid.resize(minp);

    const VecDeriv& restVector = this->restVector.getValue();

    if (!merge.getValue())
    {
        for (unsigned pid=0; pid<minp; pid++)
        {
            int tm1 = m1Indices[pid];
            int tm2 = m2Indices[pid];

            const defaulttype::Vec<3, Real> cx(1,0,0), cy(0,1,0), cz(0,0,1);

            cid[pid] = constraintId;
            constraintId += 3;

            MatrixDerivRowIterator c1_it = c1.writeLine(cid[pid]);
            c1_it.addCol(tm1, -cx);

            MatrixDerivRowIterator c2_it = c2.writeLine(cid[pid]);
            c2_it.addCol(tm2, cx);

            c1_it = c1.writeLine(cid[pid] + 1);
            c1_it.setCol(tm1, -cy);

            c2_it = c2.writeLine(cid[pid] + 1);
            c2_it.setCol(tm2, cy);

            c1_it = c1.writeLine(cid[pid] + 2);
            c1_it.setCol(tm1, -cz);

            c2_it = c2.writeLine(cid[pid] + 2);
            c2_it.setCol(tm2, cz);
        }
    }
    else
    {
        this->m_constraintIndex.setValue(constraintId);


        ///////////////// grouped constraints ///////////////
        dfree_square_total.clear();

        const DataVecCoord &d_x1 = *this->mstate1->read(core::ConstVecCoordId::position());
        const DataVecCoord &d_x2 = *this->mstate2->read(core::ConstVecCoordId::position());

        const VecCoord &x1 = d_x1.getValue();
        const VecCoord &x2 = d_x2.getValue();

        for (unsigned pid=0; pid<minp; pid++)
        {
            int tm1 = m1Indices[pid];
            int tm2 = m2Indices[pid];

            Deriv dfree_loc = x2[tm2] - x1[tm1];

            if (pid < restVector.size())
                dfree_loc -= restVector[pid];

            dfree_square_total[0]+= dfree_loc[0]*dfree_loc[0];
            dfree_square_total[1]+= dfree_loc[1]*dfree_loc[1];
            dfree_square_total[2]+= dfree_loc[2]*dfree_loc[2];
        }

        for (unsigned int i=0; i<3; i++)
        {
            if (dfree_square_total[i]>1.0e-15)
            {
                dfree_square_total[i] = sqrt(dfree_square_total[i]);
                squareXYZ[i]=derivative.getValue();
            }
            else
                squareXYZ[i]=false;
        }

        dfree.resize(minp);
        
        const DataVecCoord &d_x1free = *this->mstate1->read(core::ConstVecCoordId::freePosition());
        const DataVecCoord &d_x2free = *this->mstate2->read(core::ConstVecCoordId::freePosition());

        const VecCoord &x1free = d_x1free.getValue();
        const VecCoord &x2free = d_x2free.getValue();

        for (unsigned pid=0; pid<minp; pid++)
        {
            int tm1 = m1Indices[pid];
            int tm2 = m2Indices[pid];

            Deriv d_loc = x2[tm2] - x1[tm1];
            Deriv dfree_loc = x2free[tm2] - x1free[tm1];

            if (pid < restVector.size())
            {
                d_loc -= restVector[pid];
                dfree_loc -= restVector[pid];
            }
            dfree[pid] = dfree_loc;

            //std::cout<<" BilateralInteractionConstraint add Constraint between point "<<tm1<<" of object1 and "<< tm2<< " of object2"<<std::endl;

            const defaulttype::Vec<3, Real> cx(1.0,0,0), cy(0,1.0,0), cz(0,0,1.0);

            cid[pid] = constraintId;


            // if not grouped constraint
            // constraintId += 3;

            // contribution along x axis
            MatrixDerivRowIterator c1_it = c1.writeLine(cid[pid]);
            MatrixDerivRowIterator c2_it = c2.writeLine(cid[pid]);
            if(squareXYZ[0])
            {
                c1_it.addCol(tm1, -cx*dfree_loc[0]*2.0);
                c2_it.addCol(tm2, cx*dfree_loc[0]*2.0);
            }
            else
            {
                c1_it.addCol(tm1, -cx*sign(dfree_loc[0]) );
                c2_it.addCol(tm2, cx*sign(dfree_loc[0]));
            }


            // contribution along y axis
            c1_it = c1.writeLine(cid[pid] + 1);
            c2_it = c2.writeLine(cid[pid] + 1);
            if(squareXYZ[1])
            {

                c1_it.addCol(tm1, -cy*dfree_loc[1]*2.0);
                c2_it.addCol(tm2, cy*dfree_loc[1]*2.0);
            }
            else
            {
                c1_it.addCol(tm1, -cy*sign(dfree_loc[1]));
                c2_it.addCol(tm2, cy*sign(dfree_loc[1]));
            }

            // contribution along z axis
            c1_it = c1.writeLine(cid[pid] + 2);
            c2_it = c2.writeLine(cid[pid] + 2);
            if(squareXYZ[2])
            {
                c1_it.addCol(tm1, -cz*dfree_loc[2]*2.0);
                c2_it.addCol(tm2, cz*dfree_loc[2]*2.0);
            }
            else
            {
                c1_it.addCol(tm1, -cz*sign(dfree_loc[2]));
                c2_it.addCol(tm2, cz*sign(dfree_loc[2]));
            }
        }

        // if grouped constraint
        constraintId += 3;
    }

    c1_d.endEdit();
    c2_d.endEdit();
}


template<class DataTypes>
void BilateralInteractionConstraint<DataTypes>::getConstraintViolation(const core::ConstraintParams* cParams, defaulttype::BaseVector *v, const DataVecCoord &d_x1, const DataVecCoord &d_x2
        , const DataVecDeriv & d_v1, const DataVecDeriv & d_v2)
{
    if (!activated) return;

    const helper::vector<int> &m1Indices = m1.getValue();
    const helper::vector<int> &m2Indices = m2.getValue();

    unsigned minp = std::min(m1Indices.size(), m2Indices.size());

    const VecDeriv& restVector = this->restVector.getValue();

    if (cParams->constOrder() == core::ConstraintParams::VEL)
    {
        getVelocityViolation(v, d_x1, d_x2, d_v1, d_v2);
        return;
    }

    const VecCoord &x1 = d_x1.getValue();
    const VecCoord &x2 = d_x2.getValue();

    if (!merge.getValue())
    {
        dfree.resize(minp);

        for (unsigned pid=0; pid<minp; pid++)
        {
            dfree[pid] = x2[m2Indices[pid]] - x1[m1Indices[pid]];

            if (pid < restVector.size())
                dfree[pid] -= restVector[pid];

            v->set(cid[pid]  , dfree[pid][0]);
            v->set(cid[pid]+1, dfree[pid][1]);
            v->set(cid[pid]+2, dfree[pid][2]);
        }
    }
    else
    {
        for (unsigned pid=0; pid<minp; pid++)
        {
            dfree[pid] = x2[m2Indices[pid]] - x1[m1Indices[pid]];

            if (pid < restVector.size())
                dfree[pid] -= restVector[pid];

            for (unsigned int i=0; i<3; i++)
            {
                if(squareXYZ[i])
                    v->add(cid[pid]+i  , dfree[pid][i]*dfree[pid][i]);
                else
                {

                    v->add(cid[pid]+i  , dfree[pid][i]*sign(dfree[pid][i] ) );
                }
            }

        }
    }
}


template<class DataTypes>
void BilateralInteractionConstraint<DataTypes>::getVelocityViolation(defaulttype::BaseVector *v, const DataVecCoord &d_x1, const DataVecCoord &d_x2, const DataVecDeriv &d_v1, const DataVecDeriv &d_v2)
{
    std::cout<<"getVelocityViolation called "<<std::endl;

    const helper::vector<int> &m1Indices = m1.getValue();
    const helper::vector<int> &m2Indices = m2.getValue();

    const VecCoord &x1 = d_x1.getValue();
    const VecCoord &x2 = d_x2.getValue();
    const VecCoord &v1 = d_v1.getValue();
    const VecCoord &v2 = d_v2.getValue();

    unsigned minp = std::min(m1Indices.size(), m2Indices.size());
    const VecDeriv& restVector = this->restVector.getValue();
    std::vector<Deriv> dPrimefree;

    if (!merge.getValue())
    {
        dPrimefree.resize(minp);

        for (unsigned pid=0; pid<minp; pid++)
        {
            dPrimefree[pid] = v2[m2Indices[pid]] - v1[m1Indices[pid]];
            if (pid < restVector.size())
                dPrimefree[pid] -= restVector[pid];

            v->set(cid[pid]  , dPrimefree[pid][0]);
            v->set(cid[pid]+1, dPrimefree[pid][1]);
            v->set(cid[pid]+2, dPrimefree[pid][2]);
        }
    }
    else
    {

        dPrimefree.resize(minp);
        dfree.resize(minp);

        for (unsigned pid=0; pid<minp; pid++)
        {
            dPrimefree[pid] = v2[m2Indices[pid]] - v1[m1Indices[pid]];
            dfree[pid] = x2[m2Indices[pid]] - x1[m1Indices[pid]];

            if (pid < restVector.size())
            {
                dPrimefree[pid] -= restVector[pid];
                dfree[pid] -= restVector[pid];
            }

            std::cout<<" x2 : "<<x2[m2Indices[pid]]<<" - x1 :"<<x1[m1Indices[pid]]<<" = "<<dfree[pid]<<std::endl;
            std::cout<<" v2 : "<<v2[m2Indices[pid]]<<" - v1 :"<<v1[m1Indices[pid]]<<" = "<<dPrimefree[pid]<<std::endl;

            for (unsigned int i=0; i<3; i++)
            {
                if(squareXYZ[i])
                {
                    //std::cout<<" vel viol:"<<2*dPrimefree[pid][i]*dfree[pid][i]<<std::endl;
                    v->add(cid[pid]+i  , 2*dPrimefree[pid][i]*dfree[pid][i]);
                }
                else
                {
                    //std::cout<<" vel viol:"<<dPrimefree[pid][i]*sign(dfree[pid][i] )<<std::endl;
                    v->add(cid[pid]+i  , dPrimefree[pid][i]*sign(dfree[pid][i] ) );
                }
            }

        }
    }
}


template<class DataTypes>
void BilateralInteractionConstraint<DataTypes>::getConstraintResolution(std::vector<core::behavior::ConstraintResolution*>& resTab, unsigned int& offset)
{
    unsigned minp=std::min(m1.getValue().size(),m2.getValue().size());

    if (!merge.getValue())
    {
        prevForces.resize(minp);
        for (unsigned pid=0; pid<minp; pid++)
        {
            resTab[offset] = new BilateralConstraintResolution3Dof(&prevForces[pid]);
            offset += 3;
        }
    }
    else
    {
        prevForces.resize(1);
        resTab[offset] = new BilateralConstraintResolution3Dof(&prevForces[0]);
        offset +=3;
    }
}

template<class DataTypes>
void BilateralInteractionConstraint<DataTypes>::handleEvent(sofa::core::objectmodel::Event *event)
{
    if (sofa::core::objectmodel::KeypressedEvent* ev = dynamic_cast<sofa::core::objectmodel::KeypressedEvent*>(event))
    {
        //std::cout << "key pressed " << std::endl;
        switch(ev->getKey())
        {

        case 'A':
        case 'a':
            std::cout << "Activating constraint" << std::endl;
            activated = true;
            break;
        }
    }


    if ( /*simulation::AnimateEndEvent* ev =*/  dynamic_cast<simulation::AnimateEndEvent*>(event))
    {
        ++iteration;
        if (!activated && activateAtIteration.getValue() >= 0 && activateAtIteration.getValue() <= iteration)
        {
            std::cout << "Activating constraint" << std::endl;
            activated = true;
        }
    }
}

template<class DataTypes>
void BilateralInteractionConstraint<DataTypes>::addContact(Deriv /*norm*/, Coord P, Coord Q, Real /*contactDistance*/, int m1, int m2, Coord /*Pfree*/, Coord /*Qfree*/, long /*id*/, PersistentID /*localid*/)
{
    helper::WriteAccessor<Data<helper::vector<int> > > wm1 = this->m1;
    helper::WriteAccessor<Data<helper::vector<int> > > wm2 = this->m2;
    helper::WriteAccessor<Data<VecDeriv > > wrest = this->restVector;
    wm1.push_back(m1);
    wm2.push_back(m2);
    wrest.push_back(Q-P);
}


template<class DataTypes>
void BilateralInteractionConstraint<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
#ifndef SOFA_NO_OPENGL
    if (!vparams->displayFlags().getShowInteractionForceFields()) return;

    glDisable(GL_LIGHTING);
    glPointSize(10);
    if (activated)
        glColor4f(1,0,1,1);
    else
        glColor4f(0,1,0,1);
    glBegin(GL_POINTS);

    unsigned minp = std::min(m1.getValue().size(),m2.getValue().size());
    for (unsigned i=0; i<minp; i++)
    {
        helper::gl::glVertexT((*this->mstate1->getX())[m1.getValue()[i]]);
        helper::gl::glVertexT((*this->mstate2->getX())[m2.getValue()[i]]);
    }
    glEnd();
    glPointSize(1);
#endif /* SOFA_NO_OPENGL */
}

} // namespace constraintset

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_CONSTRAINTSET_BILATERALINTERACTIONCONSTRAINT_INL
