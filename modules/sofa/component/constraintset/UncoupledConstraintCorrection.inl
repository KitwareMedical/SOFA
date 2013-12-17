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
#ifndef SOFA_COMPONENT_CONSTRAINTSET_UNCOUPLEDCONSTRAINTCORRECTION_INL
#define SOFA_COMPONENT_CONSTRAINTSET_UNCOUPLEDCONSTRAINTCORRECTION_INL
//#define DEBUG
//#define NEW_VERSION

#include "UncoupledConstraintCorrection.h"

#include <sofa/core/behavior/ConstraintCorrection.inl>

#include <sofa/simulation/common/MechanicalVisitor.h>

#include <sofa/core/topology/BaseMeshTopology.h>

#include <sofa/core/topology/TopologyChange.h>
#include <sofa/component/topology/PointSetTopologyContainer.h>


namespace sofa
{

namespace component
{

namespace constraintset
{
using namespace sofa::core::topology;

template<class DataTypes>
UncoupledConstraintCorrection<DataTypes>::UncoupledConstraintCorrection(behavior::MechanicalState<DataTypes> *mm)
    : Inherit(mm)
    , compliance(initData(&compliance, "compliance", "compliance value on each dof"))
    , defaultCompliance(initData(&defaultCompliance, (Real)0.00001, "defaultCompliance", "Default compliance value for new dof or if all should have the same (in which case compliance vector should be empty)"))
    , f_verbose( initData(&f_verbose,false,"verbose","Dump the constraint matrix at each iteration") )
{
}

template<class DataTypes>
UncoupledConstraintCorrection<DataTypes>::~UncoupledConstraintCorrection()
{
}


template<class DataTypes>
void UncoupledConstraintCorrection<DataTypes>::init()
{
    Inherit::init();

    const VecCoord& x = *this->mstate->getX();

    if (compliance.getValue().size() == 1 && defaultCompliance.isSet() && defaultCompliance.getValue() == compliance.getValue()[0])
    {
        // the same compliance was set in the two data, only keep the default one so that the compliance vector does not need to be maintained
        compliance.setValue(VecReal());
    }

    const VecReal& comp = compliance.getValue();
    if (x.size() != comp.size() && !comp.empty())
    {
        if (!defaultCompliance.isSet() && !comp.empty())
            defaultCompliance.setValue(comp[0]); // set default compliance to first one in case it was given in the compliance vector
        if (comp.size() > 1)
            serr << "Warning compliance size ( " << comp.size() << " is not equal to the size of the mstate (" << x.size() << ")" << sendl;
        Real comp0 = (!comp.empty()) ? comp[0] : defaultCompliance.getValue();
        sout << "Using " << comp0 << " as initial compliance" << sendl;
        VecReal UsedComp;
        for (unsigned int i=0; i<x.size(); i++)
        {
            UsedComp.push_back(comp0);
        }
        // Keeps user specified compliance even if the initial MState size is null.
        if (!UsedComp.empty())
        {
            compliance.setValue(UsedComp);
        }
    }
}

template <class DataTypes>
void UncoupledConstraintCorrection<DataTypes>::reinit()
{
    Inherit::reinit();
}

template< class DataTypes >
void UncoupledConstraintCorrection< DataTypes >::handleTopologyChange()
{
    using sofa::core::topology::TopologyChange;
    using sofa::core::topology::TopologyChangeType;
    using sofa::core::topology::BaseMeshTopology;

    BaseMeshTopology *topology = this->getContext()->getMeshTopology();
    if (!topology)
        return;
    if (defaultCompliance.isSet() && compliance.getValue().empty())
        return; // uniform compliance, no need to update it

    std::list< const TopologyChange * >::const_iterator itBegin = topology->beginChange();
    std::list< const TopologyChange * >::const_iterator itEnd = topology->endChange();

    VecReal& comp = *(compliance.beginEdit());
    const Real comp0 = defaultCompliance.getValue();

    for (std::list< const TopologyChange * >::const_iterator changeIt = itBegin; changeIt != itEnd; ++changeIt)
    {
        const TopologyChangeType changeType = (*changeIt)->getChangeType();

        switch ( changeType )
        {
        case core::topology::POINTSADDED :
        {
            unsigned int nbPoints = (static_cast< const PointsAdded *> (*changeIt))->getNbAddedVertices();

            VecReal addedCompliance;

            for (unsigned int i = 0; i < nbPoints; i++)
            {
                //	std::cout << "addedCompliance --> " << comp0 << std::endl;
                addedCompliance.push_back(comp0);
            }

            comp.insert(comp.end(), addedCompliance.begin(), addedCompliance.end());

            break;
        }
        case core::topology::POINTSREMOVED :
        {
            using sofa::helper::vector;

            const vector< unsigned int > &pts = (static_cast< const PointsRemoved * >(*changeIt))->getArray();

            unsigned int lastIndexVec = comp.size() - 1;

            for (unsigned int i = 0; i < pts.size(); i++)
            {
                comp[pts[i]] = comp[lastIndexVec];
                lastIndexVec--;
            }

            comp.resize(comp.size() - pts.size());

            break;
        }
        default:
            break;
        }
    }

    compliance.endEdit();
}


template<class DataTypes>
void UncoupledConstraintCorrection<DataTypes>::getComplianceWithConstraintMerge(defaulttype::BaseMatrix* Wmerged, std::vector<int> &constraint_merge)
{
    helper::WriteAccessor<Data<MatrixDeriv> > constraintsData = *this->mstate->write(core::MatrixDerivId::holonomicC());
    MatrixDeriv& constraints = constraintsData.wref();

    MatrixDeriv constraintCopy;

    std::cout << "******\n Constraint before Merge  \n *******" << std::endl;

    MatrixDerivRowIterator rowIt = constraints.begin();
    MatrixDerivRowIterator rowItEnd = constraints.end();

    while (rowIt != rowItEnd)
    {
        constraintCopy.writeLine(rowIt.index(), rowIt.row());
        ++rowIt;
    }

    /////////// MERGE OF THE CONSTRAINTS //////////////
    constraints.clear();

    // look for the number of group;
    unsigned int numGroup = 0;
    for (unsigned int cm = 0; cm < constraint_merge.size(); cm++)
    {
        if (constraint_merge[cm] > (int) numGroup)
            numGroup = (unsigned int) constraint_merge[cm];
    }
    numGroup += 1;

    std::cout << "******\n Constraint after Merge  \n *******" << std::endl;

    for (unsigned int group = 0; group < numGroup; group++)
    {
        std::cout << "constraint[" << group << "] : " << std::endl;

        MatrixDerivRowIterator rowCopyIt = constraintCopy.begin();
        MatrixDerivRowIterator rowCopyItEnd = constraintCopy.end();

        while (rowCopyIt != rowCopyItEnd)
        {
            if (constraint_merge[rowCopyIt.index()] == (int)group)
            {
                constraints.addLine(group, rowCopyIt.row());
            }

            ++rowCopyIt;
        }
    }

    //////////// compliance computation call //////////
    this->addComplianceInConstraintSpace(ConstraintParams::defaultInstance(), Wmerged);

    /////////// BACK TO THE INITIAL CONSTRAINT SET//////////////

    constraints.clear();
    std::cout << "******\n Constraint back to initial values  \n *******" << std::endl;

    rowIt = constraintCopy.begin();
    rowItEnd = constraintCopy.end();

    while (rowIt != rowItEnd)
    {
        constraints.writeLine(rowIt.index(), rowIt.row());
        ++rowIt;
    }
}


#ifndef NEW_VERSION
template<class DataTypes>
void UncoupledConstraintCorrection<DataTypes>::addComplianceInConstraintSpace(const ConstraintParams * /*cparams*/, defaulttype::BaseMatrix *W)
{
    const MatrixDeriv& constraints = *this->mstate->getC();
    const VecReal& comp = compliance.getValue();
    const Real comp0 = defaultCompliance.getValue();

    for (MatrixDerivRowConstIterator rowIt = constraints.begin(), rowItEnd = constraints.end(); rowIt != rowItEnd; ++rowIt)
    {
        int indexCurRowConst = rowIt.index();

        if (f_verbose.getValue())
            sout << "C[" << indexCurRowConst << "]";

        for (MatrixDerivColConstIterator colIt = rowIt.begin(), colItEnd = rowIt.end(); colIt != colItEnd; ++colIt)
        {
            unsigned int dof = colIt.index();
            Deriv n = colIt.val();

            if (f_verbose.getValue())
                sout << " dof[" << dof << "]=" << n;

            int indexCurColConst;

#ifdef DEBUG
            std::cout << " [ " << dof << "]=" << n << std::endl;
#endif
            for (MatrixDerivRowConstIterator rowIt2 = rowIt; rowIt2 != rowItEnd; ++rowIt2)
            {
                indexCurColConst = rowIt2.index();

                for (MatrixDerivColConstIterator colIt2 = rowIt2.begin(), colIt2End = rowIt2.end(); colIt2 != colIt2End; ++colIt2)
                {
                    if (dof == colIt2.index())
                    {
                        double w = n * colIt2.val() * (dof < comp.size() ? comp[dof] : comp0);
                        W->add(indexCurRowConst, indexCurColConst, w);
                        if (indexCurRowConst != indexCurColConst)
                        {
                            W->add(indexCurColConst, indexCurRowConst, w);
                        }
                    }
                }
            }

            /*
            for(unsigned int curColConst = curRowConst+1; curColConst < numConstraints; curColConst++)
            {
            	indexCurColConst = mstate->getConstraintId()[curColConst];
            	W[indexCurColConst][indexCurRowConst] = W[indexCurRowConst][indexCurColConst];
            }
            */
        }
        if (f_verbose.getValue())
            sout << sendl;
    }

    /*debug : verifie qu'il n'y a pas de 0 sur la diagonale de W
    printf("\n index : ");
    for(unsigned int curRowConst = 0; curRowConst < numConstraints; curRowConst++)
    {
    	int indexCurRowConst = mstate->getConstraintId()[curRowConst];
    	printf(" %d ",indexCurRowConst);
    	if(abs(W[indexCurRowConst][indexCurRowConst]) < 0.000000001)
    		printf("\n WARNING : there is a 0 on the diagonal of matrix W");

    	if(abs(W[curRowConst][curRowConst]) <0.000000001)
    		printf("\n stop");
    }*/
}

#else

template<class DataTypes>
void UncoupledConstraintCorrection<DataTypes>::addComplianceInConstraintSpace(const ConstraintParams * /*cparams*/, defaulttype::BaseMatrix *W)
{
    const MatrixDeriv& constraints = *this->mstate->getC();
    const VecReal& comp = compliance.getValue();
    const Real comp0 = defaultCompliance.getValue();

    typedef std::list< std::pair< int, Deriv > > CIndicesAndValues;

    helper::vector< CIndicesAndValues > dofsIndexedConstraints;
    const unsigned int numDOFs = this->mstate->getSize();
    dofsIndexedConstraints.resize(numDOFs);

    for (MatrixDerivRowConstIterator rowIt = constraints.begin(), rowItEnd = constraints.end(); rowIt != rowItEnd; ++rowIt)
    {
        int indexCurRowConst = rowIt.index();

        for (MatrixDerivColConstIterator colIt = rowIt.begin(), colItEnd = rowIt.end(); colIt != colItEnd; ++colIt)
        {
            dofsIndexedConstraints[colIt.index()].push_back(std::make_pair(indexCurRowConst, colIt.val()));
        }
    }

    for (MatrixDerivRowConstIterator rowIt = constraints.begin(), rowItEnd = constraints.end(); rowIt != rowItEnd; ++rowIt)
    {
        int indexCurRowConst = rowIt.index();

        for (MatrixDerivColConstIterator colIt = rowIt.begin(), colItEnd = rowIt.end(); colIt != colItEnd; ++colIt)
        {
            unsigned int dof = colIt.index();

            CIndicesAndValues &dofsConstraint = dofsIndexedConstraints[dof];

            if (!dofsConstraint.empty())
            {
                Deriv n = colIt.val();
                const Real dofComp = dof < comp.size() ? comp[dof] : comp0;

                double w = n * dofsConstraint.front().second * dofComp;
                W->add(indexCurRowConst, indexCurRowConst, w);
                dofsConstraint.pop_front();

                for (CIndicesAndValues::const_iterator it = dofsConstraint.begin(), itEnd = dofsConstraint.end(); it != itEnd; ++it)
                {
                    w = n * it->second * dofComp;
                    W->add(indexCurRowConst, it->first, w);
                    W->add(it->first, indexCurRowConst, w);
                }
            }
        }
    }
}

#endif

template<class DataTypes>
void UncoupledConstraintCorrection<DataTypes>::getComplianceMatrix(defaulttype::BaseMatrix *m) const
{
    const VecReal& comp = compliance.getValue();
    const Real comp0 = defaultCompliance.getValue();
    const unsigned int s = this->mstate->getSize(); // comp.size();
    const unsigned int dimension = Coord::size();

    m->resize(s * dimension, s * dimension); //resize must set to zero the content of the matrix

    for (unsigned int l = 0; l < s; ++l)
    {
        for (unsigned int d = 0; d < dimension; ++d)
            m->set(dimension * l + d, dimension * l + d, (l < comp.size() ? comp[l] : comp0));
    }
}


template<class DataTypes>
void UncoupledConstraintCorrection<DataTypes>::computeDx(const Data< VecDeriv > &f_d)
{
    const VecDeriv& f = f_d.getValue();

    Data< VecDeriv > &dx_d = *this->mstate->write(core::VecDerivId::dx());
    VecDeriv& dx = *dx_d.beginEdit();

    dx.resize(f.size());
    const VecReal& comp = compliance.getValue();
    const Real comp0 = defaultCompliance.getValue();

    for (unsigned int i = 0; i < dx.size(); i++)
    {
        dx[i] = f[i] * (i < comp.size() ? comp[i] : comp0);
    }

    dx_d.endEdit();
}


template<class DataTypes>
void UncoupledConstraintCorrection<DataTypes>::computeAndApplyMotionCorrection(const core::ConstraintParams *cparams, Data< VecCoord > &x_d, Data< VecDeriv > &v_d, Data< VecDeriv > &f_d, const defaulttype::BaseVector *lambda)
{
    this->addConstraintForceInMotionSpace(f_d, lambda);

    computeDx(f_d);

    VecCoord& x = *x_d.beginEdit();
    VecDeriv& v = *v_d.beginEdit();

    const VecDeriv& dx = this->mstate->read(core::VecDerivId::dx())->getValue();

    const VecCoord& x_free = cparams->readX(this->mstate)->getValue();
    const VecDeriv& v_free = cparams->readV(this->mstate)->getValue();

    const double invDt = 1.0 / this->getContext()->getDt();

    for (unsigned int i = 0; i < dx.size(); i++)
    {
        x[i] = x_free[i];
        v[i] = v_free[i];
        x[i] += dx[i];
        v[i] += dx[i] * invDt;
    }

    x_d.endEdit();
    v_d.endEdit();
}


template<class DataTypes>
void UncoupledConstraintCorrection<DataTypes>::computeAndApplyPositionCorrection(const core::ConstraintParams *cparams, Data< VecCoord > &x_d, Data< VecDeriv > &f_d, const defaulttype::BaseVector *lambda)
{
    this->addConstraintForceInMotionSpace(f_d, lambda);

    computeDx(f_d);

    VecCoord& x = *x_d.beginEdit();

    const VecCoord& x_free = cparams->readX(this->mstate)->getValue();

    const VecDeriv& dx = this->mstate->read(core::VecDerivId::dx())->getValue();

    for (unsigned int i = 0; i < dx.size(); i++)
    {
        x[i] = x_free[i] + dx[i];
    }

    x_d.endEdit();
}


template<class DataTypes>
void UncoupledConstraintCorrection<DataTypes>::computeAndApplyVelocityCorrection(const core::ConstraintParams *cparams, Data< VecDeriv > &v_d, Data< VecDeriv > &f_d, const defaulttype::BaseVector *lambda)
{
    this->addConstraintForceInMotionSpace(f_d, lambda);

    computeDx(f_d);

    VecDeriv& v = *v_d.beginEdit();

    const VecDeriv& v_free = cparams->readV(this->mstate)->getValue();

    const VecDeriv& dx = this->mstate->read(core::VecDerivId::dx())->getValue();
//	const double invDt = 1.0 / this->getContext()->getDt();

    for (unsigned int i = 0; i < dx.size(); i++)
    {
        v[i] = v_free[i] + dx[i]/* * invDt*/;
    }

    v_d.endEdit();
}


template<class DataTypes>
void UncoupledConstraintCorrection<DataTypes>::applyContactForce(const defaulttype::BaseVector *f)
{
    helper::WriteAccessor<Data<VecDeriv> > forceData = *this->mstate->write(core::VecDerivId::externalForce());
    VecDeriv& force = forceData.wref();
    const MatrixDeriv& constraints = *this->mstate->getC();
    const VecReal& comp = compliance.getValue();
    const Real comp0 = defaultCompliance.getValue();

    force.resize((*this->mstate->getX()).size());

    MatrixDerivRowConstIterator rowItEnd = constraints.end();

    for (MatrixDerivRowConstIterator rowIt = constraints.begin(); rowIt != rowItEnd; ++rowIt)
    {
        double fC1 = f->element(rowIt.index());

        if (fC1 != 0.0)
        {
            MatrixDerivColConstIterator colItEnd = rowIt.end();

            for (MatrixDerivColConstIterator colIt = rowIt.begin(); colIt != colItEnd; ++colIt)
            {
                force[colIt.index()] += colIt.val() * fC1;
            }
        }
    }


    helper::WriteAccessor<Data<VecDeriv> > dxData = *this->mstate->write(core::VecDerivId::dx());
    VecDeriv& dx = dxData.wref();
    helper::WriteAccessor<Data<VecCoord> > xData = *this->mstate->write(core::VecCoordId::position());
    VecCoord& x = xData.wref();
    helper::WriteAccessor<Data<VecDeriv> > vData = *this->mstate->write(core::VecDerivId::velocity());
    VecDeriv& v = vData.wref();
    const VecDeriv& v_free = *this->mstate->getVfree();
    const VecCoord& x_free = *this->mstate->getXfree();
    const double invDt = 1.0/this->getContext()->getDt();

    // Euler integration... will be done in the "integrator" as soon as it exists !
    dx.resize(v.size());

    for (unsigned int i = 0; i < dx.size(); i++)
    {
        x[i] = x_free[i];
        v[i] = v_free[i];
        dx[i] = force[i] * (i<comp.size() ? comp[i] : comp0);
        x[i] += dx[i];
        v[i] += dx[i]*invDt;
    }
}


template<class DataTypes>
void UncoupledConstraintCorrection<DataTypes>::applyPredictiveConstraintForce(const core::ConstraintParams * /*cparams*/, Data< VecDeriv > &f_d, const defaulttype::BaseVector *lambda)
{
    this->setConstraintForceInMotionSpace(f_d, lambda);
}



template<class DataTypes>
void UncoupledConstraintCorrection<DataTypes>::resetContactForce()
{
    helper::WriteAccessor<Data<VecDeriv> > forceData = *this->mstate->write(core::VecDerivId::externalForce());
    VecDeriv& force = forceData.wref();

    for (unsigned i = 0; i < force.size(); ++i)
    {
        force[i] = Deriv();
    }
}


///////////////////////  new API for non building the constraint system during solving process //
template<class DataTypes>
bool UncoupledConstraintCorrection<DataTypes>::hasConstraintNumber(int index)
{
    const MatrixDeriv &constraints = *this->mstate->getC();

    return (constraints.readLine(index) != constraints.end());
}


template<class DataTypes>
void UncoupledConstraintCorrection<DataTypes>::resetForUnbuiltResolution(double * f, std::list<unsigned int>& /*renumbering*/)
{
    const MatrixDeriv& constraints = *this->mstate->getC();

    constraint_disp.clear();
    constraint_disp.resize(this->mstate->getSize());

    constraint_force.clear();
    constraint_force.resize(this->mstate->getSize());

    constraint_dofs.clear();



    for (MatrixDerivRowConstIterator rowIt = constraints.begin(); rowIt != constraints.end(); ++rowIt)
    {
        int indexC = rowIt.index();

        // buf the value of force applied on concerned dof : constraint_force
        // buf a table of indice of involved dof : constraint_dofs
        double fC = f[indexC];

        if (fC != 0.0)
        {
            for (MatrixDerivColConstIterator colIt = rowIt.begin(); colIt != rowIt.end(); ++colIt)
            {
                unsigned int dof = colIt.index();
                constraint_force[dof] += colIt.val() * fC;
                constraint_dofs.push_back(dof);
            }
        }
    }

    // constraint_dofs buff the DOF that are involved with the constraints
    constraint_dofs.unique();
}


template<class DataTypes>
void UncoupledConstraintCorrection<DataTypes>::addConstraintDisplacement(double * d, int begin, int end)
{
/// in the Vec1Types and Vec3Types case, compliance is a vector of size mstate->getSize()
/// constraint_force contains the force applied on dof involved with the contact
/// TODO : compute a constraint_disp that is updated each time a new force is provided !

    const MatrixDeriv& constraints = *this->mstate->getC();

    for (int id = begin; id <= end; id++)
    {
        MatrixDerivRowConstIterator curConstraint = constraints.readLine(id);

        if (curConstraint != constraints.end())
        {
            MatrixDerivColConstIterator colIt = curConstraint.begin();
            MatrixDerivColConstIterator colItEnd = curConstraint.end();

            while (colIt != colItEnd)
            {
                d[id] += colIt.val() * constraint_disp[colIt.index()];

                ++colIt;
            }
        }
    }
}


template<class DataTypes>
void UncoupledConstraintCorrection<DataTypes>::setConstraintDForce(double * df, int begin, int end, bool update)
{
    /// set a force difference on a set of constraints (between constraint number "begin" and constraint number "end"
    /// if update is false, do nothing
    /// if update is true, it computes the displacements due to this delta of force.
    /// As the contact are uncoupled, a displacement is obtained only on dof involved with the constraints

    const MatrixDeriv& constraints = *this->mstate->getC();
    const VecReal& comp = compliance.getValue();
    const Real comp0 = defaultCompliance.getValue();

    if (!update)
        return;

    for (int id = begin; id <= end; id++)
    {

        MatrixDerivRowConstIterator curConstraint = constraints.readLine(id);

        if (curConstraint != constraints.end())
        {
            MatrixDerivColConstIterator colIt = curConstraint.begin();
            MatrixDerivColConstIterator colItEnd = curConstraint.end();

            while (colIt != colItEnd)
            {
                const unsigned int dof = colIt.index();

                constraint_force[dof] += colIt.val() * df[id];

                Deriv dx =  constraint_force[dof] * (dof < comp.size() ? comp[dof] : comp0);

                constraint_disp[dof] = dx;

                ++colIt;
            }
        }
    }
}


template<class DataTypes>
void UncoupledConstraintCorrection<DataTypes>::getBlockDiagonalCompliance(defaulttype::BaseMatrix* W, int begin, int end)
{
    const MatrixDeriv& constraints = *this->mstate->getC();
    const VecReal& comp = compliance.getValue();
    const Real comp0 = defaultCompliance.getValue();

    for (int id1 = begin; id1 <= end; id1++)
    {

        MatrixDerivRowConstIterator curConstraint = constraints.readLine(id1);

        if (curConstraint != constraints.end())
        {
            MatrixDerivColConstIterator colIt = curConstraint.begin();
            MatrixDerivColConstIterator colItEnd = curConstraint.end();

            while (colIt != colItEnd)
            {
                Deriv n1 = colIt.val();
                unsigned int dof1 = colIt.index();

                for (int id2 = id1; id2 <= end; id2++)
                {
                    MatrixDerivRowConstIterator curConstraint2 = constraints.readLine(id2);

                    if (curConstraint2 != constraints.end())
                    {
                        MatrixDerivColConstIterator colIt2 = curConstraint2.begin();
                        MatrixDerivColConstIterator colItEnd2 = curConstraint2.end();

                        while (colIt2 != colItEnd2)
                        {
                            unsigned int dof2 = colIt2.index();

                            if (dof1 == dof2)
                            {
                                double w = n1 * colIt2.val() * (dof1 < comp.size() ? comp[dof1] : comp0);
                                W->add(id1, id2, w);
                                if (id1 != id2)
                                    W->add(id2, id1, w);
                            }

                            ++colIt2;
                        }
                    }
                }

                ++colIt;
            }
        }
    }
}

} // namespace constraintset

} // namespace component

} // namespace sofa

#endif
