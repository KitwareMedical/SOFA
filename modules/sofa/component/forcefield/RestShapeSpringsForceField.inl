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
#ifndef SOFA_COMPONENT_FORCEFIELD_RESTSHAPESPRINGFORCEFIELD_INL
#define SOFA_COMPONENT_FORCEFIELD_RESTSHAPESPRINGFORCEFIELD_INL

#include <sofa/core/behavior/ForceField.inl>
#include "RestShapeSpringsForceField.h"
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/system/config.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/helper/gl/template.h>
#include <assert.h>
#include <iostream>


namespace sofa
{

namespace component
{

namespace forcefield
{

template<class DataTypes>
RestShapeSpringsForceField<DataTypes>::RestShapeSpringsForceField()
    : points(initData(&points, "points", "points controlled by the rest shape springs"))
    , stiffness(initData(&stiffness, "stiffness", "stiffness values between the actual position and the rest shape position"))
    , angularStiffness(initData(&angularStiffness, "angularStiffness", "angularStiffness assigned when controlling the rotation of the points"))
    , pivotPoints(initData(&pivotPoints, "pivot_points", "global pivot points used when translations instead of the rigid mass centers"))
    , external_rest_shape(initData(&external_rest_shape, "external_rest_shape", "rest_shape can be defined by the position of an external Mechanical State"))
    , external_points(initData(&external_points, "external_points", "points from the external Mechancial State that define the rest shape springs"))
    , recompute_indices(initData(&recompute_indices, false, "recompute_indices", "Recompute indices (should be false for BBOX)"))
    , drawSpring(initData(&drawSpring,false,"drawSpring","draw Spring"))
    , springColor(initData(&springColor,"springColor","spring color"))
    , restMState(NULL)
//	, pp_0(NULL)
{    
}


template<class DataTypes>
void RestShapeSpringsForceField<DataTypes>::bwdInit()
{
    core::behavior::ForceField<DataTypes>::init();

    if (stiffness.getValue().empty())
    {
        std::cout << "RestShapeSpringsForceField : No stiffness is defined, assuming equal stiffness on each node, k = 100.0 " << std::endl;

        VecReal stiffs;
        stiffs.push_back(100.0);
        stiffness.setValue(stiffs);
    }

    const std::string path = external_rest_shape.getValue();

    restMState = NULL;

    if (path.size() > 0)
    {
        this->getContext()->get(restMState ,path);
    }

    if (!restMState)
    {
        useRestMState = false;

        if (path.size() > 0)
        {
            serr << "RestShapeSpringsForceField : " << external_rest_shape.getValue() << " not found" << sendl;
        }
    }
    else
    {
        useRestMState = true;

        // sout << "RestShapeSpringsForceField : Mechanical state named " << restMState->getName() << " found for RestShapeSpringFF named " << this->getName() << sendl;
    }

    this->k = stiffness.getValue();

    recomputeIndices();

#ifdef SOFA_HAVE_EIGEN2
    core::behavior::BaseMechanicalState* state = this->getContext()->getMechanicalState();
    assert(state);
    matS.resize(state->getMatrixSize(),state->getMatrixSize());
    lastUpdatedStep = -1.0;
#endif
}


template<class DataTypes>
void RestShapeSpringsForceField<DataTypes>::recomputeIndices()
{
    m_indices.clear();
    m_ext_indices.clear();

    for (unsigned int i = 0; i < points.getValue().size(); i++)
        m_indices.push_back(points.getValue()[i]);

    for (unsigned int i = 0; i < external_points.getValue().size(); i++)
        m_ext_indices.push_back(external_points.getValue()[i]);

    m_pivots = pivotPoints.getValue();

    if (m_indices.size()==0)
    {
        //	std::cout << "in RestShapeSpringsForceField no point are defined, default case: points = all points " << std::endl;

        for (unsigned int i = 0; i < (unsigned)this->mstate->getSize(); i++)
        {
            m_indices.push_back(i);
        }
    }

    if (m_ext_indices.size()==0)
    {
        // std::cout << "in RestShapeSpringsForceField no external_points are defined, default case: points = all points " << std::endl;

        if (useRestMState)
        {
            for (unsigned int i = 0; i < (unsigned)restMState->getSize(); i++)
            {
                m_ext_indices.push_back(i);
            }
        }
        else
        {
            for (unsigned int i = 0; i < (unsigned)this->mstate->getSize(); i++)
            {
                m_ext_indices.push_back(i);
            }
        }
    }

    if (m_indices.size() > m_ext_indices.size())
    {
        std::cerr << "Error : the dimention of the source and the targeted points are different " << std::endl;
        m_indices.clear();
    }
}


template<class DataTypes>
const typename RestShapeSpringsForceField<DataTypes>::DataVecCoord* RestShapeSpringsForceField<DataTypes>::getExtPosition() const
{
    return (useRestMState ? restMState->read(core::VecCoordId::position()) : this->mstate->read(core::VecCoordId::restPosition()));
}

template<class DataTypes>
void RestShapeSpringsForceField<DataTypes>::addForce(const core::MechanicalParams* /* mparams */ /* PARAMS FIRST */, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& /* v */)
{
    sofa::helper::WriteAccessor< DataVecDeriv > f1 = f;
    sofa::helper::ReadAccessor< DataVecCoord > p1 = x;
    sofa::helper::ReadAccessor< DataVecCoord > p0 = *getExtPosition();

    f1.resize(p1.size());

    if (recompute_indices.getValue())
    {
        recomputeIndices();
    }

    //Springs_dir.resize(m_indices.size() );
    if ( k.size()!= m_indices.size() )
    {
        //sout << "WARNING : stiffness is not defined on each point, first stiffness is used" << sendl;
        const Real k0 = k[0];

        for (unsigned int i=0; i<m_indices.size(); i++)
        {
            const unsigned int index = m_indices[i];

            unsigned int ext_index = m_indices[i];
            if(useRestMState)
                ext_index= m_ext_indices[i];

            Deriv dx = p1[index] - p0[ext_index];
            //Springs_dir[i] = p1[index] - p0[ext_index];
            //Springs_dir[i].normalize();
            f1[index] -=  dx * k0 ;

            //	if (dx.norm()>0.00000001)
            //		std::cout<<"force on point "<<index<<std::endl;

            //	Deriv dx = p[i] - p_0[i];
            //	f[ indices[i] ] -=  dx * k[0] ;
        }
    }
    else
    {
        for (unsigned int i=0; i<m_indices.size(); i++)
        {
            const unsigned int index = m_indices[i];
            unsigned int ext_index = m_indices[i];
            if(useRestMState)
                ext_index= m_ext_indices[i];

            Deriv dx = p1[index] - p0[ext_index];
            //Springs_dir[i] = p1[index] - p0[ext_index];
            //Springs_dir[i].normalize();
            f1[index] -=  dx * k[i];

            //	if (dx.norm()>0.00000001)
            //		std::cout<<"force on point "<<index<<std::endl;

            //	Deriv dx = p[i] - p_0[i];
            //	f[ indices[i] ] -=  dx * k[i] ;
        }
    }
}


template<class DataTypes>
void RestShapeSpringsForceField<DataTypes>::addDForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& df, const DataVecDeriv& dx)
{
    //  remove to be able to build in parallel
    // 	const VecIndex& indices = points.getValue();
    // 	const VecReal& k = stiffness.getValue();

    sofa::helper::WriteAccessor< DataVecDeriv > df1 = df;
    sofa::helper::ReadAccessor< DataVecDeriv > dx1 = dx;
    double kFactor = mparams->kFactor();

    if (k.size()!= m_indices.size() )
    {
        //sout << "WARNING : stiffness is not defined on each point, first stiffness is used" << sendl;
        const Real k0 = k[0];

        for (unsigned int i=0; i<m_indices.size(); i++)
        {
            df1[m_indices[i]] -=  dx1[m_indices[i]] * k0 * kFactor;
        }
    }
    else
    {
        for (unsigned int i=0; i<m_indices.size(); i++)
        {
            df1[m_indices[i]] -=  dx1[m_indices[i]] * k[i] * kFactor ;
        }
    }
}

template<class DataTypes>
void RestShapeSpringsForceField<DataTypes>::draw(const core::visual::VisualParams * /* vparams */ )
{

}

template<class DataTypes>
void RestShapeSpringsForceField<DataTypes>::addKToMatrix(const core::MechanicalParams* mparams /* PARAMS FIRST */, const sofa::core::behavior::MultiMatrixAccessor* matrix )
{
    //      remove to be able to build in parallel
    // 	const VecIndex& indices = points.getValue();
    // 	const VecReal& k = stiffness.getValue();   
    sofa::core::behavior::MultiMatrixAccessor::MatrixRef mref = matrix->getMatrix(this->mstate);
    sofa::defaulttype::BaseMatrix* mat = mref.matrix;
    unsigned int offset = mref.offset;
    double kFact = mparams->kFactor();   

    const int N = Coord::total_size;

    unsigned int curIndex = 0;

    if (k.size()!= m_indices.size() )
    {
        const Real k0 = k[0];
        for (unsigned int index = 0; index < m_indices.size(); index++)
        {
            curIndex = m_indices[index];

            for(int i = 0; i < N; i++)
            {

                //	for (unsigned int j = 0; j < N; j++)
                //	{
                //		mat->add(offset + N * curIndex + i, offset + N * curIndex + j, kFact * k[0]);
                //	}

                mat->add(offset + N * curIndex + i, offset + N * curIndex + i, -kFact * k0);                
            }
        }
    }
    else
    {
        for (unsigned int index = 0; index < m_indices.size(); index++)
        {
            curIndex = m_indices[index];

            for(int i = 0; i < N; i++)
            {

                //	for (unsigned int j = 0; j < N; j++)
                //	{
                //		mat->add(offset + N * curIndex + i, offset + N * curIndex + j, kFact * k[curIndex]);
                //	}

                mat->add(offset + N * curIndex + i, offset + N * curIndex + i, -kFact * k[index]);
            }
        }
    }
}

template<class DataTypes>
void RestShapeSpringsForceField<DataTypes>::addSubKToMatrix(const core::MechanicalParams* mparams /* PARAMS FIRST */, const sofa::core::behavior::MultiMatrixAccessor* matrix, const helper::vector<unsigned> & addSubIndex )
{
    //      remove to be able to build in parallel
    // 	const VecIndex& indices = points.getValue();
    // 	const VecReal& k = stiffness.getValue();
    sofa::core::behavior::MultiMatrixAccessor::MatrixRef mref = matrix->getMatrix(this->mstate);
    sofa::defaulttype::BaseMatrix* mat = mref.matrix;
    unsigned int offset = mref.offset;
    double kFact = mparams->kFactor();

    const int N = Coord::total_size;

    unsigned int curIndex = 0;

    if (k.size()!= m_indices.size() )
    {
        const Real k0 = k[0];
        for (unsigned int index = 0; index < m_indices.size(); index++)
        {
            curIndex = m_indices[index];

            bool contains=false;
            for (unsigned s=0;s<addSubIndex.size() && !contains;s++) if (curIndex==addSubIndex[s]) contains=true;
            if (!contains) continue;

            for(int i = 0; i < N; i++)
            {

                //	for (unsigned int j = 0; j < N; j++)
                //	{
                //		mat->add(offset + N * curIndex + i, offset + N * curIndex + j, kFact * k[0]);
                //	}

                mat->add(offset + N * curIndex + i, offset + N * curIndex + i, -kFact * k0);
            }
        }
    }
    else
    {
        for (unsigned int index = 0; index < m_indices.size(); index++)
        {
            curIndex = m_indices[index];

            bool contains=false;
            for (unsigned s=0;s<addSubIndex.size() && !contains;s++) if (curIndex==addSubIndex[s]) contains=true;
            if (!contains) continue;

            for(int i = 0; i < N; i++)
            {

                //	for (unsigned int j = 0; j < N; j++)
                //	{
                //		mat->add(offset + N * curIndex + i, offset + N * curIndex + j, kFact * k[curIndex]);
                //	}

                mat->add(offset + N * curIndex + i, offset + N * curIndex + i, -kFact * k[index]);
            }
        }
    }
}

#ifdef SOFA_HAVE_EIGEN2
template<class DataTypes>
const sofa::defaulttype::BaseMatrix* RestShapeSpringsForceField<DataTypes>::getStiffnessMatrix(const core::MechanicalParams* mparams) {
    sout << this->getName() << ": getStiffnessMatrix " << sendl;

    double actualTime = this->getContext()->getTime();

    sout << "Actual Time: " << actualTime << " / " << lastUpdatedStep << sendl;

    if (actualTime == lastUpdatedStep) {  //matrix already assembled in this time step
        sout << this->getName() << "Matrix already assembled in step " << actualTime << sendl;
        return &matS;
    }

    lastUpdatedStep = actualTime;
    sout << this->getName() << "Assemble matrix in step " << actualTime << sendl;

    unsigned int offset = 0;
    double kFact = mparams->kFactor();
    const int N = Coord::total_size;

    unsigned int curIndex = 0;

    typedef Eigen::Triplet<double> Triplet;
    std::vector<Triplet> tripletList;
    matS.clear();

    if (k.size()!= m_indices.size() )
    {
        const Real k0 = k[0];
        for (unsigned int index = 0; index < m_indices.size(); index++)
        {
            curIndex = m_indices[index];    
            for(int i = 0; i < N; i++)
                tripletList.push_back(Triplet(offset + N * curIndex + i, offset + N * curIndex + i, -kFact * k0));
        }        
    }
    else
    {        
        for (unsigned int index = 0; index < m_indices.size(); index++)
        {
            curIndex = m_indices[index];

            for(int i = 0; i < N; i++)
                tripletList.push_back(Triplet(offset + N * curIndex + i, offset + N * curIndex + i, -kFact * k[index]));
        }
    }        

    matS.compressedMatrix.setFromTriplets(tripletList.begin(), tripletList.end());    

    return &matS;
}
#endif

} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FORCEFIELD_RESTSHAPESPRINGFORCEFIELD_INL



