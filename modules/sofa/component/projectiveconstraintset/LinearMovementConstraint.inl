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
#ifndef SOFA_COMPONENT_PROJECTIVECONSTRAINTSET_LINEARMOVEMENTCONSTRAINT_INL
#define SOFA_COMPONENT_PROJECTIVECONSTRAINTSET_LINEARMOVEMENTCONSTRAINT_INL

#include <sofa/component/projectiveconstraintset/LinearMovementConstraint.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/behavior/ProjectiveConstraintSet.inl>
#include <sofa/simulation/common/Simulation.h>
#include <sofa/helper/gl/template.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <iostream>
#include <sofa/component/topology/TopologySubsetData.inl>



namespace sofa
{

namespace component
{

namespace projectiveconstraintset
{

using namespace core::topology;

using namespace sofa::defaulttype;
using namespace sofa::helper;
using namespace sofa::core::behavior;

// Define TestNewPointFunction
template< class DataTypes>
bool LinearMovementConstraint<DataTypes>::FCPointHandler::applyTestCreateFunction(unsigned int, const sofa::helper::vector<unsigned int> &, const sofa::helper::vector<double> &)
{
    return lc != 0;
}

// Define RemovalFunction
template< class DataTypes>
void LinearMovementConstraint<DataTypes>::FCPointHandler::applyDestroyFunction(unsigned int pointIndex, value_type &)
{
    if (lc)
    {
        lc->removeIndex((unsigned int) pointIndex);
    }
}

template <class DataTypes>
LinearMovementConstraint<DataTypes>::LinearMovementConstraint()
    : core::behavior::ProjectiveConstraintSet<DataTypes>(NULL)
    , data(new LinearMovementConstraintInternalData<DataTypes>)
    , m_indices( initData(&m_indices,"indices","Indices of the constrained points") )
    , m_keyTimes(  initData(&m_keyTimes,"keyTimes","key times for the movements") )
    , m_keyMovements(  initData(&m_keyMovements,"movements","movements corresponding to the key times") )
    , showMovement( initData(&showMovement, (bool)false, "showMovement", "Visualization of the movement to be applied to constrained dofs."))
{
    // default to indice 0
    m_indices.beginEdit()->push_back(0);
    m_indices.endEdit();

    //default valueEvent to 0
    m_keyTimes.beginEdit()->push_back( 0.0 );
    m_keyTimes.endEdit();
    m_keyMovements.beginEdit()->push_back( Deriv() );
    m_keyMovements.endEdit();

    pointHandler = new FCPointHandler(this, &m_indices);
}



template <class DataTypes>
LinearMovementConstraint<DataTypes>::~LinearMovementConstraint()
{
    if (pointHandler)
        delete pointHandler;
}

template <class DataTypes>
void LinearMovementConstraint<DataTypes>::clearIndices()
{
    m_indices.beginEdit()->clear();
    m_indices.endEdit();
}

template <class DataTypes>
void LinearMovementConstraint<DataTypes>::addIndex(unsigned int index)
{
    m_indices.beginEdit()->push_back(index);
    m_indices.endEdit();
}

template <class DataTypes>
void LinearMovementConstraint<DataTypes>::removeIndex(unsigned int index)
{
    removeValue(*m_indices.beginEdit(),index);
    m_indices.endEdit();
}

template <class DataTypes>
void LinearMovementConstraint<DataTypes>::clearKeyMovements()
{
    m_keyTimes.beginEdit()->clear();
    m_keyTimes.endEdit();
    m_keyMovements.beginEdit()->clear();
    m_keyMovements.endEdit();
}

template <class DataTypes>
void LinearMovementConstraint<DataTypes>::addKeyMovement(Real time, Deriv movement)
{
    m_keyTimes.beginEdit()->push_back( time );
    m_keyTimes.endEdit();
    m_keyMovements.beginEdit()->push_back( movement );
    m_keyMovements.endEdit();
}

// -- Constraint interface


template <class DataTypes>
void LinearMovementConstraint<DataTypes>::init()
{
    this->core::behavior::ProjectiveConstraintSet<DataTypes>::init();

    topology = this->getContext()->getMeshTopology();

    // Initialize functions and parameters
    m_indices.createTopologicalEngine(topology, pointHandler);
    m_indices.registerTopologicalData();

    x0.resize(0);
    nextM = prevM = Deriv();

    currentTime = -1.0;
    finished = false;
}


template <class DataTypes>
void LinearMovementConstraint<DataTypes>::reset()
{
    nextT = prevT = 0.0;
    nextM = prevM = Deriv();

    currentTime = -1.0;
    finished = false;
}


template <class DataTypes>
template <class DataDeriv>
void LinearMovementConstraint<DataTypes>::projectResponseT(const core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataDeriv& dx)
{
    Real cT = (Real) this->getContext()->getTime();
    if ((cT != currentTime) || !finished)
    {
        findKeyTimes();
    }

    if (finished && nextT != prevT)
    {
        const SetIndexArray & indices = m_indices.getValue();

        //set the motion to the Dofs
        for (SetIndexArray::const_iterator it = indices.begin(); it != indices.end(); ++it)
        {
            dx[*it] = Deriv();
        }
    }
}

template <class DataTypes>
void LinearMovementConstraint<DataTypes>::projectResponse(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& resData)
{
    helper::WriteAccessor<DataVecDeriv> res = resData;
    projectResponseT<VecDeriv>(mparams /* PARAMS FIRST */, res.wref());
}

template <class DataTypes>
void LinearMovementConstraint<DataTypes>::projectVelocity(const core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecDeriv& vData)
{
    helper::WriteAccessor<DataVecDeriv> dx = vData;
    Real cT = (Real) this->getContext()->getTime();
    if ((cT != currentTime) || !finished)
    {
        findKeyTimes();
    }

    if (finished && nextT != prevT)
    {
        const SetIndexArray & indices = m_indices.getValue();

        //set the motion to the Dofs
        for (SetIndexArray::const_iterator it = indices.begin(); it != indices.end(); ++it)
        {
            dx[*it] = (nextM - prevM)*(1.0 / (nextT - prevT));
        }
    }
}


template <class DataTypes>
void LinearMovementConstraint<DataTypes>::projectPosition(const core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecCoord& xData)
{
    helper::WriteAccessor<DataVecCoord> x = xData;
    Real cT = (Real) this->getContext()->getTime();

    //initialize initial Dofs positions, if it's not done
    if (x0.size() == 0)
    {
        const SetIndexArray & indices = m_indices.getValue();
        x0.resize(x.size());
        for (SetIndexArray::const_iterator it = indices.begin(); it != indices.end(); ++it)
        {
            x0[*it] = x[*it];
        }
    }

    if ((cT != currentTime) || !finished)
    {
        findKeyTimes();
    }

    //if we found 2 keyTimes, we have to interpolate a velocity (linear interpolation)
    if(finished && nextT != prevT)
    {
        interpolatePosition<Coord>(cT, x.wref());
    }
}

template <class DataTypes>
template <class MyCoord>
void LinearMovementConstraint<DataTypes>::interpolatePosition(Real cT, typename boost::disable_if<boost::is_same<MyCoord, RigidCoord<3, Real> >, VecCoord>::type& x)
{
    const SetIndexArray & indices = m_indices.getValue();
//                cerr<<"LinearMovementConstraint<DataTypes>::interpolatePosition,  current time cT = "<<cT<<endl;
//                cerr<<"LinearMovementConstraint<DataTypes>::interpolatePosition,  prevT = "<<prevT<<" ,prevM= "<<prevM<<endl;
//                cerr<<"LinearMovementConstraint<DataTypes>::interpolatePosition,  nextT = "<<nextT<<" ,nextM= "<<nextM<<endl;
    //cerr<<"LinearMovementConstraint<DataTypes>::interpolatePosition, current x = "<<x<<endl;

    Real dt = (cT - prevT) / (nextT - prevT);
//                cerr<<"LinearMovementConstraint<DataTypes>::interpolatePosition, dt = "<<dt<<endl;
    Deriv m = prevM + (nextM-prevM)*dt;

    //cerr<<"LinearMovementConstraint<DataTypes>::interpolatePosition, movement m = "<<m<<endl;

    //set the motion to the Dofs
    for (SetIndexArray::const_iterator it = indices.begin(); it != indices.end(); ++it)
    {
        x[*it] = x0[*it] + m ;
    }

    //cerr<<"LinearMovementConstraint<DataTypes>::interpolatePosition, new x = "<<x<<endl<<endl<<endl;
}

template <class DataTypes>
template <class MyCoord>
void LinearMovementConstraint<DataTypes>::interpolatePosition(Real cT, typename boost::enable_if<boost::is_same<MyCoord, RigidCoord<3, Real> >, VecCoord>::type& x)
{
    const SetIndexArray & indices = m_indices.getValue();

    Real dt = (cT - prevT) / (nextT - prevT);
    Deriv m = prevM + (nextM-prevM)*dt;
    Quater<Real> prevOrientation = Quater<Real>::createQuaterFromEuler(getVOrientation(prevM));
    Quater<Real> nextOrientation = Quater<Real>::createQuaterFromEuler(getVOrientation(nextM));

    //set the motion to the Dofs
    for (SetIndexArray::const_iterator it = indices.begin(); it != indices.end(); ++it)
    {
        x[*it].getCenter() = x0[*it].getCenter() + getVCenter(m) ;
        x[*it].getOrientation() = x0[*it].getOrientation() * prevOrientation.slerp2(nextOrientation, dt);
    }
}

template <class DataTypes>
void LinearMovementConstraint<DataTypes>::projectJacobianMatrix(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataMatrixDeriv& cData)
{
    helper::WriteAccessor<DataMatrixDeriv> c = cData;

    MatrixDerivRowIterator rowIt = c->begin();
    MatrixDerivRowIterator rowItEnd = c->end();

    while (rowIt != rowItEnd)
    {
        projectResponseT<MatrixDerivRowType>(mparams /* PARAMS FIRST */, rowIt.row());
        ++rowIt;
    }
}

template <class DataTypes>
void LinearMovementConstraint<DataTypes>::findKeyTimes()
{
    Real cT = (Real) this->getContext()->getTime();
    finished = false;

    if(m_keyTimes.getValue().size() != 0 && cT >= *m_keyTimes.getValue().begin() && cT <= *m_keyTimes.getValue().rbegin())
    {
        nextT = *m_keyTimes.getValue().begin();
        prevT = nextT;

        typename helper::vector<Real>::const_iterator it_t = m_keyTimes.getValue().begin();
        typename VecDeriv::const_iterator it_m = m_keyMovements.getValue().begin();

        //WARNING : we consider that the key-events are in chronological order
        //here we search between which keyTimes we are, to know which are the motion to interpolate
        while( it_t != m_keyTimes.getValue().end() && !finished)
        {
            if( *it_t <= cT)
            {
                prevT = *it_t;
                prevM = *it_m;
            }
            else
            {
                nextT = *it_t;
                nextM = *it_m;
                finished = true;
            }
            it_t++;
            it_m++;
        }
    }
}

// Matrix Integration interface
template <class DataTypes>
void LinearMovementConstraint<DataTypes>::applyConstraint(defaulttype::BaseMatrix *mat, unsigned int offset)
{
    //sout << "applyConstraint in Matrix with offset = " << offset << sendl;
    const unsigned int N = Deriv::size();
    const SetIndexArray & indices = m_indices.getValue();

    for (SetIndexArray::const_iterator it = indices.begin(); it != indices.end(); ++it)
    {
        // Reset Fixed Row and Col
        for (unsigned int c=0; c<N; ++c)
            mat->clearRowCol(offset + N * (*it) + c);
        // Set Fixed Vertex
        for (unsigned int c=0; c<N; ++c)
            mat->set(offset + N * (*it) + c, offset + N * (*it) + c, 1.0);
    }
}

template <class DataTypes>
void LinearMovementConstraint<DataTypes>::applyConstraint(defaulttype::BaseVector *vect, unsigned int offset)
{
    //sout << "applyConstraint in Vector with offset = " << offset << sendl;
    const unsigned int N = Deriv::size();

    const SetIndexArray & indices = m_indices.getValue();
    for (SetIndexArray::const_iterator it = indices.begin(); it != indices.end(); ++it)
    {
        for (unsigned int c=0; c<N; ++c)
            vect->clear(offset + N * (*it) + c);
    }
}

//display the path the constrained dofs will go through
template <class DataTypes>
void LinearMovementConstraint<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    if (!vparams->displayFlags().getShowBehaviorModels() || m_keyTimes.getValue().size() == 0)
        return;
    if (showMovement.getValue())
    {
        glDisable(GL_LIGHTING);
        glPointSize(10);
        glColor4f(1, 0.5, 0.5, 1);
        glBegin(GL_LINES);
        const SetIndexArray & indices = m_indices.getValue();
        for (unsigned int i = 0; i < m_keyMovements.getValue().size() - 1; i++)
        {
            for (SetIndexArray::const_iterator it = indices.begin(); it != indices.end(); ++it)
            {
                gl::glVertexT(DataTypes::getCPos(x0[*it]) + DataTypes::getDPos(m_keyMovements.getValue()[i]));
                gl::glVertexT(DataTypes::getCPos(x0[*it]) + DataTypes::getDPos(m_keyMovements.getValue()[i + 1]));
            }
        }
        glEnd();
    }
    else
    {
        const VecCoord& x = *this->mstate->getX();

        sofa::helper::vector<Vector3> points;
        Vector3 point;
        const SetIndexArray & indices = m_indices.getValue();
        for (SetIndexArray::const_iterator it = indices.begin(); it != indices.end(); ++it)
        {
            point = DataTypes::getCPos(x[*it]);
            points.push_back(point);
        }
        vparams->drawTool()->drawPoints(points, 10, Vec<4, float> (1, 0.5, 0.5, 1));
    }
}

} // namespace constraint

} // namespace component

} // namespace sofa

#endif

