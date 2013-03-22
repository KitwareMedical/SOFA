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
#ifndef SOFA_COMPONENT_PROJECTIVECONSTRAINTSET_PointConstraint_INL
#define SOFA_COMPONENT_PROJECTIVECONSTRAINTSET_PointConstraint_INL

#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/behavior/ProjectiveConstraintSet.inl>
#include <sofa/component/projectiveconstraintset/PointConstraint.h>
#include <sofa/component/linearsolver/SparseMatrix.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/simulation/common/Simulation.h>
#include <sofa/helper/gl/template.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <iostream>
using std::cerr;
using std::endl;
#include <sofa/component/topology/TopologySubsetData.inl>


#include <sofa/helper/gl/BasicShapes.h>




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



template <class DataTypes>
PointConstraint<DataTypes>::PointConstraint()
    : core::behavior::ProjectiveConstraintSet<DataTypes>(NULL)
    , f_indices( initData(&f_indices,"indices","Indices of the fixed points") )
    , _drawSize( initData(&_drawSize,0.0,"drawSize","0 -> point based rendering, >0 -> radius of spheres") )
{
}


template <class DataTypes>
PointConstraint<DataTypes>::~PointConstraint()
{
}


// -- Constraint interface


template <class DataTypes>
void PointConstraint<DataTypes>::init()
{
    this->core::behavior::ProjectiveConstraintSet<DataTypes>::init();

    //  cerr<<"PointConstraint<DataTypes>::init(), getJ = " << *getJ(0) << endl;

}

/// Update and return the jacobian. @todo update it when needed using topological engines instead of recomputing it at each call.
template <class DataTypes>
const sofa::defaulttype::BaseMatrix*  PointConstraint<DataTypes>::getJ(const core::MechanicalParams* )
{
    unsigned numBlocks = this->mstate->getSize();
    unsigned blockSize = DataTypes::deriv_total_size;
    jacobian.resize( numBlocks*blockSize,numBlocks*blockSize );
    //    cerr<<"PointConstraint<DataTypes>::getJ, numblocs = "<< numBlocks << ", block size = " << blockSize << endl;

    for(unsigned i=0; i<numBlocks*blockSize; i++ )
    {
        jacobian.set(i,i,1);
    }
    for(unsigned i=0; i<f_indices.getValue().size(); i++ )
    {
//                cerr<<"PointConstraint<DataTypes>::getJ, , set null block in " << f_indices.getValue()[i] << ", matrix before = " << jacobian << endl;
        for(unsigned j=0; j<blockSize; j++)
            jacobian.set( f_indices.getValue()[i]*blockSize+j, f_indices.getValue()[i]*blockSize+j, 0);
        cerr<<"PointConstraint<DataTypes>::getJ, , set null block in " << f_indices.getValue()[i] << ", matrix after = " << jacobian << endl;
    }

    return &jacobian;
}


template <class DataTypes>
void PointConstraint<DataTypes>::projectResponse(const core::MechanicalParams* mparams, DataVecDeriv& resData)
{
    //    cerr<<"PointConstraint<DataTypes>::projectResponse is called "<<endl;
    //    assert(false);
    helper::WriteAccessor<DataVecDeriv> res ( mparams, resData );
    const SetIndexArray & indices = f_indices.getValue(mparams);
    //serr<<"PointConstraint<DataTypes>::projectResponse, dx.size()="<<res.size()<<sendl;
    for (SetIndexArray::const_iterator it = indices.begin();
            it != indices.end();
            ++it)
    {
        res[*it] = Deriv();
    }
}

template <class DataTypes>
void PointConstraint<DataTypes>::projectJacobianMatrix(const core::MechanicalParams* mparams, DataMatrixDeriv& cData)
{
    helper::WriteAccessor<DataMatrixDeriv> c ( mparams, cData );
    const SetIndexArray & indices = f_indices.getValue(mparams);

    MatrixDerivRowIterator rowIt = c->begin();
    MatrixDerivRowIterator rowItEnd = c->end();

    while (rowIt != rowItEnd)
    {
        for (SetIndexArray::const_iterator it = indices.begin();
                it != indices.end();
                ++it)
        {
            rowIt.row().erase(*it);
        }
        ++rowIt;
    }
    //cerr<<"PointConstraint<DataTypes>::projectJacobianMatrix : helper::WriteAccessor<DataMatrixDeriv> c =  "<<endl<< c<<endl;
}

// constant velocity: do not change the velocity
template <class DataTypes>
void PointConstraint<DataTypes>::projectVelocity(const core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecDeriv& /*vData*/)
{
}

template <class DataTypes>
void PointConstraint<DataTypes>::projectPosition(const core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecCoord& /*xData*/)
{
    // nothing to do
}

// Matrix Integration interface
template <class DataTypes>
void PointConstraint<DataTypes>::applyConstraint(defaulttype::BaseMatrix *mat, unsigned int offset)
{
    //sout << "applyConstraint in Matrix with offset = " << offset << sendl;
    //cerr<<"PointConstraint<DataTypes>::applyConstraint(defaulttype::BaseMatrix *mat, unsigned int offset) is called "<<endl;
    const unsigned int N = Deriv::size();
    const SetIndexArray & indices = f_indices.getValue();

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
void PointConstraint<DataTypes>::applyConstraint(defaulttype::BaseVector *vect, unsigned int offset)
{
    //cerr<<"PointConstraint<DataTypes>::applyConstraint(defaulttype::BaseVector *vect, unsigned int offset) is called "<<endl;
    //sout << "applyConstraint in Vector with offset = " << offset << sendl;
    const unsigned int N = Deriv::size();

    const SetIndexArray & indices = f_indices.getValue();
    for (SetIndexArray::const_iterator it = indices.begin(); it != indices.end(); ++it)
    {
        for (unsigned int c=0; c<N; ++c)
            vect->clear(offset + N * (*it) + c);
    }
}




template <class DataTypes>
void PointConstraint<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
#ifndef SOFA_NO_OPENGL
    if (!vparams->displayFlags().getShowBehaviorModels()) return;
    if (!this->isActive()) return;
    const VecCoord& x = *this->mstate->getX();
    //serr<<"PointConstraint<DataTypes>::draw(), x.size() = "<<x.size()<<sendl;




    const SetIndexArray & indices = f_indices.getValue();

    if( _drawSize.getValue() == 0) // old classical drawing by points
    {
        std::vector< Vector3 > points;
        Vector3 point;
        //serr<<"PointConstraint<DataTypes>::draw(), indices = "<<indices<<sendl;
        for (SetIndexArray::const_iterator it = indices.begin();
                it != indices.end();
                ++it)
        {
            point = DataTypes::getCPos(x[*it]);
            points.push_back(point);
        }
        vparams->drawTool()->drawPoints(points, 10, Vec<4,float>(1,0.5,0.5,1));
    }
    else // new drawing by spheres
    {
        std::vector< Vector3 > points;
        Vector3 point;
        glColor4f (1.0f,0.35f,0.35f,1.0f);
        for (SetIndexArray::const_iterator it = indices.begin();
                it != indices.end();
                ++it)
        {
            point = DataTypes::getCPos(x[*it]);
            points.push_back(point);
        }
        vparams->drawTool()->drawSpheres(points, (float)_drawSize.getValue(), Vec<4,float>(1.0f,0.35f,0.35f,1.0f));
    }
#endif /* SOFA_NO_OPENGL */
}

// Specialization for rigids
#ifndef SOFA_FLOAT
template <>
void PointConstraint<Rigid3dTypes >::draw(const core::visual::VisualParams* vparams);
template <>
void PointConstraint<Rigid2dTypes >::draw(const core::visual::VisualParams* vparams);
#endif
#ifndef SOFA_DOUBLE
template <>
void PointConstraint<Rigid3fTypes >::draw(const core::visual::VisualParams* vparams);
template <>
void PointConstraint<Rigid2fTypes >::draw(const core::visual::VisualParams* vparams);
#endif



} // namespace constraint

} // namespace component

} // namespace sofa

#endif


