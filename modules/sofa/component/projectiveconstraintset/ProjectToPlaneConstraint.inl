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
#ifndef SOFA_COMPONENT_PROJECTIVECONSTRAINTSET_ProjectToPlaneConstraint_INL
#define SOFA_COMPONENT_PROJECTIVECONSTRAINTSET_ProjectToPlaneConstraint_INL

#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/behavior/ProjectiveConstraintSet.inl>
#include <sofa/component/projectiveconstraintset/ProjectToPlaneConstraint.h>
#include <sofa/component/linearsolver/SparseMatrix.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/simulation/common/Simulation.h>
#include <sofa/helper/gl/template.h>
//#include <sofa/defaulttype/RigidTypes.h>
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


// Define TestNewPointFunction
template< class DataTypes>
bool ProjectToPlaneConstraint<DataTypes>::FCPointHandler::applyTestCreateFunction(unsigned int, const sofa::helper::vector<unsigned int> &, const sofa::helper::vector<double> &)
{
    if (fc)
    {
        return false;
    }
    else
    {
        return false;
    }
}

// Define RemovalFunction
template< class DataTypes>
void ProjectToPlaneConstraint<DataTypes>::FCPointHandler::applyDestroyFunction(unsigned int pointIndex, value_type &)
{
    if (fc)
    {
        fc->removeConstraint((unsigned int) pointIndex);
    }
}

template <class DataTypes>
ProjectToPlaneConstraint<DataTypes>::ProjectToPlaneConstraint()
    : core::behavior::ProjectiveConstraintSet<DataTypes>(NULL)
    , f_indices( initData(&f_indices,"indices","Indices of the fixed points") )
    , f_drawSize( initData(&f_drawSize,0.0,"drawSize","0 -> point based rendering, >0 -> radius of spheres") )
    , f_origin( initData(&f_origin,Coord(),"origin","A point in the plane"))
    , f_normal( initData(&f_normal,Deriv(),"normal","Normal vector to the plane"))
    , data(new ProjectToPlaneConstraintInternalData<DataTypes>())
{
    // default to index 0
    f_indices.beginEdit()->push_back(0);
    f_indices.endEdit();

    pointHandler = new FCPointHandler(this, &f_indices);
}


template <class DataTypes>
ProjectToPlaneConstraint<DataTypes>::~ProjectToPlaneConstraint()
{
    if (pointHandler)
        delete pointHandler;

    delete data;
}

template <class DataTypes>
void ProjectToPlaneConstraint<DataTypes>::clearConstraints()
{
    f_indices.beginEdit()->clear();
    f_indices.endEdit();
}

template <class DataTypes>
void ProjectToPlaneConstraint<DataTypes>::addConstraint(unsigned int index)
{
    f_indices.beginEdit()->push_back(index);
    f_indices.endEdit();
}

template <class DataTypes>
void ProjectToPlaneConstraint<DataTypes>::removeConstraint(unsigned int index)
{
    removeValue(*f_indices.beginEdit(),index);
    f_indices.endEdit();
}

// -- Constraint interface


template <class DataTypes>
void ProjectToPlaneConstraint<DataTypes>::init()
{
    this->core::behavior::ProjectiveConstraintSet<DataTypes>::init();

    topology = this->getContext()->getMeshTopology();

    //  if (!topology)
    //    serr << "Can not find the topology." << sendl;

    // Initialize functions and parameters
    f_indices.createTopologicalEngine(topology, pointHandler);
    f_indices.registerTopologicalData();

    const Indices & indices = f_indices.getValue();

    unsigned int maxIndex=this->mstate->getSize();
    for (unsigned int i=0; i<indices.size(); ++i)
    {
        const unsigned int index=indices[i];
        if (index >= maxIndex)
        {
            serr << "Index " << index << " not valid!" << sendl;
            removeConstraint(index);
        }
    }

    reinit();

//  cerr<<"ProjectToPlaneConstraint<DataTypes>::init(), getJ = " << *getJ(0) << endl;

}

template <class DataTypes>
void  ProjectToPlaneConstraint<DataTypes>::reinit()
{
//    cerr<<"ProjectToPlaneConstraint<DataTypes>::getJ, numblocs = "<< numBlocks << ", block size = " << blockSize << endl;

    // normalize the normal vector
    Deriv n = f_normal.getValue();
    if( n.norm()==0 )
        n[1]=0;
    else n *= 1/n.norm();
    f_normal.setValue(n);

    // create the matrix blocks corresponding to the projection to the plane: I-nn^t or to the identity
    Block bProjection, bIdentity;
    for(unsigned i=0; i<bsize; i++)
        for(unsigned j=0; j<bsize; j++)
        {
            if(i==j)
            {
                bIdentity[i][j]   = 1;
                bProjection[i][j] = 1 - n[i]*n[j];
            }
            else
            {
                bIdentity[i][j]   = 0;
                bProjection[i][j] =    - n[i]*n[j];
            }
        }
//    cerr<<"ProjectToPlaneConstraint<DataTypes>::reinit() bIdentity[0] = " << endl << bIdentity[0] << endl;
//    cerr<<"ProjectToPlaneConstraint<DataTypes>::reinit() bProjection[0] = " << endl << bProjection[0] << endl;

    // get the indices sorted
    Indices tmp = f_indices.getValue();
    std::sort(tmp.begin(),tmp.end());

    // resize the jacobian
    unsigned numBlocks = this->mstate->getSize();
    unsigned blockSize = DataTypes::deriv_total_size;
    jacobian.resize( numBlocks*blockSize,numBlocks*blockSize );

    // fill the jacobian is ascending order
    Indices::const_iterator it= tmp.begin();
    unsigned i=0;
    for(Indices::const_iterator it= tmp.begin(); i<numBlocks; i++ )
    {
        jacobian.beginBlockRow(i);
        if( i==*it )  // constrained particle: set diagonal to projection block, and  the cursor to the next constraint
        {
            jacobian.createBlock(i,bProjection); // only one block to create
            it++;
        }
        else           // unconstrained particle: set diagonal to identity block
        {
            jacobian.createBlock(i,bIdentity); // only one block to create
        }
        jacobian.endBlockRow();
    }
    jacobian.compress();
//    cerr<<"ProjectToPlaneConstraint<DataTypes>::reinit(), jacobian = " << jacobian << endl;

}

template <class DataTypes>
void ProjectToPlaneConstraint<DataTypes>::projectMatrix( sofa::defaulttype::BaseMatrix* M, unsigned offset )
{
    J.copy(jacobian, M->colSize(), offset); // projection matrix for an assembled state
    BaseSparseMatrix* E = dynamic_cast<BaseSparseMatrix*>(M);
    assert(E);
    E->compressedMatrix = J.compressedMatrix * E->compressedMatrix * J.compressedMatrix;
}



template <class DataTypes>
void ProjectToPlaneConstraint<DataTypes>::projectResponse(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& resData)
{
    helper::WriteAccessor<DataVecDeriv> res ( mparams, resData );
    jacobian.mult(res.wref(),res.ref());
}

template <class DataTypes>
void ProjectToPlaneConstraint<DataTypes>::projectJacobianMatrix(const core::MechanicalParams* /*mparams*/ , DataMatrixDeriv& /*cData*/)
{
    serr<<"projectJacobianMatrix(const core::MechanicalParams*, DataMatrixDeriv& ) is not implemented" << sendl;
}

template <class DataTypes>
void ProjectToPlaneConstraint<DataTypes>::projectVelocity(const core::MechanicalParams* mparams, DataVecDeriv& vdata)
{
    projectResponse(mparams,vdata);
}

template <class DataTypes>
void ProjectToPlaneConstraint<DataTypes>::projectPosition(const core::MechanicalParams* mparams , DataVecCoord& xData)
{
    helper::WriteAccessor<DataVecCoord> x ( mparams, xData );
    const Deriv& n = f_normal.getValue();
    const Deriv& o = f_origin.getValue();

    const Indices& indices = f_indices.getValue();
    for(unsigned i=0; i<indices.size(); i++ )
    {
        // replace the point with its projection to the plane
        x[indices[i]] -= n * ((x[indices[i]]-o)*n);
    }

}

// Matrix Integration interface
template <class DataTypes>
void ProjectToPlaneConstraint<DataTypes>::applyConstraint(defaulttype::BaseMatrix * /*mat*/, unsigned int /*offset*/)
{
    serr << "applyConstraint is not implemented " << sendl;
}

template <class DataTypes>
void ProjectToPlaneConstraint<DataTypes>::applyConstraint(defaulttype::BaseVector * /*vect*/, unsigned int /*offset*/)
{
    serr<<"ProjectToPlaneConstraint<DataTypes>::applyConstraint(defaulttype::BaseVector *vect, unsigned int offset) is not implemented "<< sendl;
}




template <class DataTypes>
void ProjectToPlaneConstraint<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    if (!vparams->displayFlags().getShowBehaviorModels()) return;
    if (!this->isActive()) return;
    const VecCoord& x = *this->mstate->getX();

    const Indices & indices = f_indices.getValue();

    if( f_drawSize.getValue() == 0) // old classical drawing by points
    {
        std::vector< Vector3 > points;
        Vector3 point;
        //serr<<"ProjectToPlaneConstraint<DataTypes>::draw(), indices = "<<indices<<sendl;
        for (Indices::const_iterator it = indices.begin();
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
        for (Indices::const_iterator it = indices.begin();
                it != indices.end();
                ++it)
        {
            point = DataTypes::getCPos(x[*it]);
            points.push_back(point);
        }
        vparams->drawTool()->drawSpheres(points, (float)f_drawSize.getValue(), Vec<4,float>(1.0f,0.35f,0.35f,1.0f));
    }
}

//// Specialization for rigids
//#ifndef SOFA_FLOAT
//template <>
//    void ProjectToPlaneConstraint<Rigid3dTypes >::draw(const core::visual::VisualParams* vparams);
//template <>
//    void ProjectToPlaneConstraint<Rigid2dTypes >::draw(const core::visual::VisualParams* vparams);
//#endif
//#ifndef SOFA_DOUBLE
//template <>
//    void ProjectToPlaneConstraint<Rigid3fTypes >::draw(const core::visual::VisualParams* vparams);
//template <>
//    void ProjectToPlaneConstraint<Rigid2fTypes >::draw(const core::visual::VisualParams* vparams);
//#endif



} // namespace constraint

} // namespace component

} // namespace sofa

#endif


