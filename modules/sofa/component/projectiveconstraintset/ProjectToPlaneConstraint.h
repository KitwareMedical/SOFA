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
#ifndef SOFA_COMPONENT_PROJECTIVECONSTRAINTSET_ProjectToPlaneConstraint_H
#define SOFA_COMPONENT_PROJECTIVECONSTRAINTSET_ProjectToPlaneConstraint_H

#include <sofa/core/behavior/ProjectiveConstraintSet.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/defaulttype/BaseMatrix.h>
#include <sofa/defaulttype/BaseVector.h>
#include <sofa/defaulttype/VecTypes.h>
//#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/helper/vector.h>
#include <sofa/defaulttype/Mat.h>
#include <sofa/component/topology/TopologySubsetData.h>
#include <sofa/component/linearsolver/CompressedRowSparseMatrix.h>
#include <sofa/component/linearsolver/EigenSparseMatrix.h>
#include <set>

namespace sofa
{

namespace component
{

namespace projectiveconstraintset
{

using helper::vector;
using core::objectmodel::Data;
using namespace sofa::core::objectmodel;
using namespace sofa::component::topology;

/// This class can be overridden if needed for additionnal storage within template specializations.
template <class DataTypes>
class ProjectToPlaneConstraintInternalData
{

};

/** Project particles to an affine plane.
  @author Francois Faure, 2012
  @todo Optimized versions for planes parallel to the main directions
*/
template <class DataTypes>
class ProjectToPlaneConstraint : public core::behavior::ProjectiveConstraintSet<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(ProjectToPlaneConstraint,DataTypes),SOFA_TEMPLATE(sofa::core::behavior::ProjectiveConstraintSet, DataTypes));

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::MatrixDeriv MatrixDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename MatrixDeriv::RowIterator MatrixDerivRowIterator;
    typedef typename MatrixDeriv::RowType MatrixDerivRowType;
    typedef Data<VecCoord> DataVecCoord;
    typedef Data<VecDeriv> DataVecDeriv;
    typedef Data<MatrixDeriv> DataMatrixDeriv;
    typedef helper::vector<unsigned int> Indices;
    typedef sofa::component::topology::PointSubsetData< Indices > IndexSubsetData;
    typedef linearsolver::EigenBaseSparseMatrix<SReal> BaseSparseMatrix;
    typedef linearsolver::EigenSparseMatrix<DataTypes,DataTypes> SparseMatrix;
    typedef typename SparseMatrix::Block Block;                                       ///< projection matrix of a particle displacement to the plane
    enum {bsize=SparseMatrix::Nin};                                                   ///< size of a block


protected:
    ProjectToPlaneConstraint();

    virtual ~ProjectToPlaneConstraint();

public:
    IndexSubsetData f_indices;     ///< the particles to project
    Data<double> f_drawSize;
    Data<Coord> f_origin;  ///< A point in the plane
    Data<Deriv> f_normal;  ///< The normal to the plane. Normalized at init()


protected:
    ProjectToPlaneConstraintInternalData<DataTypes>* data;
    friend class ProjectToPlaneConstraintInternalData<DataTypes>;


public:
    void clearConstraints();
    void addConstraint(unsigned int index);
    void removeConstraint(unsigned int index);

    // -- Constraint interface
    virtual void init();
    virtual void reinit();

    void projectResponse(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& resData);
    void projectVelocity(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& vData);
    void projectPosition(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecCoord& xData);
    void projectJacobianMatrix(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataMatrixDeriv& cData);


    void applyConstraint(defaulttype::BaseMatrix *mat, unsigned int offset);
    void applyConstraint(defaulttype::BaseVector *vect, unsigned int offset);

    /** Project the the given matrix (Experimental API).
      Replace M with PMP, where P is the projection matrix corresponding to the projectResponse method, shifted by the given offset, i.e. P is the identity matrix with a block on the diagonal replaced by the projection matrix.
      */
    virtual void projectMatrix( sofa::defaulttype::BaseMatrix* /*M*/, unsigned /*offset*/ );


    virtual void draw(const core::visual::VisualParams* vparams);


    class FCPointHandler : public TopologySubsetDataHandler<Point, Indices >
    {
    public:
        typedef typename ProjectToPlaneConstraint<DataTypes>::Indices Indices;

        FCPointHandler(ProjectToPlaneConstraint<DataTypes>* _fc, PointSubsetData<Indices>* _data)
            : sofa::component::topology::TopologySubsetDataHandler<Point, Indices >(_data), fc(_fc) {}



        void applyDestroyFunction(unsigned int /*index*/, value_type& /*T*/);


        bool applyTestCreateFunction(unsigned int /*index*/,
                const sofa::helper::vector< unsigned int > & /*ancestors*/,
                const sofa::helper::vector< double > & /*coefs*/);
    protected:
        ProjectToPlaneConstraint<DataTypes> *fc;
    };

protected :
    /// Pointer to the current topology
    sofa::core::topology::BaseMeshTopology* topology;

    /// Handler for subset Data
    FCPointHandler* pointHandler;

    SparseMatrix jacobian; ///< projection matrix in local state
    SparseMatrix J;        ///< auxiliary variable
};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_PROJECTIVECONSTRAINTSET_ProjectToPlaneConstraint_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_BOUNDARY_CONDITION_API ProjectToPlaneConstraint<defaulttype::Vec3dTypes>;
extern template class SOFA_BOUNDARY_CONDITION_API ProjectToPlaneConstraint<defaulttype::Vec2dTypes>;
//extern template class SOFA_BOUNDARY_CONDITION_API ProjectToPlaneConstraint<defaulttype::Vec1dTypes>;
//extern template class SOFA_BOUNDARY_CONDITION_API ProjectToPlaneConstraint<defaulttype::Vec6dTypes>;
//extern template class SOFA_BOUNDARY_CONDITION_API ProjectToPlaneConstraint<defaulttype::Rigid3dTypes>;
//extern template class SOFA_BOUNDARY_CONDITION_API ProjectToPlaneConstraint<defaulttype::Rigid2dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_BOUNDARY_CONDITION_API ProjectToPlaneConstraint<defaulttype::Vec3fTypes>;
extern template class SOFA_BOUNDARY_CONDITION_API ProjectToPlaneConstraint<defaulttype::Vec2fTypes>;
//extern template class SOFA_BOUNDARY_CONDITION_API ProjectToPlaneConstraint<defaulttype::Vec1fTypes>;
//extern template class SOFA_BOUNDARY_CONDITION_API ProjectToPlaneConstraint<defaulttype::Vec6fTypes>;
//extern template class SOFA_BOUNDARY_CONDITION_API ProjectToPlaneConstraint<defaulttype::Rigid3fTypes>;
//extern template class SOFA_BOUNDARY_CONDITION_API ProjectToPlaneConstraint<defaulttype::Rigid2fTypes>;
#endif
#endif

} // namespace projectiveconstraintset

} // namespace component

} // namespace sofa


#endif
