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
#ifndef SOFA_COMPONENT_MAPPING_TriangleDeformationMapping_H
#define SOFA_COMPONENT_MAPPING_TriangleDeformationMapping_H

#include <sofa/core/Mapping.inl>
#include <sofa/component/linearsolver/EigenSparseMatrix.h>
#include <sofa/component/topology/TriangleSetTopologyContainer.h>
#include <sofa/defaulttype/Mat.h>
#include <sofa/defaulttype/Vec.h>
#include "../initFlexible.h"
#include "../types/DeformationGradientTypes.h"
#include "../shapeFunction/BaseShapeFunction.h"


namespace sofa
{
using helper::vector;

namespace component
{

namespace mapping
{

/// This class can be overridden if needed for additionnal storage within template specializations.
template<class InDataTypes, class OutDataTypes>
class TriangleDeformationMappingInternalData
{
public:
};


/** Maps triangle vertex positions to deformation gradients.

@author Francois Faure
  */
template <class TIn, class TOut>
class TriangleDeformationMapping : public core::Mapping<TIn, TOut>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(TriangleDeformationMapping,TIn,TOut), SOFA_TEMPLATE2(core::Mapping,TIn,TOut));

    typedef core::Mapping<TIn, TOut> Inherit;
    typedef TIn In;
    typedef TOut Out;
    typedef typename Out::VecCoord OutVecCoord;
    typedef typename Out::VecDeriv OutVecDeriv;
    typedef typename Out::Coord OutCoord;
    typedef typename Out::Deriv OutDeriv;
    typedef typename Out::MatrixDeriv OutMatrixDeriv;
    typedef typename Out::Real Real;
//    typedef typename Out::Frame Frame;
    typedef defaulttype::Mat<3,2,Real> Frame;
    typedef typename In::Deriv InDeriv;
    typedef typename In::MatrixDeriv InMatrixDeriv;
    typedef typename In::Coord InCoord;
    typedef typename In::VecCoord InVecCoord;
    typedef typename In::VecDeriv InVecDeriv;
    typedef linearsolver::EigenSparseMatrix<TIn,TOut>    SparseMatrixEigen;
    enum {Nin = In::deriv_total_size, Nout = Out::deriv_total_size };
    typedef defaulttype::Mat<In::deriv_total_size, In::deriv_total_size,Real>  InBlock;
    typedef defaulttype::Mat<Out::deriv_total_size, In::deriv_total_size,Real>  Block;
    typedef topology::TriangleSetTopologyContainer::SeqTriangles SeqTriangles;

    typedef core::behavior::ShapeFunctionTypes<3,2,Real> ShapeFunctionType;             // 2d shape function
    typedef core::behavior::BaseShapeFunction<ShapeFunctionType> ShapeFunction;
    typedef defaulttype::Vec<2,Real> MCoord;                                     ///< material coordinates
    typedef vector<MCoord> VMCoord;                                   ///< vector of material coordinates
    typedef defaulttype::Mat<2,2,Real> MMat;                                      ///< matrix in material coordinates
    typedef vector<MMat> VMMat;                                              ///< vector of material matrices, used to compute the deformation gradients


    Data< VMMat > f_inverseRestEdges;  ///< For each triangle, inverse matrix of edge12, edge13, normal. This is used to compute the deformation gradient based on the current edges.
    Data< SReal > f_scaleView; ///< scaling factor for the drawing of the deformation gradient

    virtual void init();

    virtual void apply(const core::MechanicalParams *mparams, Data<OutVecCoord>& out, const Data<InVecCoord>& in);

    virtual void applyJ(const core::MechanicalParams *mparams, Data<OutVecDeriv>& out, const Data<InVecDeriv>& in);

    virtual void applyJT(const core::MechanicalParams *mparams, Data<InVecDeriv>& out, const Data<OutVecDeriv>& in);

    virtual void applyJT(const core::ConstraintParams *cparams, Data<InMatrixDeriv>& out, const Data<OutMatrixDeriv>& in);

//    virtual void computeGeometricStiffness(const core::MechanicalParams *mparams);

//    virtual void applyDJT(const core::MechanicalParams* mparams, core::MultiVecDerivId parentForce, core::ConstMultiVecDerivId  childForce );

    virtual const sofa::defaulttype::BaseMatrix* getJ();
    virtual const vector<sofa::defaulttype::BaseMatrix*>* getJs();

    virtual void draw(const core::visual::VisualParams* vparams);

protected:
    TriangleDeformationMapping();
    virtual ~TriangleDeformationMapping();

    topology::TriangleSetTopologyContainer* triangleContainer;  ///< where the edges are defined
    SparseMatrixEigen jacobian;                         ///< Jacobian of the mapping
    vector<defaulttype::BaseMatrix*> baseMatrices;      ///< Jacobian of the mapping, in a vector

    Block makeBlock( Real middle, Real bottom );  ///< helper for the creation of the jacobian
};


#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_MAPPING_TriangleDeformationMapping_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_Flexible_API TriangleDeformationMapping< Vec3dTypes, F321dTypes >;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_Flexible_API TriangleDeformationMapping< Vec3fTypes, F321fTypes >;
#endif

#endif

} // namespace mapping

} // namespace component

} // namespace sofa

#endif
