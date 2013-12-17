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
#ifndef SOFA_COMPONENT_MAPPING_SUBSETMAPPING_H
#define SOFA_COMPONENT_MAPPING_SUBSETMAPPING_H


#include <sofa/core/Mapping.h>

#include <sofa/component/topology/TopologySubsetData.h>
#include <sofa/component/linearsolver/CompressedRowSparseMatrix.h>

#include <sofa/helper/vector.h>
#include <sofa/core/topology/BaseMeshTopology.h>

#include <memory>

namespace sofa
{

namespace component
{

namespace mapping
{

/// This class can be overridden if needed for additionnal storage within template specializations.
template<class InDataTypes, class OutDataTypes>
class SubsetMappingInternalData
{
public:
};

/**
 * @class SubsetMapping
 * @brief Compute a subset of input points
 */
template <class TIn, class TOut>
class SubsetMapping : public core::Mapping<TIn, TOut>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(SubsetMapping,TIn,TOut), SOFA_TEMPLATE2(core::Mapping,TIn,TOut));

    typedef core::Mapping<TIn, TOut> Inherit;
    typedef TIn In;
    typedef TOut Out;

    typedef typename In::Real         Real;
    typedef typename In::VecCoord     InVecCoord;
    typedef typename In::VecDeriv     InVecDeriv;
    typedef typename In::MatrixDeriv  InMatrixDeriv;
    typedef Data<InVecCoord>          InDataVecCoord;
    typedef Data<InVecDeriv>          InDataVecDeriv;
    typedef Data<InMatrixDeriv>       InDataMatrixDeriv;
    typedef typename In::Coord        InCoord;
    typedef typename In::Deriv        InDeriv;

    typedef typename Out::VecCoord    OutVecCoord;
    typedef typename Out::VecDeriv    OutVecDeriv;
    typedef typename Out::MatrixDeriv OutMatrixDeriv;
    typedef Data<OutVecCoord>         OutDataVecCoord;
    typedef Data<OutVecDeriv>         OutDataVecDeriv;
    typedef Data<OutMatrixDeriv>      OutDataMatrixDeriv;
    typedef typename Out::Coord       OutCoord;
    typedef typename Out::Deriv       OutDeriv;

    enum { NIn = sofa::defaulttype::DataTypeInfo<InDeriv>::Size };
    enum { NOut = sofa::defaulttype::DataTypeInfo<OutDeriv>::Size };
    typedef defaulttype::Mat<NOut, NIn, Real> MBloc;
    typedef sofa::component::linearsolver::CompressedRowSparseMatrix<MBloc> MatrixType;

    /// Correspondance array
    typedef typename InVecCoord::template rebind<unsigned int>::other IndexArray;
    typedef sofa::component::topology::PointSubsetData< IndexArray > SetIndex;
    SetIndex f_indices;

	Data < int > f_definedSubsetSize;
    Data < int > f_first;
    Data < int > f_last;
    Data < Real > f_radius;
    Data < bool > f_handleTopologyChange;
    SubsetMappingInternalData<In, Out> data;
    void postInit();
protected:
    SubsetMapping();
public:
    void clear(int reserve);

    int addPoint(int index);

    void init();

    // handle topology changes depending on the topology
    //void handleTopologyChange(core::topology::Topology* t);

    virtual ~SubsetMapping();

    virtual void apply ( const core::MechanicalParams* mparams /* PARAMS FIRST */, OutDataVecCoord& dOut, const InDataVecCoord& dIn );

    virtual void applyJ( const core::MechanicalParams* mparams /* PARAMS FIRST */, OutDataVecDeriv& dOut, const InDataVecDeriv& dIn );

    virtual void applyJT ( const core::MechanicalParams* mparams /* PARAMS FIRST */, InDataVecDeriv& dOut, const OutDataVecDeriv& dIn );

    virtual void applyJT ( const core::ConstraintParams* /*cparams*/ /* PARAMS FIRST */, InDataMatrixDeriv& dOut, const OutDataMatrixDeriv& dIn);

    const sofa::defaulttype::BaseMatrix* getJ();

protected:
    std::auto_ptr<MatrixType> matrixJ;
    bool updateJ;

    /// Pointer to the current topology
    sofa::core::topology::BaseMeshTopology* topology;
};

#ifndef SOFA_FLOAT
using sofa::defaulttype::Vec3dTypes;
using sofa::defaulttype::Vec1dTypes;
using sofa::defaulttype::Rigid3dTypes;
#endif
#ifndef SOFA_DOUBLE
using sofa::defaulttype::Vec3fTypes;
using sofa::defaulttype::Vec1fTypes;
using sofa::defaulttype::ExtVec3fTypes;
using sofa::defaulttype::Rigid3fTypes;
#endif

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_MAPPING_SUBSETMAPPING_CPP)

#ifndef SOFA_FLOAT
extern template class SOFA_BASE_MECHANICS_API SubsetMapping< Vec3dTypes, Vec3dTypes >;
extern template class SOFA_BASE_MECHANICS_API SubsetMapping< Vec1dTypes, Vec1dTypes >;
extern template class SOFA_BASE_MECHANICS_API SubsetMapping< Vec3dTypes, ExtVec3fTypes >;
extern template class SOFA_BASE_MECHANICS_API SubsetMapping< Rigid3dTypes, Rigid3dTypes >;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_BASE_MECHANICS_API SubsetMapping< Vec3fTypes, Vec3fTypes >;
extern template class SOFA_BASE_MECHANICS_API SubsetMapping< Vec1fTypes, Vec1fTypes >;
extern template class SOFA_BASE_MECHANICS_API SubsetMapping< Vec3fTypes, ExtVec3fTypes >;
extern template class SOFA_BASE_MECHANICS_API SubsetMapping< Rigid3fTypes, Rigid3fTypes >;
#endif

#ifndef SOFA_FLOAT
#ifndef SOFA_DOUBLE
extern template class SOFA_BASE_MECHANICS_API SubsetMapping< Vec3dTypes, Vec3fTypes >;
extern template class SOFA_BASE_MECHANICS_API SubsetMapping< Vec3fTypes, Vec3dTypes >;
extern template class SOFA_BASE_MECHANICS_API SubsetMapping< Vec1dTypes, Vec1fTypes >;
extern template class SOFA_BASE_MECHANICS_API SubsetMapping< Vec1fTypes, Vec1dTypes >;
extern template class SOFA_BASE_MECHANICS_API SubsetMapping< Rigid3dTypes, Rigid3fTypes >;
extern template class SOFA_BASE_MECHANICS_API SubsetMapping< Rigid3fTypes, Rigid3dTypes >;
#endif
#endif

#endif

} // namespace mapping

} // namespace component

} // namespace sofa

#endif
