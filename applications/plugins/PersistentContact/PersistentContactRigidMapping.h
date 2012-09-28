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
*                               SOFA :: Plugins                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_COMPONENT_MAPPING_PERSISTENTCONTACTRIGIDMAPPING_H
#define SOFA_COMPONENT_MAPPING_PERSISTENTCONTACTRIGIDMAPPING_H

#include <sofa/component/mapping/RigidMapping.h>

#include "PersistentContactMapping.h"
#include "PersistentContact.h"

namespace sofa
{

namespace component
{

namespace mapping
{

template <class TIn, class TOut>
class PersistentContactRigidMapping : public RigidMapping<TIn, TOut>, public PersistentContactMapping
{
public:
    SOFA_CLASS2(SOFA_TEMPLATE2(PersistentContactRigidMapping,TIn,TOut), SOFA_TEMPLATE2(RigidMapping,TIn,TOut), PersistentContactMapping);

    typedef RigidMapping<TIn, TOut> Inherit;
    typedef TIn In;
    typedef TOut Out;
    typedef Out OutDataTypes;
    typedef typename Out::VecCoord VecCoord;
    typedef typename Out::VecDeriv VecDeriv;
    typedef typename Out::Coord Coord;
    typedef typename Out::Deriv Deriv;
    typedef typename Out::MatrixDeriv OutMatrixDeriv;
    typedef typename In::Deriv InDeriv;
    typedef typename In::DRot DRot;
    typedef typename In::VecCoord InVecCoord;
    typedef typename In::VecDeriv InVecDeriv;
    typedef typename In::MatrixDeriv InMatrixDeriv;
    typedef typename Coord::value_type Real;
    enum
    {
        N = OutDataTypes::spatial_dimensions
    };
    enum
    {
        NIn = sofa::defaulttype::DataTypeInfo<InDeriv>::Size
    };
    enum
    {
        NOut = sofa::defaulttype::DataTypeInfo<Deriv>::Size
    };
    typedef defaulttype::Mat<N, N, Real> Mat;
    typedef defaulttype::Vec<N, Real> Vector;
    typedef defaulttype::Mat<NOut, NIn, Real> MBloc;
    typedef sofa::component::linearsolver::CompressedRowSparseMatrix<MBloc> MatrixType;

    PersistentContactRigidMapping();

    virtual ~PersistentContactRigidMapping() {}

    void beginAddContactPoint();

    int addContactPointFromInputMapping(const sofa::defaulttype::Vector3& pos, std::vector< std::pair<int, double> > & baryCoords);

    int keepContactPointFromInputMapping(const int index);

    void init();

    void bwdInit();

    void reset();

    void handleEvent(sofa::core::objectmodel::Event*);

    void storeFreePositionAndDx();

    void applyLinearizedPosition();

    void applyPositionAndFreePosition();

    void applyJT(const core::ConstraintParams *cparams  /* PARAMS FIRST */, Data<InMatrixDeriv>& out, const Data<OutMatrixDeriv>& in)
    {
        m_previousPosition = this->fromModel->read(core::ConstVecCoordId::position())->getValue();

//         std::cout<<"applyJT   m_previousPosition = "<<m_previousPosition<<std::endl;

        Inherit::applyJT(cparams  /* PARAMS FIRST */, out, in);
    }

protected:

    Inherit *m_inputMapping;
    bool m_init;
    VecCoord m_previousPoints;
    InVecCoord m_previousPosition;
    InVecCoord m_previousFreePosition;
    InVecDeriv m_previousDx;

    void setDefaultValues();
};

using sofa::defaulttype::Vec2dTypes;
using sofa::defaulttype::Vec3dTypes;
using sofa::defaulttype::Vec2fTypes;
using sofa::defaulttype::Vec3fTypes;
using sofa::defaulttype::Rigid2dTypes;
using sofa::defaulttype::Rigid3dTypes;
using sofa::defaulttype::Rigid2fTypes;
using sofa::defaulttype::Rigid3fTypes;

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_MAPPING_PERSISTENTCONTACTRIGIDMAPPING_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_PERSISTENTCONTACT_API PersistentContactRigidMapping< Rigid3dTypes, Vec3dTypes >;
extern template class SOFA_PERSISTENTCONTACT_API PersistentContactRigidMapping< Rigid2dTypes, Vec2dTypes >;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_PERSISTENTCONTACT_API PersistentContactRigidMapping< Rigid3fTypes, Vec3fTypes >;
extern template class SOFA_PERSISTENTCONTACT_API PersistentContactRigidMapping< Rigid2fTypes, Vec2fTypes >;
#endif

#ifndef SOFA_FLOAT
#ifndef SOFA_DOUBLE
extern template class SOFA_PERSISTENTCONTACT_API PersistentContactRigidMapping< Rigid3dTypes, Vec3fTypes >;
extern template class SOFA_PERSISTENTCONTACT_API PersistentContactRigidMapping< Rigid3fTypes, Vec3dTypes >;
extern template class SOFA_PERSISTENTCONTACT_API PersistentContactRigidMapping< Rigid2dTypes, Vec2fTypes >;
extern template class SOFA_PERSISTENTCONTACT_API PersistentContactRigidMapping< Rigid2fTypes, Vec2dTypes >;
#endif
#endif
#endif

} // namespace mapping

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_MAPPING_PERSISTENTCONTACTRIGIDMAPPING_H
