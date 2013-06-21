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
#ifndef SOFA_COMPONENT_MAPPING_BARYCENTRICMAPPINGRIGID_H
#define SOFA_COMPONENT_MAPPING_BARYCENTRICMAPPINGRIGID_H

#include <sofa/component/mapping/BarycentricMapping.h>
#include <sofa/defaulttype/RigidTypes.h>


namespace sofa
{

namespace component
{

namespace mapping
{


/// Class allowing barycentric mapping computation on a TetrahedronSetTopology in Vec3 -> Rigid case

template<class In, class Out>
class BarycentricMapperTetrahedronSetTopologyRigid : public TopologyBarycentricMapper<In,Out>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(BarycentricMapperTetrahedronSetTopologyRigid,In,Out),SOFA_TEMPLATE2(TopologyBarycentricMapper,In,Out));
    typedef TopologyBarycentricMapper<In,Out> Inherit;
    typedef typename Inherit::Real Real;
    typedef typename Inherit::OutReal OutReal;
    typedef typename Inherit::OutDeriv  OutDeriv;
    typedef typename Inherit::InDeriv  InDeriv;
    typedef typename Inherit::MappingData3D MappingData;
    typedef helper::fixed_array<MappingData,3> MappingOrientData;

    typedef typename In::VecCoord VecCoord;
    typedef typename In::VecDeriv VecDeriv;

    enum { NIn = Inherit::NIn };
    enum { NOut = Inherit::NOut };
    typedef typename Inherit::MBloc MBloc;
    typedef typename Inherit::MatrixType MatrixType;

protected:
    topology::PointData< sofa::helper::vector<MappingData > >  map;
    topology::PointData< sofa::helper::vector<MappingOrientData > >  mapOrient;

    VecCoord actualTetraPosition;

    topology::TetrahedronSetTopologyContainer*			_fromContainer;
    topology::TetrahedronSetGeometryAlgorithms<In>*	_fromGeomAlgo;
    helper::ParticleMask *maskFrom;
    helper::ParticleMask *maskTo;

    MatrixType* matrixJ;
    bool updateJ;

    /// TEMP
    VecDeriv actualOut;
    typename Out::VecCoord actualPos;

    /// TEMP

    BarycentricMapperTetrahedronSetTopologyRigid(topology::TetrahedronSetTopologyContainer* fromTopology, topology::PointSetTopologyContainer* _toTopology,
            helper::ParticleMask *_maskFrom,
            helper::ParticleMask *_maskTo)
        : TopologyBarycentricMapper<In,Out>(fromTopology, _toTopology),
          map(initData(&map,"map", "mapper data")),
          mapOrient(initData(&mapOrient,"mapOrient", "mapper data for mapped frames")),
          _fromContainer(fromTopology),
          _fromGeomAlgo(NULL),
          maskFrom(_maskFrom),
          maskTo(_maskTo),
          matrixJ(NULL),
          updateJ(true)
    {}

    virtual ~BarycentricMapperTetrahedronSetTopologyRigid() {}

public:
    void clear(int reserve=0);

    int addPointInTetra(const int index, const SReal* baryCoords);
    int addPointOrientationInTetra( const int tetraIndex, const Matrix3 baryCoorsOrient );

    void init(const typename Out::VecCoord& out, const typename In::VecCoord& in);

    void apply( typename Out::VecCoord& out, const typename In::VecCoord& in );
    void applyJ( typename Out::VecDeriv& out, const typename In::VecDeriv& in );
    void applyJT( typename In::VecDeriv& out, const typename Out::VecDeriv& in );
    void applyJT( typename In::MatrixDeriv& out, const typename Out::MatrixDeriv& in );

    virtual const sofa::defaulttype::BaseMatrix* getJ(int outSize, int inSize);

    void draw(const core::visual::VisualParams*,const typename Out::VecCoord& out, const typename In::VecCoord& in);

    //virtual int addContactPointFromInputMapping(const typename In::VecDeriv& in, const sofa::defaulttype::Vector3& /*pos*/, std::vector< std::pair<int, double> > & /*baryCoords*/);
};

template<class TInReal, class TOutReal>
class BarycentricMapperTetrahedronSetTopology< sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3,TInReal>,sofa::defaulttype::Vec<3,TInReal>,TInReal>, sofa::defaulttype::StdRigidTypes<3,TOutReal> > : public BarycentricMapperTetrahedronSetTopologyRigid< sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3,TInReal>,sofa::defaulttype::Vec<3,TInReal>,TInReal>, sofa::defaulttype::StdRigidTypes<3,TOutReal> >
{
public:
    typedef sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3,TInReal>,sofa::defaulttype::Vec<3,TInReal>,TInReal> In;
    typedef sofa::defaulttype::StdRigidTypes<3,TOutReal> Out;
    SOFA_CLASS(SOFA_TEMPLATE2(BarycentricMapperTetrahedronSetTopology,In,Out),SOFA_TEMPLATE2(BarycentricMapperTetrahedronSetTopologyRigid,In,Out));
    typedef BarycentricMapperTetrahedronSetTopologyRigid<In,Out> Inherit;

    BarycentricMapperTetrahedronSetTopology(topology::TetrahedronSetTopologyContainer* fromTopology, topology::PointSetTopologyContainer* _toTopology,
            helper::ParticleMask *_maskFrom,
            helper::ParticleMask *_maskTo)
        : Inherit(fromTopology, _toTopology, _maskFrom, _maskTo)
    {}

};

#ifndef SOFA_FLOAT
using sofa::defaulttype::Vec3dTypes;
using sofa::defaulttype::Rigid3dTypes;
#endif

#ifndef SOFA_DOUBLE
using sofa::defaulttype::Vec3fTypes;
using sofa::defaulttype::Rigid3fTypes;
#endif

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_MAPPING_BARYCENTRICMAPPINGRIGID_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_MISC_MAPPING_API BarycentricMapping< Vec3dTypes, Rigid3dTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapper< Vec3dTypes, Rigid3dTypes >;
extern template class SOFA_MISC_MAPPING_API TopologyBarycentricMapper< Vec3dTypes, Rigid3dTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperRegularGridTopology< Vec3dTypes, Rigid3dTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperSparseGridTopology< Vec3dTypes, Rigid3dTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperMeshTopology< Vec3dTypes, Rigid3dTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperEdgeSetTopology< Vec3dTypes, Rigid3dTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperTriangleSetTopology< Vec3dTypes, Rigid3dTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperQuadSetTopology< Vec3dTypes, Rigid3dTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperTetrahedronSetTopologyRigid< Vec3dTypes, Rigid3dTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperTetrahedronSetTopology< Vec3dTypes, Rigid3dTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperHexahedronSetTopology< Vec3dTypes, Rigid3dTypes >;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_MISC_MAPPING_API BarycentricMapping< Vec3fTypes, Rigid3fTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapper< Vec3fTypes, Rigid3fTypes >;
extern template class SOFA_MISC_MAPPING_API TopologyBarycentricMapper< Vec3fTypes, Rigid3fTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperRegularGridTopology< Vec3fTypes, Rigid3fTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperSparseGridTopology< Vec3fTypes, Rigid3fTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperMeshTopology< Vec3fTypes, Rigid3fTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperEdgeSetTopology< Vec3fTypes, Rigid3fTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperTriangleSetTopology< Vec3fTypes, Rigid3fTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperQuadSetTopology< Vec3fTypes, Rigid3fTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperTetrahedronSetTopologyRigid< Vec3fTypes, Rigid3fTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperTetrahedronSetTopology< Vec3fTypes, Rigid3fTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperHexahedronSetTopology< Vec3fTypes, Rigid3fTypes >;
#endif
#ifndef SOFA_FLOAT
#ifndef SOFA_DOUBLE
extern template class SOFA_MISC_MAPPING_API BarycentricMapping< Vec3dTypes, Rigid3fTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapping< Vec3fTypes, Rigid3dTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapper< Vec3dTypes, Rigid3fTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapper< Vec3fTypes, Rigid3dTypes >;
extern template class SOFA_MISC_MAPPING_API TopologyBarycentricMapper< Vec3dTypes, Rigid3fTypes >;
extern template class SOFA_MISC_MAPPING_API TopologyBarycentricMapper< Vec3fTypes, Rigid3dTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperRegularGridTopology< Vec3dTypes, Rigid3fTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperRegularGridTopology< Vec3fTypes, Rigid3dTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperSparseGridTopology< Vec3dTypes, Rigid3fTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperSparseGridTopology< Vec3fTypes, Rigid3dTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperMeshTopology< Vec3dTypes, Rigid3fTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperMeshTopology< Vec3fTypes, Rigid3dTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperEdgeSetTopology< Vec3dTypes, Rigid3fTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperEdgeSetTopology< Vec3fTypes, Rigid3dTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperTriangleSetTopology< Vec3dTypes, Rigid3fTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperTriangleSetTopology< Vec3fTypes, Rigid3dTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperQuadSetTopology< Vec3dTypes, Rigid3fTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperQuadSetTopology< Vec3fTypes, Rigid3dTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperTetrahedronSetTopologyRigid< Vec3dTypes, Rigid3fTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperTetrahedronSetTopology< Vec3dTypes, Rigid3fTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperTetrahedronSetTopologyRigid< Vec3fTypes, Rigid3dTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperTetrahedronSetTopology< Vec3fTypes, Rigid3dTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperHexahedronSetTopology< Vec3dTypes, Rigid3fTypes >;
extern template class SOFA_MISC_MAPPING_API BarycentricMapperHexahedronSetTopology< Vec3fTypes, Rigid3dTypes >;
#endif
#endif
#endif

} // namespace mapping

} // namespace component

} // namespace sofa

#endif
