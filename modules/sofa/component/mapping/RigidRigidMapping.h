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
#ifndef SOFA_COMPONENT_MAPPING_RIGIDRIGIDMAPPING_H
#define SOFA_COMPONENT_MAPPING_RIGIDRIGIDMAPPING_H

#include <sofa/component/component.h>

#include <sofa/core/Mapping.h>
#include <sofa/core/objectmodel/DataFileName.h>

#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/core/visual/VisualParams.h>
#include <vector>

namespace sofa
{

namespace component
{

namespace mapping
{

using namespace sofa::defaulttype;

template <class TIn, class TOut>
class RigidRigidMapping : public core::Mapping<TIn, TOut>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(RigidRigidMapping,TIn,TOut), SOFA_TEMPLATE2(core::Mapping,TIn,TOut));

    typedef core::Mapping<TIn, TOut> Inherit;
    typedef TIn In;
    typedef TOut Out;
    typedef Out OutDataTypes;
    typedef typename Out::VecCoord OutVecCoord;
    typedef typename Out::VecDeriv OutVecDeriv;
    typedef typename Out::Coord OutCoord;
    typedef typename Out::Deriv OutDeriv;
    typedef typename Out::MatrixDeriv OutMatrixDeriv;
    typedef typename In::Coord InCoord;
    typedef typename In::Deriv InDeriv;
    typedef typename In::VecCoord InVecCoord;
    typedef typename In::VecDeriv InVecDeriv;
    typedef typename In::MatrixDeriv InMatrixDeriv;
    typedef typename Out::Coord::value_type Real;
    enum { N=OutDataTypes::spatial_dimensions };
    typedef defaulttype::Mat<N,N,Real> Mat;
    typedef Vec<N,Real> Vector ;

protected:
    Data < OutVecCoord > points;
    OutVecCoord pointsR0;
    Mat rotation;
    class Loader;
    void load(const char* filename);
    Data< sofa::helper::vector<unsigned int> >  repartition;

public:
    Data<unsigned> index;
    sofa::core::objectmodel::DataFileName fileRigidRigidMapping;
    //axis length for display
    Data<double> axisLength;
    Data< bool > indexFromEnd;
    Data< bool > globalToLocalCoords;

    helper::ParticleMask* maskFrom;
    helper::ParticleMask* maskTo;
protected:
    RigidRigidMapping()
        : Inherit(),
          points(initData(&points, "initialPoints", "Initial position of the points")),
          repartition(initData(&repartition,"repartition","number of child frames per parent frame. If empty, all the children are attached to the parent with index given in the \"index\" attribute. If one value, each parent frame drives the given number of children frames. Otherwise, the values are the number of child frames driven by each parent frame. ")),
          index(initData(&index,(unsigned)0,"index","input frame index")),
          fileRigidRigidMapping(initData(&fileRigidRigidMapping,"fileRigidRigidMapping","Filename")),
          axisLength(initData( &axisLength, 0.7, "axisLength", "axis length for display")),
          indexFromEnd( initData ( &indexFromEnd,false,"indexFromEnd","input DOF index starts from the end of input DOFs vector") ),
          globalToLocalCoords ( initData ( &globalToLocalCoords,"globalToLocalCoords","are the output DOFs initially expressed in global coordinates" ) ),
          maskFrom(NULL),
          maskTo(NULL)
    {
        this->addAlias(&fileRigidRigidMapping,"filename");

    }

    virtual ~RigidRigidMapping()
    {
    }
public:
    virtual void init();

    virtual void apply(const core::MechanicalParams *mparams /* PARAMS FIRST */, Data<OutVecCoord>& out, const Data<InVecCoord>& in);

    virtual void applyJ(const core::MechanicalParams *mparams /* PARAMS FIRST */, Data<OutVecDeriv>& out, const Data<InVecDeriv>& in);

    virtual void applyJT(const core::MechanicalParams *mparams /* PARAMS FIRST */, Data<InVecDeriv>& out, const Data<OutVecDeriv>& in);

    virtual void applyJT(const core::ConstraintParams *cparams /* PARAMS FIRST */, Data<InMatrixDeriv>& out, const Data<OutMatrixDeriv>& in);

    virtual void computeAccFromMapping(const core::MechanicalParams *mparams /* PARAMS FIRST */, Data<OutVecDeriv>& acc_out, const Data<InVecDeriv>& v_in, const Data<InVecDeriv>& acc_in);

    virtual void applyDJT(const core::MechanicalParams* mparams /* PARAMS FIRST  = core::MechanicalParams::defaultInstance()*/, core::MultiVecDerivId parentForce, core::ConstMultiVecDerivId  childForce );

    virtual const sofa::defaulttype::BaseMatrix* getJ()
    {
        return NULL;
    }

    void draw(const core::visual::VisualParams* vparams);

    void clear();

    sofa::helper::vector<unsigned int> getRepartition() {return repartition.getValue(); }

    void setRepartition(unsigned int value);
    void setRepartition(sofa::helper::vector<unsigned int> values);

protected:

    bool getShow(const core::objectmodel::BaseObject* /*m*/, const core::visual::VisualParams* vparams) const { return vparams->displayFlags().getShowMappings(); }

    bool getShow(const core::BaseMapping* /*m*/, const core::visual::VisualParams* vparams) const { return vparams->displayFlags().getShowMechanicalMappings(); }
};

#ifndef SOFA_FLOAT
using sofa::defaulttype::Rigid2dTypes;
using sofa::defaulttype::Rigid3dTypes;
#endif

#ifndef SOFA_DOUBLE
using sofa::defaulttype::Rigid2fTypes;
using sofa::defaulttype::Rigid3fTypes;
#endif

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_MAPPING_RIGIDRIGIDMAPPING_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_RIGID_API RigidRigidMapping< Rigid3dTypes, Rigid3dTypes >;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_RIGID_API RigidRigidMapping< Rigid3fTypes, Rigid3fTypes >;
#endif

#ifndef SOFA_FLOAT
#ifndef SOFA_DOUBLE
extern template class SOFA_RIGID_API RigidRigidMapping< Rigid3dTypes, Rigid3fTypes >;
extern template class SOFA_RIGID_API RigidRigidMapping< Rigid3fTypes, Rigid3dTypes >;
#endif
#endif
#endif

} // namespace mapping

} // namespace component

} // namespace sofa

#endif
