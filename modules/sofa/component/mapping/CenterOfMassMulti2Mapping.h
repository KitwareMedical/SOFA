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
#ifndef SOFA_COMPONENT_MAPPING_CENTEROFMASSMULTI2MAPPING_H
#define SOFA_COMPONENT_MAPPING_CENTEROFMASSMULTI2MAPPING_H

#include <sofa/core/Multi2Mapping.h>
#include <sofa/core/behavior/BaseMass.h>
#include <sofa/core/VecId.h>

#include <sofa/defaulttype/Vec3Types.h>
#include <sofa/defaulttype/RigidTypes.h>

#include <sofa/component/component.h>

namespace sofa
{

namespace component
{

namespace mapping
{

template <class TIn1, class TIn2, class TOut>
class CenterOfMassMulti2Mapping : public core::Multi2Mapping<TIn1, TIn2, TOut>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE3(CenterOfMassMulti2Mapping, TIn1, TIn2, TOut), SOFA_TEMPLATE3(core::Multi2Mapping, TIn1, TIn2, TOut));

    typedef core::Multi2Mapping<TIn1, TIn2, TOut> Inherit;
    typedef TIn1 In1;
    typedef TIn2 In2;
    typedef TOut Out;

    typedef In1 In1DataTypes;
    typedef typename In1::Coord    In1Coord;
    typedef typename In1::Deriv    In1Deriv;
    typedef typename In1::VecCoord In1VecCoord;
    typedef typename In1::VecDeriv In1VecDeriv;

    typedef In2 In2DataTypes;
    typedef typename In2::Coord    In2Coord;
    typedef typename In2::Deriv    In2Deriv;
    typedef typename In2::VecCoord In2VecCoord;
    typedef typename In2::VecDeriv In2VecDeriv;

    typedef Out OutDataTypes;
    typedef typename Out::Coord   OutCoord;
    typedef typename Out::Deriv   OutDeriv;
    typedef typename Out::VecCoord OutVecCoord;
    typedef typename Out::VecDeriv OutVecDeriv;
    typedef typename OutCoord::value_type Real;

    typedef typename helper::vector<OutVecCoord*> vecOutVecCoord;
    typedef typename helper::vector<const In1VecCoord*> vecConstIn1VecCoord;
    typedef typename helper::vector<const In2VecCoord*> vecConstIn2VecCoord;

    virtual void apply(const vecOutVecCoord& outPos, const vecConstIn1VecCoord& inPos1 , const vecConstIn2VecCoord& inPos2 );
    virtual void applyJ(const helper::vector< OutVecDeriv*>& outDeriv, const helper::vector<const In1VecDeriv*>& inDeriv1, const helper::vector<const In2VecDeriv*>& inDeriv2);
    virtual void applyJT( const helper::vector<In1VecDeriv*>& outDeriv1 ,const helper::vector<In2VecDeriv*>& outDeriv2 , const helper::vector<const OutVecDeriv*>& inDeriv );
    virtual void applyDJT(const core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, core::MultiVecDerivId /*inForce*/, core::ConstMultiVecDerivId /*outForce*/) {}

    virtual void init();
    void draw(const core::visual::VisualParams* vparams);

protected:

    CenterOfMassMulti2Mapping()
    {
    }

    virtual ~CenterOfMassMulti2Mapping()
    {}

    helper::vector<const core::behavior::BaseMass*> inputBaseMass1;
    helper::vector<In1Coord> inputWeightedCOM1;
    helper::vector<In1Deriv> inputWeightedForce1;
    helper::vector<double> inputTotalMass1;

    helper::vector<const core::behavior::BaseMass*> inputBaseMass2;
    helper::vector<In2Coord> inputWeightedCOM2;
    helper::vector<In2Deriv> inputWeightedForce2;
    helper::vector<double> inputTotalMass2;

    double invTotalMass;
};

using namespace sofa::defaulttype;

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_MAPPING_CENTEROFMASSMULTI2MAPPING_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_MISC_MAPPING_API CenterOfMassMulti2Mapping< Vec3dTypes, Rigid3dTypes, Vec3dTypes >;
extern template class SOFA_MISC_MAPPING_API CenterOfMassMulti2Mapping< Vec3dTypes, Rigid3dTypes, Vec3dTypes >;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_MISC_MAPPING_API CenterOfMassMulti2Mapping< Vec3fTypes, Rigid3fTypes, Vec3fTypes >;
extern template class SOFA_MISC_MAPPING_API CenterOfMassMulti2Mapping< Vec3fTypes, Rigid3fTypes, Vec3fTypes >;
#endif
#endif

} // namespace mapping

} // namespace component

} // namespace sofa

#endif //SOFA_COMPONENT_MAPPING_CENTEROFMASSMULTI2MAPPING_H
