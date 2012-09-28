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
#ifndef SOFA_COMPONENT_PROJECTIVECONSTRAINTSET_ATTACHCONSTRAINT_H
#define SOFA_COMPONENT_PROJECTIVECONSTRAINTSET_ATTACHCONSTRAINT_H

#include <sofa/core/behavior/PairInteractionProjectiveConstraintSet.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/defaulttype/BaseMatrix.h>
#include <sofa/defaulttype/BaseVector.h>
#include <sofa/defaulttype/Vec3Types.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/helper/vector.h>
#include <sofa/component/topology/TopologySubsetData.h>
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

/// This class can be overridden if needed for additionnal storage within template specializations.
template <class DataTypes>
class AttachConstraintInternalData
{
};

/** Attach given pair of particles, projecting the positions of the second particles to the first ones.
*/
template <class DataTypes>
class AttachConstraint : public core::behavior::PairInteractionProjectiveConstraintSet<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(AttachConstraint,DataTypes),SOFA_TEMPLATE(sofa::core::behavior::PairInteractionProjectiveConstraintSet,DataTypes));

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::Real Real;

    typedef core::objectmodel::Data<VecCoord> DataVecCoord;
    typedef core::objectmodel::Data<VecDeriv> DataVecDeriv;

    typedef helper::vector<unsigned int> SetIndexArray;
    typedef sofa::component::topology::PointSubsetData< SetIndexArray > SetIndex;


protected:
    AttachConstraintInternalData<DataTypes> data;

    /// Pointer to the current topology
    sofa::core::topology::BaseMeshTopology* topology;

public:
    SetIndex f_indices1;
    SetIndex f_indices2;
    Data<Real> f_radius;
    Data<bool> f_twoWay;
    Data<bool> f_freeRotations;
    Data<bool> f_lastFreeRotation;
    Data<bool> f_restRotations;
    Data<defaulttype::Vector3> f_lastPos;
    Data<defaulttype::Vector3> f_lastDir;
    Data<bool> f_clamp;
    Data<Real> f_minDistance;

    helper::vector<bool> activeFlags;
    helper::vector<bool> constraintReleased;
    helper::vector<Real> lastDist;
    helper::vector<defaulttype::Quat> restRotations;
protected:
    AttachConstraint(core::behavior::MechanicalState<DataTypes> *mm1, core::behavior::MechanicalState<DataTypes> *mm2);
    AttachConstraint();
    virtual ~AttachConstraint();
public:
    void clearConstraints();
    void addConstraint(unsigned int index1, unsigned int index2);

    // -- Constraint interface
    void init();
    void projectResponse(const core::MechanicalParams *mparams /* PARAMS FIRST */, DataVecDeriv& dx1, DataVecDeriv& dx2);
    void projectVelocity(const core::MechanicalParams *mparams /* PARAMS FIRST */, DataVecDeriv& v1, DataVecDeriv& v2);
    void projectPosition(const core::MechanicalParams *mparams /* PARAMS FIRST */, DataVecCoord& x1, DataVecCoord& x2);

    /// Project the global Mechanical Matrix to constrained space using offset parameter
    void applyConstraint(const core::MechanicalParams *mparams /* PARAMS FIRST */, const sofa::core::behavior::MultiMatrixAccessor* matrix);

    /// Project the global Mechanical Vector to constrained space using offset parameter
    void applyConstraint(const core::MechanicalParams *mparams /* PARAMS FIRST */, defaulttype::BaseVector* vector, const sofa::core::behavior::MultiMatrixAccessor* matrix);


    virtual void draw(const core::visual::VisualParams* vparams);

protected :

    void projectPosition(Coord& x1, Coord& x2, bool /*freeRotations*/, unsigned index)
    {
        // do nothing if distance between x2 & x1 is bigger than f_minDistance
        if (f_minDistance.getValue() != -1 &&
            (x2 - x1).norm() > f_minDistance.getValue())
        {
            constraintReleased[index] = true;
            return;
        }
        constraintReleased[index] = false;

        x2 = x1;
    }

    void projectVelocity(Deriv& x1, Deriv& x2, bool /*freeRotations*/, unsigned index)
    {
        // do nothing if distance between x2 & x1 is bigger than f_minDistance
        if (constraintReleased[index]) return;

        x2 = x1;
    }

    void projectResponse(Deriv& dx1, Deriv& dx2, bool /*freeRotations*/, bool twoway, unsigned index)
    {
        // do nothing if distance between x2 & x1 is bigger than f_minDistance
        if (constraintReleased[index]) return;

        if (!twoway)
        {
            dx2 = Deriv();
        }
        else
        {
            dx1 += dx2;
            dx2 = dx1;
        }
    }

    static unsigned int DerivConstrainedSize(bool /*freeRotations*/) { return Deriv::size(); }

    void calcRestRotations();
};


#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_PROJECTIVECONSTRAINTSET_ATTACHCONSTRAINT_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_OBJECT_INTERACTION_API AttachConstraint<defaulttype::Vec3dTypes>;
extern template class SOFA_OBJECT_INTERACTION_API AttachConstraint<defaulttype::Vec2dTypes>;
extern template class SOFA_OBJECT_INTERACTION_API AttachConstraint<defaulttype::Vec1dTypes>;
extern template class SOFA_OBJECT_INTERACTION_API AttachConstraint<defaulttype::Rigid3dTypes>;
extern template class SOFA_OBJECT_INTERACTION_API AttachConstraint<defaulttype::Rigid2dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_OBJECT_INTERACTION_API AttachConstraint<defaulttype::Vec3fTypes>;
extern template class SOFA_OBJECT_INTERACTION_API AttachConstraint<defaulttype::Vec2fTypes>;
extern template class SOFA_OBJECT_INTERACTION_API AttachConstraint<defaulttype::Vec1fTypes>;
extern template class SOFA_OBJECT_INTERACTION_API AttachConstraint<defaulttype::Rigid3fTypes>;
extern template class SOFA_OBJECT_INTERACTION_API AttachConstraint<defaulttype::Rigid2fTypes>;
#endif
#endif

} // namespace projectiveconstraintset

} // namespace component

} // namespace sofa


#endif
