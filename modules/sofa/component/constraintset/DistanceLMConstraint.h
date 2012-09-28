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
#ifndef SOFA_COMPONENT_CONSTRAINTSET_DISTANCELMCONSTRAINT_H
#define SOFA_COMPONENT_CONSTRAINTSET_DISTANCELMCONSTRAINT_H

#include <sofa/core/VecId.h>
#include <sofa/core/ConstraintParams.h>
#include <sofa/core/behavior/BaseMass.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/behavior/LMConstraint.h>
#include <sofa/simulation/common/Node.h>

namespace sofa
{

namespace component
{

namespace constraintset
{

using helper::vector;
using core::objectmodel::Data;
using namespace sofa::core::objectmodel;

/// This class can be overridden if needed for additionnal storage within template specializations.
template <class DataTypes>
class DistanceLMConstraintInternalData
{
};




/** Keep two particules at an initial distance
*/
template <class DataTypes>
class DistanceLMConstraint :  public core::behavior::LMConstraint<DataTypes,DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(DistanceLMConstraint,DataTypes),SOFA_TEMPLATE2(sofa::core::behavior::LMConstraint, DataTypes, DataTypes));

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::MatrixDeriv MatrixDeriv;
    typedef typename DataTypes::MatrixDeriv::RowIterator MatrixDerivRowIterator;
    typedef typename core::behavior::MechanicalState<DataTypes> MechanicalState;

    typedef typename sofa::core::topology::BaseMeshTopology::SeqEdges SeqEdges;
    typedef typename sofa::core::topology::BaseMeshTopology::Edge Edge;

    typedef core::ConstraintParams::ConstOrder ConstOrder;

protected:
    DistanceLMConstraintInternalData<DataTypes> data;
    friend class DistanceLMConstraintInternalData<DataTypes>;


    DistanceLMConstraint( MechanicalState *dof)
        : core::behavior::LMConstraint<DataTypes,DataTypes>(dof,dof)
        , vecConstraint(Base::initData(&vecConstraint, "vecConstraint", "List of the edges to constrain"))
    {};

    DistanceLMConstraint( MechanicalState *dof1, MechanicalState * dof2)
        : core::behavior::LMConstraint<DataTypes,DataTypes>(dof1,dof2)
        , vecConstraint(Base::initData(&vecConstraint, "vecConstraint", "List of the edges to constrain"))
    {};

    DistanceLMConstraint()
        : vecConstraint(Base::initData(&vecConstraint, "vecConstraint", "List of the edges to constrain"))
    {}

    ~DistanceLMConstraint() {};
public:
    void init();
    void reinit();

    // -- LMConstraint interface

    void buildConstraintMatrix(const core::ConstraintParams* cParams /* PARAMS FIRST */, core::MultiMatrixDerivId cId, unsigned int &cIndex);
    void writeConstraintEquations(unsigned int& lineNumber, core::MultiVecId id, ConstOrder order);

    virtual void draw(const core::visual::VisualParams* vparams);

    bool isCorrectionComputedWithSimulatedDOF(core::ConstraintParams::ConstOrder /*order*/) const
    {
        simulation::Node* node1=(simulation::Node*) this->constrainedObject1->getContext();
        simulation::Node* node2=(simulation::Node*) this->constrainedObject2->getContext();
        if (node1->mechanicalMapping.empty() && node2->mechanicalMapping.empty()) return true;
        else return false;
    }

    bool useMask() const {return true;}

    std::string getTemplateName() const
    {
        return templateName(this);
    }
    static std::string templateName(const DistanceLMConstraint<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }

    //Edges involving a distance constraint
    Data< SeqEdges > vecConstraint;

protected :

    ///Compute the length of an edge given the vector of coordinates corresponding
    double lengthEdge(const Edge &e, const VecCoord &x1,const VecCoord &x2) const;
    ///Compute the direction of the constraint
    Deriv getDirection(const Edge &e, const VecCoord &x1, const VecCoord &x2) const;
    void updateRestLength();

    // Base Components of the current context
    core::topology::BaseMeshTopology *topology;

    helper::vector<  unsigned int > registeredConstraints;

    // rest length pre-computated
    sofa::helper::vector< double > l0;
};

#ifndef SOFA_FLOAT
template<>
defaulttype::Rigid3dTypes::Deriv DistanceLMConstraint<defaulttype::Rigid3dTypes>::getDirection(const Edge &e, const VecCoord &x1, const VecCoord &x2) const;
#endif
#ifndef SOFA_DOUBLE
template<>
defaulttype::Rigid3fTypes::Deriv DistanceLMConstraint<defaulttype::Rigid3fTypes>::getDirection(const Edge &e, const VecCoord &x1, const VecCoord &x2) const;
#endif

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_CONSTRAINTSET_DISTANCELMCONSTRAINT_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_CONSTRAINT_API DistanceLMConstraint<defaulttype::Vec3dTypes>;
extern template class SOFA_CONSTRAINT_API DistanceLMConstraint<defaulttype::Rigid3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_CONSTRAINT_API DistanceLMConstraint<defaulttype::Vec3fTypes>;
extern template class SOFA_CONSTRAINT_API DistanceLMConstraint<defaulttype::Rigid3fTypes>;
#endif
#endif

} // namespace constraintset

} // namespace component

} // namespace sofa

#endif
