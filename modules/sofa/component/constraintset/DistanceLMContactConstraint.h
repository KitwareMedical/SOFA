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
#ifndef SOFA_COMPONENT_CONSTRAINT_DistanceLMContactConstraint_H
#define SOFA_COMPONENT_CONSTRAINT_DistanceLMContactConstraint_H


#include <sofa/core/VecId.h>
#include <sofa/core/ConstraintParams.h>
#include <sofa/core/behavior/BaseMass.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/behavior/LMConstraint.h>
#include <sofa/component/constraintset/ContactDescription.h>
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
class DistanceLMContactConstraintInternalData
{
};




/** Keep two particules at an initial distance
*/
template <class DataTypes>
class DistanceLMContactConstraint :  public core::behavior::LMConstraint<DataTypes,DataTypes>, public ContactDescriptionHandler
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(DistanceLMContactConstraint,DataTypes),SOFA_TEMPLATE2(sofa::core::behavior::LMConstraint, DataTypes, DataTypes));

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::MatrixDeriv MatrixDeriv;
    typedef typename DataTypes::MatrixDeriv::RowIterator MatrixDerivRowIterator;

    typedef typename core::behavior::MechanicalState<DataTypes> MechanicalState;
    typedef typename sofa::core::topology::BaseMeshTopology::SeqEdges SeqEdges;
    typedef typename sofa::core::topology::BaseMeshTopology::Edge Edge;
    typedef core::ConstraintParams::ConstOrder ConstOrder;
    typedef core::behavior::ConstraintGroup ConstraintGroup;


protected:
    DistanceLMContactConstraint( MechanicalState *dof)
        : core::behavior::LMConstraint<DataTypes,DataTypes>(dof,dof)
        , pointPairs(Base::initData(&pointPairs, "pointPairs", "List of the edges to constrain"))
        , contactFriction(Base::initData(&contactFriction, "contactFriction", "Coulomb friction coefficient (same for all)"))
        , intersection(0)
    {
        initColorContactState();
    };

    DistanceLMContactConstraint( MechanicalState *dof1, MechanicalState * dof2)
        : core::behavior::LMConstraint<DataTypes,DataTypes>(dof1,dof2)
        , pointPairs(Base::initData(&pointPairs, "pointPairs", "List of the edges to constrain"))
        , contactFriction(Base::initData(&contactFriction, "contactFriction", "Coulomb friction coefficient (same for all)"))
        , intersection(0)
    {
        initColorContactState();
    };

    DistanceLMContactConstraint()
        : pointPairs(Base::initData(&pointPairs, "pointPairs", "List of the edges to constrain"))
        , contactFriction(Base::initData(&contactFriction, "contactFriction", "Coulomb friction coefficient (same for all)"))
        , intersection(0)
    {
        initColorContactState();
    }

    ~DistanceLMContactConstraint() {};
public:
    // -- LMConstraint interface
    void buildConstraintMatrix(const core::ConstraintParams* cParams /* PARAMS FIRST */, core::MultiMatrixDerivId cId, unsigned int &cIndex);

    void writeConstraintEquations(unsigned int& lineNumber, core::MultiVecId id, ConstOrder order);
    void LagrangeMultiplierEvaluation(const SReal* Wptr, const SReal* cptr, SReal* LambdaInitptr,
            core::behavior::ConstraintGroup * group);

    bool isCorrectionComputedWithSimulatedDOF(ConstOrder order) const;
    //
    void clear();
    /// register a new contact
    void addContact(unsigned m1, unsigned m2);
    virtual void draw(const core::visual::VisualParams* vparams);
    bool useMask() const {return true;}

    std::string getTemplateName() const
    {
        return templateName(this);
    }
    static std::string templateName(const DistanceLMContactConstraint<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }


protected :

    /// Each scalar constraint (up to three per contact) has an associated index
    helper::vector<  unsigned int > scalarConstraintsIndices;

public:
    /// Contacts are represented by pairs of point indices
    Data< SeqEdges > pointPairs;

    /// Friction coefficients (same for all contacts)
    Data< SReal > contactFriction;

protected:
    ///Compute the length of an edge given the vector of coordinates corresponding
    double lengthEdge(const Edge &e, const VecCoord &x1,const VecCoord &x2) const;
    /// Contact normal
    Deriv computeNormal(const Edge &e, const VecCoord &x1, const VecCoord &x2) const;
    /// Contact tangent vectors
    void computeTangentVectors( Deriv& T1, Deriv& T2, const Deriv& N );

    struct Contact
    {
        //Constrained Axis
        Deriv n,t1,t2;
        Contact() {}
        Contact( Deriv norm, Deriv tgt1, Deriv tgt2 ):n(norm),t1(tgt1),t2(tgt2),contactForce(Deriv()) {}
        Deriv contactForce;
    };

    std::map< Edge, Contact > edgeToContact;
    std::map< ConstraintGroup*, Contact* > constraintGroupToContact;
    core::collision::Intersection* intersection;
protected:
    DistanceLMContactConstraintInternalData<DataTypes> data;
    friend class DistanceLMContactConstraintInternalData<DataTypes>;


    void initColorContactState()
    {
        colorsContactState.clear();
        //Vanishing
        colorsContactState.push_back(defaulttype::Vec<4,float>(0.0f,1.0f,0.0f,1.0));
        //Sticking
        colorsContactState.push_back(defaulttype::Vec<4,float>(1.0f,0.0f,0.0f,1.0));
        //Sliding
        colorsContactState.push_back(defaulttype::Vec<4,float>(1.0f,1.0f,0.0f,1.0));
        //Sliding Direction
        colorsContactState.push_back(defaulttype::Vec<4,float>(1.0f,0.0f,1.0f,1.0));
    }

    helper::vector< defaulttype::Vec<4,float> > colorsContactState;
};


#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_CONSTRAINTSET_DistanceLMContactConstraint_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_CONSTRAINT_API DistanceLMContactConstraint<defaulttype::Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_CONSTRAINT_API DistanceLMContactConstraint<defaulttype::Vec3fTypes>;
#endif
#endif

} // namespace constraintset

} // namespace component

} // namespace sofa

#endif
