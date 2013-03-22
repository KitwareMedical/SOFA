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
#ifndef SOFA_COMPONENT_PROJECTIVECONSTRAINTSET_FIXEDTRANSLATIONCONSTRAINT_INL
#define SOFA_COMPONENT_PROJECTIVECONSTRAINTSET_FIXEDTRANSLATIONCONSTRAINT_INL

#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/component/projectiveconstraintset/FixedTranslationConstraint.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/gl/template.h>
#include <sofa/core/behavior/ProjectiveConstraintSet.inl>
#include <sofa/component/topology/TopologySubsetData.inl>

namespace sofa
{

namespace component
{

namespace projectiveconstraintset
{

using namespace sofa::core::objectmodel;
using namespace sofa::core::behavior;
using namespace sofa::helper;

using namespace sofa::core::topology;

// Define TestNewPointFunction
template< class DataTypes>
bool FixedTranslationConstraint<DataTypes>::FCPointHandler::applyTestCreateFunction(unsigned int, const sofa::helper::vector<unsigned int> &, const sofa::helper::vector<double> &)
{
    return fc != 0;
}

// Define RemovalFunction
template< class DataTypes>
void FixedTranslationConstraint<DataTypes>::FCPointHandler::applyDestroyFunction(unsigned int pointIndex, value_type &)
{
    if (fc)
    {
        fc->removeIndex((unsigned int) pointIndex);
    }
}

template< class DataTypes>
FixedTranslationConstraint<DataTypes>::FixedTranslationConstraint()
    : core::behavior::ProjectiveConstraintSet<DataTypes>(NULL)
    , f_indices( BaseObject::initData(&f_indices,"indices","Indices of the fixed points") )
    , f_fixAll( BaseObject::initData(&f_fixAll,false,"fixAll","filter all the DOF to implement a fixed object") )
    , _drawSize( BaseObject::initData(&_drawSize,0.0,"drawSize","0 -> point based rendering, >0 -> radius of spheres") )
    , f_coordinates( BaseObject::initData(&f_coordinates,"coordinates","Coordinates of the fixed points") )
{
    // default to indice 0
    f_indices.beginEdit()->push_back(0);
    f_indices.endEdit();

    pointHandler = new FCPointHandler(this, &f_indices);
}


template <class DataTypes>
FixedTranslationConstraint<DataTypes>::~FixedTranslationConstraint()
{
    if (pointHandler)
        delete pointHandler;
}

template <class DataTypes>
void FixedTranslationConstraint<DataTypes>::clearIndices()
{
    f_indices.beginEdit()->clear();
    f_indices.endEdit();
}

template <class DataTypes>
void FixedTranslationConstraint<DataTypes>::addIndex(unsigned int index)
{
    f_indices.beginEdit()->push_back(index);
    f_indices.endEdit();
}

template <class DataTypes>
void FixedTranslationConstraint<DataTypes>::removeIndex(unsigned int index)
{
    removeValue(*f_indices.beginEdit(),index);
    f_indices.endEdit();
}

// -- Constraint interface
template <class DataTypes>
void FixedTranslationConstraint<DataTypes>::init()
{
    this->core::behavior::ProjectiveConstraintSet<DataTypes>::init();

    topology = this->getContext()->getMeshTopology();

    // Initialize functions and parameters
    f_indices.createTopologicalEngine(topology, pointHandler);
    f_indices.registerTopologicalData();

    f_coordinates.createTopologicalEngine(topology);
    f_coordinates.registerTopologicalData();

}


template<int N, class T>
static inline void clearPos(defaulttype::RigidDeriv<N,T>& v)
{
    getVCenter(v).clear();
}

template<class T>
static inline void clearPos(defaulttype::Vec<6,T>& v)
{
    for (unsigned int i=0; i<3; ++i)
        v[i] = 0;
}

template <class DataTypes> template <class DataDeriv>
void FixedTranslationConstraint<DataTypes>::projectResponseT(const core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataDeriv& res)
{
    const SetIndexArray & indices = f_indices.getValue();

    if (f_fixAll.getValue() == true)
    {
        for (int i = 0; i < topology->getNbPoints(); ++i)
        {
            clearPos(res[i]);
        }
    }
    else
    {
        for (SetIndexArray::const_iterator it = indices.begin(); it
                != indices.end(); ++it)
        {
            clearPos(res[*it]);
        }
    }
}

template <class DataTypes>
void FixedTranslationConstraint<DataTypes>::projectResponse(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& resData)
{
    helper::WriteAccessor<DataVecDeriv> res = resData;
    projectResponseT(mparams /* PARAMS FIRST */, res.wref());
}

template <class DataTypes>
void FixedTranslationConstraint<DataTypes>::projectVelocity(const core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecDeriv& /*vData*/)
{

}

template <class DataTypes>
void FixedTranslationConstraint<DataTypes>::projectPosition(const core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecCoord& /*xData*/)
{

}

template <class DataTypes>
void FixedTranslationConstraint<DataTypes>::projectJacobianMatrix(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataMatrixDeriv& cData)
{
    helper::WriteAccessor<DataMatrixDeriv> c = cData;

    MatrixDerivRowIterator rowIt = c->begin();
    MatrixDerivRowIterator rowItEnd = c->end();

    while (rowIt != rowItEnd)
    {
        projectResponseT<MatrixDerivRowType>(mparams /* PARAMS FIRST */, rowIt.row());
        ++rowIt;
    }
}


template <class DataTypes>
void FixedTranslationConstraint<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
#ifndef SOFA_NO_OPENGL
    const SetIndexArray & indices = f_indices.getValue();
    if (!vparams->displayFlags().getShowBehaviorModels())
        return;
    const VecCoord& x = *this->mstate->getX();
    glDisable(GL_LIGHTING);
    glPointSize(10);
    glColor4f(1, 0.5, 0.5, 1);
    glBegin(GL_POINTS);
    if (f_fixAll.getValue() == true)
    {
        for (unsigned i = 0; i < x.size(); i++)
        {
            gl::glVertexT(x[i].getCenter());
        }
    }
    else
    {
        for (SetIndex::const_iterator it = indices.begin(); it != indices.end(); ++it)
        {
            gl::glVertexT(x[*it].getCenter());
        }
    }
    glEnd();
#endif /* SOFA_NO_OPENGL */
}

template<>
void FixedTranslationConstraint<defaulttype::Vec6dTypes>::draw(const core::visual::VisualParams* vparams);
template<>
void FixedTranslationConstraint<defaulttype::Vec6fTypes>::draw(const core::visual::VisualParams* vparams);

} // namespace constraint

} // namespace component

} // namespace sofa

#endif