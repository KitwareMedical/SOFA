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
//
// C++ Interface: EdgeSetController
//
// Description:
//
//
// Author: Pierre-Jean Bensoussan, Digital Trainers (2008)
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef SOFA_COMPONENT_CONTROLLER_EDGESETCONTROLLER_H
#define SOFA_COMPONENT_CONTROLLER_EDGESETCONTROLLER_H

#include <sofa/component/controller/MechanicalStateController.h>
#include <sofa/component/topology/EdgeSetTopologyModifier.h>


namespace sofa
{

namespace component
{
namespace topology
{
template <class T>
class EdgeSetGeometryAlgorithms;

class EdgeSetTopologyModifier;
}

namespace controller
{

/**
 * @brief EdgeSetController Class
 *
 * Provides a Mouse & Keyboard user control on an EdgeSet Topology.
 */
template<class DataTypes>
class EdgeSetController : public MechanicalStateController<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(EdgeSetController,DataTypes),SOFA_TEMPLATE(MechanicalStateController,DataTypes));
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord    Coord   ;
    typedef typename DataTypes::Deriv    Deriv   ;
    typedef typename Coord::value_type   Real    ;

    typedef MechanicalStateController<DataTypes> Inherit;
protected:
    /**
     * @brief Default Constructor.
     */
    EdgeSetController();

    /**
     * @brief Default Destructor.
     */
    virtual ~EdgeSetController() {};
public:
    /**
     * @brief SceneGraph callback initialization method.
     */
    virtual void init();

    /**
     * @name Controller Interface
     */
    //@{

    /**
     * @brief Mouse event callback.
     */
    virtual void onMouseEvent(core::objectmodel::MouseEvent *);

    /**
     * @brief Keyboard key pressed event callback.
     */
    virtual void onKeyPressedEvent(core::objectmodel::KeypressedEvent *);


    /**
     * @brief Begin Animation event callback.
     */
    virtual void onBeginAnimationStep(const double dt);

    //@}

    /**
     * @name Accessors
     */
    //@{

    virtual std::string getTemplateName() const
    {
        return templateName(this);
    }

    static std::string templateName(const EdgeSetController<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }

    //@}

    /**
     * @brief Apply the controller modifications to the controlled MechanicalState.
     */
    virtual void applyController(void);

    /**
     * @brief
     */
    virtual bool modifyTopology(void);

    /**
     * @brief
     */
    virtual void draw(const core::visual::VisualParams* vparams);

protected:
    Data<Real> step;
    Data<Real> minLength; ///< determine the minimum length of the edge set
    Data<Real> maxLength; ///< determine the maximum length of the edge set
    Data<Real> maxDepl; ///< determine the maximum deplacement in a time step
    Data<Real> speed;
    Data<bool> reversed;
    Data<int>  startingIndex;   ///< index of the edge where a topological change occurs
    Real depl;

    sofa::core::topology::BaseMeshTopology* _topology;
    sofa::component::topology::EdgeSetGeometryAlgorithms<DataTypes>* edgeGeo;
    sofa::component::topology::EdgeSetTopologyModifier* edgeMod;
    Coord refPos;
    helper::vector<Real> vertexT;

    virtual void computeVertexT();

    virtual Coord getNewRestPos(const Coord& pos, Real /*t*/, Real dt)
    {
        sofa::defaulttype::Vec<3,Real> vectrans(dt * this->mainDirection.getValue()[0], dt * this->mainDirection.getValue()[1], dt * this->mainDirection.getValue()[2]);
        vectrans = pos.getOrientation().rotate(vectrans);
        return Coord(pos.getCenter() + vectrans, pos.getOrientation());
    }

    Real edgeTLength;
};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_CONTROLLER_EDGESETCONTROLLER_CPP)
#ifndef SOFA_FLOAT
//extern template class SOFA_USER_INTERACTION_API EdgeSetController<defaulttype::Vec3dTypes>;
//extern template class SOFA_USER_INTERACTION_API EdgeSetController<defaulttype::Vec2dTypes>;
//extern template class SOFA_USER_INTERACTION_API EdgeSetController<defaulttype::Vec1dTypes>;
//extern template class SOFA_USER_INTERACTION_API EdgeSetController<defaulttype::Vec6dTypes>;
extern template class SOFA_USER_INTERACTION_API EdgeSetController<defaulttype::Rigid3dTypes>;
//extern template class SOFA_USER_INTERACTION_API EdgeSetController<defaulttype::Rigid2dTypes>;
#endif
#ifndef SOFA_DOUBLE
//extern template class SOFA_USER_INTERACTION_API EdgeSetController<defaulttype::Vec3fTypes>;
//extern template class SOFA_USER_INTERACTION_API EdgeSetController<defaulttype::Vec2fTypes>;
//extern template class SOFA_USER_INTERACTION_API EdgeSetController<defaulttype::Vec1fTypes>;
//extern template class SOFA_USER_INTERACTION_API EdgeSetController<defaulttype::Vec6fTypes>;
extern template class SOFA_USER_INTERACTION_API EdgeSetController<defaulttype::Rigid3fTypes>;
//extern template class SOFA_USER_INTERACTION_API EdgeSetController<defaulttype::Rigid2fTypes>;
#endif
#endif

} // namespace controller

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_CONTROLLER_EDGESETCONTROLLER_H
