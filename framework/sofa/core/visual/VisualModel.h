/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *
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
*                              SOFA :: Framework                              *
*                                                                             *
* Authors: The SOFA Team (see Authors.txt)                                    *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_CORE_VISUAL_VISUALMODEL_H
#define SOFA_CORE_VISUAL_VISUALMODEL_H

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/Quat.h>

namespace sofa
{

namespace core
{

namespace visual
{

class VisualParams;

/**
 *  \brief An interface which all VisualModel inherit.
 *
 *  This Interface is used for the VisualModel, which all visible objects must
 *  implement.
 *
 *  VisualModels are drawn by calling their draw method. The method update is
 *  used to recompute some internal data (such as normals) after the simulation
 *  has computed a new timestep.
 *
 *  Most VisualModel are bound by a Mapping to a BehaviorModel or
 *  MechanicalState.
 */
class SOFA_CORE_API VisualModel : public virtual objectmodel::BaseObject
{
public:
    SOFA_CLASS(VisualModel, objectmodel::BaseObject);
protected:
    /// Destructor
    virtual ~VisualModel() { }
public:
    /**
     *  \brief Initialize the textures, or other graphical resources.
     *
     *  Called once before the first frame is drawn, and if the graphical
     *  context has been recreated.
     */
    virtual void initVisual() { initTextures(); }

    /**
     *  \brief clear some graphical resources (generaly called before the deleteVisitor).
     *  \note: for more general usage you can use the cleanup visitor
     */
    virtual void clearVisual() { }

    /**
     *  \brief Initialize the textures, or other graphical resources.
     *
     *  @deprecated Use initVisual() instead.
     */
    virtual void initTextures() {}

    /**
     *  \brief Called before objects in the current branch are displayed
     */
    virtual void fwdDraw(VisualParams* /*vparams*/) {}

    /**
     *  \brief Called after objects in the current branch are displayed
     */
    virtual void bwdDraw(VisualParams* /*vparams*/) {}

    /**
     *  \brief Display the VisualModel object.
     */
    virtual void drawVisual(const VisualParams* /*vparams*/) {}
    //virtual void drawVisual() = 0;

    /**
     *  \brief Display transparent surfaces.
     *
     *  Transparent objects should use this method to get a correct display order.
     */
    virtual void drawTransparent(const VisualParams* /*vparams*/)
    {
    }

    /**
     *  \brief Display shadow-casting surfaces.
     *
     *  This method default to calling draw(). Object that do not cast any
     *  shadows, or that use a different LOD for them should reimplement it.
     */
    virtual void drawShadow(const VisualParams* vparams)
    {
        drawVisual(vparams);
    }

    /**
     *  \brief used to update the model if necessary.
     */
    virtual void updateVisual() { update(); }
    /**
    *  \brief used to update the model if necessary.
    */
    virtual void parallelUpdateVisual() { }


    /**
     *  \brief used to update the model if necessary.
     *
     *  @deprecated Use updateVisual() instead.
     */
    virtual void update() {}

    /**
     *  \brief used to add the bounding-box of this visual model to the
     *  given bounding box in order to compute the scene bounding box or
     *  cull hidden objects.
     *
     *  \return false if the visual model does not define any bounding box,
     *  which should only be the case for "debug" objects, as this lack of
     *  information might affect performances and leads to incorrect scene
     *  bounding box.
     */
    virtual bool addBBox(double* /*minBBox*/, double* /*maxBBox*/)
    {
        return false;
    }

    /// Translate the positions
    ///
    /// This method is optional, it is used when the user want to interactively change the position of an object
    virtual void applyTranslation(const double /*dx*/, const double /*dy*/, const double /*dz*/)
    {
    }

    /// Rotate the positions using Euler Angles in degree
    ///
    /// This method is optional, it is used when the user want to interactively change the position of an object
    virtual void applyRotation (const double /*rx*/, const double /*ry*/, const double /*rz*/)
    {
    }

    /// Rotate the positions
    ///
    /// This method is optional, it is used when the user want to interactively change the position of an object
    virtual void applyRotation(const defaulttype::Quat /*q*/)
    {
    }

    /// Scale the positions
    ///
    /// This method is optional, it is used when the user want to interactively change the position of an object
    virtual void applyScale(const double /*sx*/,const double /*sy*/,const double /*sz*/)
    {
    }

    /// Give the new point and its ancestors at each topological change
    virtual void updatePointAncestors(helper::vector< unsigned int >* /*vecPoint*/, helper::vector< helper::vector< unsigned int > >* /*vecAncestors*/)
    {
    }
    
    /**
     *  \brief Append this mesh to an OBJ format stream.
     *
     *  The number of vertices position, normal, and texture coordinates already written is given as parameters.
     *  This method should update them.
     */
    virtual void exportOBJ(std::string /*name*/, std::ostream* /*out*/, std::ostream* /*mtl*/, int& /*vindex*/, int& /*nindex*/, int& /*tindex*/, int& /*count*/)
    {
    }
};

} // namespace visual

} // namespace core

} // namespace sofa

#endif //SOFA_CORE_VISUAL_VISUALMODEL_H
