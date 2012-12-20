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
#ifndef SOFA_COMPONENT_VISUALMODEL_OGLVIEWPORT_H_
#define SOFA_COMPONENT_VISUALMODEL_OGLVIEWPORT_H_

#include <sofa/core/visual/VisualManager.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/BoundingBox.h>

#include <sofa/helper/gl/FrameBufferObject.h>
#include <sofa/core/visual/VisualParams.h>

namespace sofa
{

namespace component
{

namespace visualmodel
{

class OglViewport : public core::visual::VisualManager
{
public:
    typedef defaulttype::RigidCoord<3,double> RigidCoord;
    typedef core::visual::VisualParams::Viewport Viewport;

    SOFA_CLASS(OglViewport, core::visual::VisualManager);

    Data<defaulttype::Vec<2, int> > p_screenPosition;
    Data<defaulttype::Vec<2, unsigned int> > p_screenSize;
    Data<defaulttype::Vec3f> p_cameraPosition;
    Data<defaulttype::Quat> p_cameraOrientation;
    Data<RigidCoord > p_cameraRigid;
    Data<double> p_zNear, p_zFar;
    Data<double> p_fovy;
    Data<bool> p_useFBO;
    Data<bool> p_swapMainView;
    Data<bool> p_drawCamera;

    helper::gl::FrameBufferObject fbo;

protected:
    OglViewport();
    virtual ~OglViewport();
public:
    void init();
    void draw(const core::visual::VisualParams* vparams);
    void initVisual();
    void preDrawScene(core::visual::VisualParams* vp);
    bool drawScene(core::visual::VisualParams* vp);
    void postDrawScene(core::visual::VisualParams* vp);

protected:
    void renderToViewport(core::visual::VisualParams* vp);
    void renderFBOToScreen(core::visual::VisualParams* vp);

};
}

}

}

#endif /* SOFA_COMPONENT_VISUALMODEL_OGLVIEWPORT_H_ */
