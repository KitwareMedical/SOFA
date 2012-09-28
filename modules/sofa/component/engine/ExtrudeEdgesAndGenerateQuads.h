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
#ifndef SOFA_COMPONENT_ENGINE_EXTRUDEEDGESANDGENERATEQUADS_H
#define SOFA_COMPONENT_ENGINE_EXTRUDEEDGESANDGENERATEQUADS_H

#if !defined(__GNUC__) || (__GNUC__ > 3 || (_GNUC__ == 3 && __GNUC_MINOR__ > 3))
#pragma once
#endif

#include <sofa/defaulttype/Vec.h>
#include <sofa/core/DataEngine.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/component/component.h>
#include <sofa/defaulttype/Vec3Types.h>

namespace sofa
{

namespace component
{

namespace engine
{

using namespace core::behavior;
using namespace core::topology;
using namespace core::objectmodel;

/**
 * This class extrudes a quad surface into a set of hexahedra
 */
template <class DataTypes>
class ExtrudeEdgesAndGenerateQuads : public core::DataEngine
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(ExtrudeEdgesAndGenerateQuads,DataTypes),core::DataEngine);

    typedef typename DataTypes::Coord     Coord;
    typedef typename DataTypes::VecCoord  VecCoord;
    typedef typename DataTypes::Real      Real;
    typedef defaulttype::Vec<3,Real>      Vec3;

protected:

    ExtrudeEdgesAndGenerateQuads();

    ~ExtrudeEdgesAndGenerateQuads() {}
public:
    void init();

    void reinit();

    void update();

    void draw();

    virtual std::string getTemplateName() const
    {
        return templateName(this);
    }

    static std::string templateName(const ExtrudeEdgesAndGenerateQuads<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }

    bool                                             initialized;
    Data<bool>                                       isVisible;
    Data<Coord>                                      f_direction;
    Data<Real>                                       f_thickness;
    Data<Real>                                       f_thicknessIn;
    Data<Real>                                       f_thicknessOut;
    Data<int>                                        f_numberOfSections;
    Data<VecCoord>                                   f_curveVertices;
    Data< helper::vector<BaseMeshTopology::Edge> >   f_curveEdges;
    Data<VecCoord>                                   f_extrudedVertices;
    Data< helper::vector<BaseMeshTopology::Edge> >   f_extrudedEdges;
    Data< helper::vector<BaseMeshTopology::Quad> >   f_extrudedQuads;
};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_ENGINE_EXTRUDEEDGESANDGENERATEQUADS_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_ENGINE_API ExtrudeEdgesAndGenerateQuads<defaulttype::Vec3dTypes>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
extern template class SOFA_ENGINE_API ExtrudeEdgesAndGenerateQuads<defaulttype::Vec3fTypes>;
#endif //SOFA_DOUBLE
#endif

} // namespace engine

} // namespace component

} // namespace sofa

#endif
