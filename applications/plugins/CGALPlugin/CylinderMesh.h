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
*                               SOFA :: Plugins                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
/*
 * CylinderMesh.h
 *
 *  Created on: 21 mar. 2010
 *      Author: Yiyi
 */

#ifndef CGALPLUGIN_CYLINDERMESH_H
#define CGALPLUGIN_CYLINDERMESH_H

#include <sofa/defaulttype/Vec.h>
#include <sofa/core/DataEngine.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/helper/gl/template.h>

#include <math.h>
#include   <algorithm>

namespace cgal
{

template <class DataTypes>
class CylinderMesh : public sofa::core::DataEngine
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(CylinderMesh,DataTypes),sofa::core::DataEngine);

    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Point;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef sofa::helper::fixed_array<int, 3> Index;
    //    typedef sofa::helper::vector<Real> VecReal;

//        typedef sofa::core::topology::BaseMeshTopology::PointID PointID;
    //    typedef sofa::core::topology::BaseMeshTopology::Edge Edge;
    //    typedef sofa::core::topology::BaseMeshTopology::Triangle Triangle;
    //    typedef sofa::core::topology::BaseMeshTopology::Quad Quad;
    typedef sofa::core::topology::BaseMeshTopology::Tetra Tetra;
    typedef sofa::core::topology::BaseMeshTopology::Hexa Hexa;


    //    typedef sofa::core::topology::BaseMeshTopology::SeqEdges SeqEdges;
    //    typedef sofa::core::topology::BaseMeshTopology::SeqTriangles SeqTriangles;
    //    typedef sofa::core::topology::BaseMeshTopology::SeqQuads SeqQuads;
    typedef sofa::core::topology::BaseMeshTopology::SeqTetrahedra SeqTetrahedra;
    typedef sofa::core::topology::BaseMeshTopology::SeqHexahedra SeqHexahedra;


public:
    CylinderMesh();
    virtual ~CylinderMesh() { };

    void init();
    void reinit();

    void update();
    void scale();
    void orientate();
    void draw();

    virtual std::string getTemplateName() const
    {
        return templateName(this);
    }

    static std::string templateName(const CylinderMesh<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }

    //Inputs
    Data<double> m_diameter;
    Data<double> m_length;
    Data<int> m_number;
    Data<bool> m_bScale;
    Data<bool> m_viewPoints;
    Data<bool> m_viewTetras;

    //Outputs
    Data<VecCoord> m_points;
    Data<SeqTetrahedra> m_tetras;

    //Parameters
    Real m_interval;
    int m_nbVertices, m_nbCenters, m_nbBDCenters, m_nbTetras;
    int n, m, a;
    Real d, l, t;
    std::map<Index, int> m_ptID;


};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(CGALPLUGIN_CYLINDERMESH_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_CGALPLUGIN_API CylinderMesh<defaulttype::Vec3dTypes>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
extern template class SOFA_CGALPLUGIN_API CylinderMesh<defaulttype::Vec3fTypes>;
#endif //SOFA_DOUBLE
#endif

} //cgal

#endif /* CGALPLUGIN_CYLINDERMESH_H */
