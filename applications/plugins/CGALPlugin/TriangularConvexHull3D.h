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
 * TriangularConvexHull3D.h
 *
 *  Created on: 27 oct. 2009
 *      Author: froy
 */

#ifndef CGALPLUGIN_TRIANGULARCONVEXHULL3D_H
#define CGALPLUGIN_TRIANGULARCONVEXHULL3D_H

#include <sofa/defaulttype/Vec.h>
#include <sofa/core/DataEngine.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/objectmodel/Data.h>

namespace cgal
{

using namespace sofa::core::objectmodel;

template <class DataTypes>
class TriangularConvexHull3D : public sofa::core::DataEngine
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(TriangularConvexHull3D,DataTypes),sofa::core::DataEngine);

    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Point;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::VecCoord VecCoord;

    typedef sofa::core::topology::BaseMeshTopology::Triangle Triangle;
    typedef sofa::core::topology::BaseMeshTopology::Quad Quad;
    typedef sofa::core::topology::BaseMeshTopology::Tetra Tetra;

    typedef sofa::core::topology::BaseMeshTopology::SeqTriangles SeqTriangles;
    typedef sofa::core::topology::BaseMeshTopology::SeqQuads SeqQuads;
    typedef sofa::core::topology::BaseMeshTopology::SeqTetrahedra SeqTetrahedra;


public:
    TriangularConvexHull3D();
    virtual ~TriangularConvexHull3D() { };

    void init();
    void reinit();

    void update();

    virtual std::string getTemplateName() const
    {
        return templateName(this);
    }

    static std::string templateName(const TriangularConvexHull3D<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }

    //Inputs
    Data<VecCoord> f_X0;

    //Outputs
    Data<VecCoord> f_newX0;
    Data<SeqTriangles> f_triangles;
};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(CGALPLUGIN_TRIANGULARCONVEXHULL3D_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_CGALPLUGIN_API TriangularConvexHull3D<defaulttype::Vec3dTypes>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
extern template class SOFA_CGALPLUGIN_API TriangularConvexHull3D<defaulttype::Vec3fTypes>;
#endif //SOFA_DOUBLE
#endif

} //cgal

#endif /* CGALPLUGIN_TriangularConvexHull3D_H */
