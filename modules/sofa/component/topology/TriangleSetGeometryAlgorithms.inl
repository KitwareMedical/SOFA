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
#ifndef SOFA_COMPONENT_TOPOLOGY_TRIANGLESETGEOMETRYALGORITHMS_INL
#define SOFA_COMPONENT_TOPOLOGY_TRIANGLESETGEOMETRYALGORITHMS_INL

#include <sofa/component/topology/TriangleSetGeometryAlgorithms.h>
#include <sofa/component/topology/EdgeSetGeometryAlgorithms.inl>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/component/topology/CommonAlgorithms.h>

namespace sofa
{

namespace component
{

namespace topology
{
template<class DataTypes>
void TriangleSetGeometryAlgorithms< DataTypes >::init()
{
    EdgeSetGeometryAlgorithms<DataTypes>::init();

    this->reinit();
}

template<class DataTypes>
void TriangleSetGeometryAlgorithms< DataTypes >::reinit()
{
    EdgeSetGeometryAlgorithms<DataTypes>::reinit();

    if (p_recomputeTrianglesOrientation.getValue())
        this->reorderTrianglesOrientationFromNormals();
}


template< class DataTypes>
void TriangleSetGeometryAlgorithms< DataTypes >::computeTriangleAABB(const TriangleID i, Coord& minCoord, Coord& maxCoord) const
{
    const Triangle &t = this->m_topology->getTriangle(i);
    const typename DataTypes::VecCoord& p = *(this->object->getX());

    for(unsigned int i=0; i<3; ++i)
    {
        minCoord[i] = std::min(p[t[0]][i], std::min(p[t[1]][i], p[t[2]][i]));
        maxCoord[i] = std::max(p[t[0]][i], std::max(p[t[1]][i], p[t[2]][i]));
    }
}

template<class DataTypes>
typename DataTypes::Coord TriangleSetGeometryAlgorithms<DataTypes>::computeTriangleCenter(const TriangleID i) const
{
    const Triangle &t = this->m_topology->getTriangle(i);
    const typename DataTypes::VecCoord& p = *(this->object->getX());

    return (p[t[0]] + p[t[1]] + p[t[2]]) / (Real) 3.0;
}

template<class DataTypes>
typename DataTypes::Coord TriangleSetGeometryAlgorithms<DataTypes>::computeRestTriangleCenter(const TriangleID i) const
{
    const Triangle &t = this->m_topology->getTriangle(i);
    const typename DataTypes::VecCoord& p = *(this->object->getX0());

    return (p[t[0]] + p[t[1]] + p[t[2]]) / (Real) 3.0;
}

template<class DataTypes>
void TriangleSetGeometryAlgorithms<DataTypes>::computeTriangleCircumcenterBaryCoefs(Vec<3,Real> &baryCoord, const TriangleID i) const
{
    const Triangle &t = this->m_topology->getTriangle(i);
    const typename DataTypes::VecCoord& p = *(this->object->getX());
    Real a2, b2, c2; // square lengths of the 3 edges
    a2 = (p[t[1]]-p[t[0]]).norm2();
    b2 = (p[t[2]]-p[t[1]]).norm2();
    c2 = (p[t[0]]-p[t[2]]).norm2();

    Real n = a2*(-a2+b2+c2) + b2*(a2-b2+c2) + c2*(a2+b2-c2);

    baryCoord[2] = a2*(-a2+b2+c2) / n;
    baryCoord[0] = b2*(a2-b2+c2) / n;
    baryCoord[1] = c2*(a2+b2-c2) / n;

    // barycentric coordinates are defined as
    //baryCoord = sofa::defaulttype::Vec<3,Real>(a2*(-a2+b2+c2) / n, b2*(a2-b2+c2) / n, c2*(a2+b2-c2) / n);
}

template<class DataTypes>
typename DataTypes::Coord TriangleSetGeometryAlgorithms<DataTypes>::computeTriangleCircumcenter(const TriangleID i) const
{
    const Triangle &t = this->m_topology->getTriangle(i);
    const typename DataTypes::VecCoord& p = *(this->object->getX());

    sofa::defaulttype::Vec<3,Real> barycentricCoords;
    computeTriangleCircumcenterBaryCoefs(barycentricCoords, i);

    return (p[t[0]]*barycentricCoords[0] + p[t[1]]*barycentricCoords[1] + p[t[2]]*barycentricCoords[2]);
}

template< class DataTypes>
void TriangleSetGeometryAlgorithms< DataTypes >::getTriangleVertexCoordinates(const TriangleID i, Coord pnt[3]) const
{
    const Triangle &t = this->m_topology->getTriangle(i);
    const typename DataTypes::VecCoord& p = *(this->object->getX());

    for(unsigned int i=0; i<3; ++i)
    {
        pnt[i] = p[t[i]];
    }
}

template< class DataTypes>
void TriangleSetGeometryAlgorithms< DataTypes >::getRestTriangleVertexCoordinates(const TriangleID i, Coord pnt[3]) const
{
    const Triangle &t = this->m_topology->getTriangle(i);
    const typename DataTypes::VecCoord& p = *(this->object->getX0());

    for(unsigned int i=0; i<3; ++i)
    {
        pnt[i] = p[t[i]];
    }
}

template< class DataTypes>
typename DataTypes::Real TriangleSetGeometryAlgorithms< DataTypes >::computeTriangleArea( const TriangleID i) const
{
    const Triangle &t = this->m_topology->getTriangle(i);
    const typename DataTypes::VecCoord& p = *(this->object->getX());
    Real area = (Real)(areaProduct(p[t[1]]-p[t[0]], p[t[2]]-p[t[0]]) * 0.5);
    return area;
}

template< class DataTypes>
typename DataTypes::Real TriangleSetGeometryAlgorithms< DataTypes >::computeRestTriangleArea( const TriangleID i) const
{
    const Triangle &t = this->m_topology->getTriangle(i);
    const typename DataTypes::VecCoord& p = *(this->object->getX0());
    Real area = (Real) (areaProduct(p[t[1]]-p[t[0]],p[t[2]]-p[t[0]]) * 0.5);
    return area;
}

/// computes the edge length of all edges are store in the array interface
template<class DataTypes>
void TriangleSetGeometryAlgorithms<DataTypes>::computeTriangleArea( BasicArrayInterface<Real> &ai) const
{
    const sofa::helper::vector<Triangle> &ta = this->m_topology->getTriangles();
    const typename DataTypes::VecCoord& p = *(this->object->getX());

    for (unsigned int i=0; i<ta.size(); ++i)
    {
        const Triangle &t=ta[i];
        ai[i]=(Real)(areaProduct(p[t[1]]-p[t[0]],p[t[2]]-p[t[0]]) * 0.5);
    }
}

// Computes the point defined by 2 indices of vertex and 1 barycentric coordinate
template<class DataTypes>
sofa::defaulttype::Vec<3,double> TriangleSetGeometryAlgorithms< DataTypes >::computeBaryEdgePoint(unsigned int p0, unsigned int p1, double coord_p) const
{
    const typename DataTypes::VecCoord& vect_c = *(this->object->getX());

    sofa::defaulttype::Vec<3,double> c0; c0 = vect_c[p0];
    sofa::defaulttype::Vec<3,double> c1; c1 = vect_c[p1];
    return c0*(1-coord_p) + c1*coord_p;
}

template<class DataTypes>
sofa::defaulttype::Vec<3,double> TriangleSetGeometryAlgorithms< DataTypes >::computeBaryTrianglePoint(unsigned int p0, unsigned int p1, unsigned int p2, sofa::defaulttype::Vec<3,double>& coord_p) const
{
    const typename DataTypes::VecCoord& vect_c = *(this->object->getX());

    sofa::defaulttype::Vec<3,double> c0; c0 = vect_c[p0];
    sofa::defaulttype::Vec<3,double> c1; c1 = vect_c[p1];
    sofa::defaulttype::Vec<3,double> c2; c2 = vect_c[p2];
    return c0*coord_p[0] + c1*coord_p[1] + c2*coord_p[2];
}


// Computes the opposite point to ind_p
template<class DataTypes>
sofa::defaulttype::Vec<3,double> TriangleSetGeometryAlgorithms< DataTypes >::getOppositePoint(unsigned int ind_p,
        const Edge& indices,
        double coord_p) const
{
    const typename DataTypes::VecCoord& vect_c = *(this->object->getX());

    const typename DataTypes::Coord& c1 = vect_c[indices[0]];
    const typename DataTypes::Coord& c2 = vect_c[indices[1]];

    sofa::defaulttype::Vec<3,Real> p;

    if(ind_p == indices[0])
    {
        p[0]= (Real) c2[0];
        p[1]= (Real) c2[1];
        p[2]= (Real) c2[2];
    }
    else
    {
        if(ind_p == indices[1])
        {
            p[0]= (Real) c1[0];
            p[1]= (Real) c1[1];
            p[2]= (Real) c1[2];
        }
        else
        {
            p[0]= (Real) ((1.0-coord_p)*c1[0] + coord_p*c2[0]);
            p[1]= (Real) ((1.0-coord_p)*c1[1] + coord_p*c2[1]);
            p[2]= (Real) ((1.0-coord_p)*c1[2] + coord_p*c2[2]);
        }
    }

    return ((Vec<3,double>) p);
}

// Computes the normal vector of a triangle indexed by ind_t (not normed)
template<class DataTypes>
sofa::defaulttype::Vec<3,double> TriangleSetGeometryAlgorithms< DataTypes >::computeTriangleNormal(const TriangleID ind_t) const
{
    const Triangle &t = this->m_topology->getTriangle(ind_t);
    const typename DataTypes::VecCoord& vect_c = *(this->object->getX());

    const typename DataTypes::Coord& c0 = vect_c[t[0]];
    const typename DataTypes::Coord& c1 = vect_c[t[1]];
    const typename DataTypes::Coord& c2 = vect_c[t[2]];

    sofa::defaulttype::Vec<3,Real> p0;
    p0[0] = (Real) (c0[0]);
    p0[1] = (Real) (c0[1]);
    p0[2] = (Real) (c0[2]);
    sofa::defaulttype::Vec<3,Real> p1;
    p1[0] = (Real) (c1[0]);
    p1[1] = (Real) (c1[1]);
    p1[2] = (Real) (c1[2]);
    sofa::defaulttype::Vec<3,Real> p2;
    p2[0] = (Real) (c2[0]);
    p2[1] = (Real) (c2[1]);
    p2[2] = (Real) (c2[2]);

    sofa::defaulttype::Vec<3,Real> normal_t=(p1-p0).cross( p2-p0);

    return ((Vec<3,double>) normal_t);
}

// barycentric coefficients of point p in triangle (a,b,c) indexed by ind_t
template<class DataTypes>
sofa::helper::vector< double > TriangleSetGeometryAlgorithms< DataTypes >::computeTriangleBarycoefs(
    const TriangleID ind_t,
    const sofa::defaulttype::Vec<3,double> &p) const
{
    const Triangle &t=this->m_topology->getTriangle(ind_t);
    return compute3PointsBarycoefs(p, t[0], t[1], t[2]);
}

// barycentric coefficients of point p in triangle whose vertices are indexed by (ind_p1,ind_p2,ind_p3)
template<class DataTypes>
sofa::helper::vector< double > TriangleSetGeometryAlgorithms< DataTypes >::compute3PointsBarycoefs(
    const sofa::defaulttype::Vec<3,double> &p,
    unsigned int ind_p1,
    unsigned int ind_p2,
    unsigned int ind_p3,
    bool bRest) const
{
    const double ZERO = 1e-12;
    sofa::helper::vector< double > baryCoefs;

    const typename DataTypes::VecCoord& vect_c = (bRest ? *(this->object->getX0()) : *(this->object->getX()));

    const typename DataTypes::Coord& c0 = vect_c[ind_p1];
    const typename DataTypes::Coord& c1 = vect_c[ind_p2];
    const typename DataTypes::Coord& c2 = vect_c[ind_p3];

    sofa::defaulttype::Vec<3,Real> a;
    a[0] = (Real) (c0[0]);
    a[1] = (Real) (c0[1]);
    a[2] = (Real) (c0[2]);
    sofa::defaulttype::Vec<3,Real> b;
    b[0] = (Real) (c1[0]);
    b[1] = (Real) (c1[1]);
    b[2] = (Real) (c1[2]);
    sofa::defaulttype::Vec<3,Real> c;
    c[0] = (Real) (c2[0]);
    c[1] = (Real) (c2[1]);
    c[2] = (Real) (c2[2]);

    sofa::defaulttype::Vec<3,double> M = (Vec<3,double>) (b-a).cross(c-a);
    double norm2_M = M*(M);

    double coef_a, coef_b, coef_c;

    //if(norm2_M==0.0) // triangle (a,b,c) is flat
    if(norm2_M < ZERO) // triangle (a,b,c) is flat
    {
        coef_a = (double) (1.0/3.0);
        coef_b = (double) (1.0/3.0);
        coef_c = (double) (1.0 - (coef_a + coef_b));
    }
    else
    {
        sofa::defaulttype::Vec<3,Real> N =  M/norm2_M;

        coef_a = N*((b-p).cross(c-p));
        coef_b = N*((c-p).cross(a-p));
        coef_c = (double) (1.0 - (coef_a + coef_b)); //N*((a-p).cross(b-p));
    }

    baryCoefs.push_back(coef_a);
    baryCoefs.push_back(coef_b);
    baryCoefs.push_back(coef_c);

    return baryCoefs;
}

// Find the two closest points from two triangles (each of the point belonging to one triangle)
template<class DataTypes>
void TriangleSetGeometryAlgorithms< DataTypes >::computeClosestIndexPair(const TriangleID ind_ta, const TriangleID ind_tb,
        unsigned int &ind1, unsigned int &ind2) const
{
    const typename DataTypes::VecCoord& vect_c = *(this->object->getX());

    const Triangle &ta=this->m_topology->getTriangle(ind_ta);
    const Triangle &tb=this->m_topology->getTriangle(ind_tb);

    Real min_value=(Real) 0.0;
    bool is_init = false;

    for(unsigned int i=0; i<3; i++)
    {
        const typename DataTypes::Coord& ca=vect_c[ta[i]];
        sofa::defaulttype::Vec<3,Real> pa;
        pa[0] = (Real) (ca[0]);
        pa[1] = (Real) (ca[1]);
        pa[2] = (Real) (ca[2]);

        for(unsigned int j=0; j!=i && j<3; j++)
        {
            const typename DataTypes::Coord& cb=vect_c[tb[i]];
            sofa::defaulttype::Vec<3,Real> pb;
            pb[0] = (Real) (cb[0]);
            pb[1] = (Real) (cb[1]);
            pb[2] = (Real) (cb[2]);

            Real norm_v_normal = (pa-pb)*(pa-pb);
            if(!is_init)
            {
                min_value = norm_v_normal;
                ind1 = ta[i];
                ind2 = tb[j];
                is_init = true;
            }
            else
            {
                if(norm_v_normal<min_value)
                {
                    min_value = norm_v_normal;
                    ind1 = ta[i];
                    ind2 = tb[j];
                }
            }
        }
    }

    return;
}

// test if a point is inside the triangle indexed by ind_t
template<class DataTypes>
bool TriangleSetGeometryAlgorithms< DataTypes >::isPointInsideTriangle(const TriangleID ind_t,
        bool is_tested,
        const sofa::defaulttype::Vec<3,Real>& p,
        unsigned int &ind_t_test,
        bool bRest) const
{
    const double ZERO = -1e-12;
    const typename DataTypes::VecCoord& vect_c = bRest
        ? *(this->object->getX0())
        : *(this->object->getX());
    const Triangle &t=this->m_topology->getTriangle(ind_t);

    const typename DataTypes::Coord& c0=vect_c[t[0]];
    const typename DataTypes::Coord& c1=vect_c[t[1]];
    const typename DataTypes::Coord& c2=vect_c[t[2]];

    sofa::defaulttype::Vec<3,Real> ptest = p;

    sofa::defaulttype::Vec<3,Real> p0;
    p0[0] = (Real) (c0[0]);
    p0[1] = (Real) (c0[1]);
    p0[2] = (Real) (c0[2]);
    sofa::defaulttype::Vec<3,Real> p1;
    p1[0] = (Real) (c1[0]);
    p1[1] = (Real) (c1[1]);
    p1[2] = (Real) (c1[2]);
    sofa::defaulttype::Vec<3,Real> p2;
    p2[0] = (Real) (c2[0]);
    p2[1] = (Real) (c2[1]);
    p2[2] = (Real) (c2[2]);

    sofa::defaulttype::Vec<3,Real> v_normal = (p2-p0).cross(p1-p0);

    Real norm_v_normal = v_normal*(v_normal);
    if(norm_v_normal != 0.0)
    {
        sofa::defaulttype::Vec<3,Real> n_01 = (p1-p0).cross(v_normal);
        sofa::defaulttype::Vec<3,Real> n_12 = (p2-p1).cross(v_normal);
        sofa::defaulttype::Vec<3,Real> n_20 = (p0-p2).cross(v_normal);

        double v_01 = (double) ((ptest-p0)*(n_01));
        double v_12 = (double) ((ptest-p1)*(n_12));
        double v_20 = (double) ((ptest-p2)*(n_20));

        bool is_inside = (v_01 > ZERO) && (v_12 > ZERO) && (v_20 > ZERO);

        if(is_tested && (!is_inside))
        {
            sofa::helper::vector< unsigned int > shell;
            unsigned int ind_edge = 0;

            if(v_01 < 0.0)
            {
                if(v_12 < 0.0) /// vertex 1
                {
                    shell =(sofa::helper::vector< unsigned int >) (this->m_topology->getTrianglesAroundVertex(t[1]));
                }
                else
                {
                    if(v_20 < 0.0) /// vertex 0
                    {
                        shell =(sofa::helper::vector< unsigned int >) (this->m_topology->getTrianglesAroundVertex(t[0]));

                    }
                    else // v_01 < 0.0
                    {
                        ind_edge=this->m_topology->getEdgeIndex(t[0],t[1]);
                        shell =(sofa::helper::vector< unsigned int >) (this->m_topology->getTrianglesAroundEdge(ind_edge));
                    }
                }
            }
            else
            {
                if(v_12 < 0.0)
                {
                    if(v_20 < 0.0) /// vertex 2
                    {
                        shell =(sofa::helper::vector< unsigned int >) (this->m_topology->getTrianglesAroundVertex(t[2]));

                    }
                    else // v_12 < 0.0
                    {
                        ind_edge=this->m_topology->getEdgeIndex(t[1],t[2]);
                        shell =(sofa::helper::vector< unsigned int >) (this->m_topology->getTrianglesAroundEdge(ind_edge));
                    }
                }
                else // v_20 < 0.0
                {
                    ind_edge=this->m_topology->getEdgeIndex(t[2],t[0]);
                    shell =(sofa::helper::vector< unsigned int >) (this->m_topology->getTrianglesAroundEdge(ind_edge));
                }
            }

            unsigned int i =0;
            bool is_in_next_triangle=false;
            unsigned int ind_triangle=0;
            unsigned ind_t_false_init;
            unsigned int &ind_t_false = ind_t_false_init;

            if(shell.size()>1)
            {
                while(i < shell.size() && !is_in_next_triangle)
                {
                    ind_triangle=shell[i];

                    if(ind_triangle != ind_t)
                    {
                        is_in_next_triangle = isPointInTriangle(ind_triangle, false, p, ind_t_false);
                    }
                    i++;
                }

                if(is_in_next_triangle)
                {
                    ind_t_test=ind_triangle;
                    //sout << "correct to triangle indexed by " << ind_t_test << sendl;
                }
                else // not found
                {
                    //sout << "not found !!! " << sendl;
                    ind_t_test=ind_t;
                }
            }
            else
            {
                ind_t_test=ind_t;
            }
        }
        return is_inside;

    }
    else // triangle is flat
    {
        //sout << "INFO_print : triangle is flat" << sendl;
        return false;
    }
}

// test if a point is in the triangle indexed by ind_t
template<class DataTypes>
bool TriangleSetGeometryAlgorithms< DataTypes >::isPointInTriangle(const TriangleID ind_t,
        bool is_tested,
        const sofa::defaulttype::Vec<3,Real>& p,
        unsigned int &ind_t_test) const
{
    const double ZERO = 1e-12;
    const typename DataTypes::VecCoord& vect_c = *(this->object->getX());
    const Triangle &t=this->m_topology->getTriangle(ind_t);

    const typename DataTypes::Coord& c0=vect_c[t[0]];
    const typename DataTypes::Coord& c1=vect_c[t[1]];
    const typename DataTypes::Coord& c2=vect_c[t[2]];

    sofa::defaulttype::Vec<3,Real> ptest = p;

    sofa::defaulttype::Vec<3,Real> p0;
    p0[0] = (Real) (c0[0]);
    p0[1] = (Real) (c0[1]);
    p0[2] = (Real) (c0[2]);
    sofa::defaulttype::Vec<3,Real> p1;
    p1[0] = (Real) (c1[0]);
    p1[1] = (Real) (c1[1]);
    p1[2] = (Real) (c1[2]);
    sofa::defaulttype::Vec<3,Real> p2;
    p2[0] = (Real) (c2[0]);
    p2[1] = (Real) (c2[1]);
    p2[2] = (Real) (c2[2]);

    sofa::defaulttype::Vec<3,Real> v_normal = (p2-p0).cross(p1-p0);

    Real norm_v_normal = v_normal*(v_normal);
    //if(norm_v_normal != 0.0)
    if(norm_v_normal > ZERO)
    {
        sofa::defaulttype::Vec<3,Real> n_01 = (p1-p0).cross(v_normal);
        sofa::defaulttype::Vec<3,Real> n_12 = (p2-p1).cross(v_normal);
        sofa::defaulttype::Vec<3,Real> n_20 = (p0-p2).cross(v_normal);

        double v_01 = (double) ((ptest-p0)*(n_01));
        double v_12 = (double) ((ptest-p1)*(n_12));
        double v_20 = (double) ((ptest-p2)*(n_20));

        //bool is_inside = (v_01 > 0.0) && (v_12 > 0.0) && (v_20 > 0.0);
        bool is_inside = (v_01 > -ZERO) && (v_12 > -ZERO) && (v_20 >= -ZERO);

        if(is_tested && (!is_inside))
        {
            sofa::helper::vector< unsigned int > shell;
            unsigned int ind_edge = 0;

            //if(v_01 < 0.0)
            if(v_01 < -ZERO)
            {
                //if(v_12 < 0.0) /// vertex 1
                if(v_12 < -ZERO) /// vertex 1
                {
                    shell =(sofa::helper::vector< unsigned int >) (this->m_topology->getTrianglesAroundVertex(t[1]));
                }
                else
                {
                    //if(v_20 < 0.0) /// vertex 0
                    if(v_20 < -ZERO) /// vertex 0
                    {
                        shell =(sofa::helper::vector< unsigned int >) (this->m_topology->getTrianglesAroundVertex(t[0]));

                    }
                    else // v_01 < 0.0
                    {
                        ind_edge=this->m_topology->getEdgeIndex(t[0],t[1]);
                        shell =(sofa::helper::vector< unsigned int >) (this->m_topology->getTrianglesAroundEdge(ind_edge));
                    }
                }
            }
            else
            {
                //if(v_12 < 0.0)
                if(v_12 < -ZERO)
                {
                    //if(v_20 < 0.0) /// vertex 2
                    if(v_20 < -ZERO) /// vertex 2
                    {
                        shell =(sofa::helper::vector< unsigned int >) (this->m_topology->getTrianglesAroundVertex(t[2]));

                    }
                    else // v_12 < 0.0
                    {
                        ind_edge=this->m_topology->getEdgeIndex(t[1],t[2]);
                        shell =(sofa::helper::vector< unsigned int >) (this->m_topology->getTrianglesAroundEdge(ind_edge));
                    }
                }
                else // v_20 < 0.0
                {
                    ind_edge=this->m_topology->getEdgeIndex(t[2],t[0]);
                    shell =(sofa::helper::vector< unsigned int >) (this->m_topology->getTrianglesAroundEdge(ind_edge));
                }
            }

            unsigned int i =0;
            bool is_in_next_triangle=false;
            unsigned int ind_triangle=0;
            unsigned ind_t_false_init;
            unsigned int &ind_t_false = ind_t_false_init;

            if(shell.size()>1)
            {
                while(i < shell.size() && !is_in_next_triangle)
                {
                    ind_triangle=shell[i];

                    if(ind_triangle != ind_t)
                    {
                        is_in_next_triangle = isPointInTriangle(ind_triangle, false, p, ind_t_false);
                    }
                    i++;
                }

                if(is_in_next_triangle)
                {
                    ind_t_test=ind_triangle;
                    //sout << "correct to triangle indexed by " << ind_t_test << sendl;
                }
                else // not found
                {
                    //sout << "not found !!! " << sendl;
                    ind_t_test=ind_t;
                }
            }
            else
            {
                ind_t_test=ind_t;
            }
        }
        return is_inside;

    }
    else // triangle is flat
    {
        //sout << "INFO_print : triangle is flat" << sendl;
        return false;
    }
}

// Tests how to triangularize a quad whose vertices are defined by (p_q1, p_q2, ind_q3, ind_q4) according to the Delaunay criterion
template<class DataTypes>
bool TriangleSetGeometryAlgorithms< DataTypes >::isQuadDeulaunayOriented(const typename DataTypes::Coord& p_q1,
        const typename DataTypes::Coord& p_q2,
        unsigned int ind_q3,
        unsigned int ind_q4)
{
    sofa::helper::vector< double > baryCoefs;

    const typename DataTypes::VecCoord& vect_c = *(this->object->getX());

    const typename DataTypes::Coord& c3 = vect_c[ind_q3];
    const typename DataTypes::Coord& c4 = vect_c[ind_q4];

    return isQuadDeulaunayOriented(p_q1, p_q2, c3, c4);
}

/** \brief Tests how to triangularize a quad whose vertices are defined by (p1, p2, p3, p4) according to the Delaunay criterion
 *
 */
template<class DataTypes>
bool TriangleSetGeometryAlgorithms< DataTypes >::isQuadDeulaunayOriented(const typename DataTypes::Coord& p1,
        const typename DataTypes::Coord& p2,
        const typename DataTypes::Coord& p3,
        const typename DataTypes::Coord& p4)
{
    Coord tri1[3], tri2[3];

    tri1[0] = p1; tri1[1] = p2; tri1[2] = p3;
    tri2[0] = p3; tri2[1] = p4; tri2[2] = p1;


    //Test if one vertex is inside the triangle fromed by the the 3 others
    Coord CommonEdge[2], oppositeVertices[2];

    oppositeVertices[0] = p1; sofa::defaulttype::Vec<3,double> A; A = p1;
    CommonEdge[0] = p2;       sofa::defaulttype::Vec<3,double> C; C = p2;
    CommonEdge[1] = p4;       sofa::defaulttype::Vec<3,double> B; B = p3;
    oppositeVertices[1] = p3; sofa::defaulttype::Vec<3,double> D; D = p4;

    bool intersected = false;

    Coord inter = this->compute2EdgesIntersection (CommonEdge, oppositeVertices, intersected);

    if (intersected)
    {

        sofa::defaulttype::Vec<3,double> X; DataTypes::get(X[0], X[1], X[2], inter);

        double ABAX = (A - B)*(A - X);
        double CDCX = (C - D)*(C - X);

        if ( (ABAX < 0) || ((A - X).norm2() > (A - B).norm2()) )
            return true;
        else if (	(CDCX < 0) || ((C - X).norm2() > (C - D).norm2()) )
            return false;
    }

    sofa::defaulttype::Vec<3,double> G = (A+B+C)/3.0;

    if((G-C)*(G-C) <= (G-D)*(G-D))
        return true;
    else
        return false;
}


template<class DataTypes>
bool TriangleSetGeometryAlgorithms< DataTypes >::isDiagonalsIntersectionInQuad (const typename DataTypes::Coord triangle1[3],const  typename DataTypes::Coord triangle2[3])
{

    Coord CommonEdge[2], oppositeVertices[2];
    unsigned int cpt = 0;
    bool test = false;

    for (unsigned int i = 0; i<3; i++)
    {
        test = false;
        for (unsigned int j = 0; j<3; j++)
            if(triangle1[i] == triangle2[j])
            {
                test = true;
                break;
            }

        if(test)
        {
            CommonEdge[cpt] = triangle1[i];
            cpt++;
        }
        else
            oppositeVertices[0] = triangle1[i];
    }


    for (unsigned int i = 0; i<3; i++)
    {
        test = false;
        for (unsigned int j = 0; j<2; j++)
            if (triangle2[i] == CommonEdge[j])
            {
                test = true;
                break;
            }

        if (!test)
        {
            oppositeVertices[1] = triangle2[i];
            break;
        }
    }

    bool intersected = false;

    Coord inter = this->compute2EdgesIntersection (CommonEdge, oppositeVertices, intersected);

    if (intersected)
    {
        sofa::defaulttype::Vec<3,double> A; DataTypes::get(A[0], A[1], A[2], CommonEdge[0]);
        sofa::defaulttype::Vec<3,double> B; DataTypes::get(B[0], B[1], B[2], CommonEdge[1]);

        sofa::defaulttype::Vec<3,double> C; DataTypes::get(C[0], C[1], C[2], oppositeVertices[0]);
        sofa::defaulttype::Vec<3,double> D; DataTypes::get(D[0], D[1], D[2], oppositeVertices[1]);

        sofa::defaulttype::Vec<3,double> X; DataTypes::get(X[0], X[1], X[2], inter);

        double ABAX = (A - B)*(A - X);
        double CDCX = (C - D)*(C - X);

        if ( (ABAX < 0) || (CDCX < 0) || ((A - X).norm2() > (A - B).norm2()) || ((C - X).norm2() > (C - D).norm2()) )
            return false;
        else
            return true;
    }

    return false;
}


// Test if a triangle indexed by ind_triangle (and incident to the vertex indexed by ind_p) is included or not in the plane defined by (ind_p, plane_vect)
template<class DataTypes>
bool TriangleSetGeometryAlgorithms< DataTypes >::isTriangleInPlane(const TriangleID ind_t,
        const unsigned int ind_p,
        const sofa::defaulttype::Vec<3,Real>&plane_vect) const
{
    const Triangle &t=this->m_topology->getTriangle(ind_t);

    // HYP : ind_p==t[0] or ind_p==t[1] or ind_p==t[2]

    const typename DataTypes::VecCoord& vect_c = *(this->object->getX());

    unsigned int ind_1;
    unsigned int ind_2;

    if(ind_p==t[0])
    {
        ind_1=t[1];
        ind_2=t[2];
    }
    else
    {
        if(ind_p==t[1])
        {
            ind_1=t[2];
            ind_2=t[0];
        }
        else // ind_p==t[2]
        {
            ind_1=t[0];
            ind_2=t[1];
        }
    }

    const typename DataTypes::Coord& c0=vect_c[ind_p];
    const typename DataTypes::Coord& c1=vect_c[ind_1];
    const typename DataTypes::Coord& c2=vect_c[ind_2];

    sofa::defaulttype::Vec<3,Real> p0;
    p0[0] = (Real) (c0[0]);
    p0[1] = (Real) (c0[1]);
    p0[2] = (Real) (c0[2]);
    sofa::defaulttype::Vec<3,Real> p1;
    p1[0] = (Real) (c1[0]);
    p1[1] = (Real) (c1[1]);
    p1[2] = (Real) (c1[2]);
    sofa::defaulttype::Vec<3,Real> p2;
    p2[0] = (Real) (c2[0]);
    p2[1] = (Real) (c2[1]);
    p2[2] = (Real) (c2[2]);

    return((p1-p0)*( plane_vect)>=0.0 && (p1-p0)*( plane_vect)>=0.0);
}

// Prepares the duplication of a vertex
template<class DataTypes>
void TriangleSetGeometryAlgorithms< DataTypes >::prepareVertexDuplication(const unsigned int ind_p,
        const TriangleID ind_t_from,
        const TriangleID ind_t_to,
        const Edge& indices_from,
        const double &coord_from,
        const Edge& indices_to,
        const double &coord_to,
        sofa::helper::vector< unsigned int > &triangles_list_1,
        sofa::helper::vector< unsigned int > &triangles_list_2) const
{
    //HYP : if coord_from or coord_to == 0.0 or 1.0, ind_p is distinct from ind_from and from ind_to

    const typename DataTypes::VecCoord& vect_c = *(this->object->getX());

    const typename DataTypes::Coord& c_p = vect_c[ind_p];
    sofa::defaulttype::Vec<3,Real> point_p;
    point_p[0]= (Real) c_p[0];
    point_p[1]= (Real) c_p[1];
    point_p[2]= (Real) c_p[2];

    sofa::defaulttype::Vec<3,Real> point_from=(Vec<3,Real>) getOppositePoint(ind_p, indices_from, coord_from);
    sofa::defaulttype::Vec<3,Real> point_to=(Vec<3,Real>) getOppositePoint(ind_p, indices_to, coord_to);

    //Vec<3,Real> point_from=(Vec<3,Real>) computeBaryEdgePoint((sofa::helper::vector< unsigned int>&) indices_from, coord_from);
    //Vec<3,Real> point_to=(Vec<3,Real>) computeBaryEdgePoint((sofa::helper::vector< unsigned int>&) indices_to, coord_to);

    sofa::defaulttype::Vec<3,Real> vect_from = point_from - point_p;
    sofa::defaulttype::Vec<3,Real> vect_to = point_p - point_to;

    //sout << "INFO_print : vect_from = " << vect_from <<  sendl;
    //sout << "INFO_print : vect_to = " << vect_to <<  sendl;

    sofa::defaulttype::Vec<3,Real> normal_from;
    sofa::defaulttype::Vec<3,Real> normal_to;

    sofa::defaulttype::Vec<3,Real> plane_from;
    sofa::defaulttype::Vec<3,Real> plane_to;

    if((coord_from!=0.0) && (coord_from!=1.0))
    {
        normal_from=(Vec<3,Real>) computeTriangleNormal(ind_t_from);
        plane_from=vect_from.cross( normal_from); // inverse ??
    }
    else
    {
        // HYP : only 2 edges maximum are adjacent to the same triangle (otherwise : compute the one which minimizes the normed dotProduct and which gives the positive cross)

        unsigned int ind_edge;

        if(coord_from==0.0)
        {
            ind_edge=this->m_topology->getEdgeIndex(indices_from[0], ind_p);
        }
        else // coord_from==1.0
        {
            ind_edge=this->m_topology->getEdgeIndex(indices_from[1], ind_p);
        }

        if (this->m_topology->getNbEdges()>0)
        {
            sofa::helper::vector< unsigned int > shell =(sofa::helper::vector< unsigned int >) (this->m_topology->getTrianglesAroundEdge(ind_edge));
            unsigned int ind_triangle=shell[0];
            unsigned int i=0;
            bool is_in_next_triangle=false;

            if(shell.size()>1)
            {
                while(i < shell.size() || !is_in_next_triangle)
                {
                    if(shell[i] != ind_t_from)
                    {
                        ind_triangle=shell[i];
                        is_in_next_triangle=true;
                    }
                    i++;
                }
            }
            else
            {
                return;
            }

            if(is_in_next_triangle)
            {
                sofa::defaulttype::Vec<3,Real> normal_from_1=(Vec<3,Real>) computeTriangleNormal(ind_triangle);
                sofa::defaulttype::Vec<3,Real> normal_from_2=(Vec<3,Real>) computeTriangleNormal(ind_t_from);

                normal_from=(normal_from_1+normal_from_2)/2.0;
                plane_from=vect_from.cross( normal_from);
            }
            else
            {
                return;
            }
        }
    }

    if((coord_to!=0.0) && (coord_to!=1.0))
    {
        normal_to=(Vec<3,Real>) computeTriangleNormal(ind_t_to);

        plane_to=vect_to.cross( normal_to);
    }
    else
    {
        // HYP : only 2 edges maximum are adjacent to the same triangle (otherwise : compute the one which minimizes the normed dotProduct and which gives the positive cross)

        unsigned int ind_edge;

        if(coord_to==0.0)
        {
            ind_edge=this->m_topology->getEdgeIndex(indices_to[0], ind_p);
        }
        else // coord_to==1.0
        {
            ind_edge=this->m_topology->getEdgeIndex(indices_to[1], ind_p);
        }

        if (this->m_topology->getNbEdges()>0)
        {
            sofa::helper::vector< unsigned int > shell =(sofa::helper::vector< unsigned int >) (this->m_topology->getTrianglesAroundEdge(ind_edge));
            unsigned int ind_triangle=shell[0];
            unsigned int i=0;
            bool is_in_next_triangle=false;

            if(shell.size()>1)
            {
                while(i < shell.size() || !is_in_next_triangle)
                {
                    if(shell[i] != ind_t_to)
                    {
                        ind_triangle=shell[i];
                        is_in_next_triangle=true;
                    }
                    i++;
                }
            }
            else
            {
                return;
            }

            if(is_in_next_triangle)
            {
                sofa::defaulttype::Vec<3,Real> normal_to_1=(Vec<3,Real>) computeTriangleNormal(ind_triangle);
                sofa::defaulttype::Vec<3,Real> normal_to_2=(Vec<3,Real>) computeTriangleNormal(ind_t_to);

                normal_to=(normal_to_1+normal_to_2)/2.0;
                plane_to=vect_to.cross( normal_to);
            }
            else
            {
                return;
            }
        }
    }

    if (this->m_topology->getNbPoints()>0)
    {
        sofa::helper::vector< unsigned int > shell =(sofa::helper::vector< unsigned int >) (this->m_topology->getTrianglesAroundVertex(ind_p));
        unsigned int ind_triangle=shell[0];
        unsigned int i=0;

        bool is_in_plane_from;
        bool is_in_plane_to;

        if(shell.size()>1)
        {
            sofa::defaulttype::Vec<3,Real> normal_test = plane_from.cross( plane_to);
            Real value_test =   normal_test*(normal_from+normal_to);

            if(value_test<=0.0)
            {
                //sout << "INFO_print : CONVEXE, value_test = " << value_test <<  sendl;
                //sout << "INFO_print : shell.size() = " << shell.size() << ", ind_t_from = " << ind_t_from << ", ind_t_to = " << ind_t_to <<  sendl;

                while(i < shell.size())
                {
                    ind_triangle=shell[i];

                    is_in_plane_from=isTriangleInPlane(ind_triangle,ind_p, (const sofa::defaulttype::Vec<3,double>&) plane_from);
                    is_in_plane_to=isTriangleInPlane(ind_triangle,ind_p, (const sofa::defaulttype::Vec<3,double>&) plane_to);

                    if((ind_triangle != ind_t_from) && (ind_triangle != ind_t_to))
                    {
                        if(is_in_plane_from || is_in_plane_to)
                        {
                            triangles_list_1.push_back(ind_triangle);
                        }
                        else
                        {
                            triangles_list_2.push_back(ind_triangle);
                        }
                    }
                    i++;
                }
            }
            else // value_test>0.0
            {
                //sout << "INFO_print : CONCAVE, value_test = " << value_test <<  sendl;
                //sout << "INFO_print : shell.size() = " << shell.size() << ", ind_t_from = " << ind_t_from << ", ind_t_to = " << ind_t_to <<  sendl;

                while(i < shell.size())
                {
                    ind_triangle=shell[i];

                    is_in_plane_from=isTriangleInPlane(ind_triangle,ind_p, (const sofa::defaulttype::Vec<3,double>&) plane_from);
                    is_in_plane_to=isTriangleInPlane(ind_triangle,ind_p, (const sofa::defaulttype::Vec<3,double>&) plane_to);

                    if((ind_triangle != ind_t_from) && (ind_triangle != ind_t_to))
                    {
                        if(is_in_plane_from && is_in_plane_to)
                        {
                            triangles_list_1.push_back(ind_triangle);
                        }
                        else
                        {
                            triangles_list_2.push_back(ind_triangle);
                        }
                    }
                    i++;
                }
            }
        }
        else
        {
            return;
        }
    }
    else
    {
        return;
    }
}

// Computes the intersection of the segment from point a to point b and the triangle indexed by t
template<class DataTypes>
bool TriangleSetGeometryAlgorithms< DataTypes >::computeSegmentTriangleIntersection(bool is_entered,
        const sofa::defaulttype::Vec<3,double>& a,
        const sofa::defaulttype::Vec<3,double>& b,
        const TriangleID ind_t,
        sofa::helper::vector<unsigned int> &indices,
        double &baryCoef, double& coord_kmin) const
{
    // HYP : point a is in triangle indexed by t
    // is_entered == true => indices.size() == 2





    unsigned int ind_first=0;
    unsigned int ind_second=0;

    if(indices.size()>1)
    {
        ind_first=indices[0];
        ind_second=indices[1];
    }

    indices.clear();

    bool is_validated = false;
    bool is_intersected = false;

    const Triangle &t=this->m_topology->getTriangle(ind_t);
    const typename DataTypes::VecCoord& vect_c = *(this->object->getX());

    bool is_full_01=(is_entered && ((t[0] == ind_first && t[1] == ind_second) || (t[1] == ind_first && t[0] == ind_second)));
    bool is_full_12=(is_entered && ((t[1] == ind_first && t[2] == ind_second) || (t[2] == ind_first && t[1] == ind_second)));
    bool is_full_20=(is_entered && ((t[2] == ind_first && t[0] == ind_second) || (t[0] == ind_first && t[2] == ind_second)));

    const typename DataTypes::Coord& c0=vect_c[t[0]];
    const typename DataTypes::Coord& c1=vect_c[t[1]];
    const typename DataTypes::Coord& c2=vect_c[t[2]];

    sofa::defaulttype::Vec<3,Real> p0;
    p0[0] = (Real) (c0[0]);
    p0[1] = (Real) (c0[1]);
    p0[2] = (Real) (c0[2]);
    sofa::defaulttype::Vec<3,Real> p1;
    p1[0] = (Real) (c1[0]);
    p1[1] = (Real) (c1[1]);
    p1[2] = (Real) (c1[2]);
    sofa::defaulttype::Vec<3,Real> p2;
    p2[0] = (Real) (c2[0]);
    p2[1] = (Real) (c2[1]);
    p2[2] = (Real) (c2[2]);

    sofa::defaulttype::Vec<3,Real> pa;
    pa[0] = (Real) (a[0]);
    pa[1] = (Real) (a[1]);
    pa[2] = (Real) (a[2]);
    sofa::defaulttype::Vec<3,Real> pb;
    pb[0] = (Real) (b[0]);
    pb[1] = (Real) (b[1]);
    pb[2] = (Real) (b[2]);

    sofa::defaulttype::Vec<3,Real> v_normal = (p2-p0).cross(p1-p0);
    //Vec<3,Real> v_normal = (Vec<3,Real>) computeTriangleNormal(ind_t);

    Real norm_v_normal = v_normal.norm(); // WARN : square root COST

    if(norm_v_normal != 0.0)
    {

        v_normal/=norm_v_normal;

        sofa::defaulttype::Vec<3,Real> v_ab = pb-pa;
        sofa::defaulttype::Vec<3,Real> v_ab_proj = v_ab - v_normal * dot(v_ab,v_normal); // projection (same values if incision in the plan)
        sofa::defaulttype::Vec<3,Real> pb_proj = v_ab_proj + pa;

        sofa::defaulttype::Vec<3,Real> v_01 = p1-p0;
        sofa::defaulttype::Vec<3,Real> v_12 = p2-p1;
        sofa::defaulttype::Vec<3,Real> v_20 = p0-p2;

        sofa::defaulttype::Vec<3,Real> n_proj =v_ab_proj.cross(v_normal);

        sofa::defaulttype::Vec<3,Real> n_01 = v_01.cross(v_normal);
        sofa::defaulttype::Vec<3,Real> n_12 = v_12.cross(v_normal);
        sofa::defaulttype::Vec<3,Real> n_20 = v_20.cross(v_normal);

        Real norm2_v_ab_proj = v_ab_proj*(v_ab_proj); //dot product WARNING

        if(norm2_v_ab_proj != 0.0) // pb_proj != pa
        {
            double coord_t=0.0;
            double coord_k=0.0;

            double is_initialized=false;
            coord_kmin=0.0;

            double coord_test1;
            double coord_test2;

            double s_t;
            double s_k;

            if(!is_full_01)
            {
                /// Test of edge (p0,p1) :
                s_t = (p0-p1)*n_proj;
                s_k = (pa-pb_proj)*n_01;

                // s_t == 0.0 iff s_k == 0.0

                if(s_t==0.0) // (pa,pb_proj) and (p0,p1) are parallel
                {
                    if((p0-pa)*(n_proj)==0.0) // (pa,pb_proj) and (p0,p1) are on the same line
                    {
                        coord_test1 = (pa-p0)*(pa-pb_proj)/norm2_v_ab_proj; // HYP : pb_proj != pa
                        coord_test2 = (pa-p1)*(pa-pb_proj)/norm2_v_ab_proj; // HYP : pb_proj != pa

                        if(coord_test1>=0)
                        {
                            coord_k=coord_test1;
                            coord_t=0.0;
                        }
                        else
                        {
                            coord_k=coord_test2;
                            coord_t=1.0;
                        }

                        is_intersected = (coord_k > 0.0 && (coord_t >= 0.0 && coord_t <= 1.0));

                    }
                    else  // (pa,pb_proj) and (p0,p1) are parallel and disjoint
                    {
                        is_intersected=false;
                    }
                }
                else // s_t != 0.0 and s_k != 0.0
                {
                    coord_k=double((pa-p0)*(n_01))*1.0/double(s_k);
                    coord_t=double((p0-pa)*(n_proj))*1.0/double(s_t);

                    is_intersected = ((coord_k > 0.0) && (coord_t >= 0.0 && coord_t <= 1.0));
                }

                if(is_intersected)
                {
                    if((!is_initialized) || (coord_k > coord_kmin))
                    {
                        indices.clear();
                        indices.push_back(t[0]);
                        indices.push_back(t[1]);
                        baryCoef=coord_t;
                        coord_kmin=coord_k;
                    }

                    is_initialized=true;
                }

                is_validated = is_validated || is_initialized;
            }



            if(!is_full_12)
            {
                /// Test of edge (p1,p2) :

                s_t = (p1-p2)*(n_proj);
                s_k = (pa-pb_proj)*(n_12);

                // s_t == 0.0 iff s_k == 0.0

                if(s_t==0.0) // (pa,pb_proj) and (p1,p2) are parallel
                {
                    if((p1-pa)*(n_proj)==0.0) // (pa,pb_proj) and (p1,p2) are on the same line
                    {
                        coord_test1 = (pa-p1)*(pa-pb_proj)/norm2_v_ab_proj; // HYP : pb_proj != pa
                        coord_test2 = (pa-p2)*(pa-pb_proj)/norm2_v_ab_proj; // HYP : pb_proj != pa

                        if(coord_test1>=0)
                        {
                            coord_k=coord_test1;
                            coord_t=0.0;
                        }
                        else
                        {
                            coord_k=coord_test2;
                            coord_t=1.0;
                        }

                        is_intersected = (coord_k > 0.0 && (coord_t >= 0.0 && coord_t <= 1.0));
                    }
                    else // (pa,pb_proj) and (p1,p2) are parallel and disjoint
                    {
                        is_intersected=false;
                    }

                }
                else   // s_t != 0.0 and s_k != 0.0
                {

                    coord_k=double((pa-p1)*(n_12))*1.0/double(s_k);
                    coord_t=double((p1-pa)*(n_proj))*1.0/double(s_t);

                    is_intersected = ((coord_k > 0.0) && (coord_t >= 0.0 && coord_t <= 1.0));
                }

                if(is_intersected)
                {
                    if((!is_initialized) || (coord_k > coord_kmin))
                    {
                        indices.clear();
                        indices.push_back(t[1]);
                        indices.push_back(t[2]);
                        baryCoef=coord_t;
                        coord_kmin=coord_k;
                    }

                    is_initialized=true;
                }

                is_validated = is_validated || is_initialized;
            }



            if(!is_full_20)
            {
                /// Test of edge (p2,p0) :

                s_t = (p2-p0)*(n_proj);
                s_k = (pa-pb_proj)*(n_20);

                // s_t == 0.0 iff s_k == 0.0

                if(s_t==0.0) // (pa,pb_proj) and (p2,p0) are parallel
                {
                    if((p2-pa)*(n_proj)==0.0) // (pa,pb_proj) and (p2,p0) are on the same line
                    {
                        coord_test1 = (pa-p2)*(pa-pb_proj)/norm2_v_ab_proj; // HYP : pb_proj != pa
                        coord_test2 = (pa-p0)*(pa-pb_proj)/norm2_v_ab_proj; // HYP : pb_proj != pa

                        if(coord_test1>=0)
                        {
                            coord_k=coord_test1;
                            coord_t=0.0;
                        }
                        else
                        {
                            coord_k=coord_test2;
                            coord_t=1.0;
                        }

                        is_intersected = (coord_k > 0.0 && (coord_t >= 0.0 && coord_t <= 1.0));

                    }
                    else // (pa,pb_proj) and (p2,p0) are parallel and disjoint
                    {
                        is_intersected = false;
                    }
                }
                else // s_t != 0.0 and s_k != 0.0
                {
                    coord_k=double((pa-p2)*(n_20))*1.0/double(s_k);
                    coord_t=double((p2-pa)*(n_proj))*1.0/double(s_t);

                    is_intersected = ((coord_k > 0.0) && (coord_t >= 0.0 && coord_t <= 1.0));
                }

                if(is_intersected)
                {
                    if((!is_initialized) || (coord_k > coord_kmin))
                    {
                        indices.clear();
                        indices.push_back(t[2]);
                        indices.push_back(t[0]);
                        baryCoef=coord_t;
                        coord_kmin=coord_k;
                    }

                    is_initialized = true;
                }
                is_validated = is_validated || is_initialized;
            }


        }
        else
        {
            //std::cout << "points a and b are projected to the same point on triangle t" << std::endl;
            is_validated = false; // points a and b are projected to the same point on triangle t
        }
    }
    else
    {
        //std::cout << "triangle t is flat" << std::endl;
        is_validated = false; // triangle t is flat
    }

    return is_validated;
}





// Computes the list of points (edge,coord) intersected by the segment from point a to point b
// and the triangular mesh
template<class DataTypes>
bool TriangleSetGeometryAlgorithms< DataTypes >::computeIntersectedPointsList(const PointID last_point,
        const sofa::defaulttype::Vec<3,double>& a,
        const sofa::defaulttype::Vec<3,double>& b,
        unsigned int& ind_ta,
        unsigned int& ind_tb,
        sofa::helper::vector< unsigned int > &triangles_list,
        sofa::helper::vector<unsigned int> &edges_list,
        sofa::helper::vector< double >& coords_list,
        bool& is_on_boundary) const
{

    bool is_validated=true;
    bool is_intersected=true;

    sofa::defaulttype::Vec<3,double> c_t_test = a;

    is_on_boundary = false;

    sofa::helper::vector<unsigned int> indices;

    double coord_t=0.0;
    double coord_k=0.0;
    double coord_k_test=0.0;
    double dist_min=0.0;

    sofa::defaulttype::Vec<3,double> p_current=a;

    TriangleID ind_t_current=ind_ta;
    EdgeID ind_edge;
    PointID ind_index;
    TriangleID ind_triangle = ind_ta;
    is_intersected=computeSegmentTriangleIntersection(false, p_current, b, (const unsigned int) ind_t_current, indices, coord_t, coord_k);


    // In case the ind_t is not the good one.
    if ( (!is_intersected || indices[0] == last_point || indices[1] == last_point) && (last_point != BaseMeshTopology::InvalidID))
    {

        const sofa::helper::vector< unsigned int >& shell = this->m_topology->getTrianglesAroundVertex (last_point);

        for (unsigned int i = 0; i<shell.size(); i++)
        {
            if (shell [i] != ind_t_current)
                is_intersected=computeSegmentTriangleIntersection(false, p_current, b, shell[i], indices, coord_t, coord_k);

            if (is_intersected && indices[0] != last_point && indices[1] != last_point)
            {
                ind_t_current = shell[i];
                ind_ta = ind_t_current;
                break;
            }
        }
    }

    if (ind_ta == ind_tb)
    {
#ifndef NDEBUG
        std::cout << "INFO_print - TriangleSetTopology.inl : Cut is not reached because inputs elements are the same element." << std::endl;
#endif
        return false;
    }


#ifndef NDEBUG
    std::cout << "*********************************" << std::endl;
    std::cout << "ind_t_current: " << ind_t_current << std::endl;
    std::cout << "p_current: " << p_current << std::endl;
    std::cout << "coord_t: " << coord_t << std::endl;
    std::cout << "coord_k: " << coord_k << std::endl;
    std::cout << "indices: " << indices << std::endl;
    std::cout << "last_point: " << last_point << std::endl;
    std::cout << "a: " << a << std::endl;
    std::cout << "b: " << b << std::endl;
    std::cout << "is_intersected: "<< is_intersected << std::endl;
    std::cout << "*********************************" << std::endl;
#endif

    coord_k_test=coord_k;
    dist_min=(b-a)*(b-a);

    while((coord_k_test<1.0 && is_validated) && is_intersected)
    {
        ind_edge=this->m_topology->getEdgeIndex(indices[0],indices[1]);
        edges_list.push_back(ind_edge);
        triangles_list.push_back(ind_t_current);
        if (this->m_topology->getEdge(ind_edge)[0] == indices[0])
            coords_list.push_back(coord_t);
        else
            coords_list.push_back(1.0-coord_t);

        const typename DataTypes::VecCoord& vect_c = *(this->object->getX());

        sofa::defaulttype::Vec<3,double> c_t_current; // WARNING : conversion from 'double' to 'float', possible loss of data ! // typename DataTypes::Coord
        c_t_current[0]=(double) ((1.0-coord_t)*((double) (vect_c[indices[0]][0]))+coord_t*((double) (vect_c[indices[1]][0])));
        c_t_current[1]=(double) ((1.0-coord_t)*((double) (vect_c[indices[0]][1]))+coord_t*((double) (vect_c[indices[1]][1])));
        c_t_current[2]=(double) ((1.0-coord_t)*((double) (vect_c[indices[0]][2]))+coord_t*((double) (vect_c[indices[1]][2])));

        p_current=c_t_current;

        sofa::defaulttype::Vec<3,Real> p_t_aux;
        p_t_aux[0] = (Real) (c_t_current[0]);
        p_t_aux[1] = (Real) (c_t_current[1]);
        p_t_aux[2] = (Real) (c_t_current[2]);




        if(coord_t==0.0 || coord_t==1.0) // current point indexed by ind_t_current is on a vertex
        {
            //std::cout << "INFO_print : INPUT ON A VERTEX !!!" <<  std::endl;
            if(coord_t==0.0)
            {
                ind_index=indices[0];
            }
            else // coord_t==1.0
            {
                ind_index=indices[1];
            }

            if (this->m_topology->getNbPoints() >0)
            {
                sofa::helper::vector< unsigned int > shell =(sofa::helper::vector< unsigned int >) (this->m_topology->getTrianglesAroundVertex(ind_index));
                ind_triangle=shell[0];
                unsigned int i=0;
                //bool is_in_next_triangle=false;
                bool is_test_init=false;

                unsigned int ind_from = ind_t_current;

                if(shell.size()>1) // at leat one neighbor triangle which is not indexed by ind_t_current
                {
                    is_on_boundary=false;

                    while(i < shell.size())
                    {
                        if(shell[i] != ind_from)
                        {
                            ind_triangle=shell[i];

                            const Triangle &t=this->m_topology->getTriangle(ind_triangle);

                            const typename DataTypes::Coord& c0=vect_c[t[0]];
                            const typename DataTypes::Coord& c1=vect_c[t[1]];
                            const typename DataTypes::Coord& c2=vect_c[t[2]];

                            sofa::defaulttype::Vec<3,Real> p0_aux;
                            p0_aux[0] = (Real) (c0[0]);
                            p0_aux[1] = (Real) (c0[1]);
                            p0_aux[2] = (Real) (c0[2]);
                            sofa::defaulttype::Vec<3,Real> p1_aux;
                            p1_aux[0] = (Real) (c1[0]);
                            p1_aux[1] = (Real) (c1[1]);
                            p1_aux[2] = (Real) (c1[2]);
                            sofa::defaulttype::Vec<3,Real> p2_aux;
                            p2_aux[0] = (Real) (c2[0]);
                            p2_aux[1] = (Real) (c2[1]);
                            p2_aux[2] = (Real) (c2[2]);

                            is_intersected=computeSegmentTriangleIntersection(true, p_current, b, ind_triangle, indices, coord_t, coord_k);

                            if(is_intersected)
                            {
                                sofa::defaulttype::Vec<3,double> c_t_test; // WARNING : conversion from 'double' to 'float', possible loss of data ! // typename DataTypes::Coord
                                c_t_test[0]=(double) ((1.0-coord_t)*((double) (vect_c[indices[0]][0]))+coord_t*((double) (vect_c[indices[1]][0])));
                                c_t_test[1]=(double) ((1.0-coord_t)*((double) (vect_c[indices[0]][1]))+coord_t*((double) (vect_c[indices[1]][1])));
                                c_t_test[2]=(double) ((1.0-coord_t)*((double) (vect_c[indices[0]][2]))+coord_t*((double) (vect_c[indices[1]][2])));

                                double dist_test=(b-c_t_test)*(b-c_t_test);

                                if(is_test_init)
                                {
                                    if(dist_test<dist_min && coord_k<=1) //dist_test<dist_min
                                    {
                                        coord_k_test=coord_k;
                                        dist_min=dist_test;
                                        ind_t_current=ind_triangle;
                                    }
                                }
                                else
                                {
                                    is_test_init=true;
                                    coord_k_test=coord_k;
                                    dist_min=dist_test;
                                    ind_t_current=ind_triangle;
                                }
                            }
                        }

                        i=i+1;
                    }

                    is_intersected=is_test_init;

                }
                else
                {
                    is_on_boundary=true;
                    is_validated=false;
                }
            }
            else
            {
                is_validated=false;
            }
        }
        else // current point indexed by ind_t_current is on an edge, but not on a vertex
        {
            ind_edge=this->m_topology->getEdgeIndex(indices[0],indices[1]);

            if (this->m_topology->getNbEdges()>0)
            {
                sofa::helper::vector< unsigned int > shell =(sofa::helper::vector< unsigned int >) (this->m_topology->getTrianglesAroundEdge(ind_edge));

                ind_triangle=shell[0];
                unsigned int i=0;

                bool is_test_init=false;

                unsigned int ind_from = ind_t_current;

                if(shell.size()>0) // at leat one neighbor triangle which is not indexed by ind_t_current
                {
                    is_on_boundary=false;

                    while(i < shell.size())
                    {
                        if(shell[i] != ind_from)
                        {
                            ind_triangle=shell[i];

                            const Triangle &t=this->m_topology->getTriangle(ind_triangle);

                            const typename DataTypes::Coord& c0=vect_c[t[0]];
                            const typename DataTypes::Coord& c1=vect_c[t[1]];
                            const typename DataTypes::Coord& c2=vect_c[t[2]];

                            sofa::defaulttype::Vec<3,Real> p0_aux;
                            p0_aux[0] = (Real) (c0[0]);
                            p0_aux[1] = (Real) (c0[1]);
                            p0_aux[2] = (Real) (c0[2]);
                            sofa::defaulttype::Vec<3,Real> p1_aux;
                            p1_aux[0] = (Real) (c1[0]);
                            p1_aux[1] = (Real) (c1[1]);
                            p1_aux[2] = (Real) (c1[2]);
                            sofa::defaulttype::Vec<3,Real> p2_aux;
                            p2_aux[0] = (Real) (c2[0]);
                            p2_aux[1] = (Real) (c2[1]);
                            p2_aux[2] = (Real) (c2[2]);

                            is_intersected=computeSegmentTriangleIntersection(true, p_current, b, ind_triangle, indices, coord_t, coord_k);

                            if(is_intersected)
                            {
                                //Vec<3,double> c_t_test; // WARNING : conversion from 'double' to 'float', possible loss of data ! // typename DataTypes::Coord
                                c_t_test[0]=(double) ((1.0-coord_t)*((double) (vect_c[indices[0]][0]))+coord_t*((double) (vect_c[indices[1]][0])));
                                c_t_test[1]=(double) ((1.0-coord_t)*((double) (vect_c[indices[0]][1]))+coord_t*((double) (vect_c[indices[1]][1])));
                                c_t_test[2]=(double) ((1.0-coord_t)*((double) (vect_c[indices[0]][2]))+coord_t*((double) (vect_c[indices[1]][2])));

                                double dist_test=(b-c_t_test)*(b-c_t_test);

                                if(is_test_init)
                                {
                                    if(dist_test<dist_min && coord_k<=1) //dist_test<dist_min
                                    {
                                        coord_k_test=coord_k;
                                        dist_min=dist_test;
                                        ind_t_current=ind_triangle;
                                    }
                                }
                                else
                                {
                                    is_test_init=true;
                                    coord_k_test=coord_k;
                                    dist_min=dist_test;
                                    ind_t_current=ind_triangle;
                                }
                            }
                        }

                        i=i+1;
                    }


                    is_intersected=is_test_init;
                }
                else
                {
                    is_on_boundary=true;
                    is_validated=false;
                }
            }
            else
            {
                is_validated=false;
            }
        }
    }

    if (ind_tb == BaseMeshTopology::InvalidID)
        ind_tb = ind_triangle;

    bool is_reached = (ind_tb==ind_triangle && coord_k_test>=1.0);

#ifndef NDEBUG
    if(is_reached)
    {
        std::cout << "INFO_print - TriangleSetTopology.inl : Cut is reached" << std::endl;
    }

    if(is_on_boundary)
    {

        std::cout << "INFO_print - TriangleSetTopology.inl : Cut meets a mesh boundary" << std::endl;
    }
#endif


    if(!is_reached && !is_on_boundary)
    {
#ifndef NDEBUG
        std::cout << "INFO_print - TriangleSetTopology.inl : Cut is not reached" << std::endl;
#endif
//      ind_tb=ind_triangle;
//      return true;
    }

    return (is_reached && is_validated && is_intersected); // b is in triangle indexed by ind_t_current

}



template <typename DataTypes>
bool TriangleSetGeometryAlgorithms<DataTypes>::computeIntersectedObjectsList (const PointID last_point,
        const sofa::defaulttype::Vec<3,double>& a, const sofa::defaulttype::Vec<3,double>& b,
        unsigned int& ind_ta, unsigned int& ind_tb,// A verifier pourquoi la ref!
        sofa::helper::vector< sofa::core::topology::TopologyObjectType>& topoPath_list,
        sofa::helper::vector<unsigned int>& indices_list,
        sofa::helper::vector< sofa::defaulttype::Vec<3, double> >& coords_list) const
{
    //// QUICK FIX TO USE THE NEW PATH DECLARATION (WITH ONLY EDGES COMING FROM PREVIOUS FUNCTION)
    //// ** TODO: create the real function handle different objects intersection **
    // QUICK FIX for fracture: border points (a and b) can be a point.

    // Output declarations
    sofa::helper::vector<unsigned int> triangles_list;
    sofa::helper::vector<unsigned int> edges_list;
    sofa::helper::vector< double > coordsEdge_list;
    bool pathOK;
    bool isOnPoint = false;
    bool is_on_boundary = false;

    // using old function:
    pathOK = this->computeIntersectedPointsList (last_point, a, b, ind_ta, ind_tb, triangles_list, edges_list, coordsEdge_list, is_on_boundary);

#ifndef NDEBUG
    std::cout << "*********************************" << std::endl;
    std::cout << "last_point: " << last_point << std::endl;
    std::cout << "a: " << a << std::endl;
    std::cout << "b: " << b << std::endl;
    std::cout << "triangles_list: "<< triangles_list << std::endl;
    std::cout << "edges_list: "<< edges_list << std::endl;
    std::cout << "coordsEdge_list: "<< coordsEdge_list << std::endl;
    std::cout << "*********************************" << std::endl;
#endif

    if (pathOK)
    {
        // creating new declaration path:
        sofa::defaulttype::Vec<3,double> baryCoords;

        // 1 - First point a (for the moment: always a point in a triangle)
        if (last_point != BaseMeshTopology::InvalidID)
        {
            topoPath_list.push_back (core::topology::POINT);
            indices_list.push_back (last_point);
            const typename DataTypes::VecCoord& realC = *(this->object->getX());
            for (unsigned int i = 0; i<3; i++)
                baryCoords[i]=realC[last_point][i];
        }
        else
        {
            sofa::helper::vector< double > coefs_a = computeTriangleBarycoefs (ind_ta, a);
            topoPath_list.push_back (core::topology::TRIANGLE);
            indices_list.push_back (ind_ta);
            for (unsigned int i = 0; i<3; i++)
                baryCoords[i]=coefs_a[i];
        }

        coords_list.push_back (baryCoords);


        // 2 - All edges intersected (only edges for now)
        for (unsigned int i = 0; i< edges_list.size(); i++)
        {
            topoPath_list.push_back (core::topology::EDGE);
            indices_list.push_back (edges_list[i]);

            baryCoords[0] = coordsEdge_list[i];
            baryCoords[1] = 0.0; // or 1 - coordsEdge_list[i] ??
            baryCoords[2] = 0.0;

            coords_list.push_back (baryCoords);
        }

        // 3 - Last point b (for the moment: always a point in a triangle)
        sofa::helper::vector< double > coefs_b = computeTriangleBarycoefs (ind_tb, b);

        for (unsigned int i = 0; i<3; i++)
            if (coefs_b[i] > 0.9999 )
            {
                topoPath_list.push_back (core::topology::POINT);
                indices_list.push_back (this->m_topology->getTriangle (ind_tb)[i]);
                isOnPoint = true;
                break;
            }

        if (!isOnPoint)
        {
            topoPath_list.push_back (core::topology::TRIANGLE);
            indices_list.push_back (ind_tb);
        }
        for (unsigned int i = 0; i<3; i++)
            baryCoords[i]=coefs_b[i];

        coords_list.push_back (baryCoords);
    }

    return pathOK;
}





/// Get the triangle in a given direction from a point.
template <typename DataTypes>
int TriangleSetGeometryAlgorithms<DataTypes>::getTriangleInDirection(PointID p, const sofa::defaulttype::Vec<3,double>& dir) const
{
    const typename DataTypes::VecCoord& vect_c = *(this->object->getX());
    const sofa::helper::vector<TriangleID> &shell=this->m_topology->getTrianglesAroundVertex(p);
    sofa::defaulttype::Vec<3,Real> dtest = dir;
    for (unsigned int i=0; i<shell.size(); ++i)
    {
        unsigned int ind_t = shell[i];
        const Triangle &t=this->m_topology->getTriangle(ind_t);

        const typename DataTypes::Coord& c0=vect_c[t[0]];
        const typename DataTypes::Coord& c1=vect_c[t[1]];
        const typename DataTypes::Coord& c2=vect_c[t[2]];

        sofa::defaulttype::Vec<3,Real> p0;
        p0[0] = (Real) (c0[0]);
        p0[1] = (Real) (c0[1]);
        p0[2] = (Real) (c0[2]);
        sofa::defaulttype::Vec<3,Real> p1;
        p1[0] = (Real) (c1[0]);
        p1[1] = (Real) (c1[1]);
        p1[2] = (Real) (c1[2]);
        sofa::defaulttype::Vec<3,Real> p2;
        p2[0] = (Real) (c2[0]);
        p2[1] = (Real) (c2[1]);
        p2[2] = (Real) (c2[2]);

        sofa::defaulttype::Vec<3,Real> e1, e2;
        if (t[0] == p) { e1 = p1-p0; e2 = p2-p0; }
        else if (t[1] == p) { e1 = p2-p1; e2 = p0-p1; }
        else { e1 = p0-p2; e2 = p1-p2; }

        sofa::defaulttype::Vec<3,Real> v_normal = (e2).cross(e1);

        if(v_normal.norm2() > 1e-20)
        {
            sofa::defaulttype::Vec<3,Real> n_01 = e1.cross(v_normal);
            sofa::defaulttype::Vec<3,Real> n_02 = e2.cross(v_normal);

            double v_01 = (double) ((dtest)*(n_01));
            double v_02 = (double) ((dtest)*(n_02));

            bool is_inside = (v_01 >= 0.0) && (v_02 < 0.0);
            if (is_inside) return ind_t;
        }
    }
    return (TriangleID)-1;
}

/// Write the current mesh into a msh file
template <typename DataTypes>
void TriangleSetGeometryAlgorithms<DataTypes>::writeMSHfile(const char *filename) const
{
    std::ofstream myfile;
    myfile.open (filename);

    const typename DataTypes::VecCoord& vect_c = *(this->object->getX());

    const unsigned int numVertices = vect_c.size();

    myfile << "$NOD\n";
    myfile << numVertices <<"\n";

    for (unsigned int i=0; i<numVertices; ++i)
    {
        double x = (double) vect_c[i][0];
        double y = (double) vect_c[i][1];
        double z = (double) vect_c[i][2];

        myfile << i+1 << " " << x << " " << y << " " << z <<"\n";
    }

    myfile << "$ENDNOD\n";
    myfile << "$ELM\n";

    const sofa::helper::vector<Triangle> &ta=this->m_topology->getTriangles();

    myfile << ta.size() <<"\n";

    for (unsigned int i=0; i<ta.size(); ++i)
    {
        myfile << i+1 << " 2 6 6 3 " << ta[i][0]+1 << " " << ta[i][1]+1 << " " << ta[i][2]+1 <<"\n";
    }

    myfile << "$ENDELM\n";

    myfile.close();
}


template <typename DataTypes>
void TriangleSetGeometryAlgorithms<DataTypes>::reorderTrianglesOrientationFromNormals()
{
    sofa::defaulttype::Vec<3,Real> firstNormal = computeTriangleNormal(0);

    if (p_flipNormals.getValue())
        firstNormal = -firstNormal;

    sofa::helper::vector<TriangleID> _neighTri = this->m_topology->getElementAroundElement(0);
    sofa::helper::vector<TriangleID> _neighTri2, buffK, buffKK;
    unsigned int cpt_secu = 0, max = this->m_topology->getNbTriangles();
    sofa::defaulttype::Vec<3,Real> triNormal;
    bool pair = true;


    while (!_neighTri.empty() && cpt_secu < max)
    {
        for (unsigned int i=0; i<_neighTri.size(); ++i)
        {
            unsigned int triId = _neighTri[i];
            triNormal = this->computeTriangleNormal(triId);
            //std::cout << "triNormal: "<< triNormal << std::endl;
            double prod = (firstNormal*triNormal)/(firstNormal.norm()*triNormal.norm());
            //std::cout << "prod: "<< prod << std::endl;
            if (prod < 0.15) //change orientation
                this->m_topology->reOrientateTriangle(triId);
        }

        _neighTri2 = this->m_topology->getElementAroundElements(_neighTri);

        if (pair)
        {
            buffK = _neighTri;
            pair = false;

            _neighTri.clear();
            for (unsigned int i=0; i<_neighTri2.size(); ++i)
            {
                bool find = false;
                unsigned int id = _neighTri2[i];
                for (unsigned int j=0; j<buffKK.size(); ++j)
                    if (id == buffKK[j])
                    {
                        find = true;
                        break;
                    }

                if (!find)
                    _neighTri.push_back(id);
            }
        }
        else
        {
            buffKK = _neighTri;
            pair = true;

            _neighTri.clear();
            for (unsigned int i=0; i<_neighTri2.size(); ++i)
            {
                bool find = false;
                unsigned int id = _neighTri2[i];
                for (unsigned int j=0; j<buffK.size(); ++j)
                    if (id == buffK[j])
                    {
                        find = true;
                        break;
                    }

                if (!find)
                    _neighTri.push_back(id);
            }
        }

#ifndef NDEBUG
        std::cout << "_neighTri: "<< _neighTri << std::endl;
        std::cout << "_neighTri2: "<< _neighTri2 << std::endl;
        std::cout << "buffk: "<< buffK << std::endl;
        std::cout << "buffkk: "<< buffKK << std::endl;
#endif
        //_neighTri = _neighTri2;
        cpt_secu++;
    }

    if(cpt_secu == max)
        std::cerr << "WARNING: TriangleSetGeometryAlgorithms: reorder triangle orientation reach security end of loop." << std::endl;

    return;
}



template<class Real>
bool is_point_in_triangle(const sofa::defaulttype::Vec<3,Real>& p, const sofa::defaulttype::Vec<3,Real>& a, const sofa::defaulttype::Vec<3,Real>& b, const sofa::defaulttype::Vec<3,Real>& c)
{
    const double ZERO = 1e-6;

    sofa::defaulttype::Vec<3,Real> ptest = p;
    sofa::defaulttype::Vec<3,Real> p0 = a;
    sofa::defaulttype::Vec<3,Real> p1 = b;
    sofa::defaulttype::Vec<3,Real> p2 = c;

    sofa::defaulttype::Vec<3,Real> v_normal = (p2-p0).cross(p1-p0);

    Real norm_v_normal = v_normal*(v_normal);
    if(norm_v_normal > ZERO)
    {
        if(fabs((ptest-p0)*(v_normal)) < ZERO) // p is in the plane defined by the triangle (p0,p1,p2)
        {

            sofa::defaulttype::Vec<3,Real> n_01 = (p1-p0).cross(v_normal);
            sofa::defaulttype::Vec<3,Real> n_12 = (p2-p1).cross(v_normal);
            sofa::defaulttype::Vec<3,Real> n_20 = (p0-p2).cross(v_normal);

            return (((ptest-p0)*(n_01) > -ZERO) && ((ptest-p1)*(n_12) > -ZERO) && ((ptest-p2)*(n_20) > -ZERO));

        }
        else // p is not in the plane defined by the triangle (p0,p1,p2)
        {
            //sout << "INFO_print : p is not in the plane defined by the triangle (p0,p1,p2)" << sendl;
            return false;
        }

    }
    else // triangle is flat
    {
        //sout << "INFO_print : triangle is flat" << sendl;
        return false;
    }
}


/// Test if a point p is in the right halfplane

template<class Real>
bool is_point_in_halfplane(const sofa::defaulttype::Vec<3,Real>& p, unsigned int e0, unsigned int e1,
        const sofa::defaulttype::Vec<3,Real>& a, const sofa::defaulttype::Vec<3,Real>& b, const sofa::defaulttype::Vec<3,Real>& c,
        unsigned int ind_p0, unsigned int ind_p1, unsigned int ind_p2)
{
    sofa::defaulttype::Vec<3,Real> ptest = p;

    sofa::defaulttype::Vec<3,Real> p0 = a;
    sofa::defaulttype::Vec<3,Real> p1 = b;
    sofa::defaulttype::Vec<3,Real> p2 = c;

    sofa::defaulttype::Vec<3,Real> v_normal = (p2-p0).cross(p1-p0);

    Real norm_v_normal = (v_normal)*(v_normal);

    if(norm_v_normal != 0.0)
    {
        //v_normal/=norm_v_normal;

        if(ind_p0==e0 || ind_p0==e1)
        {
            sofa::defaulttype::Vec<3,Real> n_12 = (p2-p1).cross(v_normal);
            return ((ptest-p1)*(n_12) >= 0.0);
        }
        else
        {
            if(ind_p1==e0 || ind_p1==e1)
            {
                sofa::defaulttype::Vec<3,Real> n_20 = (p0-p2).cross(v_normal);
                return ((ptest-p2)*(n_20) >= 0.0);
            }
            else
            {
                if(ind_p2==e0 || ind_p2==e1)
                {
                    sofa::defaulttype::Vec<3,Real> n_01 = (p1-p0).cross(v_normal);
                    return ((ptest-p0)*(n_01) >= 0.0);
                }
                else
                {
                    return false; // not expected
                }
            }
        }
    }
    else // triangle is flat
    {
        //sout << "INFO_print : triangle is flat" << sendl;
        return false;
    }
}


template<class DataTypes>
void TriangleSetGeometryAlgorithms<DataTypes>::initPointAdded(unsigned int index, const core::topology::AncestorElem &ancestorElem
        , const helper::vector< VecCoord* >& coordVecs, const helper::vector< VecDeriv* >& derivVecs)
{
    using namespace sofa::core::topology;

    if (ancestorElem.type != TRIANGLE)
    {
        EdgeSetGeometryAlgorithms< DataTypes >::initPointAdded(index, ancestorElem, coordVecs, derivVecs);
    }
    else
    {
        const Triangle &t = this->m_topology->getTriangle(ancestorElem.index);

        for (unsigned int i = 0; i < coordVecs.size(); i++)
        {
            VecCoord &curVecCoord = *coordVecs[i];
            Coord& curCoord = curVecCoord[index];

            const Coord &c0 = curVecCoord[t[0]];
            const Coord &c1 = curVecCoord[t[1]];
            const Coord &c2 = curVecCoord[t[2]];

            // Compute normal (ugly but doesn't require template specialization...)
            defaulttype::Vec<3,Real> p0;
            p0[0] = (Real) (c0[0]);
            p0[1] = (Real) (c0[1]);
            p0[2] = (Real) (c0[2]);
            defaulttype::Vec<3,Real> p1;
            p1[0] = (Real) (c1[0]);
            p1[1] = (Real) (c1[1]);
            p1[2] = (Real) (c1[2]);
            defaulttype::Vec<3,Real> p2;
            p2[0] = (Real) (c2[0]);
            p2[1] = (Real) (c2[1]);
            p2[2] = (Real) (c2[2]);

            defaulttype::Vec<3,Real> p0p1 = p1 - p0;
            defaulttype::Vec<3,Real> p0p2 = p2 - p0;
            p0p1.normalize();
            p0p2.normalize();

            defaulttype::Vec<3,Real> n = p0p1.cross(p0p2);
            n.normalize();

            defaulttype::Vec<3,Real> newCurCoord = p0 + p0p1 * ancestorElem.localCoords[0] + p0p2 * ancestorElem.localCoords[1] + n * ancestorElem.localCoords[2];
            curCoord = Coord();
     
            DataTypes::add(curCoord, newCurCoord[0], newCurCoord[1], newCurCoord[2]);
        }
    }
}


template<class DataTypes>
void TriangleSetGeometryAlgorithms<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
#ifndef SOFA_NO_OPENGL
    EdgeSetGeometryAlgorithms<DataTypes>::draw(vparams);

    // Draw Triangles indices
    if (showTriangleIndices.getValue())
    {
        Mat<4,4, GLfloat> modelviewM;
        const VecCoord& coords = *(this->object->getX());
        const sofa::defaulttype::Vec3f& color = _drawColor.getValue();
        glColor3f(color[0]-0.2f, color[1]-0.2f, color[2]-0.2f);
        glDisable(GL_LIGHTING);
        float scale = this->getIndicesScale();

        //for triangles:
        scale = scale/2;

        const sofa::helper::vector<Triangle> &triangleArray = this->m_topology->getTriangles();

        for (unsigned int i =0; i<triangleArray.size(); i++)
        {

            Triangle the_tri = triangleArray[i];
            Coord vertex1 = coords[ the_tri[0] ];
            Coord vertex2 = coords[ the_tri[1] ];
            Coord vertex3 = coords[ the_tri[2] ];
            sofa::defaulttype::Vec3f center; center = (DataTypes::getCPos(vertex1)+DataTypes::getCPos(vertex2)+DataTypes::getCPos(vertex3))/3;

            std::ostringstream oss;
            oss << i;
            std::string tmp = oss.str();
            const char* s = tmp.c_str();
            glPushMatrix();

            glTranslatef(center[0], center[1], center[2]);
            glScalef(scale,scale,scale);

            // Makes text always face the viewer by removing the scene rotation
            // get the current modelview matrix
            glGetFloatv(GL_MODELVIEW_MATRIX , modelviewM.ptr() );
            modelviewM.transpose();

            sofa::defaulttype::Vec3f temp = modelviewM.transform(center);

            //glLoadMatrixf(modelview);
            glLoadIdentity();

            glTranslatef(temp[0], temp[1], temp[2]);
            glScalef(scale,scale,scale);

            while(*s)
            {
                glutStrokeCharacter(GLUT_STROKE_ROMAN, *s);
                s++;
            }

            glPopMatrix();

        }
    }


    // Draw Triangles
    if (_draw.getValue())
    {
        const sofa::helper::vector<Triangle> &triangleArray = this->m_topology->getTriangles();

        if (!triangleArray.empty()) // Draw triangle surfaces
        {
            const VecCoord& coords = *(this->object->getX());

            glDisable(GL_LIGHTING);
            const sofa::defaulttype::Vec3f& color = _drawColor.getValue();
            glColor3f(color[0], color[1], color[2]);
            glBegin(GL_TRIANGLES);
            for (unsigned int i = 0; i<triangleArray.size(); i++)
            {
                const Triangle& t = triangleArray[i];

                for (unsigned int j = 0; j<3; j++)
                {
                    sofa::defaulttype::Vec3f coordP; coordP = DataTypes::getCPos(coords[t[j]]);
                    glVertex3f(coordP[0], coordP[1], coordP[2]);
                }
            }
            glEnd();

            glColor3f(color[0]-0.2f, color[1]-0.2f, color[2]-0.2f);
            glBegin(GL_LINES);
            const sofa::helper::vector<Edge> &edgeArray = this->m_topology->getEdges();

            if (!edgeArray.empty()) // Draw triangle edges for better display
            {
                for (unsigned int i = 0; i<edgeArray.size(); i++)
                {
                    const Edge& e = edgeArray[i];
                    sofa::defaulttype::Vec3f coordP1; coordP1 = DataTypes::getCPos(coords[e[0]]);
                    sofa::defaulttype::Vec3f coordP2; coordP2 = DataTypes::getCPos(coords[e[1]]);
                    glVertex3f(coordP1[0], coordP1[1], coordP1[2]);
                    glVertex3f(coordP2[0], coordP2[1], coordP2[2]);
                }
            }
            else
            {
                for (unsigned int i = 0; i<triangleArray.size(); i++)
                {
                    const Triangle& t = triangleArray[i];
                    sofa::helper::vector <sofa::defaulttype::Vec3f> triCoord;

                    for (unsigned int j = 0; j<3; j++)
                    {
                        sofa::defaulttype::Vec3f p; p = DataTypes::getCPos(coords[t[j]]);
                        triCoord.push_back(p);
                    }

                    for (unsigned int j = 0; j<3; j++)
                    {
                        glVertex3f(triCoord[j][0], triCoord[j][1], triCoord[j][2]);
                        glVertex3f(triCoord[(j+1)%3][0], triCoord[(j+1)%3][1], triCoord[(j+1)%3][2]);
                    }
                }
            }
            glEnd();
        }
    }


    if (_drawNormals.getValue())
    {
        const VecCoord& coords = *(this->object->getX());
        const sofa::helper::vector<Triangle> &triangleArray = this->m_topology->getTriangles();
        unsigned int nbrTtri = triangleArray.size();

        Coord point2;
        sofa::defaulttype::Vec<3,double> colors;
        SReal normalLength = _drawNormalLength.getValue();

        glDisable(GL_LIGHTING);
        glBegin(GL_LINES);

        for (unsigned int i =0; i<nbrTtri; i++)
        {
            Triangle _tri = triangleArray[i];
            sofa::defaulttype::Vec<3,double> normal = this->computeTriangleNormal(i);
            normal.normalize();

            // compute bary triangle
            Coord vertex1 = coords[ _tri[0] ];
            Coord vertex2 = coords[ _tri[1] ];
            Coord vertex3 = coords[ _tri[2] ];
            sofa::defaulttype::Vec3d center; center = (DataTypes::getCPos(vertex1)+DataTypes::getCPos(vertex2)+DataTypes::getCPos(vertex3))/3;
            sofa::defaulttype::Vec3d point2 = center + normal*normalLength;

            for(unsigned int j=0; j<3; j++)
                colors[j] = fabs (normal[j]);

            glColor3f ((float)colors[0], (float)colors[1], (float)colors[2]);

            glVertex3d(center[0], center[1], center[2]);
            glVertex3d(point2[0], point2[1], point2[2]);
        }
        glEnd();
    }

#endif /* SOFA_NO_OPENGL */
}





} // namespace topology

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENTS_TriangleSetTOPOLOGY_INL
