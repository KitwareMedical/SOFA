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
#include <sofa/core/ObjectFactory.h>
#include <sofa/component/loader/GridMeshCreator.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/system/SetDirectory.h>

namespace sofa
{

namespace component
{

namespace loader
{

using namespace sofa::defaulttype;
using namespace sofa::core::loader;

SOFA_DECL_CLASS(GridMeshCreator)

int GridMeshCreatorClass = core::RegisterObject("Procedural creation of a two-dimensional mesh.")
        .add< GridMeshCreator >()
        ;



GridMeshCreator::GridMeshCreator(): MeshLoader()
    , resolution( initData(&resolution,Vec2i(2,2),"resolution","Number of vertices in each direction"))
    , trianglePattern( initData(&trianglePattern,2,"trianglePattern","0: no triangles, 1: alternate triangles, 2: upward triangles, 3: downward triangles"))
{
}


void GridMeshCreator::insertUniqueEdge( unsigned a, unsigned b )
{
    if( uniqueEdges.find(Edge(b,a))==uniqueEdges.end() ) // symmetric not found
        uniqueEdges.insert(Edge(a,b));                   // redundant elements are automatically pruned
}

void GridMeshCreator::insertTriangle(unsigned a, unsigned b, unsigned c)
{
    helper::vector<Triangle >& my_triangles = *(triangles.beginEdit());

    my_triangles.push_back(Triangle( a,b,c ) );
    insertUniqueEdge(a,b);
    insertUniqueEdge(b,c);
    insertUniqueEdge(c,a);

    triangles.endEdit();
}

void GridMeshCreator::insertQuad(unsigned a, unsigned b, unsigned c, unsigned d)
{
    helper::vector<Quad >& my_quads = *(quads.beginEdit());

    my_quads.push_back( Quad( a,b,c,d ) );
    insertUniqueEdge(a,b);
    insertUniqueEdge(b,c);
    insertUniqueEdge(c,d);
    insertUniqueEdge(d,a);

    quads.endEdit();
}


bool GridMeshCreator::load()
{
    helper::WriteAccessor<Data<vector<sofa::defaulttype::Vector3> > > my_positions (positions);
    unsigned numX = resolution.getValue()[0], numY=resolution.getValue()[1];

    // Warning: Vertex creation order must be consistent with method vert.
    for(unsigned y=0; y<numY; y++)
    {
        for(unsigned x=0; x<numX; x++)
        {
            my_positions.push_back( Vector3(x * 1./(numX-1), y * 1./(numY-1), 0) );
            //            cerr<<"GridMeshCreator::load, add point " << Vector3(i * 1./(numX-1), j * 1./(numY-1), 0) << endl;
        }
    }

    uniqueEdges.clear();

    if( trianglePattern.getValue()==0 ) // quads
        for(unsigned y=0; y<numY-1; y++ )
        {
            for(unsigned x=0; x<numX-1; x++ )
            {
                insertQuad( vert(x,y), vert(x+1,y), vert(x+1,y+1), vert(x,y+1) );
            }
        }
    else if( trianglePattern.getValue()==1 ) // alternate
        for(unsigned y=0; y<numY-1; y++ )
        {
            for(unsigned x=0; x<numX-1; x++ )
            {
                if( (x+y)%2 == 0 )
                {
                    insertTriangle( vert(x,y), vert(x+1,y), vert(x+1,y+1) );
                    insertTriangle( vert(x,y), vert(x+1,y+1), vert(x,y+1)   ) ;
                }
                else
                {
                    insertTriangle( vert(x  ,y), vert(x+1,y)  , vert(x,y+1) ) ;
                    insertTriangle( vert(x+1,y), vert(x+1,y+1), vert(x,y+1)    ) ;
                }
            }
        }
    else if( trianglePattern.getValue()==2 ) // upward
        for(unsigned y=0; y<numY-1; y++ )
        {
            for(unsigned x=0; x<numX-1; x++ )
            {
                insertTriangle( vert(x,y), vert(x+1,y), vert(x+1,y+1) ) ;
                insertTriangle( vert(x,y), vert(x+1,y+1), vert(x,y+1)   ) ;
            }
        }
    else if( trianglePattern.getValue()==3 ) // downward
        for(unsigned y=0; y<numY-1; y++ )
        {
            for(unsigned x=0; x<numX-1; x++ )
            {
                insertTriangle( vert(x  ,y), vert(x+1,y)  , vert(x,y+1) ) ;
                insertTriangle( vert(x+1,y), vert(x+1,y+1), vert(x,y+1)    ) ;
            }
        }

    helper::vector<Edge >& my_edges = *(edges.beginEdit());
    for( std::set<Edge>::const_iterator i=uniqueEdges.begin(),iend=uniqueEdges.end(); i!=iend; i++ )
        my_edges.push_back( *i );
    edges.endEdit();

    return true;

}


}
}
}

