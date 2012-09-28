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

#include <map>
#include <sofa/helper/gl/template.h>
#include <sofa/core/ObjectFactory.h>

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/behavior/BaseMechanicalState.h>
#include <sofa/core/topology/BaseMeshTopology.h>

#include <sofa/component/topology/SparseGridTopology.h>

#include <sofa/core/loader/VoxelLoader.h>

#include <sofa/component/visualmodel/SlicedVolumetricModel.h>
#include <sofa/core/visual/VisualParams.h>


#define GETCOORD(i) Coord((Real)_mstate->getPX(i), (Real)_mstate->getPY(i), (Real)_mstate->getPZ(i) )



namespace sofa
{

namespace component
{

namespace visualmodel
{

SOFA_DECL_CLASS(SlicedVolumetricModel)

int SlicedVolumetricModelClass = core::RegisterObject("A simple visualization for a cloud of points.")
        .add< SlicedVolumetricModel >()
        ;

using namespace sofa::defaulttype;
using namespace sofa::core::topology;



const int SlicedVolumetricModel::__edges__[12][2] = {{ 0,1 }, { 3,2 }, { 4,5 }, { 7,6 }, { 0,3 }, { 1,2 }, { 4,7 }, { 5,6 }, { 0,4 }, { 1,5 }, { 2,6 }, { 3,7 }};



SlicedVolumetricModel::SlicedVolumetricModel() //const std::string &name, std::string filename, std::string loader, std::string textureName)
    :
    alpha(initData(&alpha, 0.2f, "alpha", "Opacity of the billboards. 1.0 is 100% opaque.")),
    color(initData(&color, std::string("white"), "color", "Billboard color.")),
    _nbPlanes(initData(&_nbPlanes, 100, "nbSlices", "Number of billboards.")),
    _topology(NULL),
    _mstate(NULL),
    texture_data(NULL),
    _first(1)
{
}

SlicedVolumetricModel::~SlicedVolumetricModel()
{
    if(texture_data != NULL)
        delete [] texture_data;
}

void SlicedVolumetricModel::init()
{
    getContext()->get(_topology);
    if(_topology)
        _mstate = dynamic_cast<core::behavior::BaseMechanicalState*>(_topology->getContext()->getMechanicalState());
    else
        getContext()->get(_mstate);

// 	_topology->init();
    _mstate->init();

    VisualModel::init();

    core::loader::VoxelLoader *loader;
    getContext()->get(loader);
    if(loader)
    {
        loader->createSegmentation3DTexture( &texture_data, _width, _height, _depth );
    }


    if( topology::SparseGridTopology* sparseGrid = dynamic_cast<topology::SparseGridTopology*>(_topology ) )
    {
        _minBBox[0] = sparseGrid->getXmin();
        _minBBox[1] = sparseGrid->getYmin();
        _minBBox[2] = sparseGrid->getZmin();
        _maxBBox[0] = sparseGrid->getXmax();
        _maxBBox[1] = sparseGrid->getYmax();
        _maxBBox[2] = sparseGrid->getZmax();
    }
    else
    {
        _minBBox[0]=_minBBox[1]=_minBBox[2]=999999999;
        _maxBBox[0]=_maxBBox[1]=_maxBBox[2]=-999999999;

    }

    _nbPlanesOld = _nbPlanes.getValue();

    const Coord& p0 = GETCOORD(_topology->getHexahedron(0)[0]);
    const Coord& p7 = GETCOORD(_topology->getHexahedron(0)[6]);
    _radius = (p7-p0).norm() / 2;




    _textureCoordinates.resize( _mstate->getSize() );
    for( int i=0; i<_mstate->getSize(); ++i)
    {
        const Coord& p = GETCOORD( i );
        _textureCoordinates[i][0] = (Real)((p[0]- _minBBox[0]) / (_maxBBox[0] - _minBBox[0]));
        _textureCoordinates[i][1] = (Real)((p[1]- _minBBox[1]) / (_maxBBox[1] - _minBBox[1]));
        _textureCoordinates[i][2] = (Real)((p[2]- _minBBox[2]) / (_maxBBox[2] - _minBBox[2]));
    }

    reinit();

    updateVisual();
}

void SlicedVolumetricModel::reinit()
{
    setColor(color.getValue());

    if( _nbPlanesOld != _nbPlanes.getValue() || _first )
    {
        // 	if( _nbPlanes.getValue()>2048)_nbPlanes.setValue(2048);
        alpha.setValue((alpha.getValue()*Real(_nbPlanesOld))/Real(_nbPlanes.getValue()));
        _planeSeparations = (Real)((_maxBBox[0]-_minBBox[0]) / (Real)_nbPlanes.getValue());
// 		cerr<<"_planeSeparations : "<<_planeSeparations<<endl;
        _nbPlanesOld = _nbPlanes.getValue();
    }

}

void SlicedVolumetricModel::drawTransparent(const core::visual::VisualParams* vparams)
{
    if(!vparams->displayFlags().getShowVisualModels()) return;

    if( _first )
    {

        _first = false;

#ifdef SOFA_HAVE_GLEW
        glewInit();
#endif

//   	// set up our OpenGL state
// 	glDisable(GL_DEPTH_TEST);
// 		glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE); // our texture colors will replace the untextured colors
        glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);





        // request 1 texture name from OpenGL
// 		glGenTextures(numInstance, &_texname);
        glGenTextures(1, &_texname);
        // tell OpenGL we're going to be setting up the texture name it gave us
        glBindTexture(GL_TEXTURE_3D, _texname);
        // when this texture needs to be shrunk to fit on small polygons, use linear interpolation of the texels to determine the color
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        // when this texture needs to be magnified to fit on a big polygon, use linear interpolation of the texels to determine the color
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        // we want the texture to repeat over the S axis, so if we specify coordinates out of range we still get textured.
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP);
        // same as above for T axis
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP);
        // same as above for R axis
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP);
        // this is a 3d texture, level 0 (max detail), GL should store it in RGB8 format, its WIDTHxHEIGHTxDEPTH in size,
        // it doesnt have a border, we're giving it to GL in RGB format as a series of unsigned bytes, and texels is where the texel data is.
        glTexImage3D(GL_TEXTURE_3D, 0, GL_ALPHA, _width, _height, _depth, 0, GL_ALPHA, GL_UNSIGNED_BYTE, texture_data);



        delete [] texture_data;
        texture_data = NULL;
        return;
    }


// 	glPushAttrib(GL_ALL_ATTRIB_BITS);


    glDisable(GL_LIGHTING);
    glPolygonMode (GL_FRONT,GL_FILL );


    glEnable(GL_BLEND);
// 	glDisable(GL_DEPTH_TEST);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);







    float mat[16];
    glGetFloatv( GL_MODELVIEW_MATRIX, mat );
    vRight=Coord( mat[0], mat[4], mat[8] );
    vUp=Coord( mat[1], mat[5], mat[9] );
    _planeNormal = vRight.cross( vUp);
    _planeNormal.normalize();

    glBindTexture(GL_TEXTURE_3D, _texname);


    glColor4f( r,g,b, alpha.getValue());

    glEnable(GL_TEXTURE_3D);


    glBegin(GL_TRIANGLES);
    findAndDrawTriangles();
    glEnd();

    glDisable(GL_TEXTURE_3D);


}





void SlicedVolumetricModel::findAndDrawTriangles()
{
    int actualPlane=0;


    Coord lastPoint;
    Real maxLastPoint = (Real)999999999;



    for(int i = 0 ; i < _mstate->getSize(); ++i )
    {
        Coord p = GETCOORD( i );


        Real actualLastPoint = _planeNormal * p;
        if( actualLastPoint < maxLastPoint )
        {
            maxLastPoint = actualLastPoint;
            lastPoint = p;
        }
    }


    lastPoint += _planeNormal * .1;

    std::list<int>positiveCubes;
    for(int i=0; i<_topology->getNbHexahedra(); ++i)
        positiveCubes.push_back( i );


    int nbintersections;
    do
    {

        nbintersections = 0;

        // trouver le centre du plan de coupe
        Coord planeCenter = lastPoint + _planeNormal * (actualPlane * _planeSeparations);
        Real planeConstant = _planeNormal * planeCenter;

        EdgesMap _edgesMap;

// 		for(Octree::CellPtrList::iterator itcell=_octree->getLeafLists(animal::octree::GEOMETRY).begin();   itcell!=_octree->getLeafLists(animal::octree::GEOMETRY).end();itcell++)

        // seulement les nouveaux cubes potentiellement intersectable, ie proches ou devant le plan
        for(std::list<int>::iterator itcell=positiveCubes.begin(); itcell!=positiveCubes.end(); /*++itcell*/)
        {
            const BaseMeshTopology::Hexa& cell = _topology->getHexahedron( *itcell );

            Coord cubebarycenter = GETCOORD( cell[0] ) + (GETCOORD( cell[6] ) - GETCOORD( cell[0] ) ) / 2.0;

            Real dist = (_planeNormal * cubebarycenter) - planeConstant; //distance du centre du cube au plan


            if( fabs(dist) >= _radius)
            {
                if( dist>0 ) // du bon cote mais plus loin, on garde pour plus tard
                {
                    ++itcell;
                    continue;
                }
                else // pas du bon cote, on oublie le cube
                {
                    std::list<int>::iterator it = positiveCubes.erase( itcell );
                    itcell = it;
                    continue;
                }
            }


            // find intersections
            helper::vector<Intersection> intersections;
            for(int i=0; i<12; ++i)
            {
                int e0 = __edges__[i][0];
                int e1 = __edges__[i][1];
                Coord s0 = GETCOORD( cell[e0] );
                Coord s1 = GETCOORD( cell[e1] );
                Edge e(cell[e0],cell[e1]);



                EdgesMap::iterator em = _edgesMap.find( e );
                if(  em != _edgesMap.end() )
                {
                    intersections.push_back( (*em).second );
                }
                else
                {
                    Coord dir = s1-s0;
                    Coord dirnormalized=dir;
                    dirnormalized.normalize();

                    Real where;
                    int howmany = intersectionSegmentPlane( s0,s1, dirnormalized, _planeNormal, planeConstant, where );
                    if( howmany == 1 )
                    {

                        Coord w = s0 + dir * where;


                        Coord st0 ( _textureCoordinates[cell[e0]] );
                        Coord dir2 = _textureCoordinates[cell[e1]]-st0;


                        Intersection inter( w, st0 + dir2 * where);
                        intersections.push_back( inter );

                        _edgesMap[e]=inter;
                    }
                    else if(howmany==2)
                    {
// 						cerr<<"intersect une ligne entiere"<<sendl;

                        Intersection inter( s0, _textureCoordinates[cell[e0]]);
                        intersections.push_back( inter );
                        inter = Intersection( s1, _textureCoordinates[cell[e1]]);
                        intersections.push_back( inter );
                    }


                }
            }

// 			cerr<<"intersections.size() : "<<intersections.size()<<endl;

            if( intersections.size() <2 )
            {
                ++itcell;
                continue;
            }
// 			else cerr<<"pas assez inter"<<sendl;



            nbintersections += intersections.size();












            // trier les intersections

            helper::vector<std::pair<Real,int> > neg; // angle + indice
            helper::vector<std::pair<Real,int> > pos;
            helper::vector<int> nul;



// 			Coord middle(0,0,0);
// 			for(unsigned int i=0;i<intersections.size();++i)
// 				middle+=intersections[i].first;
// 			animal::v_teq(middle,1.0/intersections.size());
// 			Coord referenceLine = intersections[0].first - middle;


            Coord referenceLine = intersections[1].first - intersections[0].first;



            Coord referenceLine2( referenceLine[1],- referenceLine[0], 0);
// 		// 	Coord referenceLine2( referenceLine[1], referenceLine[0], (-2*referenceLine[0]*referenceLine[1])/referenceLine[2]);
// 		// 	Coord tmp,referenceLine2;
// 		// 	animal::v_eq_cross( tmp,referenceLine, intersections[2] - intersections[0]);
// 		// 	animal::v_eq_cross( referenceLine2,referenceLine, tmp); // est-ce que la line2 a besoin d'etre dans le plan ???
// 		// 	cerr<<"---\n"<<sendl;



            for(unsigned int i=2; i<intersections.size(); ++i) // les cas 0 et 1 sont trait�s � la mano
            {
                Coord actualline = intersections[i].first-intersections[0].first;

                Real angle1 = referenceLine * actualline;
                Real angle2 = referenceLine2 * actualline ;

                if( angle2<0.0)
                    neg.push_back( std::pair<Real,int>(angle1, i) );
                else
                    pos.push_back( std::pair<Real,int>(angle1, i) );

                // 		cerr<<i<<" : "<<angle1<<" "<<angle2<<endl;
            }


// 			for(unsigned int i=1;i<intersections.size();++i) // les cas 0 et 1 sont trait�s � la mano
// 			{
// 				Coord actualline = intersections[i].first-middle;
            //
// 				Real angle1 = animal::v_dot( referenceLine,  actualline);
// 				Real angle2 = animal::v_dot( referenceLine2, actualline );
            //
// 				if( angle2<0.0)
// 					neg.push_back( std::pair<Real,int>(angle1, i) );
// 				else
// 					pos.push_back( std::pair<Real,int>(angle1, i) );
            //
// 		// 		cerr<<i<<" : "<<angle1<<" "<<angle2<<endl;
// 			}




            stable_sort( pos.begin(),pos.end());
            stable_sort( neg.begin(),neg.end());

            helper::vector<int> tripoints;
// 			tripoints.push_back(0);

            glPointSize(30.0);
            for(unsigned  int i=0; i<pos.size(); ++i)
            {
                tripoints.push_back(pos[i].second);
            }

            tripoints.push_back(1);

            for( int i=neg.size()-1; i>=0; --i)
            {
                tripoints.push_back(neg[i].second);
            }





            for( unsigned int i=0; i<tripoints.size()-1; ++i)
            {
                helper::gl::glTexCoordT(intersections[0].second);
                helper::gl::glVertexT(intersections[0].first);
                helper::gl::glTexCoordT(intersections[tripoints[i]].second);
                helper::gl::glVertexT(intersections[tripoints[i]].first);
                helper::gl::glTexCoordT(intersections[tripoints[i+1]].second);
                helper::gl::glVertexT(intersections[tripoints[i+1]].first);
            }







// 			for( unsigned int i=0;i<tripoints.size()-1;++i)
// 			{
// 				glTexCoord3fv(middle);
// 				glVertex3fv(middle);
// 				glTexCoord3fv(intersections[tripoints[i]].second);
// 				glVertex3fv(intersections[tripoints[i]].first);
// 				glTexCoord3fv(intersections[tripoints[i+1]].second);
// 				glVertex3fv(intersections[tripoints[i+1]].first);
// 			}


            ++itcell;
        }

        if(actualPlane>_nbPlanes.getValue()*.9 && !nbintersections)break;

        ++actualPlane;

    }
    while( true );

}



/// return 0->no intersection, 1->1 intersection, 2->line on plane
int SlicedVolumetricModel::intersectionSegmentPlane( const Coord&s0,const Coord&s1,  const Coord &segmentDirection, const Coord& planeNormal, const Real&planeConstant,Real & m_fLineT /*where is the intersection on the segment*/)
{
// 	++_debugNbInteresctionComputations;

    Real fDdN = segmentDirection * planeNormal;
    Real fSDistance = (planeNormal * s0) - planeConstant;

    if (fabs(fDdN) > 1.0e-5)
    {
        // The line is not parallel to the plane, so they must intersect.
        m_fLineT = -fSDistance/fDdN;

        // The line intersects the plane, but possibly at a point that is
        // not on the segment.
        Real norm = (s1-s0).norm();
        if( m_fLineT>0 && fabs(m_fLineT) <= norm  )
        {
            m_fLineT /= norm;
            return 1;
        }
        else
            return 0;

    }

    // The Line and plane are parallel.  Determine if they are numerically
    // close enough to be coincident.
    else if (fabs(fSDistance) <= 1.0e-5)
    {
        // The line is coincident with the plane, so choose t = 0 for the
        // parameter.
        m_fLineT = (Real)0.0;
        return 2;
    }

    else return 0;
}




void SlicedVolumetricModel::setColor(float r, float g, float b)
{
    this->r = r;
    this->g = g;
    this->b = b;
}

static int hexval(char c)
{
    if (c>='0' && c<='9') return c-'0';
    else if (c>='a' && c<='f') return (c-'a')+10;
    else if (c>='A' && c<='F') return (c-'A')+10;
    else return 0;
}

void SlicedVolumetricModel::setColor(std::string color)
{
    if (color.empty()) return;
    float r = 1.0f;
    float g = 1.0f;
    float b = 1.0f;
    if (color[0]>='0' && color[0]<='9')
    {
        sscanf(color.c_str(),"%f %f %f", &r, &g, &b);
    }
    else if (color[0]=='#' && color.length()>=7)
    {
        r = (hexval(color[1])*16+hexval(color[2]))/255.0f;
        g = (hexval(color[3])*16+hexval(color[4]))/255.0f;
        b = (hexval(color[5])*16+hexval(color[6]))/255.0f;
    }
    else if (color[0]=='#' && color.length()>=4)
    {
        r = (hexval(color[1])*17)/255.0f;
        g = (hexval(color[2])*17)/255.0f;
        b = (hexval(color[3])*17)/255.0f;
    }
    else if (color == "white")    { r = 1.0f; g = 1.0f; b = 1.0f; }
    else if (color == "black")    { r = 0.0f; g = 0.0f; b = 0.0f; }
    else if (color == "red")      { r = 1.0f; g = 0.0f; b = 0.0f; }
    else if (color == "green")    { r = 0.0f; g = 1.0f; b = 0.0f; }
    else if (color == "blue")     { r = 0.0f; g = 0.0f; b = 1.0f; }
    else if (color == "cyan")     { r = 0.0f; g = 1.0f; b = 1.0f; }
    else if (color == "magenta")  { r = 1.0f; g = 0.0f; b = 1.0f; }
    else if (color == "yellow")   { r = 1.0f; g = 1.0f; b = 0.0f; }
    else if (color == "gray")     { r = 0.5f; g = 0.5f; b = 0.5f; }
    else
    {
        serr << "Unknown color "<<color<<sendl;
        return;
    }
    setColor(r,g,b);
}

} // namespace visualmodel

} // namespace component

} // namespace sofa

