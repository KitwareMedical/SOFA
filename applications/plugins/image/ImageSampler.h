/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 4      *
*                (c) 2006-2009 MGH, INRIA, USTL, UJF, CNRS                    *
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
#ifndef SOFA_IMAGE_IMAGESAMPLER_H
#define SOFA_IMAGE_IMAGESAMPLER_H

#include "initImage.h"
#include "ImageTypes.h"
#include "ImageAlgorithms.h"
#include <sofa/core/DataEngine.h>
#include <sofa/component/component.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/visual/VisualParams.h>

#include <sofa/core/objectmodel/Event.h>
#include <sofa/simulation/common/AnimateEndEvent.h>

#include <sofa/defaulttype/Vec.h>
#include <sofa/helper/OptionsGroup.h>

#include "BranchingImage.h"

#define REGULAR 0
#define LLOYD 1


namespace sofa
{

namespace component
{

namespace engine
{

using helper::vector;
using defaulttype::Vec;
using defaulttype::Mat;
using cimg_library::CImg;




/// Default implementation does not compile
template <int imageTypeLabel>
struct ImageSamplerSpecialization
{
};


/// Specialization for regular Image
template <>
struct ImageSamplerSpecialization<defaulttype::IMAGELABEL_IMAGE>
{

    template<class ImageSampler>
    static void init( ImageSampler* )
    {
    }

    template<class ImageSampler>
    static void regularSampling( ImageSampler* sampler, const bool atcorners=false, const bool recursive=false )
    {
        typedef typename ImageSampler::Real Real;
        typedef typename ImageSampler::Coord Coord;
        typedef typename ImageSampler::Edge Edge;
        typedef typename ImageSampler::Hexa Hexa;
        typedef typename ImageSampler::T T;


        // get tranform and image at time t
        typename ImageSampler::raImage in(sampler->image);
        typename ImageSampler::raTransform inT(sampler->transform);
        const CImg<T>& inimg = in->getCImg(sampler->time);

        // data access
        typename ImageSampler::waPositions pos(sampler->position);       pos.clear();
        typename ImageSampler::waEdges e(sampler->edges);                e.clear();
        typename ImageSampler::waEdges g(sampler->graphEdges);           g.clear();
        typename ImageSampler::waHexa h(sampler->hexahedra);             h.clear();

        // convert to single channel boolean image
        CImg<bool> img(inimg.width(),inimg.height(),inimg.depth(),1,false);
        if(atcorners)
        {
            CImg_2x2x2(I,bool);
            cimg_for2x2x2(inimg,x,y,z,0,I,bool) if(Iccc || Iccn || Icnc || Incc || Innc || Incn || Icnn || Innn) img(x,y,z)=true;
        }
        else cimg_forXYZC(inimg,x,y,z,c) if(inimg(x,y,z,c)) img(x,y,z)=true;

        // count non empty voxels
        unsigned int nb=0;
        cimg_foroff(img,off) if(img[off]) nb++;
        pos.resize(nb);
        // record indices of previous y line and z plane for connectivity
        CImg<unsigned int> pLine(inimg.width()),nLine(inimg.width());
        CImg<unsigned int> pPlane(inimg.width(),inimg.height()),nPlane(inimg.width(),inimg.height());
        // fill pos and edges
        nb=0;
        cimg_forZ(img,z)
        {
            cimg_forY(img,y)
            {
                cimg_forX(img,x)
                {
                    if(img(x,y,z))
                    {
                        // pos
                        if(atcorners) pos[nb]=Coord(x+0.5,y+0.5,z+0.5);
                        else pos[nb]=Coord(x,y,z);
                        // edges
                        if(x) if(img(x-1,y,z)) e.push_back(Edge(nb-1,nb));
                        if(y) if(img(x,y-1,z)) e.push_back(Edge(pLine(x),nb));
                        if(z) if(img(x,y,z-1)) e.push_back(Edge(pPlane(x,y),nb));
                        // hexa
                        if(x && y && z) if(img(x-1,y,z) && img(x,y-1,z) && img(x,y,z-1) && img(x-1,y-1,z) && img(x-1,y,z-1)  && img(x,y-1,z-1)   && img(x-1,y-1,z-1) )
                                h.push_back(Hexa(nb,pLine(x),pLine(x-1),nb-1,pPlane(x,y),pPlane(x,y-1),pPlane(x-1,y-1),pPlane(x-1,y) ));

                        nLine(x)=nb; nPlane(x,y)=nb;
                        nb++;
                    }
                }
                nLine.swap(pLine);
            }
            nPlane.swap(pPlane);
        }

        if(recursive)
        {
            vector<unsigned int> indices; indices.resize(pos.size()); for(unsigned int i=0; i<pos.size(); i++) indices[i]=i;
            sampler->subdivide(indices);
        }

        for(unsigned int i=0; i<pos.size(); i++) pos[i]=inT->fromImage(pos[i]);
    }


    template<class ImageSampler>
    static void uniformSampling( ImageSampler* sampler,const unsigned int nb=0,  const bool bias=false, const unsigned int lloydIt=100,const bool useDijkstra=false )
    {
        typedef typename ImageSampler::Real Real;
        typedef typename ImageSampler::Coord Coord;
        typedef typename ImageSampler::Edge Edge;
        typedef typename ImageSampler::Hexa Hexa;
        typedef typename ImageSampler::T T;

        clock_t timer = clock();

        // get tranform and image at time t
        typename ImageSampler::raImage in(sampler->image);
        typename ImageSampler::raTransform inT(sampler->transform);
        const CImg<T>& inimg = in->getCImg(sampler->time);
        const CImg<T>* biasFactor=bias?&inimg:NULL;

        // data access
        typename ImageSampler::raPositions fpos(sampler->fixedPosition);
        std::vector<Vec<3,Real> >& pos = *sampler->position.beginEdit();    pos.clear();
        typename ImageSampler::waEdges e(sampler->edges);                e.clear();
        typename ImageSampler::waEdges g(sampler->graphEdges);           g.clear();
        typename ImageSampler::waHexa h(sampler->hexahedra);             h.clear();

        // init voronoi and distances
        CImg<unsigned int>  voronoi(inimg.width(),inimg.height(),inimg.depth(),1,0);
        typename ImageSampler::waDist distData(sampler->distances);
        typename ImageSampler::imCoord dim = in->getDimensions(); dim[3]=dim[4]=1; distData->setDimensions(dim);
        CImg<Real>& dist = distData->getCImg(); dist.fill(-1);
        cimg_forXYZC(inimg,x,y,z,c) if(inimg(x,y,z,c)) dist(x,y,z)=cimg::type<Real>::max();

        // list of seed points
        std::set<std::pair<Real,sofa::defaulttype::Vec<3,int> > > trial;

        // farthest point sampling using geodesic distances
        vector<unsigned int> fpos_voronoiIndex;
        vector<unsigned int> pos_voronoiIndex;
        for(unsigned int i=0; i<fpos.size(); i++)
        {
            fpos_voronoiIndex.push_back(i+1);
            AddSeedPoint<Real>(trial,dist,voronoi, sampler->transform.getValue(), fpos[i],fpos_voronoiIndex[i]);
        }
        while(pos.size()<nb)
        {
            Real dmax=0;  Coord pmax;
            cimg_forXYZ(dist,x,y,z) if(dist(x,y,z)>dmax) { dmax=dist(x,y,z); pmax =Coord(x,y,z); }
            if(dmax)
            {
                pos_voronoiIndex.push_back(fpos.size()+pos.size()+1);
                pos.push_back(inT->fromImage(pmax));
                AddSeedPoint<Real>(trial,dist,voronoi, sampler->transform.getValue(), pos.back(),pos_voronoiIndex.back());
                if(useDijkstra) dijkstra<Real,T>(trial,dist, voronoi, sampler->transform.getValue().getScale(), biasFactor);
                else fastMarching<Real,T>(trial,dist, voronoi, sampler->transform.getValue().getScale(),biasFactor );
            }
            else break;
        }
        //voronoi.display();

        unsigned int it=0;
        bool converged =(it>=lloydIt)?true:false;

        while(!converged)
        {
            if(Lloyd<Real,T>(pos,pos_voronoiIndex,dist,voronoi,sampler->transform.getValue(),NULL)) // one lloyd iteration
            {
                // recompute distance from scratch
                cimg_foroff(dist,off) if(dist[off]!=-1) dist[off]=cimg::type<Real>::max();
                for(unsigned int i=0; i<fpos.size(); i++) AddSeedPoint<Real>(trial,dist,voronoi, sampler->transform.getValue(), fpos[i], fpos_voronoiIndex[i]);
                for(unsigned int i=0; i<pos.size(); i++) AddSeedPoint<Real>(trial,dist,voronoi, sampler->transform.getValue(), pos[i], pos_voronoiIndex[i]);
                if(useDijkstra) dijkstra<Real,T>(trial,dist, voronoi,  sampler->transform.getValue().getScale(), biasFactor);
                else fastMarching<Real,T>(trial,dist, voronoi,  sampler->transform.getValue().getScale(), biasFactor);
                it++; if(it>=lloydIt) converged=true;
            }
            else converged=true;
        }

        if(sampler->f_printLog.getValue())
        {
            std::cout<<"ImageSampler: Completed in "<< it <<" Lloyd iterations ("<< (clock() - timer) / (float)CLOCKS_PER_SEC <<"s )"<<std::endl;
        }

        sampler->position.endEdit();
    }


    template<class ImageSampler>
    static void recursiveUniformSampling( ImageSampler* sampler,const unsigned int nb=0,  const bool bias=false, const unsigned int lloydIt=100,const bool useDijkstra=false,  const unsigned int N=1 )
    {
        typedef typename ImageSampler::Real Real;
        typedef typename ImageSampler::Coord Coord;
        typedef typename ImageSampler::Edge Edge;
        typedef typename ImageSampler::Hexa Hexa;
        typedef typename ImageSampler::T T;

        clock_t timer = clock();

        // get tranform and image at time t
        typename ImageSampler::raImage in(sampler->image);
        typename ImageSampler::raTransform inT(sampler->transform);
        const CImg<T>& inimg = in->getCImg(sampler->time);
        const CImg<T>* biasFactor=bias?&inimg:NULL;

        // data access
        typename ImageSampler::raPositions fpos(sampler->fixedPosition);
        std::vector<Vec<3,Real> >& pos = *sampler->position.beginEdit();    pos.clear();
        typename ImageSampler::waEdges e(sampler->edges);                e.clear();
        typename ImageSampler::waEdges g(sampler->graphEdges);           g.clear();
        typename ImageSampler::waHexa h(sampler->hexahedra);             h.clear();

        // init voronoi and distances
        CImg<unsigned int>  voronoi(inimg.width(),inimg.height(),inimg.depth(),1,0);
        typename ImageSampler::waDist distData(sampler->distances);
        typename ImageSampler::imCoord dim = in->getDimensions(); dim[3]=dim[4]=1; distData->setDimensions(dim);
        CImg<Real>& dist = distData->getCImg(); dist.fill(-1);
        cimg_forXYZC(inimg,x,y,z,c) if(inimg(x,y,z,c)) dist(x,y,z)=cimg::type<Real>::max();

        // list of seed points
        std::set<std::pair<Real,sofa::defaulttype::Vec<3,int> > > trial;

        // fixed points
        vector<unsigned int> fpos_voronoiIndex;
        for(unsigned int i=0; i<fpos.size(); i++)
        {
            fpos_voronoiIndex.push_back(i+1);
            AddSeedPoint<Real>(trial,dist,voronoi, sampler->transform.getValue(), fpos[i],fpos_voronoiIndex[i]);
        }

        // new points
        vector<unsigned int> pos_voronoiIndex;
        while(pos.size()<nb)
        {
            std::vector<Vec<3,Real> > newpos;
            vector<unsigned int> newpos_voronoiIndex;

            // farthest sampling of N points
            unsigned int currentN = N;
            if(!pos.size()) currentN = 1; // special case at the beginning: we start by adding just one point
            else if(pos.size()+N>nb) currentN = nb-pos.size();  // when trying to add more vertices than necessary
            while(newpos.size()<currentN)
            {
                Real dmax=0;  Coord pmax;
                cimg_forXYZ(dist,x,y,z) if(dist(x,y,z)>dmax) { dmax=dist(x,y,z); pmax =Coord(x,y,z); }
                if(!dmax) break;

                newpos_voronoiIndex.push_back(fpos.size()+pos.size()+newpos.size()+1);
                newpos.push_back(inT->fromImage(pmax));
                AddSeedPoint<Real>(trial,dist,voronoi, sampler->transform.getValue(), newpos.back(),newpos_voronoiIndex.back());
                if(useDijkstra) dijkstra<Real,T>(trial,dist, voronoi, sampler->transform.getValue().getScale(), biasFactor);
                else fastMarching<Real,T>(trial,dist, voronoi, sampler->transform.getValue().getScale(),biasFactor );
            }

            // lloyd iterations for the N points
            unsigned int it=0;
            bool converged =(it>=lloydIt)?true:false;

            while(!converged)
            {
                if(Lloyd<Real,T>(newpos,newpos_voronoiIndex,dist,voronoi,sampler->transform.getValue(),NULL))
                {
                    // recompute distance from scratch
                    cimg_foroff(dist,off) if(dist[off]!=-1) dist[off]=cimg::type<Real>::max();
                    for(unsigned int i=0; i<fpos.size(); i++) AddSeedPoint<Real>(trial,dist,voronoi, sampler->transform.getValue(), fpos[i], fpos_voronoiIndex[i]);
                    for(unsigned int i=0; i<pos.size(); i++) AddSeedPoint<Real>(trial,dist,voronoi, sampler->transform.getValue(), pos[i], pos_voronoiIndex[i]);
                    for(unsigned int i=0; i<newpos.size(); i++) AddSeedPoint<Real>(trial,dist,voronoi, sampler->transform.getValue(), newpos[i], newpos_voronoiIndex[i]);
                    if(useDijkstra) dijkstra<Real,T>(trial,dist, voronoi,  sampler->transform.getValue().getScale(), biasFactor);
                    else fastMarching<Real,T>(trial,dist, voronoi,  sampler->transform.getValue().getScale(), biasFactor);
                    it++; if(it>=lloydIt) converged=true;
                }
                else converged=true;
            }

            // check neighbors of the new voronoi cell and add graph edges
            unsigned int nbold = fpos.size()+pos.size();
            for(unsigned int i=0; i<newpos.size() && pos.size()<nb; i++)
            {
                std::set<unsigned int> neighb;
                CImg_3x3x3(I,unsigned int);
                cimg_for3x3x3(voronoi,x,y,z,0,I,unsigned int)
                if(Iccc==newpos_voronoiIndex[i])
                {
                    if(Incc && Incc<=nbold) neighb.insert(Incc);
                    if(Icnc && Icnc<=nbold) neighb.insert(Icnc);
                    if(Iccn && Iccn<=nbold) neighb.insert(Iccn);
                    if(Ipcc && Ipcc<=nbold) neighb.insert(Ipcc);
                    if(Icpc && Icpc<=nbold) neighb.insert(Icpc);
                    if(Iccp && Iccp<=nbold) neighb.insert(Iccp);
                }
                for(typename std::set<unsigned int>::iterator itr=neighb.begin(); itr!=neighb.end(); itr++)
                {
                    g.push_back(Edge(*itr-1,newpos_voronoiIndex[i]-1));
                    //if(*itr>fpos.size()) g.push_back(Edge(*itr-fpos.size()-1,newpos_voronoiIndex[i]-1));
                }
                pos.push_back(newpos[i]);
                pos_voronoiIndex.push_back(newpos_voronoiIndex[i]);
            }

            if(newpos.size()<currentN) break; // check possible failure in point insertion (not enough voxels)
        }

        if(sampler->f_printLog.getValue())
        {
            std::cout<<"ImageSampler: Completed in "<< (clock() - timer) / (float)CLOCKS_PER_SEC <<"s "<<std::endl;
        }

        sampler->position.endEdit();
    }
};



/// Specialization for BranchingImage
template <>
struct ImageSamplerSpecialization<defaulttype::IMAGELABEL_BRANCHINGIMAGE>
{
    template<class ImageSampler>
    static void init( ImageSampler* sampler )
    {
        sampler->addAlias( &sampler->image, "branchingImage" );
    }

    template<class ImageSampler>
    static void regularSampling( ImageSampler* sampler, const bool atcorners=false, const bool recursive=false )
    {
        if( !atcorners )
        {
            sampler->serr<<"ImageSampler::regularSampling - only at corner is implemented\n";
        }


        typedef typename ImageSampler::Real Real;
        typedef typename ImageSampler::Coord Coord;
        typedef typename ImageSampler::Edge Edge;
        typedef typename ImageSampler::Hexa Hexa;
        typedef typename ImageSampler::T T;


        // get tranform and image at time t
        typename ImageSampler::raImage in(sampler->image);
        typename ImageSampler::raTransform inT(sampler->transform);
        const typename ImageSampler::ImageTypes::BranchingImage3D& inimg = in->imgList[sampler->time];

        // data access
        typename ImageSampler::waPositions pos(sampler->position);       pos.clear();
        typename ImageSampler::waEdges e(sampler->edges);                e.clear();
        typename ImageSampler::waEdges g(sampler->graphEdges);           g.clear();
        typename ImageSampler::waHexa h(sampler->hexahedra);             h.clear();

        const typename ImageSampler::ImageTypes::Dimension& dim = in->getDimension();


        unsigned index1d = 0;

        {
        std::map< unsigned, std::map<unsigned, unsigned> > hindices; // for each index1d, for each superimposed offset -> hexa index

        // first add hexa with linked vertex
        unsigned indexVertex = 0;
        for( unsigned z=0 ; z<dim[ImageSampler::ImageTypes::DIMENSION_Z] ; ++z )
        for( unsigned y=0 ; y<dim[ImageSampler::ImageTypes::DIMENSION_Y] ; ++y )
        for( unsigned x=0 ; x<dim[ImageSampler::ImageTypes::DIMENSION_X] ; ++x )
        {
            for( unsigned v=0 ; v<inimg[index1d].size() ; ++v )
            {
                h.push_back( Hexa( indexVertex, indexVertex+1, indexVertex+2, indexVertex+3, indexVertex+4, indexVertex+5, indexVertex+6, indexVertex+7 ) );
                indexVertex += 8;
                Hexa& hexa = h[h.size()-1];
                hindices[index1d][v] = h.size()-1;

                if( x>0 )
                {
                    for( unsigned n=0 ; n<inimg[index1d][v].neighbours[ImageSampler::ImageTypes::LEFT].size() ; ++n )
                    {
                        const unsigned neighbourIndex = in->getNeighbourIndex(ImageSampler::ImageTypes::LEFT,index1d);
                        const unsigned neighbourOffset = inimg[index1d][v].neighbours[ImageSampler::ImageTypes::LEFT][n];
                        Hexa& neighbor = h[hindices[neighbourIndex][neighbourOffset]];
                        mergeVertexIndex( h, hexa[0], neighbor[1] );
                        mergeVertexIndex( h, hexa[4], neighbor[5] );
                        mergeVertexIndex( h, hexa[7], neighbor[6] );
                        mergeVertexIndex( h, hexa[3], neighbor[2] );
                    }
                }
                if( y>0 )
                {
                    for( unsigned n=0 ; n<inimg[index1d][v].neighbours[ImageSampler::ImageTypes::BOTTOM].size() ; ++n )
                    {
                        const unsigned neighbourIndex = in->getNeighbourIndex(ImageSampler::ImageTypes::BOTTOM,index1d);
                        const unsigned neighbourOffset = inimg[index1d][v].neighbours[ImageSampler::ImageTypes::BOTTOM][n];
                        Hexa& neighbor = h[hindices[neighbourIndex][neighbourOffset]];
                        mergeVertexIndex( h, hexa[0], neighbor[3] );
                        mergeVertexIndex( h, hexa[1], neighbor[2] );
                        mergeVertexIndex( h, hexa[4], neighbor[7] );
                        mergeVertexIndex( h, hexa[5], neighbor[6] );
                    }
                }
                if( z>0 )
                {
                    for( unsigned n=0 ; n<inimg[index1d][v].neighbours[ImageSampler::ImageTypes::BACK].size() ; ++n )
                    {
                        const unsigned neighbourIndex = in->getNeighbourIndex(ImageSampler::ImageTypes::BACK,index1d);
                        const unsigned neighbourOffset = inimg[index1d][v].neighbours[ImageSampler::ImageTypes::BACK][n];
                        Hexa& neighbor = h[hindices[neighbourIndex][neighbourOffset]];
                        mergeVertexIndex( h, hexa[0], neighbor[4] );
                        mergeVertexIndex( h, hexa[1], neighbor[5] );
                        mergeVertexIndex( h, hexa[2], neighbor[6] );
                        mergeVertexIndex( h, hexa[3], neighbor[7] );
                    }
                }
            }
            index1d++;
        }
        }


        {
            // give a continue index from 0 to max (without hole)
            unsigned continueIndex = 0;
            std::map<unsigned,unsigned> continueMap;

            for( unsigned i=0 ; i<h.size() ; ++i )
            for( unsigned j=0 ; j<8 ; ++j )
            {
                if( continueMap.find(h[i][j])==continueMap.end() ) continueMap[h[i][j]]=continueIndex++;
            }
            for( unsigned i=0 ; i<h.size() ; ++i )
            for( unsigned j=0 ; j<8 ; ++j )
            {
                h[i][j] = continueMap[h[i][j]];
            }
            pos.resize( continueIndex );
        }

        {
            std::map<unsigned,bool> alreadyAddedPos;
            std::map<Edge,bool> alreadyAddedEdge;
            index1d = 0;
            unsigned indexHexa = 0;
            for( unsigned z=0 ; z<dim[ImageSampler::ImageTypes::DIMENSION_Z] ; ++z )
            for( unsigned y=0 ; y<dim[ImageSampler::ImageTypes::DIMENSION_Y] ; ++y )
            for( unsigned x=0 ; x<dim[ImageSampler::ImageTypes::DIMENSION_X] ; ++x )
            {
                for( unsigned v=0 ; v<inimg[index1d].size() ; ++v )
                {

                    static const int hexaIndices[8][3] = { {0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}  };

                    Hexa& hexa = h[indexHexa];
                    for( unsigned j=0 ; j<8 ; ++j )
                    {
                        if( alreadyAddedPos.find(hexa[j])==alreadyAddedPos.end() )
                        {
                            alreadyAddedPos[hexa[j]] = true;
                            pos[hexa[j]] = Coord(x+hexaIndices[j][0],y+hexaIndices[j][1],z+hexaIndices[j][2]);
                        }
                    }

                    Edge edge( hexa[0], hexa[1] );
                    if( alreadyAddedEdge.find( edge )==alreadyAddedEdge.end() )
                    {
                        alreadyAddedEdge[edge] = true;
                        e.push_back( edge );
                    }
                    edge = Edge( hexa[1], hexa[2] );
                    if( alreadyAddedEdge.find( edge )==alreadyAddedEdge.end() )
                    {
                        alreadyAddedEdge[edge] = true;
                        e.push_back( edge );
                    }
                    edge = Edge( hexa[2], hexa[3] );
                    if( alreadyAddedEdge.find( edge )==alreadyAddedEdge.end() )
                    {
                        alreadyAddedEdge[edge] = true;
                        e.push_back( edge );
                    }
                    edge = Edge( hexa[3], hexa[0] );
                    if( alreadyAddedEdge.find( edge )==alreadyAddedEdge.end() )
                    {
                        alreadyAddedEdge[edge] = true;
                        e.push_back( edge );
                    }
                    edge = Edge( hexa[4], hexa[5] );
                    if( alreadyAddedEdge.find( edge )==alreadyAddedEdge.end() )
                    {
                        alreadyAddedEdge[edge] = true;
                        e.push_back( edge );
                    }
                    edge = Edge( hexa[5], hexa[6] );
                    if( alreadyAddedEdge.find( edge )==alreadyAddedEdge.end() )
                    {
                        alreadyAddedEdge[edge] = true;
                        e.push_back( edge );
                    }
                    edge = Edge( hexa[6], hexa[7] );
                    if( alreadyAddedEdge.find( edge )==alreadyAddedEdge.end() )
                    {
                        alreadyAddedEdge[edge] = true;
                        e.push_back( edge );
                    }
                    edge = Edge( hexa[7], hexa[4] );
                    if( alreadyAddedEdge.find( edge )==alreadyAddedEdge.end() )
                    {
                        alreadyAddedEdge[edge] = true;
                        e.push_back( edge );
                    }
                    edge = Edge( hexa[0], hexa[4] );
                    if( alreadyAddedEdge.find( edge )==alreadyAddedEdge.end() )
                    {
                        alreadyAddedEdge[edge] = true;
                        e.push_back( edge );
                    }
                    edge = Edge( hexa[1], hexa[5] );
                    if( alreadyAddedEdge.find( edge )==alreadyAddedEdge.end() )
                    {
                        alreadyAddedEdge[edge] = true;
                        e.push_back( edge );
                    }
                    edge = Edge( hexa[2], hexa[6] );
                    if( alreadyAddedEdge.find( edge )==alreadyAddedEdge.end() )
                    {
                        alreadyAddedEdge[edge] = true;
                        e.push_back( edge );
                    }
                    edge = Edge( hexa[3], hexa[7] );
                    if( alreadyAddedEdge.find( edge )==alreadyAddedEdge.end() )
                    {
                        alreadyAddedEdge[edge] = true;
                        e.push_back( edge );
                    }

                    ++indexHexa;
                }

                ++index1d;
            }
        }


        if(recursive)
        {
            vector<unsigned int> indices; indices.resize(pos.size()); for(unsigned int i=0; i<pos.size(); i++) indices[i]=i;
            sampler->subdivide(indices);
        }

        for(unsigned int i=0; i<pos.size(); i++) pos[i]=inT->fromImage(pos[i]);
    }

    template<class Hexas>
    static void mergeVertexIndex( Hexas& h, unsigned index0, unsigned index1 )
    {
        for( unsigned i=0 ; i<h.size() ; ++i )
        for( unsigned j=0 ; j<8 ; ++j )
        {
            if( h[i][j]==index1 ) h[i][j]=index0;
        }
    }

    template<class ImageSampler>
    static void uniformSampling( ImageSampler* sampler,const unsigned int /*nb*/=0,  const bool /*bias*/=false, const unsigned int /*lloydIt*/=100,const bool /*useDijkstra*/=false )
    {
        sampler->serr<<"ImageSampler::uniformSampling is not yet immplemented for BranchingImage\n";
    }

    template<class ImageSampler>
    static void recursiveUniformSampling( ImageSampler* sampler,const unsigned int /*nb*/=0,  const bool /*bias*/=false, const unsigned int /*lloydIt*/=100,const bool /*useDijkstra*/=false,  const unsigned int /*N*/=1 )
    {
        sampler->serr<<"ImageSampler::recursiveUniformSampling is not yet immplemented for BranchingImage\n";
    }
};





/**
 * This class samples an object represented by an image
 */


template <class _ImageTypes>
class ImageSampler : public core::DataEngine
{
    friend struct ImageSamplerSpecialization<defaulttype::IMAGELABEL_IMAGE>;
    friend struct ImageSamplerSpecialization<defaulttype::IMAGELABEL_BRANCHINGIMAGE>;

public:

    typedef core::DataEngine Inherited;
    SOFA_CLASS(SOFA_TEMPLATE(ImageSampler,_ImageTypes),Inherited);

    typedef SReal Real;

    //@name Image data
    /**@{*/
    typedef _ImageTypes ImageTypes;
    typedef typename ImageTypes::T T;
    typedef typename ImageTypes::imCoord imCoord;
    typedef helper::ReadAccessor<Data< ImageTypes > > raImage;
    Data< ImageTypes > image;
    /**@}*/

    //@name Transform data
    /**@{*/
    typedef defaulttype::ImageLPTransform<Real> TransformType;
    typedef typename TransformType::Coord Coord;
    typedef helper::ReadAccessor<Data< TransformType > > raTransform;
    Data< TransformType > transform;
    /**@}*/

    //@name option data
    /**@{*/
    typedef vector<double> ParamTypes;
    typedef helper::ReadAccessor<Data< ParamTypes > > raParam;

    Data<helper::OptionsGroup> method;
    Data< bool > computeRecursive;
    Data< ParamTypes > param;
    /**@}*/

    //@name sample data (points+connectivity)
    /**@{*/
    typedef vector<Vec<3,Real> > SeqPositions;
    typedef helper::ReadAccessor<Data< SeqPositions > > raPositions;
    typedef helper::WriteAccessor<Data< SeqPositions > > waPositions;
    Data< SeqPositions > position;
    Data< SeqPositions > fixedPosition;

    typedef typename core::topology::BaseMeshTopology::Edge Edge;
    typedef typename core::topology::BaseMeshTopology::SeqEdges SeqEdges;
    typedef helper::ReadAccessor<Data< SeqEdges > > raEdges;
    typedef helper::WriteAccessor<Data< SeqEdges > > waEdges;
    Data< SeqEdges > edges;
    Data< SeqEdges > graphEdges;

    typedef typename core::topology::BaseMeshTopology::Hexa Hexa;
    typedef typename core::topology::BaseMeshTopology::SeqHexahedra SeqHexahedra;
    typedef helper::ReadAccessor<Data< SeqHexahedra > > raHexa;
    typedef helper::WriteAccessor<Data< SeqHexahedra > > waHexa;
    Data< SeqHexahedra > hexahedra;
    /**@}*/

    //@name distances (may be used for shape function computation)
    /**@{*/
    typedef defaulttype::Image<Real> DistTypes;
    typedef helper::ReadAccessor<Data< DistTypes > > raDist;
    typedef helper::WriteAccessor<Data< DistTypes > > waDist;
    Data< DistTypes > distances;
    /**@}*/


    //@name visu data
    /**@{*/
    Data< bool > showSamples;
    Data< bool > showEdges;
    Data< bool > showGraph;
    /**@}*/

    virtual std::string getTemplateName() const    { return templateName(this);    }
    static std::string templateName(const ImageSampler<ImageTypes>* = NULL) { return ImageTypes::Name();    }
    ImageSampler()    :   Inherited()
        , image(initData(&image,ImageTypes(),"image",""))
        , transform(initData(&transform,TransformType(),"transform",""))
        , method ( initData ( &method,"method","method (param)" ) )
        , computeRecursive(initData(&computeRecursive,false,"computeRecursive","if true: insert nodes recursively and build the graph"))
        , param ( initData ( &param,"param","Parameters" ) )
        , position(initData(&position,SeqPositions(),"position","output positions"))
        , fixedPosition(initData(&fixedPosition,SeqPositions(),"fixedPosition","user defined sample positions"))
        , edges(initData(&edges,SeqEdges(),"edges","edges connecting neighboring nodes"))
        , graphEdges(initData(&graphEdges,SeqEdges(),"graphEdges","oriented graph connecting parent to child nodes"))
        , hexahedra(initData(&hexahedra,SeqHexahedra(),"hexahedra","output hexahedra"))
        , distances(initData(&distances,DistTypes(),"distances",""))
        , showSamples(initData(&showSamples,false,"showSamples","show samples"))
        , showEdges(initData(&showEdges,false,"showEdges","show edges"))
        , showGraph(initData(&showGraph,false,"showGraph","show graph"))
        , time((unsigned int)0)
    {
        image.setReadOnly(true);
        transform.setReadOnly(true);
        f_listening.setValue(true);

        helper::OptionsGroup methodOptions(2,"0 - Regular sampling (at voxel center(0) or corners (1)) "
                ,"1 - Uniform sampling using Fast Marching and Lloyd relaxation (nbSamples | bias distances=false | nbiterations=100  | FastMarching(0)/Dijkstra(1)=1)"
                                          );
        methodOptions.setSelectedItem(REGULAR);
        method.setValue(methodOptions);

        ImageSamplerSpecialization<ImageTypes::label>::init( this );
    }

    virtual void init()
    {
        addInput(&image);
        addInput(&transform);
        addInput(&fixedPosition);
        addOutput(&position);
        addOutput(&edges);
        addOutput(&graphEdges);
        addOutput(&hexahedra);
        addOutput(&distances);
        setDirtyValue();
    }

    virtual void reinit() { update(); }

protected:

    unsigned int time;

    virtual void update()
    {
        cleanDirty();

        raParam params(this->param);

        if(this->method.getValue().getSelectedId() == REGULAR)
        {
            // get params
            bool atcorners=false; if(params.size())   atcorners=(bool)params[0];

            // sampling
            regularSampling(atcorners, computeRecursive.getValue());
        }
        else if(this->method.getValue().getSelectedId() == LLOYD)
        {
            // get params
            unsigned int nb=0;        if(params.size())       nb=(unsigned int)params[0];
            bool bias=false;          if(params.size()>1)     bias=(bool)params[1];
            unsigned int lloydIt=100; if(params.size()>2)     lloydIt=(unsigned int)params[2];
            bool Dij=true;            if(params.size()>3)     Dij=(bool)params[3];
            unsigned int N=1;        //curently does not work// if(params.size()>4)     N=(unsigned int)params[4];

            // sampling
            if(!computeRecursive.getValue()) uniformSampling(nb,bias,lloydIt,Dij);
            else recursiveUniformSampling(nb,bias,lloydIt,Dij,N);
        }

        if(this->f_printLog.getValue())
        {
            if(this->position.getValue().size())    std::cout<<"ImageSampler: "<< this->position.getValue().size() <<" generated samples"<<std::endl;
            if(this->edges.getValue().size())       std::cout<<"ImageSampler: "<< this->edges.getValue().size() <<" generated edges"<<std::endl;
            if(this->hexahedra.getValue().size())   std::cout<<"ImageSampler: "<< this->hexahedra.getValue().size() <<" generated hexahedra"<<std::endl;
            if(this->graphEdges.getValue().size())       std::cout<<"ImageSampler: "<< this->graphEdges.getValue().size() <<" generated dependencies"<<std::endl;
        }
    }

    void handleEvent(sofa::core::objectmodel::Event *event)
    {
        if ( dynamic_cast<simulation::AnimateEndEvent*>(event))
        {
            raImage in(this->image);
            raTransform inT(this->transform);

            // get current time modulo dimt
            const unsigned int dimt=in->getDimensions()[4];
            if(!dimt) return;
            Real t=inT->toImage(this->getContext()->getTime()) ;
            t-=(Real)((int)((int)t/dimt)*dimt);
            t=(t-floor(t)>0.5)?ceil(t):floor(t); // nearest
            if(t<0) t=0.0; else if(t>=(Real)dimt) t=(Real)dimt-1.0; // clamp

            if(this->time!=(unsigned int)t) { this->time=(unsigned int)t; update(); }
        }
    }

    virtual void draw(const core::visual::VisualParams* vparams)
    {
        if (!vparams->displayFlags().getShowVisualModels()) return;

        raPositions pos(this->position);
        raPositions fpos(this->fixedPosition);
        raEdges e(this->edges);
        raEdges g(this->graphEdges);

        if (this->showSamples.getValue()) vparams->drawTool()->drawPoints(this->position.getValue(),5.0,defaulttype::Vec4f(0.2,1,0.2,1));
        if (this->showSamples.getValue()) vparams->drawTool()->drawPoints(this->fixedPosition.getValue(),7.0,defaulttype::Vec4f(1,0.2,0.2,1));
        if (this->showEdges.getValue())
        {
            std::vector<defaulttype::Vector3> points;
            points.resize(2*e.size());
            for (unsigned int i=0; i<e.size(); ++i)
            {
                points[2*i][0]=pos[e[i][0]][0];            points[2*i][1]=pos[e[i][0]][1];            points[2*i][2]=pos[e[i][0]][2];
                points[2*i+1][0]=pos[e[i][1]][0];          points[2*i+1][1]=pos[e[i][1]][1];          points[2*i+1][2]=pos[e[i][1]][2];
            }
            vparams->drawTool()->drawLines(points,2.0,defaulttype::Vec4f(0.7,1,0.7,1));
        }
        if (this->showGraph.getValue())
        {
            std::vector<defaulttype::Vector3> points;
            points.resize(2*g.size());
            for (unsigned int i=0; i<g.size(); ++i)
                for (unsigned int j=0; j<2; ++j)
                {
                    if(g[i][j]<fpos.size()) {points[2*i+j][0]=fpos[g[i][j]][0];            points[2*i+j][1]=fpos[g[i][j]][1];            points[2*i+j][2]=fpos[g[i][j]][2];}
                    else {points[2*i+j][0]=pos[g[i][j]-fpos.size()][0];            points[2*i+j][1]=pos[g[i][j]-fpos.size()][1];            points[2*i+j][2]=pos[g[i][j]-fpos.size()][2];}

                }
            vparams->drawTool()->drawLines(points,2.0,defaulttype::Vec4f(1,1,0.5,1));
        }
    }


    /**
    * put regularly spaced samples at each non empty voxel center or corners
    * generated topology: edges + hexahedra
    * @param atcorners : put samples at voxel corners instead of centers ?
    * if @param buildgraph = true, several resolutions are recursively built; a graph is generated relating each higher resolution node to its parent nodes
    */
    void regularSampling ( const bool atcorners=false , const bool recursive=false )
    {
        ImageSamplerSpecialization<ImageTypes::label>::regularSampling( this, atcorners, recursive );
    }


    /// subdivide positions indexed in indices in eight sub-lists, add new points in this->position and run recursively
    void subdivide(vector<unsigned int> &indices)
    {
        waPositions pos(this->position);
        waEdges g(this->graphEdges);
        unsigned int nb=indices.size();

        // detect leaf
        if(nb<=(unsigned int)8) return;

        // computes center and bounding box
        typedef std::pair<Real,unsigned int> distanceToPoint;
        typedef std::set<distanceToPoint> distanceSet;

        Coord C; Coord BB[2];
        for(unsigned int dir=0; dir<3; dir++)
        {
            distanceSet q;
            for(unsigned int i=0; i<nb; i++) {unsigned int index=indices[i]; q.insert(distanceToPoint(pos[index][dir],i));}
            typename distanceSet::iterator it=q.begin();
            BB[0][dir]=q.begin()->first; BB[1][dir]=q.rbegin()->first;
            C[dir]=(BB[1][dir]+BB[0][dir])*0.5;   while(it->first<C[dir]) it++;       // mean
            // for(unsigned int count=0; count<nb/2; count++) it++;   // median
            Real c=it->first;
            it--; if(C[dir]-it->first<c-C[dir]) c=it->first;
            C[dir]=c;
            //            std::cout<<"dir="<<dir<<":"; for( it=q.begin(); it!=q.end(); it++)  std::cout<<it->first <<" "; std::cout<<std::endl; std::cout<<"C="<<C[dir]<<std::endl;
        }
        //        for(unsigned int i=0;i<nb;i++) std::cout<<"("<<pos[indices[i]]<<") ";  std::cout<<std::endl;
        Coord p;
        typename vector<Coord>::iterator it;
        // add corners
        unsigned int corners[8]= {addPoint(Coord(BB[0][0],BB[0][1],BB[0][2]),pos,indices),addPoint(Coord(BB[1][0],BB[0][1],BB[0][2]),pos,indices),addPoint(Coord(BB[0][0],BB[1][1],BB[0][2]),pos,indices),addPoint(Coord(BB[1][0],BB[1][1],BB[0][2]),pos,indices),addPoint(Coord(BB[0][0],BB[0][1],BB[1][2]),pos,indices),addPoint(Coord(BB[1][0],BB[0][1],BB[1][2]),pos,indices),addPoint(Coord(BB[0][0],BB[1][1],BB[1][2]),pos,indices),addPoint(Coord(BB[1][0],BB[1][1],BB[1][2]),pos,indices)};
        // add cell center
        unsigned int center=addPoint(Coord(C[0],C[1],C[2]),pos,indices);
        // add face centers
        unsigned int faces[6]= {addPoint(Coord(BB[0][0],C[1],C[2]),pos,indices),addPoint(Coord(BB[1][0],C[1],C[2]),pos,indices),addPoint(Coord(C[0],BB[0][1],C[2]),pos,indices),addPoint(Coord(C[0],BB[1][1],C[2]),pos,indices),addPoint(Coord(C[0],C[1],BB[0][2]),pos,indices),addPoint(Coord(C[0],C[1],BB[1][2]),pos,indices)};
        // add edge centers
        unsigned int edgs[12]= {addPoint(Coord(C[0],BB[0][1],BB[0][2]),pos,indices),addPoint(Coord(C[0],BB[1][1],BB[0][2]),pos,indices),addPoint(Coord(C[0],BB[0][1],BB[1][2]),pos,indices),addPoint(Coord(C[0],BB[1][1],BB[1][2]),pos,indices),addPoint(Coord(BB[0][0],C[1],BB[0][2]),pos,indices),addPoint(Coord(BB[1][0],C[1],BB[0][2]),pos,indices),addPoint(Coord(BB[0][0],C[1],BB[1][2]),pos,indices),addPoint(Coord(BB[1][0],C[1],BB[1][2]),pos,indices),addPoint(Coord(BB[0][0],BB[0][1],C[2]),pos,indices),addPoint(Coord(BB[1][0],BB[0][1],C[2]),pos,indices),addPoint(Coord(BB[0][0],BB[1][1],C[2]),pos,indices),addPoint(Coord(BB[1][0],BB[1][1],C[2]),pos,indices)};
        // connect
        bool connect=true;
        for(unsigned int i=0; i<6; i++) if(center==faces[i]) connect=false; for(unsigned int i=0; i<12; i++) if(center==edgs[i]) connect=false;
        if(connect) for(unsigned int i=0; i<8; i++) addEdge(Edge(corners[i],center),g);
        connect=true; for(unsigned int i=0; i<12; i++) if(faces[0]==edgs[i]) connect=false; if(connect) { addEdge(Edge(corners[0],faces[0]),g); addEdge(Edge(corners[2],faces[0]),g); addEdge(Edge(corners[4],faces[0]),g); addEdge(Edge(corners[6],faces[0]),g); }
        connect=true; for(unsigned int i=0; i<12; i++) if(faces[1]==edgs[i]) connect=false; if(connect) { addEdge(Edge(corners[1],faces[1]),g); addEdge(Edge(corners[3],faces[1]),g); addEdge(Edge(corners[5],faces[1]),g); addEdge(Edge(corners[7],faces[1]),g); }
        connect=true; for(unsigned int i=0; i<12; i++) if(faces[2]==edgs[i]) connect=false; if(connect) { addEdge(Edge(corners[0],faces[2]),g); addEdge(Edge(corners[1],faces[2]),g); addEdge(Edge(corners[4],faces[2]),g); addEdge(Edge(corners[5],faces[2]),g); }
        connect=true; for(unsigned int i=0; i<12; i++) if(faces[3]==edgs[i]) connect=false; if(connect) { addEdge(Edge(corners[2],faces[3]),g); addEdge(Edge(corners[3],faces[3]),g); addEdge(Edge(corners[6],faces[3]),g); addEdge(Edge(corners[7],faces[3]),g); }
        connect=true; for(unsigned int i=0; i<12; i++) if(faces[4]==edgs[i]) connect=false; if(connect) { addEdge(Edge(corners[0],faces[4]),g); addEdge(Edge(corners[1],faces[4]),g); addEdge(Edge(corners[2],faces[4]),g); addEdge(Edge(corners[3],faces[4]),g); }
        connect=true; for(unsigned int i=0; i<12; i++) if(faces[5]==edgs[i]) connect=false; if(connect) { addEdge(Edge(corners[4],faces[5]),g); addEdge(Edge(corners[5],faces[5]),g); addEdge(Edge(corners[6],faces[5]),g); addEdge(Edge(corners[7],faces[5]),g); }
        if(edgs[0]!=corners[0] && edgs[0]!=corners[1]) {addEdge(Edge(corners[0],edgs[0]),g); addEdge(Edge(corners[1],edgs[0]),g);}
        if(edgs[1]!=corners[2] && edgs[1]!=corners[3]) {addEdge(Edge(corners[2],edgs[1]),g); addEdge(Edge(corners[3],edgs[1]),g);}
        if(edgs[2]!=corners[4] && edgs[2]!=corners[5]) {addEdge(Edge(corners[4],edgs[2]),g); addEdge(Edge(corners[5],edgs[2]),g);}
        if(edgs[3]!=corners[6] && edgs[3]!=corners[7]) {addEdge(Edge(corners[6],edgs[3]),g); addEdge(Edge(corners[7],edgs[3]),g);}
        if(edgs[4]!=corners[0] && edgs[4]!=corners[2]) {addEdge(Edge(corners[0],edgs[4]),g); addEdge(Edge(corners[2],edgs[4]),g);}
        if(edgs[5]!=corners[1] && edgs[5]!=corners[3]) {addEdge(Edge(corners[1],edgs[5]),g); addEdge(Edge(corners[3],edgs[5]),g);}
        if(edgs[6]!=corners[4] && edgs[6]!=corners[6]) {addEdge(Edge(corners[4],edgs[6]),g); addEdge(Edge(corners[6],edgs[6]),g);}
        if(edgs[7]!=corners[5] && edgs[7]!=corners[7]) {addEdge(Edge(corners[5],edgs[7]),g); addEdge(Edge(corners[7],edgs[7]),g);}
        if(edgs[8]!=corners[0] && edgs[8]!=corners[4]) {addEdge(Edge(corners[0],edgs[8]),g); addEdge(Edge(corners[4],edgs[8]),g);}
        if(edgs[9]!=corners[1] && edgs[9]!=corners[5]) {addEdge(Edge(corners[1],edgs[9]),g); addEdge(Edge(corners[5],edgs[9]),g);}
        if(edgs[10]!=corners[2] && edgs[10]!=corners[6]) {addEdge(Edge(corners[2],edgs[10]),g); addEdge(Edge(corners[6],edgs[10]),g);}
        if(edgs[11]!=corners[3] && edgs[11]!=corners[7]) {addEdge(Edge(corners[3],edgs[11]),g); addEdge(Edge(corners[7],edgs[11]),g);}

        // check in which octant lies each point
        vector<vector<unsigned int> > octant(8); for(unsigned int i=0; i<8; i++)  octant[i].reserve(nb);
        for(unsigned int i=0; i<indices.size(); i++)
        {
            unsigned int index=indices[i];
            if(pos[index][0]<=C[0] && pos[index][1]<=C[1] && pos[index][2]<=C[2]) octant[0].push_back(index);
            if(pos[index][0]>=C[0] && pos[index][1]<=C[1] && pos[index][2]<=C[2]) octant[1].push_back(index);
            if(pos[index][0]<=C[0] && pos[index][1]>=C[1] && pos[index][2]<=C[2]) octant[2].push_back(index);
            if(pos[index][0]>=C[0] && pos[index][1]>=C[1] && pos[index][2]<=C[2]) octant[3].push_back(index);
            if(pos[index][0]<=C[0] && pos[index][1]<=C[1] && pos[index][2]>=C[2]) octant[4].push_back(index);
            if(pos[index][0]>=C[0] && pos[index][1]<=C[1] && pos[index][2]>=C[2]) octant[5].push_back(index);
            if(pos[index][0]<=C[0] && pos[index][1]>=C[1] && pos[index][2]>=C[2]) octant[6].push_back(index);
            if(pos[index][0]>=C[0] && pos[index][1]>=C[1] && pos[index][2]>=C[2]) octant[7].push_back(index);
        }
        for(unsigned int i=0; i<8; i++)
        {
            // std::cout<<i<<" : "<<octant[i]<<std::endl;
            subdivide(octant[i]);
        }
    }

    // add point p in pos if not already there are return its index
    unsigned int addPoint(const Coord p, waPositions& pos, vector<unsigned int> &indices)
    {
        unsigned int ret ;
        typename vector<Coord>::iterator it=std::find(pos.begin(),pos.end(),p);
        if(it==pos.end()) {ret=pos.size(); indices.push_back(ret); pos.push_back(p); }
        else ret=it-pos.begin();
        return ret;
    }

    // add edge e in edg if not already there
    void addEdge(const Edge e, waEdges& edg)
    {
        if(e[0]==e[1]) return;
        typename vector<Edge>::iterator it=edg.begin();
        while(it!=edg.end() && ((*it)[0]!=e[0] || (*it)[1]!=e[1])) it++; // to replace std::find that does not compile here for some reasons..
        if(it==edg.end()) edg.push_back(e);
    }

    /**
    * @brief computes a uniform sample distribution (=point located at the center of their Voronoi cell) based on farthest point sampling + Lloyd (=kmeans) relaxation
    * @param nb : target number of samples
    * @param bias : bias distances using the input image ?
    * @param lloydIt : maximum number of Lloyd iterations.
    */

    void uniformSampling (const unsigned int nb=0,  const bool bias=false, const unsigned int lloydIt=100,const bool useDijkstra=false)
    {
        ImageSamplerSpecialization<ImageTypes::label>::uniformSampling( this, nb, bias, lloydIt, useDijkstra );
    }



    /**
    * same as above except that relaxation is done at each insertion of N samples
    * a graph is generated relating the new samples to its neighbors at the instant of insertion
    */

    void recursiveUniformSampling ( const unsigned int nb=0,  const bool bias=false, const unsigned int lloydIt=100,const bool useDijkstra=false,  const unsigned int N=1)
    {
        ImageSamplerSpecialization<ImageTypes::label>::recursiveUniformSampling( this, nb, bias, lloydIt, useDijkstra, N );
    }



};



} // namespace engine

} // namespace component

} // namespace sofa

#endif // SOFA_IMAGE_IMAGESAMPLER_H
