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
#ifndef SOFA_ImageGaussPointSAMPLER_H
#define SOFA_ImageGaussPointSAMPLER_H

#include "../initFlexible.h"
#include "../quadrature/BaseGaussPointSampler.h"
#include "../deformationMapping/BaseDeformationMapping.h"

#include "../types/PolynomialBasis.h"

#include <image/ImageTypes.h>
#include <image/BranchingImage.h>
#include <image/ImageAlgorithms.h>

#include <sofa/helper/rmath.h>

#include <set>
#include <map>


namespace sofa
{
namespace component
{
namespace engine
{

using helper::vector;

/**
 * This class samples an object represented by an image with gauss points
 */

/// Default implementation does not compile
template <int imageTypeLabel>
struct ImageGaussPointSamplerSpecialization
{
};


/// Specialization for regular Image
template <>
struct ImageGaussPointSamplerSpecialization<defaulttype::IMAGELABEL_IMAGE>
{
    typedef unsigned int IndT;
    typedef defaulttype::Image<IndT> IndTypes;
    typedef defaulttype::Image<bool> MaskTypes;

    template<class ImageGaussPointSampler>
    static void init(ImageGaussPointSampler* This)
    {
        typedef typename ImageGaussPointSampler::IndTypes IndTypes;
        typedef typename ImageGaussPointSampler::waInd waInd;
        typedef typename ImageGaussPointSampler::DistTypes DistTypes;
        typedef typename ImageGaussPointSampler::raDist raDist;
        typedef typename ImageGaussPointSampler::waDist waDist;
        typedef typename ImageGaussPointSampler::waPositions waPositions;
        typedef typename ImageGaussPointSampler::waVolume waVolume;

        // retrieve data
        raDist rweights(This->f_w);             if(rweights->isEmpty())  { This->serr<<"Weights not found"<<This->sendl; return; }

        // init pos, vol, reg data; voronoi (=region data) and distances (=error image)
        typename ImageGaussPointSampler::imCoord dim = rweights->getDimensions();
        dim[DistTypes::DIMENSION_S]=dim[DistTypes::DIMENSION_T]=1;

        waPositions pos(This->f_position);          pos.clear();                // pos is cleared since it is always initialized with one point, so user placed points are not allowed for now..
        waVolume vol(This->f_volume);   vol.clear();

        waInd wreg(This->f_region);        wreg->setDimensions(dim);
        typename IndTypes::CImgT& regimg = wreg->getCImg();        regimg.fill(0);

        waDist werr(This->f_error);        werr->setDimensions(dim);
        typename DistTypes::CImgT& dist = werr->getCImg();        dist.fill(-1.0);

        This->Reg.clear();
    }


    /// midpoint integration : put samples uniformly and weight them by their volume

    template<class ImageGaussPointSampler>
    static void midpoint(ImageGaussPointSampler* This)
    {
        typedef typename ImageGaussPointSampler::Real Real;
        typedef typename ImageGaussPointSampler::IndTypes IndTypes;
        typedef typename ImageGaussPointSampler::raInd raInd;
        typedef typename ImageGaussPointSampler::waInd waInd;
        typedef typename ImageGaussPointSampler::DistTypes DistTypes;
        typedef typename ImageGaussPointSampler::DistT DistT;
        typedef typename ImageGaussPointSampler::waDist waDist;
        typedef typename ImageGaussPointSampler::SeqPositions SeqPositions;
        typedef typename ImageGaussPointSampler::Coord Coord;
        typedef typename ImageGaussPointSampler::waPositions waPositions;
        typedef typename ImageGaussPointSampler::waVolume waVolume;
        typedef typename ImageGaussPointSampler::indList indList;
        typedef typename ImageGaussPointSampler::raTransform raTransform;
        typedef typename ImageGaussPointSampler::factType factType;

        typedef Vec<3,int> iCoord;
        typedef std::pair<DistT,iCoord > DistanceToPoint;

        // retrieve data
        raInd rindices(This->f_index);          if(rindices->isEmpty())  { This->serr<<"Indices not found"<<This->sendl; return; }        const typename IndTypes::CImgT& indices = rindices->getCImg();
        raTransform transform(This->f_transform);
        const Coord voxelsize(transform->getScale());

        waPositions pos(This->f_position);
        waInd wreg(This->f_region);        typename IndTypes::CImgT& regimg = wreg->getCImg();
        waDist werr(This->f_error);        typename DistTypes::CImgT& dist = werr->getCImg();   // use error image as a container for distances

        // init soft regions (=more than one parent) where uniform sampling will be done
        // rigid regions (one parent) are kept in the list of region (and dist remains=-1 so they will not be sampled)
        for(unsigned int i=0; i<This->Reg.size();i++)
        {
            if(This->Reg[i].parentsToNodeIndex.size()>1)
            {
                cimg_forXYZ(regimg,x,y,z) if(regimg(x,y,z)==*(This->Reg[i].voronoiIndices.begin()) )
                {
                    dist(x,y,z)=cimg::type<DistT>::max();
                    regimg(x,y,z)=0;
                }
                This->Reg.erase (This->Reg.begin()+i); i--;  //  erase region (soft regions will be generated after uniform sampling)
            }
        }
        unsigned int nbrigid = This->Reg.size();

        // fixed points = points set by the user in soft regions.
        // Disabled for now since pos is cleared
        SeqPositions fpos_voxelIndex;
        vector<unsigned int> fpos_voronoiIndex;
        for(unsigned int i=0; i<pos.size(); i++)
        {
            Coord p = transform->toImageInt(pos[i]);
            if(indices.containsXYZC(p[0],p[1],p[2]))
            {
                indList l;
                cimg_forC(indices,v) if(indices(p[0],p[1],p[2],v)) l.insert(indices(p[0],p[1],p[2],v)-1);
                if(l.size()>1) { fpos_voxelIndex.push_back(p); fpos_voronoiIndex.push_back(i+1+nbrigid); }
            }
        }

        // target nb of points
        unsigned int nb = (fpos_voxelIndex.size()+nbrigid>This->targetNumber.getValue())?fpos_voxelIndex.size()+nbrigid:This->targetNumber.getValue();
        unsigned int nbsoft = nb-nbrigid;
        if(This->f_printLog.getValue()) std::cout<<This->getName()<<": Number of rigid/soft regions : "<<nbrigid<<"/"<<nbsoft<< std::endl;

        // init seeds for uniform sampling
        std::set<DistanceToPoint> trial;

        // farthest point sampling using geodesic distances
        SeqPositions newpos_voxelIndex;
        vector<unsigned int> newpos_voronoiIndex;

        for(unsigned int i=0; i<fpos_voxelIndex.size(); i++) AddSeedPoint<DistT>(trial,dist,regimg, fpos_voxelIndex[i],fpos_voronoiIndex[i]);
        while(newpos_voxelIndex.size()+fpos_voxelIndex.size()<nbsoft)
        {
            DistT dmax=0;  Coord pmax;
            cimg_forXYZ(dist,x,y,z) if(dist(x,y,z)>dmax) { dmax=dist(x,y,z); pmax =Coord(x,y,z); }
            if(dmax)
            {
                newpos_voxelIndex.push_back(pmax);
                newpos_voronoiIndex.push_back(fpos_voxelIndex.size()+nbrigid+newpos_voxelIndex.size());
                AddSeedPoint<DistT>(trial,dist,regimg, newpos_voxelIndex.back(),newpos_voronoiIndex.back());
                if(This->useDijkstra.getValue()) dijkstra<DistT,DistT>(trial,dist, regimg,voxelsize);
                else fastMarching<DistT,DistT>(trial,dist, regimg,voxelsize);
            }
            else break;
        }

        // Loyd
        unsigned int it=0;
        bool converged =(it>=This->iterations.getValue())?true:false;
        while(!converged)
        {
            converged=!(Lloyd<DistT>(newpos_voxelIndex,newpos_voronoiIndex,regimg));
            // recompute voronoi
            cimg_foroff(dist,off) if(dist[off]!=-1) dist[off]=cimg::type<DistT>::max();
            for(unsigned int i=0; i<fpos_voxelIndex.size(); i++) AddSeedPoint<DistT>(trial,dist,regimg, fpos_voxelIndex[i],fpos_voronoiIndex[i]);
            for(unsigned int i=0; i<newpos_voxelIndex.size(); i++) AddSeedPoint<DistT>(trial,dist,regimg, newpos_voxelIndex[i],newpos_voronoiIndex[i]);
            if(This->useDijkstra.getValue()) dijkstra<DistT,DistT>(trial,dist, regimg,voxelsize); else fastMarching<DistT,DistT>(trial,dist, regimg,voxelsize);
            it++; if(it>=This->iterations.getValue()) converged=true;
        }
        if(This->f_printLog.getValue()) std::cout<<This->getName()<<": Completed in "<< it <<" Lloyd iterations"<<std::endl;

        // create soft regions and update teir data
        for(unsigned int i=0; i<fpos_voxelIndex.size(); i++)           // Disabled for now since fpos is empty
        {
            indList l; cimg_forXYZ(regimg,x,y,z) if(regimg(x,y,z)==fpos_voronoiIndex[i]) { cimg_forC(indices,v) if(indices(x,y,z,v)) l.insert(indices(x,y,z,v)-1); }   // collect indices over the region
            if(l.size())
            {
                factType reg(l,fpos_voronoiIndex[i]); reg.center=transform->fromImage(fpos_voxelIndex[i]);
                This->Reg.push_back(reg);
            }
        }
        for(unsigned int i=0; i<newpos_voxelIndex.size(); i++)
        {
            indList l; cimg_forXYZ(regimg,x,y,z) if(regimg(x,y,z)==newpos_voronoiIndex[i]) { cimg_forC(indices,v) if(indices(x,y,z,v)) l.insert(indices(x,y,z,v)-1); }   // collect indices over the region
            if(l.size())
            {
                factType reg(l,newpos_voronoiIndex[i]); reg.center=transform->fromImage(newpos_voxelIndex[i]);
                This->Reg.push_back(reg);
            }
        }
        // update rigid regions (might contain soft material due to voronoi proximity)
        for(unsigned int i=0; i<nbrigid; i++)
        {
            indList l; cimg_forXYZ(regimg,x,y,z) if(regimg(x,y,z)==*(This->Reg[i].voronoiIndices.begin()) ) { cimg_forC(indices,v) if(indices(x,y,z,v)) l.insert(indices(x,y,z,v)-1); }   // collect indices over the region
            This->Reg[i].setParents(l);
        }

        // update nb voxels in each region (used later in weight fitting)
        for(unsigned int i=0; i<This->Reg.size(); i++)
        {
            This->Reg[i].nb=0; cimg_foroff(regimg,off)  if(regimg(off) == *(This->Reg[i].voronoiIndices.begin())) This->Reg[i].nb++;
        }
    }


    /// Identify regions sharing similar parents
    /// returns a list of region containing the parents, the number of voxels and center; and fill the voronoi image
    template<class ImageGaussPointSampler>
    static void Cluster_SimilarIndices(ImageGaussPointSampler* This)
    {
        typedef typename ImageGaussPointSampler::Real Real;
        typedef typename ImageGaussPointSampler::IndTypes IndTypes;
        typedef typename ImageGaussPointSampler::raInd raInd;
        typedef typename ImageGaussPointSampler::waInd waInd;
        typedef typename ImageGaussPointSampler::MaskTypes MaskTypes;
        typedef typename ImageGaussPointSampler::raMask raMask;
        typedef typename ImageGaussPointSampler::indList indList;
        typedef typename ImageGaussPointSampler::raTransform raTransform;
        typedef typename ImageGaussPointSampler::Coord Coord;
        typedef typename ImageGaussPointSampler::raPositions raPositions;
        typedef typename ImageGaussPointSampler::factType factType;

        // retrieve data
        raInd rindices(This->f_index);          if(rindices->isEmpty())  { This->serr<<"Indices not found"<<This->sendl; return; }        const typename IndTypes::CImgT& indices = rindices->getCImg();
        raMask rmask(This->f_mask);        const typename MaskTypes::CImgT* mask = rmask->isEmpty()?NULL:&rmask->getCImg();
        waInd wreg(This->f_region);        typename IndTypes::CImgT& regimg = wreg->getCImg();
        raTransform transform(This->f_transform);

        // map to find repartitions-> region index
        typedef std::map<indList, unsigned int> indMap;
        indMap List;

        // allows user to fix points. Currently disabled since pos is cleared
        raPositions pos(This->f_position);
        const unsigned int initialPosSize=pos.size();
        for(unsigned int i=0; i<initialPosSize; i++)
        {
            Coord p = transform->toImageInt(pos[i]);
            if(indices.containsXYZC(p[0],p[1],p[2]))
            {
                indList l;
                cimg_forC(indices,v) if(indices(p[0],p[1],p[2],v)) l.insert(indices(p[0],p[1],p[2],v)-1);
                List[l]=i;
                This->Reg.push_back(factType(l,i+1));
                regimg(p[0],p[1],p[2])=i+1;
            }
        }

        // traverse index image to identify regions with unique indices
        cimg_forXYZ(indices,x,y,z)
                if(indices(x,y,z))
                if(!mask || (*mask)(x,y,z))
        {
            indList l;
            cimg_forC(indices,v) if(indices(x,y,z,v)) l.insert(indices(x,y,z,v)-1);
            typename indMap::iterator it=List.find(l);
            unsigned int index;
            if(it==List.end()) { index=List.size(); List[l]=index;  This->Reg.push_back(factType(l,index+1)); This->Reg.back().nb=1; }
            else { index=it->second; This->Reg[index].nb++;}

            This->Reg[index].center+=transform->fromImage(Coord(x,y,z));
            regimg(x,y,z)=*(This->Reg[index].voronoiIndices.begin());
        }

        // average to get centroid (may not be inside the region if not convex)
        for(unsigned int i=0; i<This->Reg.size(); i++) This->Reg[i].center/=(Real)This->Reg[i].nb;
    }

    /// subdivide region[index] in two regions
    template<class ImageGaussPointSampler>
    static void subdivideRegion(ImageGaussPointSampler* This,const unsigned int index)
    {
        typedef typename ImageGaussPointSampler::Real Real;
        typedef typename ImageGaussPointSampler::IndTypes IndTypes;
        typedef typename ImageGaussPointSampler::waInd waInd;
        typedef typename ImageGaussPointSampler::DistTypes DistTypes;
        typedef typename ImageGaussPointSampler::DistT DistT;
        typedef typename ImageGaussPointSampler::waDist waDist;
        typedef typename ImageGaussPointSampler::raTransform raTransform;
        typedef typename ImageGaussPointSampler::Coord Coord;
        typedef typename ImageGaussPointSampler::factType factType;

        typedef Vec<3,int> iCoord;
        typedef std::pair<DistT,iCoord > DistanceToPoint;

        // retrieve data
        raTransform transform(This->f_transform);
        const Coord voxelsize(transform->getScale());

        waInd wreg(This->f_region);        typename IndTypes::CImgT& regimg = wreg->getCImg();
        waDist werr(This->f_error);        typename DistTypes::CImgT& dist = werr->getCImg();

        // compute
        vector<Coord> pos(2);
        vector<unsigned int> vorindex;
        vorindex.push_back(*(This->Reg[index].voronoiIndices.begin()));
        vorindex.push_back(This->Reg.size()+1);
        for(unsigned int i=0; i<This->Reg.size(); i++) if(vorindex[1]==*(This->Reg[i].voronoiIndices.begin())) vorindex[1]++; // check that the voronoi index is unique. not necessary in principle

        // get closest/farthest point from c and init distance image
        Real dmin=cimg::type<Real>::max(),dmax=0;
        cimg_forXYZ(regimg,x,y,z)
                if(regimg(x,y,z)==vorindex[0])
        {
            dist(x,y,z)=cimg::type<DistT>::max();
            Coord p = Coord(x,y,z);
            Real d = (transform->fromImage(p)-This->Reg[index].center).norm2();
            if(dmin>d) {dmin=d; pos[0]=p;}
            if(dmax<d) {dmax=d; pos[1]=p;}
        }
        else dist(x,y,z)=(DistT)(-1);

        // Loyd relaxation
        std::set<DistanceToPoint> trial;
        unsigned int it=0;
        bool converged =(it>=This->iterations.getValue())?true:false;

        for(unsigned int i=0; i<2; i++) AddSeedPoint<DistT>(trial,dist,regimg, pos[i],vorindex[i]);
        if(This->useDijkstra.getValue()) dijkstra<DistT,DistT>(trial,dist, regimg,voxelsize); else fastMarching<DistT,DistT>(trial,dist, regimg,voxelsize);
        //dist.display();
        //regimg.display();
        while(!converged)
        {
            converged=!(Lloyd<DistT>(pos,vorindex,regimg));
            // recompute voronoi
            cimg_foroff(dist,off) if(dist[off]!=-1) dist[off]=cimg::type<DistT>::max();
            for(unsigned int i=0; i<2; i++) AddSeedPoint<DistT>(trial,dist,regimg, pos[i],vorindex[i]);
            if(This->useDijkstra.getValue()) dijkstra<DistT,DistT>(trial,dist, regimg,voxelsize); else fastMarching<DistT,DistT>(trial,dist, regimg,voxelsize);
            it++; if(it>=This->iterations.getValue()) converged=true;
        }

        // add region
        factType reg;
        reg.parentsToNodeIndex=This->Reg[index].parentsToNodeIndex;
        reg.voronoiIndices.insert(vorindex[1]);
        reg.center=transform->fromImage(pos[1]);
        reg.nb=0; cimg_foroff(regimg,off)  if(regimg(off) == vorindex[1]) reg.nb++;
        This->Reg.push_back(reg);

        // update old region data
        This->Reg[index].center=transform->fromImage(pos[0]);
        This->Reg[index].nb=0; cimg_foroff(regimg,off)  if(regimg(off) == vorindex[0]) This->Reg[index].nb++;
    }



    /// update Polynomial Factors from the voxel map
    template<class ImageGaussPointSampler>
    static void fillPolynomialFactors(ImageGaussPointSampler* This,const unsigned int factIndex, const bool writeErrorImg=false)
    {
        typedef typename ImageGaussPointSampler::Real Real;
        typedef typename ImageGaussPointSampler::IndTypes IndTypes;
        typedef typename ImageGaussPointSampler::raInd raInd;
        typedef typename ImageGaussPointSampler::DistTypes DistTypes;
        typedef typename ImageGaussPointSampler::DistT DistT;
        typedef typename ImageGaussPointSampler::raDist raDist;
        typedef typename ImageGaussPointSampler::waDist waDist;
        typedef typename ImageGaussPointSampler::Coord Coord;
        typedef typename ImageGaussPointSampler::indList indList;
        typedef typename ImageGaussPointSampler::indListIt indListIt;
        typedef typename ImageGaussPointSampler::raTransform raTransform;
        typedef typename ImageGaussPointSampler::factType factType;

        // retrieve data
        raDist rweights(This->f_w);             if(rweights->isEmpty())  { This->serr<<"Weights not found"<<This->sendl; return; }  const typename DistTypes::CImgT& weights = rweights->getCImg();
        raInd rindices(This->f_index);          if(rindices->isEmpty())  { This->serr<<"Indices not found"<<This->sendl; return; }  const typename IndTypes::CImgT& indices = rindices->getCImg();
        raInd rreg(This->f_region);        const typename IndTypes::CImgT& regimg = rreg->getCImg();
        raTransform transform(This->f_transform);
        const Coord voxelsize(transform->getScale());

        // list of absolute coords
        factType &fact = This->Reg[factIndex];
        vector<Coord> pi(fact.nb);

        // weights (one line for each parent)
        typename ImageGaussPointSampler::Matrix wi(fact.parentsToNodeIndex.size(),fact.nb); wi.setZero();

        // get them from images
        unsigned int count=0;
        cimg_forXYZ(regimg,x,y,z)
        {
            indListIt it=fact.voronoiIndices.find(regimg(x,y,z));
            if(it!=fact.voronoiIndices.end())
            {
                cimg_forC(indices,v) if(indices(x,y,z,v))
                {
                    std::map<unsigned int,unsigned int>::iterator pit=fact.parentsToNodeIndex.find(indices(x,y,z,v)-1);
                    if(pit!=fact.parentsToNodeIndex.end())  wi(pit->second,count)= (Real)weights(x,y,z,v);
                }
                pi[count]= transform->fromImage(Coord(x,y,z));
                count++;
            }
        }

        fact.fill(wi,pi,This->fillOrder(),voxelsize,This->volOrder());

        //  std::cout<<"pt "<<*(fact.voronoiIndices.begin())-1<<" : "<<fact.center<<std::endl<<std::endl<<std::endl<<pi<<std::endl<<std::endl<<wi<<std::endl;
        //test: fact.directSolve(wi,pi); std::cout<<"Jacobi err="<<fact.getError()<<std::endl;

        // write error into output image
        if(writeErrorImg)
        {
            waDist werr(This->f_error); typename DistTypes::CImgT& outimg = werr->getCImg();
            count=0;
            cimg_forXYZ(regimg,x,y,z)
            {
                indListIt it=fact.voronoiIndices.find(regimg(x,y,z));
                if(it!=fact.voronoiIndices.end()) { outimg(x,y,z)=fact.getError(pi[count],wi.col(count)); count++; }
            }
        }
    }

};

/// Specialization for branching Image
template <>
struct ImageGaussPointSamplerSpecialization<defaulttype::IMAGELABEL_BRANCHINGIMAGE>
{
    typedef unsigned int IndT;
    typedef defaulttype::BranchingImage<IndT> IndTypes;
    typedef defaulttype::BranchingImage<bool> MaskTypes;

    template<class ImageGaussPointSampler>
    static void init(ImageGaussPointSampler* This)
    {
        typedef typename ImageGaussPointSampler::IndTypes IndTypes;
        typedef typename ImageGaussPointSampler::waInd waInd;
        typedef typename ImageGaussPointSampler::DistTypes DistTypes;
        typedef typename ImageGaussPointSampler::raDist raDist;
        typedef typename ImageGaussPointSampler::waDist waDist;
        typedef typename ImageGaussPointSampler::waPositions waPositions;
        typedef typename ImageGaussPointSampler::waVolume waVolume;

        // retrieve data
        raDist rweights(This->f_w);         const DistTypes& weights = rweights.ref();
        if(weights.isEmpty())  { This->serr<<"Weights not found"<<This->sendl; return; }

        // init pos, vol, reg data; voronoi (=region data) and distances (=error image)
        typename ImageGaussPointSampler::imCoord dim = weights.getDimensions();
        dim[DistTypes::DIMENSION_S]=dim[DistTypes::DIMENSION_T]=1;

        waPositions pos(This->f_position);          pos.clear();                // pos is cleared since it is always initialized with one point, so user placed points are not allowed for now..
        waVolume vol(This->f_volume);   vol.clear();

        waInd wreg(This->f_region);        IndTypes& regimg = wreg.wref();
        regimg.setDimensions(dim);
        regimg.cloneTopology(weights,0);

        waDist werr(This->f_error);       DistTypes& dist = werr.wref();
        dist.setDimensions(dim);
        dist.cloneTopology(weights,-1.0);

        This->Reg.clear();
    }


    /// midpoint integration : put samples uniformly and weight them by their volume

    template<class ImageGaussPointSampler>
    static void midpoint(ImageGaussPointSampler* This)
    {
        typedef typename ImageGaussPointSampler::Real Real;
        typedef typename ImageGaussPointSampler::IndTypes IndTypes;
        typedef typename ImageGaussPointSampler::raInd raInd;
        typedef typename ImageGaussPointSampler::waInd waInd;
        typedef typename ImageGaussPointSampler::DistTypes DistTypes;
        typedef typename ImageGaussPointSampler::DistT DistT;
        typedef typename ImageGaussPointSampler::waDist waDist;
        typedef typename ImageGaussPointSampler::Coord Coord;
        typedef typename ImageGaussPointSampler::waPositions waPositions;
        typedef typename ImageGaussPointSampler::waVolume waVolume;
        typedef typename ImageGaussPointSampler::indList indList;
        typedef typename ImageGaussPointSampler::raTransform raTransform;
        typedef typename ImageGaussPointSampler::factType factType;

        typedef typename DistTypes::VoxelIndex VoxelIndex;
        typedef std::pair<DistT,VoxelIndex > DistanceToPoint;
        typedef vector<VoxelIndex> SeqPositions;

        // retrieve data
        raInd rindices(This->f_index);          const IndTypes& indices = rindices.ref();
        if(indices.isEmpty())  { This->serr<<"Indices not found"<<This->sendl; return; }
        raTransform transform(This->f_transform);
        const Coord voxelsize(transform->getScale());

        waPositions pos(This->f_position);
        waInd wreg(This->f_region);        IndTypes& regimg = wreg.wref();
        waDist werr(This->f_error);        DistTypes& dist = werr.wref();   // use error image as a container for distances

        // init soft regions (=more than one parent) where uniform sampling will be done
        // rigid regions (one parent) are kept in the list of region (and dist remains=-1 so they will not be sampled)
        for(unsigned int i=0; i<This->Reg.size();i++)
        {
            if(This->Reg[i].parentsToNodeIndex.size()>1)
            {
                bimg_forCVoffT(regimg,c,v,off1D,t) if(regimg(off1D,v,c,t)==*(This->Reg[i].voronoiIndices.begin()) )
                {
                    dist(off1D,v,c,t)=cimg::type<DistT>::max();
                    regimg(off1D,v,c,t)=0;
                }
                This->Reg.erase (This->Reg.begin()+i); i--;  //  erase region (soft regions will be generated after uniform sampling)
            }
        }
        unsigned int nbrigid = This->Reg.size();

        // fixed points = points set by the user in soft regions.
        // Disabled for now since pos is cleared
        SeqPositions fpos_voxelIndex;
        vector<unsigned int> fpos_voronoiIndex;
        for(unsigned int i=0; i<pos.size(); i++)
        {
            Coord p = transform->toImageInt(pos[i]);
            if(indices.isInside(p[0],p[1],p[2]))
            {
                VoxelIndex vi(indices.index3Dto1D(p[0],p[1],p[2]),0);
                indList l;
                bimg_forC(indices,v) if(indices(vi,v)) l.insert(indices(vi,v)-1);
                if(l.size()>1) { fpos_voxelIndex.push_back(vi); fpos_voronoiIndex.push_back(i+1+nbrigid); }
            }
        }

        // target nb of points
        unsigned int nb = (fpos_voxelIndex.size()+nbrigid>This->targetNumber.getValue())?fpos_voxelIndex.size()+nbrigid:This->targetNumber.getValue();
        unsigned int nbsoft = nb-nbrigid;
        if(This->f_printLog.getValue()) std::cout<<This->getName()<<": Number of rigid/soft regions : "<<nbrigid<<"/"<<nbsoft<< std::endl;

        // init seeds for uniform sampling
        std::set<DistanceToPoint> trial;

        // farthest point sampling using geodesic distances
        SeqPositions newpos_voxelIndex;
        vector<unsigned int> newpos_voronoiIndex;

        for(unsigned int i=0; i<fpos_voxelIndex.size(); i++) AddSeedPoint<DistT>(trial,dist,regimg, fpos_voxelIndex[i],fpos_voronoiIndex[i]);
        while(newpos_voxelIndex.size()+fpos_voxelIndex.size()<nbsoft)
        {
            DistT dmax=0;  VoxelIndex indMax;
            bimg_forCVoffT(dist,c,v,off1D,t) if(dist(off1D,v,c,t)>dmax) { dmax=dist(off1D,v,c,t); indMax = VoxelIndex(off1D,v); }
            if(dmax)
            {
                newpos_voxelIndex.push_back(indMax);
                newpos_voronoiIndex.push_back(fpos_voxelIndex.size()+nbrigid+newpos_voxelIndex.size());
                AddSeedPoint<DistT>(trial,dist,regimg, newpos_voxelIndex.back(),newpos_voronoiIndex.back());
                if(This->useDijkstra.getValue()) dijkstra<DistT,DistT>(trial,dist, regimg,voxelsize);
                else fastMarching<DistT,DistT>(trial,dist, regimg,voxelsize);
            }
            else break;
        }

        // Loyd
        unsigned int it=0;
        bool converged =(it>=This->iterations.getValue())?true:false;
        while(!converged)
        {
            converged=!(Lloyd<DistT>(newpos_voxelIndex,newpos_voronoiIndex,regimg));
            // recompute voronoi
            bimg_forCVoffT(dist,c,v,off1D,t) if(dist(off1D,v,c,t)!=-1) dist(off1D,v,c,t)=cimg::type<Real>::max();
            for(unsigned int i=0; i<fpos_voxelIndex.size(); i++) AddSeedPoint<DistT>(trial,dist,regimg, fpos_voxelIndex[i],fpos_voronoiIndex[i]);
            for(unsigned int i=0; i<newpos_voxelIndex.size(); i++) AddSeedPoint<DistT>(trial,dist,regimg, newpos_voxelIndex[i],newpos_voronoiIndex[i]);
            if(This->useDijkstra.getValue()) dijkstra<DistT,DistT>(trial,dist, regimg,voxelsize); else fastMarching<DistT,DistT>(trial,dist, regimg,voxelsize);
            it++; if(it>=This->iterations.getValue()) converged=true;
        }
        if(This->f_printLog.getValue()) std::cout<<This->getName()<<": Completed in "<< it <<" Lloyd iterations"<<std::endl;

        // create soft regions and update teir data
        for(unsigned int i=0; i<fpos_voxelIndex.size(); i++)           // Disabled for now since fpos is empty
        {
            indList l;
            bimg_forVoffT(regimg,v,off1D,t0) if(regimg(off1D,v)==fpos_voronoiIndex[i])  {    bimg_forC(indices,c) if(indices(off1D,v,c,t0)) l.insert(indices(off1D,v,c,t0)-1); }   // collect indices over the region
            if(l.size())
            {
                factType reg(l,fpos_voronoiIndex[i]);
                unsigned x,y,z; regimg.index1Dto3D(fpos_voxelIndex[i].index1d, x,y,z);
                reg.center=transform->fromImage(Coord(x,y,z));
                // todo: memorize fpos_voxelIndex[i].offset in cell
                This->Reg.push_back(reg);
            }
        }
        for(unsigned int i=0; i<newpos_voxelIndex.size(); i++)
        {
            indList l;
            bimg_forVoffT(regimg,v,off1D,t0) if(regimg(off1D,v)==newpos_voronoiIndex[i])  {    bimg_forC(indices,c) if(indices(off1D,v,c,t0)) l.insert(indices(off1D,v,c,t0)-1); }   // collect indices over the region
            if(l.size())
            {
                factType reg(l,newpos_voronoiIndex[i]);
                unsigned x,y,z; regimg.index1Dto3D(newpos_voxelIndex[i].index1d, x,y,z);
                reg.center=transform->fromImage(Coord(x,y,z));
                // todo: memorize newpos_voxelIndex[i].offset in cell
                This->Reg.push_back(reg);
            }
        }
        // update rigid regions (might contain soft material due to voronoi proximity)
        for(unsigned int i=0; i<nbrigid; i++)
        {
            indList l;
            bimg_forVoffT(regimg,v,off1D,t0) if(regimg(off1D,v)==*(This->Reg[i].voronoiIndices.begin()))  {    bimg_forC(indices,c) if(indices(off1D,v,c,t0)) l.insert(indices(off1D,v,c,t0)-1); }   // collect indices over the region
            This->Reg[i].setParents(l);
        }

        // update nb voxels in each region (used later in weight fitting)
        for(unsigned int i=0; i<This->Reg.size(); i++)
        {
            This->Reg[i].nb=0;
            bimg_forVoffT(regimg,v,off1D,t0) if(regimg(off1D,v)==*(This->Reg[i].voronoiIndices.begin()))  This->Reg[i].nb++;
        }
    }



    /// Identify regions sharing similar parents
    /// returns a list of region containing the parents, the number of voxels and center; and fill the voronoi image
    template<class ImageGaussPointSampler>
    static void Cluster_SimilarIndices(ImageGaussPointSampler* This)
    {
        typedef typename ImageGaussPointSampler::Real Real;
        typedef typename ImageGaussPointSampler::IndTypes IndTypes;
        typedef typename ImageGaussPointSampler::raInd raInd;
        typedef typename ImageGaussPointSampler::waInd waInd;
        typedef typename ImageGaussPointSampler::MaskTypes MaskTypes;
        typedef typename ImageGaussPointSampler::raMask raMask;
        typedef typename ImageGaussPointSampler::indList indList;
        typedef typename ImageGaussPointSampler::raTransform raTransform;
        typedef typename ImageGaussPointSampler::Coord Coord;
        typedef typename ImageGaussPointSampler::raPositions raPositions;
        typedef typename ImageGaussPointSampler::factType factType;

        typedef typename IndTypes::VoxelIndex VoxelIndex;

        // retrieve data
        raInd rindices(This->f_index);     if(rindices->isEmpty())  { This->serr<<"Indices not found"<<This->sendl; return; }        const IndTypes& indices = rindices.ref();
        raMask rmask(This->f_mask);        const MaskTypes* mask = rmask->isEmpty()?NULL:&rmask.ref();
        waInd wreg(This->f_region);        IndTypes& regimg = wreg.wref();
        raTransform transform(This->f_transform);

        // map to find repartitions-> region index
        typedef std::map<indList, unsigned int> indMap;
        indMap List;

        // allows user to fix points. Currently disabled since pos is cleared
        raPositions pos(This->f_position);
        const unsigned int initialPosSize=pos.size();
        for(unsigned int i=0; i<initialPosSize; i++)
        {
            Coord p = transform->toImageInt(pos[i]);
            if(indices.isInside(p[0],p[1],p[2]))
            {
                VoxelIndex vi(indices.index3Dto1D(p[0],p[1],p[2]),0);
                indList l;
                bimg_forC(indices,v) if(indices(vi,v)) l.insert(indices(vi,v)-1);
                List[l]=i;
                This->Reg.push_back(factType(l,i+1));
                regimg(vi)=i+1;
            }
        }

        // traverse index image to identify regions with unique indices
        bimg_forVoffT(indices,v,off1D,t)
                if(indices(off1D,v))
                if(!mask || (*mask)(off1D,v) )
        {
            indList l;
            bimg_forC(indices,c) if(indices(off1D,v,c)) l.insert(indices(off1D,v,c)-1);
            typename indMap::iterator it=List.find(l);
            unsigned int index;
            if(it==List.end()) { index=List.size(); List[l]=index;  This->Reg.push_back(factType(l,index+1)); This->Reg.back().nb=1; }
            else { index=it->second; This->Reg[index].nb++;}
            unsigned x,y,z; regimg.index1Dto3D(off1D, x,y,z);
            This->Reg[index].center+=transform->fromImage(Coord(x,y,z));
            // todo: memorize v in cell
            regimg(off1D,v)=*(This->Reg[index].voronoiIndices.begin());
        }

        // average to get centroid (may not be inside the region if not convex)
        for(unsigned int i=0; i<This->Reg.size(); i++) This->Reg[i].center/=(Real)This->Reg[i].nb;
    }

    /// subdivide region[index] in two regions
    template<class ImageGaussPointSampler>
    static void subdivideRegion(ImageGaussPointSampler* This,const unsigned int index)
    {
        typedef typename ImageGaussPointSampler::Real Real;
        typedef typename ImageGaussPointSampler::IndTypes IndTypes;
        typedef typename ImageGaussPointSampler::waInd waInd;
        typedef typename ImageGaussPointSampler::DistTypes DistTypes;
        typedef typename ImageGaussPointSampler::DistT DistT;
        typedef typename ImageGaussPointSampler::waDist waDist;
        typedef typename ImageGaussPointSampler::raTransform raTransform;
        typedef typename ImageGaussPointSampler::Coord Coord;
        typedef typename ImageGaussPointSampler::factType factType;

        typedef typename DistTypes::VoxelIndex VoxelIndex;
        typedef std::pair<DistT,VoxelIndex > DistanceToPoint;

        // retrieve data
        raTransform transform(This->f_transform);
        const Coord voxelsize(transform->getScale());

        waInd wreg(This->f_region);        IndTypes& regimg = wreg.wref();
        waDist werr(This->f_error);        DistTypes& dist = werr.wref();

        // compute
        vector<VoxelIndex> pos(2);
        vector<unsigned int> vorindex;
        vorindex.push_back(*(This->Reg[index].voronoiIndices.begin()));
        vorindex.push_back(This->Reg.size()+1);
        for(unsigned int i=0; i<This->Reg.size(); i++) if(vorindex[1]==*(This->Reg[i].voronoiIndices.begin())) vorindex[1]++; // check that the voronoi index is unique. not necessary in principle

        // get closest/farthest point from c and init distance image
        Real dmin=cimg::type<Real>::max(),dmax=0;
        bimg_forVoffT(regimg,v,off1D,t0) if(regimg(off1D,v)==vorindex[0])
        {
            dist(off1D,v)=cimg::type<DistT>::max();
            unsigned x,y,z; regimg.index1Dto3D(off1D, x,y,z);
            Coord p = Coord(x,y,z);
            Real d = (transform->fromImage(p)-This->Reg[index].center).norm2();
            if(dmin>d) {dmin=d; pos[0]=VoxelIndex(off1D,v);}
            if(dmax<d) {dmax=d; pos[1]=VoxelIndex(off1D,v);}
        }
        else dist(off1D,v)=(DistT)(-1);

        // Loyd relaxation
        std::set<DistanceToPoint> trial;
        unsigned int it=0;
        bool converged =(it>=This->iterations.getValue())?true:false;

        for(unsigned int i=0; i<2; i++) AddSeedPoint<DistT>(trial,dist,regimg, pos[i],vorindex[i]);
        if(This->useDijkstra.getValue()) dijkstra<DistT,DistT>(trial,dist, regimg,voxelsize); else fastMarching<DistT,DistT>(trial,dist, regimg,voxelsize);
        //dist.display();
        //regimg.display();
        while(!converged)
        {
            converged=!(Lloyd<DistT>(pos,vorindex,regimg));
            // recompute voronoi
            bimg_forVoffT(dist,v,off1D,t) if(dist(off1D,v)!=-1) dist(off1D,v)=cimg::type<DistT>::max();
            for(unsigned int i=0; i<2; i++) AddSeedPoint<DistT>(trial,dist,regimg, pos[i],vorindex[i]);
            if(This->useDijkstra.getValue()) dijkstra<DistT,DistT>(trial,dist, regimg,voxelsize); else fastMarching<DistT,DistT>(trial,dist, regimg,voxelsize);
            it++; if(it>=This->iterations.getValue()) converged=true;
        }

        // add region
        factType reg;
        reg.parentsToNodeIndex=This->Reg[index].parentsToNodeIndex;
        reg.voronoiIndices.insert(vorindex[1]);
        unsigned x,y,z; regimg.index1Dto3D(pos[1].index1d, x,y,z);
        reg.center=transform->fromImage(Coord(x,y,z));
        // todo: memorize pos[1].offset in cell
        reg.nb=0;
        bimg_forVoffT(regimg,v,off1D,t0) if(regimg(off1D,v)==vorindex[1]) reg.nb++;
        This->Reg.push_back(reg);

        // update old region data
        regimg.index1Dto3D(pos[0].index1d, x,y,z);
        This->Reg[index].center=transform->fromImage(Coord(x,y,z));
        // todo: memorize pos[0].offset in cell
        This->Reg[index].nb=0;
        bimg_forVoffT(regimg,v,off1D,t0) if(regimg(off1D,v)==vorindex[0])  This->Reg[index].nb++;
    }



    /// update Polynomial Factors from the voxel map
    template<class ImageGaussPointSampler>
    static void fillPolynomialFactors(ImageGaussPointSampler* This,const unsigned int factIndex, const bool writeErrorImg=false)
    {
        typedef typename ImageGaussPointSampler::Real Real;
        typedef typename ImageGaussPointSampler::IndTypes IndTypes;
        typedef typename ImageGaussPointSampler::raInd raInd;
        typedef typename ImageGaussPointSampler::DistTypes DistTypes;
        typedef typename ImageGaussPointSampler::DistT DistT;
        typedef typename ImageGaussPointSampler::raDist raDist;
        typedef typename ImageGaussPointSampler::waDist waDist;
        typedef typename ImageGaussPointSampler::Coord Coord;
        typedef typename ImageGaussPointSampler::indList indList;
        typedef typename ImageGaussPointSampler::indListIt indListIt;
        typedef typename ImageGaussPointSampler::raTransform raTransform;
        typedef typename ImageGaussPointSampler::factType factType;

        // retrieve data
        raDist rweights(This->f_w);             if(rweights->isEmpty())  { This->serr<<"Weights not found"<<This->sendl; return; }  const DistTypes& weights = rweights.ref();
        raInd rindices(This->f_index);          if(rindices->isEmpty())  { This->serr<<"Indices not found"<<This->sendl; return; }  const IndTypes& indices = rindices.ref();
        raInd rreg(This->f_region);        const IndTypes& regimg = rreg.ref();
        raTransform transform(This->f_transform);
        const Coord voxelsize(transform->getScale());

        // list of absolute coords
        factType &fact = This->Reg[factIndex];
        vector<Coord> pi(fact.nb);

        // weights (one line for each parent)
        typename ImageGaussPointSampler::Matrix wi(fact.parentsToNodeIndex.size(),fact.nb); wi.setZero();

        // get them from images
        unsigned int count=0;
        bimg_forVoffT(regimg,v,off1D,t0)
        {
            indListIt it=fact.voronoiIndices.find(regimg(off1D,v));
            if(it!=fact.voronoiIndices.end())
            {
                bimg_forC(indices,c) if(indices(off1D,v,c))
                {
                    std::map<unsigned int,unsigned int>::iterator pit=fact.parentsToNodeIndex.find(indices(off1D,v,c)-1);
                    if(pit!=fact.parentsToNodeIndex.end())  wi(pit->second,count)= (Real)weights(off1D,v,c);
                }
                unsigned x,y,z; regimg.index1Dto3D(off1D, x,y,z);
                pi[count]= transform->fromImage(Coord(x,y,z));
                count++;
            }
        }

        fact.fill(wi,pi,This->fillOrder(),voxelsize,This->volOrder());

        //  std::cout<<"pt "<<*(fact.voronoiIndices.begin())-1<<" : "<<fact.center<<std::endl<<std::endl<<std::endl<<pi<<std::endl<<std::endl<<wi<<std::endl;
        //test: fact.directSolve(wi,pi); std::cout<<"Jacobi err="<<fact.getError()<<std::endl;

        // write error into output image
        if(writeErrorImg)
        {
            waDist werr(This->f_error); DistTypes& outimg = werr.wref();
            count=0;
            bimg_forVoffT(regimg,v,off1D,t0)
            {
                indListIt it=fact.voronoiIndices.find(regimg(off1D,v));
                if(it!=fact.voronoiIndices.end()) { outimg(off1D,v)=fact.getError(pi[count],wi.col(count)); count++; }
            }
        }
    }

};



template <class ImageTypes_>
class ImageGaussPointSampler : public BaseGaussPointSampler
{
    friend struct ImageGaussPointSamplerSpecialization<defaulttype::IMAGELABEL_IMAGE>;
    friend struct ImageGaussPointSamplerSpecialization<defaulttype::IMAGELABEL_BRANCHINGIMAGE>;
    typedef ImageGaussPointSamplerSpecialization<ImageTypes_::label> ImageGaussPointSamplerSpec;

public:
    typedef BaseGaussPointSampler Inherit;
    SOFA_CLASS(SOFA_TEMPLATE(ImageGaussPointSampler,ImageTypes_),Inherit);

    /** @name  GaussPointSampler types */
    //@{
    typedef Inherit::Real Real;
    typedef Inherit::Coord Coord;
    typedef Inherit::SeqPositions SeqPositions;
    typedef Inherit::raPositions raPositions;
    typedef Inherit::waPositions waPositions;
    //@}

    /** @name  Image data */
    //@{
    typedef typename ImageGaussPointSamplerSpec::IndT IndT;
    typedef typename ImageGaussPointSamplerSpec::IndTypes IndTypes;
    typedef helper::ReadAccessor<Data< IndTypes > > raInd;
    typedef helper::WriteAccessor<Data< IndTypes > > waInd;
    Data< IndTypes > f_index;

    typedef ImageTypes_ DistTypes;
    typedef typename DistTypes::T DistT;
    typedef typename DistTypes::imCoord imCoord;
    typedef helper::ReadAccessor<Data< DistTypes > > raDist;
    typedef helper::WriteAccessor<Data< DistTypes > > waDist;
    Data< DistTypes > f_w;

    typedef typename ImageGaussPointSamplerSpec::MaskTypes MaskTypes;
    typedef helper::ReadAccessor<Data< MaskTypes > > raMask;
    Data< MaskTypes > f_mask;

    typedef defaulttype::ImageLPTransform<Real> TransformType;
    typedef helper::ReadAccessor<Data< TransformType > > raTransform;
    Data< TransformType > f_transform;
    //@}

    /** @name  region data */
    //@{
    Data< IndTypes > f_region;
    Data< DistTypes > f_error;
    //@}

    /** @name  Options */
    //@{
    Data<bool> f_clearData;
    Data<unsigned int> targetNumber;
    Data<bool> useDijkstra;
    Data<unsigned int> iterations;
    //@}

    virtual std::string getTemplateName() const    { return templateName(this); }
    static std::string templateName(const ImageGaussPointSampler<ImageTypes_>* = NULL) { return ImageTypes_::Name(); }

    virtual void init()
    {
        Inherit::init();

        addInput(&f_index);
        addInput(&f_w);
        addInput(&f_transform);
        addInput(&f_mask);
        addOutput(&f_region);
        addOutput(&f_error);
        setDirtyValue();


        this->getContext()->get( deformationMapping, core::objectmodel::BaseContext::Local);
    }

    virtual void reinit() { update(); }

    virtual void bwdInit() {  updateMapping(); }

protected:
    ImageGaussPointSampler()    :   Inherit()
      , f_index(initData(&f_index,IndTypes(),"indices","image of dof indices"))
      , f_w(initData(&f_w,DistTypes(),"weights","weight image"))
      , f_mask(initData(&f_mask,MaskTypes(),"mask","optional mask to restrict the sampling region"))
      , f_transform(initData(&f_transform,TransformType(),"transform",""))
      , f_region(initData(&f_region,IndTypes(),"region","sample region : labeled image with sample indices"))
      , f_error(initData(&f_error,DistTypes(),"error","weigth fitting error"))
      , f_clearData(initData(&f_clearData,true,"clearData","clear region and error images after computation"))
      , targetNumber(initData(&targetNumber,(unsigned int)0,"targetNumber","target number of samples"))
      , useDijkstra(initData(&useDijkstra,true,"useDijkstra","Use Dijkstra for geodesic distance computation (use fastmarching otherwise)"))
      , iterations(initData(&iterations,(unsigned int)100,"iterations","maximum number of Lloyd iterations"))
      , deformationMapping(NULL)
    {
    }

    virtual ~ImageGaussPointSampler()
    {

    }

    /** @name  region types */
    //@{
    typedef Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>  Matrix;
    typedef Eigen::Matrix<Real,Eigen::Dynamic,1>  Vector;
    typedef typename std::set<unsigned int> indList;  ///< list of parent indices
    typedef typename indList::iterator indListIt;
    typedef typename defaulttype::PolynomialFitFactors<Real> factType;
    vector<factType> Reg;  ///< data related to each voronoi region
    //@}

    // polynomial orders
    inline unsigned int fillOrder() const {return 1;}     // For the mapping, we use first order fit (not 2nd order, to have translation invariance of elastons)
    inline unsigned int fitOrder() const {return (this->f_order.getValue()==1)?0:1;} // for elastons, we measure the quality of the integration using first order least squares fit
    inline unsigned int volOrder() const {return (this->f_order.getValue()==1)?0:4;} // for elastons, we generate volume moments up to order 4

    static const int spatial_dimensions=3;
    mapping::BasePointMapper<spatial_dimensions,Real>* deformationMapping; ///< link to local deformation mapping for weights update

    virtual void update()
    {
        cleanDirty();

        ImageGaussPointSamplerSpec::init(this);
        ImageGaussPointSamplerSpec::Cluster_SimilarIndices(this);

        if(this->f_order.getValue()==1)                                     ImageGaussPointSamplerSpec::midpoint(this);
        else if(this->f_method.getValue().getSelectedId() == GAUSSLEGENDRE) serr<<"GAUSSLEGENDRE quadrature not yet implemented"<<sendl;
        else if(this->f_method.getValue().getSelectedId() == NEWTONCOTES)   serr<<"NEWTONCOTES quadrature not yet implemented"<<sendl;
        else if(this->f_method.getValue().getSelectedId() == ELASTON)       this->elaston();

        this->fitWeights();

        if(this->f_clearData.getValue())
        {
            waDist err(this->f_error); err->clear();
            waInd reg(this->f_region); reg->clear();
        }

        this->updateMapping();

        if(this->f_printLog.getValue()) if(this->f_position.getValue().size())    std::cout<<this->getName()<<": "<< this->f_position.getValue().size() <<" generated samples"<<std::endl;
    }



    /// elaston integration : put samples so as to maximize weight linearity

    void elaston()
    {
        // retrieve data
        waPositions pos(this->f_position);

        // fit weights
        for(unsigned int i=pos.size(); i<this->Reg.size(); i++)  {ImageGaussPointSamplerSpec::fillPolynomialFactors(this,i); this->Reg[i].solve(this->fitOrder());}

        // subdivide region with largest error until target number is reached
        while(this->Reg.size()<this->targetNumber.getValue())
        {
            Real maxerr=-1;
            unsigned int maxindex=0;
            for(unsigned int i=0; i<this->Reg.size(); i++) { Real e=this->Reg[i].getError(); if(maxerr<e) {maxerr=e; maxindex=i;} }
            if(maxerr==0) break;
            ImageGaussPointSamplerSpec::subdivideRegion(this,maxindex);

            ImageGaussPointSamplerSpec::fillPolynomialFactors(this,maxindex); this->Reg[maxindex].solve(this->fitOrder());
            ImageGaussPointSamplerSpec::fillPolynomialFactors(this,this->Reg.size()-1);    this->Reg.back().solve(this->fitOrder());
        }
    }


    /// fit weights to obtain final weights and derivatives
    /// optionaly write error image
    void fitWeights()
    {
        Real err=0;
        for(unsigned int i=0; i<this->Reg.size(); i++)
        {
            ImageGaussPointSamplerSpec::fillPolynomialFactors(this,i,!this->f_clearData.getValue());
            this->Reg[i].solve(this->fitOrder());
            err+=this->Reg[i].getError();
            //if(this->f_printLog.getValue()) std::cout<<this->getName()<<"GaussPointSampler: weight fitting error on sample "<<i<<" = "<<this->Reg[i].getError()<< std::endl;
        }
//        waDist werr(this->f_error);        typename DistTypes& errimg = werr.wref();
//        cimg_forXYZ(errimg,x,y,z) if(errimg(x,y,z)==-1) errimg(x,y,z)=0; // clean error output image (used as a container for distances)
        if(this->f_printLog.getValue()) std::cout<<this->getName()<<": total error = "<<err<<std::endl;
    }

    /// update mapping with weights fitted over a region (order 2)
    /// typically done in bkwinit (to overwrite weights computed in the mapping using shape function interpolation)
    virtual void updateMapping()
    {
        if(!deformationMapping) {serr<<"deformationMapping not found -> cannot map Gauss points"<<sendl; return;}

        unsigned int nb = Reg.size();

        waPositions pos(this->f_position);
        waVolume vol(this->f_volume);

        pos.resize ( nb );
        vol.resize ( nb );

        vector<vector<unsigned int> > index(nb);
        vector<vector<Real> > w(nb);
        vector<vector<Vec<spatial_dimensions,Real> > > dw(nb);
        vector<vector<Mat<spatial_dimensions,spatial_dimensions,Real> > > ddw(nb);

        Mat<spatial_dimensions,spatial_dimensions,Real> I=Mat<spatial_dimensions,spatial_dimensions,Real>::Identity(); // could be image orientation
        vector<Mat<spatial_dimensions,spatial_dimensions,Real> > F0((int)nb,I);

        for(unsigned int i=0; i<nb; i++)
        {
            factType* reg=&Reg[i];

            reg->solve(fillOrder());
            pos[i]=reg->center;
            vol[i].resize(reg->vol.rows());  for(unsigned int j=0; j<vol[i].size(); j++) vol[i][j]=reg->vol(j);
            reg->getMapping(index[i],w[i],dw[i],ddw[i]);
        }

        // test
        /*for(unsigned int i=0; i<nb; i++)
        {
            Real sumw=0; for(unsigned int j=0; j<w[i].size(); j++) { sumw+=w[i][j]; }
            Vec<spatial_dimensions,Real>  sumdw; for(unsigned int j=0; j<dw[i].size(); j++) sumdw+=dw[i][j];
            if(sumdw.norm()>1E-2 || fabs(sumw-1)>1E-2) std::cout<<"error on "<<i<<" : "<<sumw<<","<<sumdw<<std::endl;
        }*/

        if(this->f_printLog.getValue())  std::cout<<this->getName()<<" : "<< nb <<" gauss points exported"<<std::endl;
        deformationMapping->resizeOut(pos.ref(),index,w,dw,ddw,F0);
    }



};

}
}
}

#endif
