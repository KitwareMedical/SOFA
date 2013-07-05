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
#ifndef FLEXIBLE_BaseImageShapeFunction_H
#define FLEXIBLE_BaseImageShapeFunction_H

#include "../initFlexible.h"
#include "../shapeFunction/BaseShapeFunction.h"
#include "../types/PolynomialBasis.h"

#include <image/ImageTypes.h>
#include <image/BranchingImage.h>
#include <image/ImageAlgorithms.h>

#include <sofa/helper/rmath.h>
#include <sofa/helper/OptionsGroup.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <string>

namespace sofa
{
namespace component
{
namespace shapefunction
{

using core::behavior::BaseShapeFunction;
using defaulttype::Mat;
using defaulttype::Vec;





/// Default implementation does not compile
template <int imageTypeLabel>
struct BaseImageShapeFunctionSpecialization
{
};


/// Specialization for regular Image
template <>
struct BaseImageShapeFunctionSpecialization<defaulttype::IMAGELABEL_IMAGE>
{
    typedef SReal DistT;
    typedef defaulttype::Image<DistT> DistTypes;
    typedef unsigned int IndT;
    typedef defaulttype::Image<IndT> IndTypes;

    template<class BaseImageShapeFunction>
    static void constructor( BaseImageShapeFunction* This )
    {
        This->f_cell.setDisplayed( false );
    }

    /// interpolate weights and their derivatives at a spatial position
    template<class BaseImageShapeFunction>
    static void computeShapeFunction( BaseImageShapeFunction* This, const typename BaseImageShapeFunction::Coord& childPosition, typename BaseImageShapeFunction::MaterialToSpatial& /*M*/, typename BaseImageShapeFunction::VRef& ref, typename BaseImageShapeFunction::VReal& w, typename BaseImageShapeFunction::VGradient* dw=NULL, typename BaseImageShapeFunction::VHessian* ddw=NULL, const int /*cell*/=-1 )
    {
        typedef typename BaseImageShapeFunction::Real Real;
        typedef typename BaseImageShapeFunction::IndT IndT;
        typedef typename BaseImageShapeFunction::DistT DistT;
        typedef typename BaseImageShapeFunction::Coord Coord;

        // get transform
        typename BaseImageShapeFunction::raTransform inT(This->transform);

        // get precomputed indices and weights
        typename BaseImageShapeFunction::raInd indData(This->f_index);
        typename BaseImageShapeFunction::raDist weightData(This->f_w);
        if(indData->isEmpty() || weightData->isEmpty()) { This->serr<<"Weights not available"<<This->sendl; return; }

        const typename BaseImageShapeFunction::IndTypes::CImgT& indices = indData->getCImg();
        const typename BaseImageShapeFunction::DistTypes::CImgT& weights = weightData->getCImg();

        // interpolate weights in neighborhood
        Coord p = inT->toImage(childPosition);
        Coord P;  for (unsigned int j=0; j<3; j++)  P[j]=sofa::helper::round(p[j]);
        unsigned int order=0;
        /*if(ddw) order=2; else */  // do not use order 2 for local weight interpolation. Order two is used only in weight fitting over regions
        if(dw) order=1;

        //get closest voxel with non zero weights
        if(P[0]<=0 || P[1]<=0 || P[2]<=0 || P[0]>=indices.width()-1 || P[1]>=indices.height()-1 || P[2]>=indices.depth()-1 ||
                indices(P[0],P[1],P[2],0)==0)
        {
            Real dmin=cimg::type<Real>::max();
            cimg_for_insideXYZ(indices,x,y,z,1) if(indices(x,y,z,0)) {Real d=(Coord(x,y,z)-p).norm2(); if(d<dmin) { P.set(x,y,z); dmin=d; } }
            if(dmin==cimg::type<Real>::max()) return;
        }

        // prepare neighborood
        sofa::defaulttype::Vec<27,  Coord > lpos;      // precomputed local positions
        int count=0;
        for (int k=-1; k<=1; k++) for (int j=-1; j<=1; j++) for (int i=-1; i<=1; i++) lpos[count++]= inT->fromImage(P+Coord(i,j,k)) - childPosition;

        // get indices at P
        int index=0;
        unsigned int nbRef=This->f_nbRef.getValue();
        for (unsigned int r=0; r<nbRef; r++)
        {
            IndT ind=indices(P[0],P[1],P[2],r);
            if(ind>0)
            {
                vector<DistT> val; val.reserve(27);
                vector<Coord> pos; pos.reserve(27);
                // add neighbors with same index
                count=0;
                for (int k=-1; k<=1; k++) for (int j=-1; j<=1; j++) for (int i=-1; i<=1; i++)
                {
                    for (unsigned int r2=0;r2<nbRef;r2++)
                        if(indices(P[0]+i,P[1]+j,P[2]+k,r2)==ind)
                        {
                            val.push_back(weights(P[0]+i,P[1]+j,P[2]+k,r2));
                            pos.push_back(lpos[count]);
                        }
                    count++;
                }
                // fit weights
                vector<Real> coeff;
                defaulttype::PolynomialFit(coeff,val,pos, order);
                //std::cout<<ind<<":"<<coeff[0]<<", err= "<<getPolynomialFit_Error(coeff,val,pos)<< std::endl;
                if(!dw) defaulttype::getPolynomialFit_differential(coeff,w[index]);
                else if(!ddw) defaulttype::getPolynomialFit_differential(coeff,w[index],&(*dw)[index]);
                else defaulttype::getPolynomialFit_differential(coeff,w[index],&(*dw)[index],&(*ddw)[index]);
                ref[index]=ind-1; // remove offset from indices image
                if(w[index]<=0) // clamp negative weights
                {
                    w[index]=0;
                    if(dw) (*dw)[index].fill(0);
                    if(ddw) (*ddw)[index].fill(0);
                    index--;
                }
                index++;

            }
        }
        // remove unecessary weights
        ref.resize(index);
        w.resize(index);
        if(dw)  dw->resize(index);
        if(ddw) ddw->resize(index);
    }
};


/// Specialization for Branching Image
template <>
struct BaseImageShapeFunctionSpecialization<defaulttype::IMAGELABEL_BRANCHINGIMAGE>
{
    typedef SReal DistT;
    typedef defaulttype::BranchingImage<DistT> DistTypes;
    typedef unsigned int IndT;
    typedef defaulttype::BranchingImage<IndT> IndTypes;

    typedef IndTypes::ConnectionVoxel ConnectionVoxel;
    typedef IndTypes::VoxelIndex VoxelIndex;
    typedef IndTypes::Neighbours Neighbours;

    template<class BaseImageShapeFunction>
    static void constructor( BaseImageShapeFunction* This )
    {
        This->addAlias( &This->image, "branchingImage" );
    }


    /// interpolate weights and their derivatives at a spatial position
    template<class BaseImageShapeFunction>
    static void computeShapeFunction( BaseImageShapeFunction* This, const typename BaseImageShapeFunction::Coord& childPosition, typename BaseImageShapeFunction::MaterialToSpatial& /*M*/, typename BaseImageShapeFunction::VRef& ref, typename BaseImageShapeFunction::VReal& w, typename BaseImageShapeFunction::VGradient* dw=NULL, typename BaseImageShapeFunction::VHessian* ddw=NULL, const int cell=-1)
    {
        typedef typename BaseImageShapeFunction::Real Real;
        typedef typename BaseImageShapeFunction::IndT IndT;
        typedef typename BaseImageShapeFunction::DistT DistT;
        typedef typename BaseImageShapeFunction::Coord Coord;

        // get transform
        typename BaseImageShapeFunction::raTransform inT(This->transform);

        // get precomputed indices and weights
        typename BaseImageShapeFunction::raInd indData(This->f_index);
        typename BaseImageShapeFunction::raDist weightData(This->f_w);
        if(indData->isEmpty() || weightData->isEmpty()) { This->serr<<"Weights not available"<<This->sendl; return; }

        const typename BaseImageShapeFunction::IndTypes::BranchingImage3D& indices = indData->imgList[0];
        const typename BaseImageShapeFunction::DistTypes::BranchingImage3D& weights = weightData->imgList[0];

        // interpolate weights in neighborhood
        Coord p = inT->toImage( childPosition );
        Coord P;  for (unsigned int j=0; j<3; j++)  P[j]=sofa::helper::round(p[j]);
        VoxelIndex voxelIndex (indData->index3Dto1D( P[0], P[1], P[2] ) , 0);
        if(cell>0) voxelIndex.offset=cell-1;
        unsigned int order=0;
        /*if(ddw) order=2; else */  // do not use order 2 for local weight interpolation. Order two is used only in weight fitting over regions
        if(dw) order=1;

        // get closest voxel with non zero weights
        if( !indData.ref().isInside((int)P[0],(int)p[1],(int)P[2]) ||
                indices[voxelIndex.index1d].size()>voxelIndex.offset )
        {
            Real dmin=cimg::type<Real>::max();
            bimg_forXYZ(indData.ref(),x,y,z) if(indices[indData->index3Dto1D(x,y,z)].size()) {Real d=(Coord(x,y,z)-p).norm2(); if(d<dmin) { P.set(x,y,z); dmin=d;  voxelIndex.index1d=indData->index3Dto1D( P[0], P[1], P[2] );   } }
            if(dmin==cimg::type<Real>::max()) return;
        }
        if(voxelIndex.offset>=indices[voxelIndex.index1d].size()) voxelIndex.offset=0;

        // prepare neighborood
        const ConnectionVoxel& indV = indices[voxelIndex.index1d][voxelIndex.offset];
        const Neighbours& indN = indV.neighbours;
        vector< Coord > lpos;
        lpos.push_back ( inT->fromImage(P) - childPosition ); // add central voxel
        for (unsigned int n=0; n<indN.size(); n++) // add neighbors
        {
            Vec<3,unsigned int> P2;  indData->index1Dto3D( indN[n].index1d, P2[0], P2[1], P2[2] );
            lpos.push_back ( inT->fromImage(P2) - childPosition );
        }

        // get indices at P
        int index=0;
        const unsigned int nbRef = This->f_nbRef.getValue();
        for (unsigned int r=0; r<nbRef; r++)
        {
            IndT ind=indV[r];
            if(ind>0)
            {
                vector<DistT> val; val.reserve(lpos.size());
                vector<Coord> pos; pos.reserve(lpos.size());

                val.push_back(weights[voxelIndex.index1d][voxelIndex.offset][r]);
                pos.push_back(lpos[0]);

                // add neighbors with same index
                for (unsigned int n=0; n<indN.size(); n++)
                    for (unsigned int r2=0;r2<nbRef;r2++)
                        if(indices[indN[n].index1d][indN[n].offset][r2]==ind)
                        {
                            val.push_back(weights[indN[n].index1d][indN[n].offset][r2]);
                            pos.push_back(lpos[n+1]);
                        }
                // fit weights
                vector<Real> coeff;
                defaulttype::PolynomialFit(coeff,val,pos, order);
                //std::cout<<ind<<":"<<coeff[0]<<", err= "<<getPolynomialFit_Error(coeff,val,pos)<< std::endl;
                if(!dw) defaulttype::getPolynomialFit_differential(coeff,w[index]);
                else if(!ddw) defaulttype::getPolynomialFit_differential(coeff,w[index],&(*dw)[index]);
                else defaulttype::getPolynomialFit_differential(coeff,w[index],&(*dw)[index],&(*ddw)[index]);
                ref[index]=ind-1;  // remove offset from indices image
                if(w[index]<=0) // clamp negative weights
                {
                    w[index]=0;
                    if(dw) (*dw)[index].fill(0);
                    if(ddw) (*ddw)[index].fill(0);
                    index--;
                }
                index++;
            }
        }
        // remove unecessary weights
        ref.resize(index);
        w.resize(index);
        if(dw)  dw->resize(index);
        if(ddw) ddw->resize(index);
    }
};





/**
abstract class for shape functions computed from a set of images (typically rasterized objects)
  */
template <class ShapeFunctionTypes_,class ImageTypes_>
class BaseImageShapeFunction : public BaseShapeFunction<ShapeFunctionTypes_>
{
    friend struct BaseImageShapeFunctionSpecialization<defaulttype::IMAGELABEL_IMAGE>;
    friend struct BaseImageShapeFunctionSpecialization<defaulttype::IMAGELABEL_BRANCHINGIMAGE>;

public:
    SOFA_ABSTRACT_CLASS(SOFA_TEMPLATE2(BaseImageShapeFunction, ShapeFunctionTypes_,ImageTypes_) , SOFA_TEMPLATE(BaseShapeFunction, ShapeFunctionTypes_));
    typedef BaseShapeFunction<ShapeFunctionTypes_> Inherit;

    /** @name  Shape function types */
    //@{
    typedef typename Inherit::Real Real;
    typedef typename Inherit::Coord Coord;
    typedef typename Inherit::VCoord VCoord;
    enum {material_dimensions=Inherit::material_dimensions};
    typedef typename Inherit::VReal VReal;
    typedef typename Inherit::VGradient VGradient;
    typedef typename Inherit::VHessian VHessian;
    typedef typename Inherit::VRef VRef;

    typedef typename Inherit::Gradient Gradient;
    typedef typename Inherit::Hessian Hessian;
    enum {spatial_dimensions=Inherit::spatial_dimensions};
    typedef typename Inherit::MaterialToSpatial MaterialToSpatial;
    typedef typename Inherit::VMaterialToSpatial VMaterialToSpatial;
    //@}

    /** @name  Image data */
    //@{
    typedef ImageTypes_ ImageTypes;
    typedef typename ImageTypes::T T;
    typedef typename ImageTypes::imCoord imCoord;
    typedef helper::ReadAccessor<Data< ImageTypes > > raImage;
    Data< ImageTypes > image;

    typedef defaulttype::ImageLPTransform<Real> TransformType;
    typedef helper::ReadAccessor<Data< TransformType > > raTransform;
    Data< TransformType > transform;

    typedef typename BaseImageShapeFunctionSpecialization<ImageTypes::label>::DistT DistT;
    typedef typename BaseImageShapeFunctionSpecialization<ImageTypes::label>::DistTypes DistTypes;
    typedef helper::ReadAccessor<Data< DistTypes > > raDist;
    typedef helper::WriteAccessor<Data< DistTypes > > waDist;
    Data< DistTypes > f_w;

    typedef typename BaseImageShapeFunctionSpecialization<ImageTypes::label>::IndT IndT;
    typedef typename BaseImageShapeFunctionSpecialization<ImageTypes::label>::IndTypes IndTypes;
    typedef helper::ReadAccessor<Data< IndTypes > > raInd;
    typedef helper::WriteAccessor<Data< IndTypes > > waInd;
    Data< IndTypes > f_index;

    // only used for branching image
     Data< vector<int> > f_cell;    ///< indices required by shape function in case of overlapping elements
    //@}

    virtual std::string getTemplateName() const    { return templateName(this); }
    static std::string templateName(const BaseImageShapeFunction<ShapeFunctionTypes_,ImageTypes_>* = NULL) { return ShapeFunctionTypes_::Name()+std::string(",")+ImageTypes_::Name(); }


    /// interpolate weights and their derivatives at a spatial position
    void computeShapeFunction(const Coord& childPosition, MaterialToSpatial& M, VRef& ref, VReal& w, VGradient* dw=NULL,VHessian* ddw=NULL, const int cell=-1)
    {
        // resize input
        unsigned int nbRef=this->f_nbRef.getValue();
        ref.resize(nbRef); ref.fill(0);
        w.resize(nbRef); w.fill(0);
        if(dw) { dw->resize(nbRef); for (unsigned int j=0; j<nbRef; j++ ) (*dw)[j].fill(0); }
        if(ddw) { ddw->resize(nbRef); for (unsigned int j=0; j<nbRef; j++ ) (*ddw)[j].fill(0); }

        // material to world transformation = image orientation
        helper::Quater<Real> q = helper::Quater< Real >::createQuaterFromEuler(this->transform.getValue().getRotation() * (Real)M_PI / (Real)180.0);
        Mat<3,3,Real> R; q.toMatrix(R);
        for ( unsigned int i = 0; i < BaseImageShapeFunction::spatial_dimensions; i++ )  for ( unsigned int j = 0; j < BaseImageShapeFunction::material_dimensions; j++ ) M[i][j]=R[i][j];

        BaseImageShapeFunctionSpecialization<ImageTypes::label>::computeShapeFunction( this, childPosition, M, ref, w, dw, ddw, cell );

        // normalize
        this->normalize(w,dw,ddw);
    }

    virtual void init()
    {
        Inherit::init();
    }

protected:
    BaseImageShapeFunction()
        :Inherit()
        , image(initData(&image,ImageTypes(),"image",""))
        , transform(initData(&transform,TransformType(),"transform",""))
        , f_w(initData(&f_w,DistTypes(),"weights",""))
        , f_index(initData(&f_index,IndTypes(),"indices",""))
        , f_cell ( initData ( &f_cell,"cell","indices of surimposed voxels required in case of overlapping elements" ) )
    {
        image.setReadOnly(true);
        transform.setReadOnly(true);

        image.setGroup("input");
        transform.setGroup("input");
        f_w.setGroup("output");
        f_index.setGroup("output");

        BaseImageShapeFunctionSpecialization<ImageTypes::label>::constructor( this );
    }

    virtual ~BaseImageShapeFunction()
    {

    }


};


}
}
}


#endif
