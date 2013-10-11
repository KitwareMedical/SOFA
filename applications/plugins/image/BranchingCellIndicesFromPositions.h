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
#ifndef SOFA_IMAGE_BranchingCellIndicesFromPositions_H
#define SOFA_IMAGE_BranchingCellIndicesFromPositions_H

#include "initImage.h"
#include "ImageTypes.h"
#include "BranchingImage.h"
#include <sofa/component/component.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/simulation/common/AnimateEndEvent.h>
#include <sofa/defaulttype/Vec.h>

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

/**
 * Returns global index of branching image voxels at sample locations, given a fine image of superimposed offsets
 */


template <class _ImageTypes,class _BranchingImageTypes>
class BranchingCellIndicesFromPositions : public core::DataEngine
{
public:
    typedef core::DataEngine Inherited;
    SOFA_CLASS(SOFA_TEMPLATE2(BranchingCellIndicesFromPositions,_ImageTypes,_BranchingImageTypes),Inherited);

    typedef SReal Real;

    typedef _ImageTypes ImageTypes;
    typedef typename ImageTypes::T T;
    typedef typename ImageTypes::imCoord imCoord;
    typedef helper::ReadAccessor<Data< ImageTypes > > raImage;
    Data< ImageTypes > image;

    typedef _BranchingImageTypes BranchingImageTypes;
    typedef typename BranchingImageTypes::T bT;
    typedef helper::ReadAccessor<Data< BranchingImageTypes > > raBranchingImage;
    Data< BranchingImageTypes > branchingImage;

    typedef defaulttype::ImageLPTransform<Real> TransformType;
    typedef typename TransformType::Coord Coord;
    typedef helper::ReadAccessor<Data< TransformType > > raTransform;
    Data< TransformType > transform;
    Data< TransformType > branchingImageTransform;

    typedef vector<Vec<3,Real> > SeqPositions;
    typedef helper::ReadAccessor<Data< SeqPositions > > raPositions;
    Data< SeqPositions > position;

    typedef vector<unsigned int> valuesType;
    typedef helper::WriteAccessor<Data< valuesType > > waValues;
    Data< valuesType > cell;  ///< output interpolated values

    virtual std::string getTemplateName() const    { return templateName(this);    }
    static std::string templateName(const BranchingCellIndicesFromPositions<ImageTypes,BranchingImageTypes>* = NULL) { return ImageTypes::Name()+std::string(",")+BranchingImageTypes::Name();    }

    BranchingCellIndicesFromPositions()    :   Inherited()
        , image(initData(&image,ImageTypes(),"image",""))
        , branchingImage(initData(&branchingImage,BranchingImageTypes(),"branchingImage",""))
        , transform(initData(&transform,TransformType(),"transform",""))
        , branchingImageTransform(initData(&branchingImageTransform,TransformType(),"branchingImageTransform",""))
        , position(initData(&position,SeqPositions(),"position","input positions"))
        , cell( initData ( &cell,"cell","cell indices." ) )
        , time((unsigned int)0)
    {
        this->addAlias(&branchingImageTransform, "branchingTransform");
        image.setReadOnly(true);
        branchingImage.setReadOnly(true);
        transform.setReadOnly(true);
        branchingImageTransform.setReadOnly(true);
        f_listening.setValue(true);
    }

    virtual void init()
    {
        addInput(&image);
        addInput(&branchingImage);
        addInput(&transform);
        addInput(&branchingImageTransform);
        addInput(&position);
        addOutput(&cell);
        setDirtyValue();
    }

    virtual void reinit() { update(); }

protected:

    unsigned int time;

    virtual void update()
    {
        cleanDirty();

        raImage in(this->image);
        raBranchingImage bin(this->branchingImage);
        raTransform inT(this->transform);
        raTransform binT(this->branchingImageTransform);
        raPositions pos(this->position);

        // get images at time t
        const typename ImageTypes::CImgT& img = in->getCImg(this->time);
        const typename BranchingImageTypes::BranchingImage3D& bimg = bin->imgList[this->time];

        waValues val(this->cell);
        unsigned int outval=0;
        val.resize(pos.size());

            for(unsigned int i=0; i<pos.size(); i++)
            {
                Coord Tp = inT->toImageInt(pos[i]);
                if(!in->isInside((int)Tp[0],(int)Tp[1],(int)Tp[2]))  val[i] = outval;
                else
                {
                    typename BranchingImageTypes::VoxelIndex vi;
                    Coord bTp = binT->toImageInt(pos[i]);
                    vi.index1d=bin->index3Dto1D(bTp[0],bTp[1],bTp[2]);
                    vi.offset=img.atXYZ(Tp[0],Tp[1],Tp[2]);
                    if(vi.offset<=0) val[i] = outval;
                    else if(vi.offset>bimg[vi.index1d].size()) val[i] = outval;
                    else
                    {
                        // count to retrieve global index.
                        val[i]=0;
                        for(unsigned int i1d=0;i1d<vi.index1d;i1d++) val[i]+=bimg[i1d].size();
                        val[i]+=vi.offset-1;
                    }
                }
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

};


} // namespace engine
} // namespace component
} // namespace sofa

#endif // SOFA_IMAGE_BranchingCellIndicesFromPositions_H
