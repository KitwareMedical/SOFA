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
#ifndef SOFA_IMAGE_IMAGEOPERATION_H
#define SOFA_IMAGE_IMAGEOPERATION_H

#include "initImage.h"
#include "ImageTypes.h"
#include <sofa/core/DataEngine.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/helper/rmath.h>
#include <sofa/helper/OptionsGroup.h>

#include <sofa/component/component.h>

#define ADDITION 0
#define SUBTRACTION 1
#define MULTIPLICATION 2
#define DIVISION 3


namespace sofa
{

namespace component
{

namespace engine
{

using helper::vector;
using cimg_library::CImg;
using cimg_library::CImgList;

/**
 * This class computes an image as an operation between two images
 */


template <class _ImageTypes>
class ImageOperation : public core::DataEngine
{
public:
    typedef core::DataEngine Inherited;
    SOFA_CLASS(SOFA_TEMPLATE(ImageOperation,_ImageTypes),Inherited);

    typedef _ImageTypes ImageTypes;
    typedef typename ImageTypes::T T;
    typedef typename ImageTypes::imCoord imCoord;
    typedef helper::WriteAccessor<Data< ImageTypes > > waImage;
    typedef helper::ReadAccessor<Data< ImageTypes > > raImage;

    Data<helper::OptionsGroup> operation;

    Data< ImageTypes > inputImage1;
    Data< ImageTypes > inputImage2;

    Data< ImageTypes > outputImage;

    virtual std::string getTemplateName() const    { return templateName(this);    }
    static std::string templateName(const ImageOperation<ImageTypes>* = NULL) { return ImageTypes::Name(); }

    ImageOperation()    :   Inherited()
      , operation ( initData ( &operation,"operation","operation" ) )
      , inputImage1(initData(&inputImage1,ImageTypes(),"inputImage1",""))
      , inputImage2(initData(&inputImage2,ImageTypes(),"inputImage2",""))
      , outputImage(initData(&outputImage,ImageTypes(),"outputImage",""))
    {
        inputImage1.setReadOnly(true);  this->addAlias(&inputImage1, "image1");
        inputImage2.setReadOnly(true);  this->addAlias(&inputImage2, "image2");
        outputImage.setReadOnly(true);  this->addAlias(&outputImage, "image");
        helper::OptionsGroup operationOptions(4	,"0 - Addition"
                                           ,"1 - Subtraction"
                                           ,"2 - Multiplication"
                                           ,"3 - Division"
                                           );
        operationOptions.setSelectedItem(SUBTRACTION);
        operation.setValue(operationOptions);
    }

    virtual ~ImageOperation() {}

    virtual void init()
    {
        addInput(&inputImage1);
        addInput(&inputImage2);
        addOutput(&outputImage);
        setDirtyValue();
    }

    virtual void reinit() { update(); }

protected:

    virtual void update()
    {
        cleanDirty();

        raImage in1(this->inputImage1);
        raImage in2(this->inputImage2);
        waImage out(this->outputImage);

        if(in1->isEmpty() || in2->isEmpty()) return;

        const CImgList<T>& inimg1 = in1->getCImgList() , inimg2 = in2->getCImgList();
        CImgList<T>& img = out->getCImgList();
        img.assign(inimg1);	// copy

        switch(this->operation.getValue().getSelectedId())
        {
        case ADDITION:            cimglist_for(img,l) cimg_forXYZC(img(l),x,y,z,c) img(l)(x,y,z,c)+=inimg2(l)(x,y,z,c);            break;
        case SUBTRACTION:         cimglist_for(img,l) cimg_forXYZC(img(l),x,y,z,c) img(l)(x,y,z,c)-=inimg2(l)(x,y,z,c);            break;
        case MULTIPLICATION:      cimglist_for(img,l) cimg_forXYZC(img(l),x,y,z,c) img(l)(x,y,z,c)*=inimg2(l)(x,y,z,c);            break;
        case DIVISION:            cimglist_for(img,l) cimg_forXYZC(img(l),x,y,z,c) img(l)(x,y,z,c)/=inimg2(l)(x,y,z,c);            break;
        default:            break;
        }
    }

};


} // namespace engine

} // namespace component

} // namespace sofa

#endif // SOFA_IMAGE_IMAGEOPERATION_H
