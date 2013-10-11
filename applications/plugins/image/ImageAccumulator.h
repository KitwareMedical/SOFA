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
#ifndef SOFA_IMAGE_IMAGEACCUMULATOR_H
#define SOFA_IMAGE_IMAGEACCUMULATOR_H

#include "initImage.h"
#include "ImageTypes.h"
#include <sofa/core/DataEngine.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/simulation/common/AnimateBeginEvent.h>
#include <sofa/simulation/common/AnimateEndEvent.h>
#include <sofa/component/component.h>
#include <sofa/helper/system/thread/CTime.h>


namespace sofa
{

namespace component
{

namespace engine
{

using helper::vector;
using cimg_library::CImg;

/**
 * This class wraps images from a video stream into a single image
 */


template <class _ImageTypes>
class ImageAccumulator : public core::DataEngine
{
public:
    typedef core::DataEngine Inherited;
    SOFA_CLASS(SOFA_TEMPLATE(ImageAccumulator,_ImageTypes),Inherited);

    typedef _ImageTypes ImageTypes;
    typedef typename ImageTypes::T T;
    typedef typename ImageTypes::imCoord imCoord;
    typedef helper::WriteAccessor<Data< ImageTypes > > waImage;
    typedef helper::ReadAccessor<Data< ImageTypes > > raImage;

    typedef SReal Real;
    typedef defaulttype::ImageLPTransform<Real> TransformType;
    typedef typename TransformType::Coord Coord;
    typedef helper::WriteAccessor<Data< TransformType > > waTransform;
    typedef helper::ReadAccessor<Data< TransformType > > raTransform;

    typedef sofa::helper::system::thread::ctime_t ctime_t;
    typedef sofa::helper::system::thread::CTime CTime;

    Data< bool > accumulate;
    Data< ImageTypes > inputImage;
    Data< TransformType > inputTransform;
    Data< ImageTypes > outputImage;
    Data< TransformType > outputTransform;

    virtual std::string getTemplateName() const    { return templateName(this);    }
    static std::string templateName(const ImageAccumulator<ImageTypes>* = NULL) { return ImageTypes::Name(); }

    ImageAccumulator()    :   Inherited()
        , accumulate(initData(&accumulate,false,"accumulate","accumulate ?"))
        , inputImage(initData(&inputImage,ImageTypes(),"inputImage",""))
        , inputTransform(initData(&inputTransform,TransformType(),"inputTransform",""))
        , outputImage(initData(&outputImage,ImageTypes(),"outputImage",""))
        , outputTransform(initData(&outputTransform,TransformType(),"outputTransform",""))
        , SimuTime(0.0)
        , count(0.0)
    {
        inputImage.setReadOnly(true);
        inputTransform.setReadOnly(true);
        outputImage.setReadOnly(true);
        outputTransform.setReadOnly(true);
        f_listening.setValue(true);
    }

    virtual ~ImageAccumulator() {}

    virtual void init()
    {
        addInput(&inputImage);
        addInput(&inputTransform);
        addOutput(&outputImage);
        addOutput(&outputTransform);
        setDirtyValue();
    }

    virtual void reinit() { update(); }

protected:
    double SimuTime;
    ctime_t t0,t;
    int count;

    virtual void update()
    {
        cleanDirty();
        if(SimuTime==this->getContext()->getTime()) return; // check if simutime has changed
        SimuTime=this->getContext()->getTime();

        if(!accumulate.getValue()) return;

        raImage in(this->inputImage);
        if(in->isEmpty()) return;
        waImage out(this->outputImage);
        raTransform inT(this->inputTransform);
        waTransform outT(this->outputTransform);

        if(out->isEmpty()) {t0=CTime::getTime(); outT->operator=(inT);}
        else { count++; t=CTime::getTime(); outT->getScaleT()=0.000001*(t-t0)/(Real)count; } // update time scale to fit acquisition rate

        out->getCImgList().push_back(in->getCImg(0));
    }

    void handleEvent(sofa::core::objectmodel::Event *event)
    {
        if ( /*simulation::AnimateEndEvent* ev =*/  dynamic_cast<simulation::AnimateEndEvent*>(event)) update();
    }
};


} // namespace engine

} // namespace component

} // namespace sofa

#endif // SOFA_IMAGE_IMAGEACCUMULATOR_H
