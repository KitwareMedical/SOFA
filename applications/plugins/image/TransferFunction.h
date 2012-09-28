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
#ifndef SOFA_IMAGE_TRANSFERFUNCTION_H
#define SOFA_IMAGE_TRANSFERFUNCTION_H

#include "initImage.h"
#include "ImageTypes.h"
#include <sofa/core/DataEngine.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/helper/rmath.h>
#include <sofa/helper/OptionsGroup.h>

#include <sofa/component/component.h>
#include <map>

#define LINEAR 0


namespace sofa
{
namespace component
{
namespace engine
{

using helper::vector;
using namespace cimg_library;

/**
 * This class transforms pixel intensities
 */


template <class _InImageTypes,class _OutImageTypes>
class TransferFunction : public core::DataEngine
{
public:
    typedef core::DataEngine Inherited;
    SOFA_CLASS(SOFA_TEMPLATE2(TransferFunction,_InImageTypes,_OutImageTypes),Inherited);

    typedef _InImageTypes InImageTypes;
    typedef typename InImageTypes::T Ti;
    typedef typename InImageTypes::imCoord imCoord;
    typedef helper::WriteAccessor<Data< InImageTypes > > waImagei;
    typedef helper::ReadAccessor<Data< InImageTypes > > raImagei;

    typedef _OutImageTypes OutImageTypes;
    typedef typename OutImageTypes::T To;
    typedef helper::WriteAccessor<Data< OutImageTypes > > waImageo;
    typedef helper::ReadAccessor<Data< OutImageTypes > > raImageo;

    typedef std::map<Ti,To> iomap;
    typedef typename iomap::const_iterator iomapit;


    typedef vector<double> ParamTypes;
    typedef helper::ReadAccessor<Data< ParamTypes > > raParam;

    Data<helper::OptionsGroup> filter;
    Data< ParamTypes > param;

    Data< InImageTypes > inputImage;

    Data< OutImageTypes > outputImage;

    virtual std::string getTemplateName() const    { return templateName(this);    }
    static std::string templateName(const TransferFunction<InImageTypes,OutImageTypes>* = NULL) { return InImageTypes::Name()+std::string(",")+OutImageTypes::Name(); }

    TransferFunction()    :   Inherited()
        , filter ( initData ( &filter,"filter","Filter" ) )
        , param ( initData ( &param,"param","Parameters" ) )
        , inputImage(initData(&inputImage,InImageTypes(),"inputImage",""))
        , outputImage(initData(&outputImage,OutImageTypes(),"outputImage",""))
    {
        inputImage.setReadOnly(true);
        outputImage.setReadOnly(true);
        helper::OptionsGroup filterOptions(1	,"0 - Piecewise Linear ( i1, o1, i2, o2 ...)"
                                          );
        filterOptions.setSelectedItem(LINEAR);
        filter.setValue(filterOptions);
    }

    virtual ~TransferFunction() {}

    virtual void init()
    {
        addInput(&inputImage);
        addOutput(&outputImage);
        setDirtyValue();
    }

    virtual void reinit() { update(); }

protected:

    virtual void update()
    {
        cleanDirty();
        raParam p(this->param);
        raImagei in(this->inputImage);
        const CImgList<Ti>& inimg = in->getCImgList();

        waImageo out(this->outputImage);
        imCoord dim=in->getDimensions();
        out->setDimensions(dim);
        CImgList<To>& img = out->getCImgList();

        switch(this->filter.getValue().getSelectedId())
        {
        case LINEAR:
        {
            iomap mp; for(unsigned int i=0; i<p.size(); i+=2) mp[(Ti)p[i]]=(To)p[i+1];
            cimglist_for(inimg,l) cimg_forXYZC(inimg(l),x,y,z,c) img(l)(x,y,z,c)=Linear_TransferFunction(inimg(l)(x,y,z,c),mp);
        }
        break;

        default:
            img.assign(in->getCImgList());	// copy

            break;
        }

    }


    inline To Linear_TransferFunction(const Ti& vi, const iomap & mp) const
    {
        To vo=mp.begin()->second;
        iomapit mit;
        for (iomapit it=mp.begin(); it!=mp.end(); it++)
        {
            if (it->first>vi && it!=mp.begin())
            {
                double alpha=(((double)it->first-(double)vi)/((double)it->first-(double)mit->first));
                double v= alpha*(double)mit->second + (1.-alpha)*(double)it->second;
                return (To)v;
            }
            else vo=it->second;
            mit=it;
        }
        return vo;
    }


};





} // namespace engine
} // namespace component
} // namespace sofa

#endif // SOFA_IMAGE_TRANSFERFUNCTION_H
