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
#ifndef FLEXIBLE_ShapeFunctionDiscretizer_H
#define FLEXIBLE_ShapeFunctionDiscretizer_H

#include "../initFlexible.h"
#include "../shapeFunction/BaseShapeFunction.h"
#include <image/ImageTypes.h>

#include <sofa/core/DataEngine.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/defaulttype/Vec.h>

#include <sofa/helper/rmath.h>

#include <sofa/component/component.h>


namespace sofa
{
namespace component
{
namespace engine
{

using helper::vector;
using cimg_library::CImg;

/**
 * This class discretize shape functions in an image
 */


template <class ImageTypes_>
class ShapeFunctionDiscretizer : public core::DataEngine
{
public:
    typedef core::DataEngine Inherited;
    SOFA_CLASS(SOFA_TEMPLATE(ShapeFunctionDiscretizer,ImageTypes_),Inherited);

    typedef SReal Real;

    /** @name  Image data */
    //@{
    typedef ImageTypes_ ImageTypes;
    typedef typename ImageTypes::T T;
    typedef typename ImageTypes::imCoord imCoord;
    typedef helper::ReadAccessor<Data< ImageTypes > > raImage;
    Data< ImageTypes > f_image;

    typedef defaulttype::ImageLPTransform<Real> TransformType;
    typedef helper::ReadAccessor<Data< TransformType > > raTransform;
    Data< TransformType > f_transform;

    typedef unsigned int IndT;
    typedef defaulttype::Image<IndT> IndTypes;
    typedef helper::WriteAccessor<Data< IndTypes > > waInd;
    Data< IndTypes > f_index;

    typedef double DistT;
    typedef defaulttype::Image<DistT> DistTypes;
    typedef helper::WriteAccessor<Data< DistTypes > > waDist;
    Data< DistTypes > f_w;
    //@}

    /** @name  Shape Function types    */
    //@{
    enum { spatial_dimensions = 3 };
    enum { material_dimensions = 3 };
    typedef core::behavior::ShapeFunctionTypes<spatial_dimensions,material_dimensions,Real> ShapeFunctionType;
    typedef core::behavior::BaseShapeFunction<ShapeFunctionType> BaseShapeFunction;
    typedef typename BaseShapeFunction::VReal VReal;
    typedef typename BaseShapeFunction::VRef VRef;
    typedef typename BaseShapeFunction::MaterialToSpatial MaterialToSpatial ;
    typedef typename BaseShapeFunction::Coord Coord;
    BaseShapeFunction* _shapeFunction;        ///< where the weights are computed
    //@}


    virtual std::string getTemplateName() const    { return templateName(this);    }
    static std::string templateName(const ShapeFunctionDiscretizer<ImageTypes>* = NULL) { return ImageTypes::Name(); }

    ShapeFunctionDiscretizer()    :   Inherited()
      , f_image(initData(&f_image,ImageTypes(),"image",""))
      , f_transform(initData(&f_transform,TransformType(),"transform",""))
      , f_index(initData(&f_index,IndTypes(),"indices",""))
      , f_w(initData(&f_w,DistTypes(),"weights",""))
      , _shapeFunction(NULL)
    {
        f_image.setReadOnly(true);
        f_transform.setReadOnly(true);
    }

    virtual ~ShapeFunctionDiscretizer() {}

    virtual void init()
    {
        if( !_shapeFunction ) this->getContext()->get(_shapeFunction,core::objectmodel::BaseContext::SearchUp);
        if ( !_shapeFunction ) serr << "ShapeFunction<"<<ShapeFunctionType::Name()<<"> component not found" << sendl;

        addInput(&f_image);
        addInput(&f_transform);
        addOutput(&f_w);
        addOutput(&f_index);
        setDirtyValue();
    }

    virtual void reinit() { update(); }

protected:

    virtual void update()
    {
        cleanDirty();

        if( !_shapeFunction ) return;

        // read input image and transform
        raImage in(this->f_image);
        raTransform inT(this->f_transform);
        if(in->isEmpty())  { serr<<"Image not found"<<sendl; return; }
        const CImg<T>& inimg = in->getCImg(0);  // suppose time=0

        // init indices and weights images
        const unsigned int nbref=_shapeFunction->f_nbRef.getValue();
        imCoord dim = in->getDimensions();
        dim[3]=nbref;
        dim[4]=1;

        waInd indData(this->f_index); indData->setDimensions(dim);
        CImg<IndT>& indices = indData->getCImg(); indices.fill(0);
        waDist weightData(this->f_w);         weightData->setDimensions(dim);
        CImg<DistT>& weights = weightData->getCImg(); weights.fill(0);

        // fill indices and weights images
#ifdef USING_OMP_PRAGMAS
#pragma omp parallel for
#endif
        for(int z=0; z<inimg.depth(); z++)
            for(int y=0; y<inimg.height(); y++)
                for(int x=0; x<inimg.width(); x++)
                    if(inimg(x,y,z))
                    {
                        Coord p=inT->fromImage(Coord(x,y,z));
                        VReal w; VRef ref; MaterialToSpatial M;
                        _shapeFunction->computeShapeFunction(p,M,ref,w);
                        for(unsigned int i=0;i<nbref;i++)
                        {
                            indices(x,y,z,i)=ref[i]+1;
                            weights(x,y,z,i)=w[i];
                        }
                    }

    }

};


} // namespace engine

} // namespace component

} // namespace sofa

#endif // SOFA_IMAGE_ShapeFunctionDiscretizer_H
