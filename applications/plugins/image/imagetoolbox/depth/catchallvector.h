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
#ifndef SOFA_IMAGE_CATCHALLVECTOR_H
#define SOFA_IMAGE_CATCHALLVECTOR_H

#include "initImage.h"
#include <sofa/core/DataEngine.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/defaulttype/Vec.h>

#include <sofa/component/component.h>
#include <sofa/helper/OptionsGroup.h>

#include "ImageContainer.h"

namespace sofa
{
namespace component
{
namespace engine
{

template <class _Type>
class CatchAllVector : public core::DataEngine
{
public:
    typedef core::DataEngine Inherited;
    SOFA_CLASS(SOFA_TEMPLATE(CatchAllVector,_Type),Inherited);

    typedef _Type Type;

    CatchAllVector()    :   Inherited()
        ,_data(initData(&_data,"data","data"))
    {

    }

    ~CatchAllVector() {}

    void init()
    {
        reinit();
    }

    void reinit()
    {
    }

    Data< vector<Type> > _data;

protected:

    virtual void update()
    {
    }

public:

   /* virtual void draw(const core::visual::VisualParams*)
    {
    }

    void handleEvent(core::objectmodel::Event *event)
    {
    }*/
};





} // namespace engine
} // namespace component
} // namespace sofa

#endif // SOFA_IMAGE_CatchAllVector_H
