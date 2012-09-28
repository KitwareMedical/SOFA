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
*                               SOFA :: Plugins                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef OGRESHADERENTRYPOINT_H
#define OGRESHADERENTRYPOINT_H

#include <sofa/helper/fixed_array.h>
#include <sofa/core/visual/VisualModel.h>

namespace sofa
{

namespace component
{

namespace visualmodel
{

class OgreShaderEntryPoint: public core::visual::VisualModel
{
public:
    SOFA_CLASS(OgreShaderEntryPoint, core::visual::VisualModel);

    OgreShaderEntryPoint():
        techniqueIndex(initData(&techniqueIndex, 0, "techniqueIndex", "Index of the technique where we have to add the Texture Unit"))
        , passIndex(initData(&passIndex, 0, "passIndex", "Index of the pass where we have to add the Texture Unit"))
    {}

    void setTechniqueIndex(int entry) {techniqueIndex.setValue(entry);}
    int getTechniqueIndex() const {return techniqueIndex.getValue();}

    void setPassIndex(int entry) {passIndex.setValue(entry);}
    int getPassIndex() const {return passIndex.getValue();}

protected:
    Data<int> techniqueIndex;
    Data<int> passIndex;

};
}
}
}
#endif
