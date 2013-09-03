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
#ifndef OGRESHADERTEXTUREUNIT_H
#define OGRESHADERTEXTUREUNIT_H

#include <sofa/core/objectmodel/DataFileName.h>
#include "OgreShaderEntryPoint.h"

namespace sofa
{

namespace component
{

namespace visualmodel
{

class OgreShaderTextureUnit : public OgreShaderEntryPoint
{
public:
    SOFA_CLASS(OgreShaderTextureUnit, OgreShaderEntryPoint);

    OgreShaderTextureUnit():
        textureIndex(initData(&textureIndex, 0, "textureIndex", "Index of the texture in the pass"))
        , textureName(initData(&textureName,"textureName", "File to use for the Texture"))
    {
    };
    virtual ~OgreShaderTextureUnit() {};

    void setTextureIndex(int entry) {textureIndex.setValue(entry);}
    int getTextureIndex() const {return textureIndex.getValue();}

    void setTextureName(const std::string &filename) {textureName.setValue(filename);}
    const std::string &getTextureName() const {return textureName.getValue();}

protected:
    Data<int> textureIndex;
    core::objectmodel::DataFileName textureName;

};
}
}
}

#endif

