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
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_COMPONENT_ENGINE_TRANSFORMPOSITION_INL
#define SOFA_COMPONENT_ENGINE_TRANSFORMPOSITION_INL

#if !defined(__GNUC__) || (__GNUC__ > 3 || (_GNUC__ == 3 && __GNUC_MINOR__ > 3))
#pragma once
#endif

#include <sofa/component/engine/TransformPosition.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/gl/template.h>
#include <math.h>
#include <sofa/helper/RandomGenerator.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <boost/algorithm/string.hpp>

namespace sofa
{

namespace component
{

namespace engine
{

using namespace sofa::helper;
using namespace sofa::defaulttype;
using namespace core::objectmodel;

template <class DataTypes>
TransformPosition<DataTypes>::TransformPosition()
    : f_origin( initData(&f_origin, "origin", "a 3d point on the plane") )
    , f_inputX( initData (&f_inputX, "input_position", "input array of 3d points") )
    , f_outputX( initData (&f_outputX, "output_position", "output array of 3d points projected on a plane") )
    , f_normal(initData(&f_normal, "normal", "plane normal") )
    , f_translation(initData(&f_translation, "translation", "translation vector ") )
    , f_rotation(initData(&f_rotation, "rotation", "rotation vector ") )
    , f_scale(initData(&f_scale, Coord(1.0,1.0,1.0), "scale", "scale factor") )
    , f_affineMatrix(initData(&f_affineMatrix, AffineMatrix(), "matrix", "4x4 affine matrix") )
    , f_method(initData(&f_method, "method", "transformation method either translation or scale or rotation or random or projectOnPlane") )
    , f_seed(initData(&f_seed, (long) 0, "seedValue", "the seed value for the random generator") )
    , f_maxRandomDisplacement(initData(&f_maxRandomDisplacement, (Real) 1.0, "maxRandomDisplacement", "the maximum displacement around initial position for the random transformation") )
    , f_fixedIndices( initData(&f_fixedIndices,"fixedIndices","Indices of the entries that are not transformed") )
    , f_filename(initData(&f_filename, "filename", "filename of an affine matrix. Supported extensions are: .trm, .tfm, .xfm and .txt(read as .xfm)") )
    , mstate(NULL), x0(NULL)
{
    f_method.beginEdit()->setNames(9,
        "projectOnPlane",
        "translation",
        "rotation",
        "random",
        "scale",
        "scaleTranslation",
        "scaleRotationTranslation",
        "affine",
        "fromFile");
    f_method.endEdit();
}

template <class DataTypes>
void TransformPosition<DataTypes>::selectTransformationMethod()
{
    if (f_method.getValue().getSelectedItem()=="projectOnPlane")
    {
        transformationMethod=PROJECT_ON_PLANE;
    }
    else if (f_method.getValue().getSelectedItem()=="translation")
    {
        transformationMethod=TRANSLATION;
    }
    else if (f_method.getValue().getSelectedItem()=="rotation")
    {
        transformationMethod=ROTATION;
    }
    else if (f_method.getValue().getSelectedItem()=="random")
    {
        transformationMethod=RANDOM;
    }
    else if (f_method.getValue().getSelectedItem()=="scale")
    {
        transformationMethod=SCALE;
    }
    else if (f_method.getValue().getSelectedItem()=="scaleTranslation")
    {
        transformationMethod=SCALE_TRANSLATION;
    }
    else if (f_method.getValue().getSelectedItem()=="scaleRotationTranslation")
    {
        transformationMethod=SCALE_ROTATION_TRANSLATION;
    }
    else if (f_method.getValue().getSelectedItem()=="affine")
    {
        transformationMethod=AFFINE;
    }
    else if (f_method.getValue().getSelectedItem()=="fromFile")
    {
        transformationMethod=AFFINE;
        if (f_filename.isSet())
        {
            std::string fname = f_filename.getValue();
            if (fname.size()>=4 && fname.substr(fname.size()-4)==".trm")
                getTransfoFromTrm();
            else if (fname.size()>=4 && (fname.substr(fname.size()-4)==".txt" || fname.substr(fname.size()-4)==".xfm"))
                getTransfoFromTxt();
            if (fname.size()>=4 && fname.substr(fname.size()-4)==".tfm")
                getTransfoFromTfm();
            else
                serr << "Unknown extension. Will use affine instead." << sendl;
        }
        else
        {
            serr << "Filename not set. Will use affine instead" <<sendl;
        }
    }
    else
    {
        transformationMethod=TRANSLATION;
        serr << "Error : Method " << f_method.getValue().getSelectedItem() << " is unknown. Wil use translation instead." <<sendl;
    }
}

template <class DataTypes>
void TransformPosition<DataTypes>::init()
{
    Coord& normal = *(f_normal.beginEdit());

    /// check if the normal is of norm 1
    if (fabs((normal.norm2()-1.0))>1e-10)
        normal/=normal.norm();

    f_normal.endEdit();

    addInput(&f_inputX);
    addOutput(&f_outputX);

    setDirtyValue();
}

template <class DataTypes>
void TransformPosition<DataTypes>::reinit()
{
    update();
}

/**************************************************
 * .tfm spec:
 * 12 values in the lines begining by "Parameters"
 **************************************************/
template <class DataTypes>
void TransformPosition<DataTypes>::getTransfoFromTfm()
{
    std::string fname(this->f_filename.getFullPath());
    sout << "Loading .tfm file " << fname << sendl;

    std::ifstream stream(fname.c_str());
    if (stream)
    {
        std::string line;
        AffineMatrix mat;

        bool found = false;
        while (getline(stream,line) && !found)
        {
            if (line.find("Parameters")!=std::string::npos)
            {
                found=true;

                typedef std::vector<std::string> vecString;
                vecString vLine;
                boost::split(vLine, line, boost::is_any_of(" "), boost::token_compress_on);

                std::vector<Real> values;
                for (vecString::iterator it = vLine.begin(); it < vLine.end(); it++)
                {
                    std::string c = *it;
                    if ( c.find_first_of("1234567890.-") != std::string::npos)
                        values.push_back(atof(c.c_str()));
                }

                if (values.size() != 12)
                    serr << "Error in file " << fname << sendl;
                else
                {
                    for(unsigned int i = 0 ; i < 3; i++)
                    {
                        for (unsigned int j = 0 ; j < 3; j++)
                        {
                            mat[i][j] = values[i*3+j];//rotation matrix
                        }
                        mat[i][3] = values[values.size()-1-i];//translation
                    }
                }
            }
        }

        if (!found) serr << "Transformation not found in " << fname << sendl;
        f_affineMatrix.setValue(mat);
    }
    else
    {
        serr << "Could not open file " << fname << sendl << "Matrix set to identity" << sendl;
    }
}

/**************************************************
 * .trm spec:
 * 1st line: 3 values for translation
 * then 3 lines of 3 values for the rotation matrix
 **************************************************/
template <class DataTypes>
void TransformPosition<DataTypes>::getTransfoFromTrm()
{
    std::string fname(this->f_filename.getFullPath());
    sout << "Loading .trm file " << fname << sendl;

    std::ifstream stream(fname.c_str());
    if (stream)
    {
        std::string line;
        unsigned int nbLines = 0;
        AffineMatrix mat;

        while (getline(stream,line))
        {
            if (line == "") continue;
            nbLines++;

            if (nbLines > 4)
            {
                serr << "File with more than 4 lines" << sendl;
                break;
            }

            std::vector<std::string> vLine;
            boost::split(vLine, line, boost::is_any_of(" "), boost::token_compress_on);

            if (vLine.size()>3 )
            {
                for (unsigned int i = 3; i < vLine.size();i++)
                {
                    if (vLine[i]!="")
                    {
                        serr << "Should be a line of 3 values" << sendl;
                        break;
                    }
                }
            }
            else if (vLine.size()<3) {serr << "Should be a line of 3 values" << sendl;continue;}

            if (nbLines == 1)
            {
                //translation vector
                Coord tr;
                for ( unsigned int i = 0; i < std::min((unsigned int)vLine.size(),(unsigned int)3); i++)
                {
                    tr[i] = mat[i][3] = atof(vLine[i].c_str());
                }
                f_translation.setValue(tr);

            }
            else
            {
                //rotation matrix
                for ( unsigned int i = 0; i < std::min((unsigned int)vLine.size(),(unsigned int)3); i++)
                    mat[nbLines-2][i] = atof(vLine[i].c_str());
            }

        }
        f_affineMatrix.setValue(mat);
    }
    else
    {
        serr << "Could not open file " << fname << sendl << "Matrix set to identity" << sendl;
    }

}

/**************************************************
 * .txt and .xfm spec:
 * 4 lines of 4 values for an affine matrix
 **************************************************/
template <class DataTypes>
void TransformPosition<DataTypes>::getTransfoFromTxt()
{
    std::string fname(this->f_filename.getFullPath());
    sout << "Loading matrix file " << fname << sendl;

    std::ifstream stream(fname.c_str());
    if (stream)
    {
        std::string line;
        unsigned int nbLines = 0;
        AffineMatrix mat;

        while (getline(stream,line))
        {
            if (line == "") continue;
            nbLines++;

            if (nbLines > 4)
            {
                serr << "Matrix is not 4x4" << sendl;
                break;
            }

            std::vector<std::string> vLine;
            boost::split(vLine, line, boost::is_any_of(" "), boost::token_compress_on);

            if (vLine.size()>4 )
            {
                for (unsigned int i = 4; i < vLine.size();i++)
                {
                    if (vLine[i]!="")
                    {
                        serr << "Matrix is not 4x4." << sendl;
                        break;
                    }
                }
            }
            else if (vLine.size()<4) {serr << "Matrix is not 4x4." << sendl;continue;}

            for ( unsigned int i = 0; i < std::min((unsigned int)vLine.size(),(unsigned int)4); i++)
                mat[nbLines-1][i] = atof(vLine[i].c_str());
        }
        f_affineMatrix.setValue(mat);

    }
    else
    {
        serr << "Could not open file " << fname << sendl << "Matrix set to identity" << sendl;
    }
}


template <class DataTypes>
void TransformPosition<DataTypes>::update()
{
    cleanDirty();

    selectTransformationMethod();

    helper::ReadAccessor< Data<VecCoord> > in = f_inputX;
    helper::WriteAccessor< Data<VecCoord> > out = f_outputX;
    helper::ReadAccessor< Data<Coord> > normal = f_normal;
    helper::ReadAccessor< Data<Coord> > origin = f_origin;
    helper::ReadAccessor< Data<Coord> > translation = f_translation;
    helper::ReadAccessor< Data<Coord> > scale = f_scale;
    helper::ReadAccessor< Data<Coord> > rotation = f_rotation;
    helper::ReadAccessor< Data<Real> > maxDisplacement = f_maxRandomDisplacement;
    helper::ReadAccessor< Data<long> > seed = f_seed;
    helper::ReadAccessor< Data<SetIndex> > fixedIndices = f_fixedIndices;

    out.resize(in.size());
    unsigned int i;
    switch(transformationMethod)
    {
    case PROJECT_ON_PLANE :
        for (i=0; i< in.size(); ++i)
        {
            out[i]=in[i]+normal.ref()*dot((origin.ref()-in[i]),normal.ref());
        }
        break;
    case RANDOM :
    {
        sofa::helper::RandomGenerator rg;
        double dis=(double) maxDisplacement.ref();
        if (seed.ref()!=0)
            rg.initSeed(seed.ref());
        for (i=0; i< in.size(); ++i)
        {

            out[i]=in[i]+Coord((Real)rg.random<double>(-dis,dis),(Real)rg.random<double>(-dis,dis),(Real)rg.random<double>(-dis,dis));
        }
    }
    break;
    case TRANSLATION :
        for (i=0; i< in.size(); ++i)
        {
            out[i]=in[i]+translation.ref();
        }
        break;
    case SCALE :
        for (i=0; i< in.size(); ++i)
        {
            out[i]=in[i].linearProduct(scale.ref());
        }
        break;
    case SCALE_TRANSLATION :
        for (i=0; i< in.size(); ++i)
        {
            out[i]=in[i].linearProduct(scale.ref()) +translation.ref();
        }
        break;
    case ROTATION :
    {
        Quaternion q=helper::Quater<Real>::createQuaterFromEuler( rotation.ref()*M_PI/180.0);

        for (i=0; i< in.size(); ++i)
        {
            out[i]=q.rotate(in[i]);
        }
    }
    break;
    case SCALE_ROTATION_TRANSLATION :
    {
        Quaternion q=helper::Quater<Real>::createQuaterFromEuler( rotation.ref()*M_PI/180.0);

        for (i=0; i< in.size(); ++i)
        {
            out[i]=q.rotate(in[i].linearProduct(scale.ref())) +translation.ref();
        }
        break;
    }
    case AFFINE:
        for (i=0; i< in.size(); ++i)
        {
            Vec4 coord = f_affineMatrix.getValue()*Vec4(in[i], 1);
            if ( fabs(coord[3]) > 1e-10)
                out[i]=coord/coord[3];
        }
        break;
    }
    /// assumes the set of fixed indices is small compared to the whole set
    SetIndex::const_iterator it=fixedIndices.ref().begin();
    for (; it!=fixedIndices.ref().end(); ++it)
    {
        out[*it]=in[*it];
    }

}

template <class DataTypes>
void TransformPosition<DataTypes>::draw()
{

}


} // namespace engine

} // namespace component

} // namespace sofa

#endif
