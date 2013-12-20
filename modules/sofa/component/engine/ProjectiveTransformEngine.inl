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
#ifndef SOFA_COMPONENT_ENGINE_PROJECTIVETRANSFORMENGINE_INL
#define SOFA_COMPONENT_ENGINE_PROJECTIVETRANSFORMENGINE_INL

#include <sofa/component/engine/ProjectiveTransformEngine.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/rmath.h> //M_PI

namespace sofa
{

namespace component
{

namespace engine
{

template <class DataTypes>
ProjectiveTransformEngine<DataTypes>::ProjectiveTransformEngine()
    : f_inputX ( initData (&f_inputX, "input_position", "input array of 3d points") )
    , f_outputX( initData (&f_outputX, "output_position", "output array of projected 3d points") )
    , proj_mat(initData(&proj_mat, "proj_mat", "projection matrix ") )
    , focal_distance(initData(&focal_distance, (Real)1,"focal_distance", "focal distance ") )
{
}


template <class DataTypes>
void ProjectiveTransformEngine<DataTypes>::init()
{
    addInput(&f_inputX);
    addInput(&focal_distance);
    addOutput(&f_outputX);
    setDirtyValue();
}

template <class DataTypes>
void ProjectiveTransformEngine<DataTypes>::reinit()
{
    update();
}

template <class DataTypes>
void ProjectiveTransformEngine<DataTypes>::update()
{
    cleanDirty();

    helper::ReadAccessor< Data<VecCoord> > in = f_inputX;
    helper::WriteAccessor< Data<VecCoord> > out = f_outputX;
    helper::ReadAccessor< Data<Real> > fdist = focal_distance;
    helper::ReadAccessor< Data<ProjMat> > pmat = proj_mat;

    out.resize(in.size());
    double f=double(fdist.ref());
    double r,s;    
    ProjMat P=ProjMat(pmat.ref());

    unsigned j,k;
    for (j=0;j<3; j++)
    	{
	for (k=0; k<4; k++)
		std::cout << P[j][k] << " ";
	std::cout << std::endl;
	}

    unsigned int i;
    for (i=0; i< in.size(); ++i)
        {
	out[i]=Vec3(in[i][0],in[i][1],in[i][2]);
       	    out[i]=P*Vec4(out[i],1);
	    s = out[i][2];
	    if (fabs(s) < 1e-10) s=s<0 ? -1e-10 : 1e-10; // handle undefined case where out[i][2] == 0 -> set it to 1e-10 (and keep its sign)
	    r = f/out[i][2]; 
	    out[i] *= r;
	std::cout << in[i] << " <-> " << out[i] << std::endl;
        }
}



} // namespace engine

} // namespace component

} // namespace sofa

#endif
