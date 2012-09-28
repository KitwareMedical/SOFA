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
#ifndef SOFA_COMPONENT_CONTAINER_INTERPOLATEDIMPLICITSURFACE_H
#define SOFA_COMPONENT_CONTAINER_INTERPOLATEDIMPLICITSURFACE_H


#include "ImplicitSurfaceContainer.h"
#include <sofa/core/objectmodel/DataFileName.h>


namespace sofa
{
namespace component
{
namespace container
{
using namespace sofa::defaulttype;



class DomainCache
{
public:
    bool insideImg; // shows if the domain lies inside the valid image region or outside
    Vec3d bbMin, bbMax; // bounding box (min and max) of the domain
    double val[8]; // corner values of the domain
};



class SOFA_VOLUMETRIC_DATA_API InterpolatedImplicitSurface : public virtual ImplicitSurface
{

public:

    SOFA_CLASS(InterpolatedImplicitSurface,ImplicitSurface);
protected:
    InterpolatedImplicitSurface();
    ~InterpolatedImplicitSurface();
public:
    virtual void init();
    virtual double getValue( defaulttype::Vec3d &transformedPos );
    virtual double getValue( defaulttype::Vec3d &transformedPos, int &domain );
    virtual unsigned int getDomain( sofa::defaulttype::Vec3d &pos, int ref_domain ) { (void)pos; return ref_domain; }



protected:

    bool loadImage( const char *mhdFile );
    void updateCache( DomainCache *cache, defaulttype::Vec3d& pos );
    int getNextDomain();


protected:

    sofa::core::objectmodel::DataFileName distanceMapHeader;
    Data< int > maxDomains;
    Data< double > dx, dy, dz;    // translation of original image

    int usedDomains;              // number of domains already given out
    unsigned int imgSize[3];      // number of voxels
    double spacing[3];            // physical distance between two neighboring voxels
    double scale[3];              // (1/spacing)
    double imgMin[3], imgMax[3];  // physical locations of the centers of both corner voxels
    float *imgData;               // raw data
    unsigned int deltaOfs[8];     // offsets to define 8 corners of cube for interpolation
    std::vector<DomainCache> domainCache;

};


} // namespace container
} // namespace component
} // namespace sofa

#endif
