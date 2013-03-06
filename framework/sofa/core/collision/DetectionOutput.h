/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *
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
*                              SOFA :: Framework                              *
*                                                                             *
* Authors: The SOFA Team (see Authors.txt)                                    *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_CORE_COLLISION_DETECTIONOUTPUT_H
#define SOFA_CORE_COLLISION_DETECTIONOUTPUT_H

#include <sofa/helper/system/config.h>
#include <sofa/core/CollisionElement.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/helper/vector.h>
#include <iostream>

namespace sofa
{

namespace core
{

namespace collision
{

using sofa::defaulttype::Vector3;

// uncomment if you want to use the freePoint information
// #define DETECTIONOUTPUT_FREEMOTION
#define DETECTIONOUTPUT_BARYCENTRICINFO



/**
 *  \brief Abstract description of a set of contact point.
 */

class DetectionOutputVector
{
protected:
    virtual ~DetectionOutputVector() {}
public:
    /// Clear the content of this vector
    virtual void clear() = 0;
    /// Current size (number of detected contacts
    virtual unsigned int size() const = 0;
    /// Test if the vector is empty
    bool empty() const { return size()==0; }
    /// Delete this vector from memory once the contact pair is no longer active
    virtual void release() { delete this; }
};


/**
 *  \brief Generic description of a contact point, used for most collision models except special cases such as GPU-based collisions.
 *
 *  Each contact point is described by :
 *
 *  \item elem: pair of colliding elements.
 *  \item id: unique id of the contact for the given pair of collision models.
 *  \item point: contact points on the surface of each model.
 *  \item normal: normal of the contact, pointing outward from the first model.
 *  \item value: signed distance (negative if objects are interpenetrating).
 *  \item deltaT: estimated of time of contact.
 *
 *  The contact id is used to filter redundant contacts (only the contact with
 *  the smallest distance is kept), and to store persistant data over time for
 *  the response.
 *
 */

class DetectionOutput
{
public:
    /// Pair of colliding elements.
    std::pair<core::CollisionElementIterator, core::CollisionElementIterator> elem;
    typedef int64_t ContactId;
    /// Unique id of the contact for the given pair of collision models.
    ContactId id;
    /// Contact points on the surface of each model. They are expressed in the local coordinate system of the model if any is defined..
    Vector3 point[2];
#ifdef DETECTIONOUTPUT_FREEMOTION
    Vector3 freePoint[2]; ///< free Point in contact on each element
#endif
#ifdef DETECTIONOUTPUT_BARYCENTRICINFO
    Vector3 baryCoords[2]; ///< provides the barycentric Coordinates (alpha, beta, gamma) of each contact points over the element
    ///< (alpha is used for lines / alpha and beta for triangles / alpha, beta and gamma for tetrahedrons)
#endif

    /// Normal of the contact, pointing outward from the first model
    Vector3 normal;
    /*
    /// Signed distance (negative if objects are interpenetrating). If using a proximity-based detection, this is the actual distance between the objets minus the specified contact distance.
    */
    /// Store information for the collision Response. Depending on the kind of contact, can be a distance, or a pression, ...
    double value;
    /// If using a continuous collision detection, estimated of time of contact.
    double deltaT;
    DetectionOutput()
        : elem( (sofa::core::CollisionModel* )NULL,
                (sofa::core::CollisionModel* ) NULL), id(0), value(0.0), deltaT(0.0)
    {
    }
};



/**
 *  \brief Generic description of a set of contact point between two specific collision models
 */

template<class CM1, class CM2>
class TDetectionOutputVector : public DetectionOutputVector, public sofa::helper::vector<DetectionOutput>
{
public:
    typedef sofa::helper::vector<DetectionOutput> Vector;
    virtual ~TDetectionOutputVector() {}
    /// Clear the content of this vector
    virtual void clear()
    {
        return this->Vector::clear();
    }
    /// Current size (number of detected contacts)
    virtual unsigned int size() const
    {
        return this->Vector::size();
    }
};

} // namespace collision

} // namespace core

} // namespace sofa

#endif
