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
#ifndef SOFA_CORE_COLLISION_CONTACTMANAGER_H
#define SOFA_CORE_COLLISION_CONTACTMANAGER_H

#include <sofa/core/collision/CollisionAlgorithm.h>
#include <sofa/core/collision/Contact.h>
#include <sofa/core/collision/NarrowPhaseDetection.h>

#include <vector>

namespace sofa
{

namespace core
{

namespace collision
{
/**
 * @brief Given a set of detected contact points, create contact response components
 */

class ContactManager : public virtual CollisionAlgorithm
{
public:
    SOFA_ABSTRACT_CLASS(ContactManager, CollisionAlgorithm);

    typedef NarrowPhaseDetection::DetectionOutputMap DetectionOutputMap;
    typedef sofa::helper::vector<Contact::SPtr> ContactVector;
protected:
    /// Constructor
    ContactManager()
        : intersectionMethod(NULL)
    {
    }

    /// Destructor
    virtual ~ContactManager() { }
public:
    virtual void createContacts(DetectionOutputMap& outputs) = 0;

    virtual const ContactVector& getContacts() { return contacts; }

    /// virtual methods used for cleaning the pipeline after a dynamic graph node deletion.
    /**
     * Contacts can be attached to a deleted node and their deletion is a problem for the pipeline.
     * @param c is the list of deleted contacts.
     */
    virtual void removeContacts(const ContactVector &/*c*/) { }

    /// virtual because subclasses might do precomputations based on intersection algorithms
    virtual void setIntersectionMethod(Intersection* v) { intersectionMethod = v;    }
    Intersection* getIntersectionMethod() const         { return intersectionMethod; }

    virtual std::string getContactResponse(core::CollisionModel* model1, core::CollisionModel* model2)=0;

protected:
    /// Current intersection method
    Intersection* intersectionMethod;

    ContactVector contacts;


    /// All intersection methods
    std::map<Instance,Intersection*> storedIntersectionMethod;

    std::map<Instance,ContactVector> storedContacts;

    virtual void changeInstance(Instance inst)
    {
        storedIntersectionMethod[instance] = intersectionMethod;
        intersectionMethod = storedIntersectionMethod[inst];
        storedContacts[instance].swap(contacts);
        contacts.swap(storedContacts[inst]);
    }
};

} // namespace collision

} // namespace core

} // namespace sofa

#endif
