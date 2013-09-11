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
#ifndef SOFA_CORE_OBJECTMODEL_BASECONTEXT_H
#define SOFA_CORE_OBJECTMODEL_BASECONTEXT_H

#include <sofa/core/objectmodel/Base.h>
#include <sofa/core/objectmodel/BaseLink.h>
#include <sofa/core/objectmodel/Tag.h>
#include <sofa/core/objectmodel/ClassInfo.h>
#include <sofa/core/ExecParams.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <set>
#ifdef SOFA_SMP
#include <IterativePartition.h>
#endif

namespace sofa
{

namespace simulation
{
class Visitor;
}

namespace core
{

// forward declaration of classes accessible from the context

class BaseState;

namespace topology
{
class Topology;
class BaseMeshTopology;
} // namespace topology

namespace behavior
{
class BaseMechanicalState;
class BaseMass;
} // namespace behavior

namespace visual
{
class Shader;
} // namespace visual

namespace objectmodel
{
class BaseObject;
class Event;

/**
 *  \brief Base class for Context classes, storing shared variables and parameters.
 *
 *  A Context contains values or pointers to variables and parameters shared
 *  by a group of objects, typically refering to the same simulated body.
 *  Derived classes can defined simple isolated contexts or more powerful
 *  hierarchical representations (scene-graphs), in which case the context also
 *  implements the BaseNode interface.
 *
 * \author Jeremie Allard
 */
class SOFA_CORE_API BaseContext : public virtual Base
{
public:
    SOFA_CLASS(BaseContext, Base);

    /// @name Types defined for local coordinate system handling
    /// @{

    typedef defaulttype::SolidTypes<SReal> SolidTypes;

    typedef SolidTypes::Transform Frame;
    typedef SolidTypes::Vec Vec3;
    typedef SolidTypes::Rot Quat;
    typedef SolidTypes::Mat Mat33;
    typedef SolidTypes::SpatialVector SpatialVector;
    /// @}
protected:
    BaseContext();
    virtual ~BaseContext();
public:
    /// Get the default Context object, that contains the default values for
    /// all parameters and can be used when no local context is defined.
    static BaseContext* getDefault();

    /// Specification of where to search for queried objects
    enum SearchDirection { SearchUp = -1, Local = 0, SearchDown = 1, SearchRoot = 2, SearchParents = 3 };

    /// @name Parameters
    /// @{

    /// The Context is active
    virtual bool isActive() const;
#ifdef SOFA_SMP
    virtual bool is_partition() const;
#endif

    /// State of the context
    virtual void setActive(bool) {};

    /// Simulation time
    virtual double getTime() const;

    /// Simulation timestep
    virtual double getDt() const;

    /// Animation flag
    virtual bool getAnimate() const;

#ifdef SOFA_SMP
    virtual int  getProcessor() const;
    virtual void  setProcessor(int) {}
    virtual Iterative::IterativePartition*  getPartition() const;
#endif

#ifdef SOFA_SUPPORT_MULTIRESOLUTION
    /// Multiresolution support (UNSTABLE)
    virtual int getCurrentLevel() const;

    /// Multiresolution support (UNSTABLE)
    virtual int getCoarsestLevel() const;

    /// Multiresolution support (UNSTABLE)
    virtual int getFinestLevel() const;

    /// Multiresolution support (UNSTABLE)
    //     virtual unsigned int nbLevels() const;
#endif

    /// @}

#ifdef SOFA_SUPPORT_MOVING_FRAMES
    /// @name Local Coordinate System
    /// @{
    /// Projection from the local coordinate system to the world coordinate system.
    virtual const Frame& getPositionInWorld() const;
    /// Projection from the local coordinate system to the world coordinate system.
    virtual void setPositionInWorld(const Frame&)
    {}

    /// Spatial velocity (linear, angular) of the local frame with respect to the world
    virtual const SpatialVector& getVelocityInWorld() const;
    /// Spatial velocity (linear, angular) of the local frame with respect to the world
    virtual void setVelocityInWorld(const SpatialVector&)
    {}

    /// Linear acceleration of the origin induced by the angular velocity of the ancestors
    virtual const Vec3& getVelocityBasedLinearAccelerationInWorld() const;
    /// Linear acceleration of the origin induced by the angular velocity of the ancestors
    virtual void setVelocityBasedLinearAccelerationInWorld(const Vec3& )
    {}
    /// Gravity in local coordinates  TODO: replace with world coordinates
    virtual Vec3 getLocalGravity() const;
    ///// Gravity in local coordinates
    //virtual void setGravity( const Vec3& ) { }
    /// @}
#endif

    /// Gravity in local coordinates
    virtual const Vec3& getGravity() const;
    /// Gravity in local coordinates
    virtual void setGravity( const Vec3& )
    { }

    /// Get the root context of the graph
    virtual BaseContext* getRootContext() const;

    /// @name Containers
    /// @{

    /// Mechanical Degrees-of-Freedom
    virtual core::BaseState* getState() const;

    /// Mechanical Degrees-of-Freedom
    virtual behavior::BaseMechanicalState* getMechanicalState() const;

    /// Topology
    virtual core::topology::Topology* getTopology() const;

    /// Mesh Topology (unified interface for both static and dynamic topologies)
    virtual core::topology::BaseMeshTopology* getMeshTopology() const;

    /// Mass
    virtual core::behavior::BaseMass* getMass() const;

    /// Global Shader
    virtual core::visual::Shader* getShader() const;

    /// Generic object access, possibly searching up or down from the current context
    ///
    /// Note that the template wrapper method should generally be used to have the correct return type,
    virtual void* getObject(const ClassInfo& class_info, SearchDirection dir = SearchUp) const;

    /// Generic object access, given a set of required tags, possibly searching up or down from the current context
    ///
    /// Note that the template wrapper method should generally be used to have the correct return type,
    virtual void* getObject(const ClassInfo& class_info, const TagSet& tags, SearchDirection dir = SearchUp) const;

    /// Generic object access, given a path from the current context
    ///
    /// Note that the template wrapper method should generally be used to have the correct return type,
    virtual void* getObject(const ClassInfo& class_info, const std::string& path) const;

    class GetObjectsCallBack
    {
    public:
        virtual ~GetObjectsCallBack() {}
        virtual void operator()(void* ptr) = 0;
    };

    /// Generic list of objects access, possibly searching up or down from the current context
    ///
    /// Note that the template wrapper method should generally be used to have the correct return type,
    virtual void getObjects(const ClassInfo& class_info, GetObjectsCallBack& container, SearchDirection dir = SearchUp) const;

    /// Generic list of objects access, given a set of required tags, possibly searching up or down from the current context
    ///
    /// Note that the template wrapper method should generally be used to have the correct return type,
    virtual void getObjects(const ClassInfo& class_info, GetObjectsCallBack& container, const TagSet& tags, SearchDirection dir = SearchUp) const;


    /// Generic object access template wrapper, possibly searching up or down from the current context
    template<class T>
    T* get(SearchDirection dir = SearchUp) const
    {
        return reinterpret_cast<T*>(this->getObject(classid(T), dir));
    }


    /// Generic object access template wrapper, possibly searching up or down from the current context
    template<class T>
    void get(T*& ptr, SearchDirection dir = SearchUp) const
    {
        ptr = this->get<T>(dir);
    }

    /// Generic object access template wrapper, possibly searching up or down from the current context
    template<class T>
    void get(boost::intrusive_ptr<T>& ptr, SearchDirection dir = SearchUp) const
    {
        ptr = this->get<T>(dir);
    }

    /// Generic object access template wrapper, given a required tag, possibly searching up or down from the current context
    template<class T>
    T* get(const Tag& tag, SearchDirection dir = SearchUp) const
    {
        return reinterpret_cast<T*>(this->getObject(classid(T), TagSet(tag), dir));
    }

    /// Generic object access template wrapper, given a required tag, possibly searching up or down from the current context
    template<class T>
    void get(T*& ptr, const Tag& tag, SearchDirection dir = SearchUp) const
    {
        ptr = this->get<T>(tag, dir);
    }

    /// Generic object access template wrapper, given a required tag, possibly searching up or down from the current context
    template<class T>
    void get(boost::intrusive_ptr<T>& ptr, const Tag& tag, SearchDirection dir = SearchUp) const
    {
        ptr = this->get<T>(tag, dir);
    }

    /// Generic object access template wrapper, given a set of required tags, possibly searching up or down from the current context
    template<class T>
    T* get(const TagSet& tags, SearchDirection dir = SearchUp) const
    {
        return reinterpret_cast<T*>(this->getObject(classid(T), tags, dir));
    }

    /// Generic object access template wrapper, given a set of required tags, possibly searching up or down from the current context
    template<class T>
    void get(T*& ptr, const TagSet& tags, SearchDirection dir = SearchUp) const
    {
        ptr = this->get<T>(tags, dir);
    }

    /// Generic object access template wrapper, given a set of required tags, possibly searching up or down from the current context
    template<class T>
    void get(boost::intrusive_ptr<T>& ptr, const TagSet& tags, SearchDirection dir = SearchUp) const
    {
        ptr = this->get<T>(tags, dir);
    }

    /// Generic object access template wrapper, given a path from the current context
    template<class T>
    T* get(const std::string& path) const
    {
        return reinterpret_cast<T*>(this->getObject(classid(T), path));
    }

    /// Generic object access template wrapper, given a path from the current context
    template<class T>
    void get(T*& ptr, const std::string& path) const
    {
        ptr = this->get<T>(path);
    }

    /// Generic object access template wrapper, given a path from the current context
    template<class T>
    void get(boost::intrusive_ptr<T>& ptr, const std::string& path) const
    {
        ptr = this->get<T>(path);
    }

    template<class T, class Container>
    class GetObjectsCallBackT : public GetObjectsCallBack
    {
    public:
        Container* dest;
        GetObjectsCallBackT(Container* d) : dest(d) {}
        virtual void operator()(void* ptr)
        {
            dest->push_back(reinterpret_cast<T*>(ptr));
        }
    };

    /// Generic list of objects access template wrapper, possibly searching up or down from the current context
    template<class T, class Container>
    void get(Container* list, SearchDirection dir = SearchUp) const
    {
        GetObjectsCallBackT<T,Container> cb(list);
        this->getObjects(classid(T), cb, dir);
    }

    /// Generic list of objects access template wrapper, given a required tag, possibly searching up or down from the current context
    template<class T, class Container>
    void get(Container* list, const Tag& tag, SearchDirection dir = SearchUp) const
    {
        GetObjectsCallBackT<T,Container> cb(list);
        this->getObjects(classid(T), cb, TagSet(tag), dir);
    }

    /// Generic list of objects access template wrapper, given a set of required tags, possibly searching up or down from the current context
    template<class T, class Container>
    void get(Container* list, const TagSet& tags, SearchDirection dir = SearchUp) const
    {
        GetObjectsCallBackT<T,Container> cb(list);
        this->getObjects(classid(T), cb, tags, dir);
    }

    /// @}

    /// @name Parameters Setters
    /// @{


    /// Simulation timestep
    virtual void setDt( double /*dt*/ )
    { }

    /// Animation flag
    virtual void setAnimate(bool /*val*/)
    { }

#ifdef SOFA_SUPPORT_MULTIRESOLUTION
    /// Multiresolution support (UNSTABLE) : Set the current level, return false if l >= coarsestLevel
    virtual bool setCurrentLevel(int )
    {
        return false;
    }

    /// Multiresolution support (UNSTABLE)
    virtual void setCoarsestLevel(int ) {}

    /// Multiresolution support (UNSTABLE)
    virtual void setFinestLevel(int ) {}
#endif

    /// @}

    /// @name Variables Setters
    /// @{

    /// Mechanical Degrees-of-Freedom
    virtual void setMechanicalState( BaseObject* )
    { }

    /// Topology
    virtual void setTopology( BaseObject* )
    { }

    /// @}


    /// Test if the given context is an ancestor of this context.
    /// An ancestor is a parent or (recursively) the parent of an ancestor.
    ///
    /// This method is an alias to BaseNode::hasAncestor, so that dynamic
    /// casts are not required to test relationships between contexts.
    virtual bool hasAncestor(const BaseContext* /*context*/) const
    {
        return false;
    }

    /// @name Adding/Removing objects. Note that these methods can fail if the context doesn't support attached objects
    /// @{

    /// Add an object, or return false if not supported
    virtual bool addObject( boost::intrusive_ptr<BaseObject> /*obj*/ )
    {
        return false;
    }

    /// Remove an object, or return false if not supported
    virtual bool removeObject( boost::intrusive_ptr<BaseObject> /*obj*/ )
    {
        return false;
    }

    /// @}

    /// @name Visitors.
    /// @{

    /// apply an action
    virtual void executeVisitor( simulation::Visitor* );

    /// Propagate an event
    virtual void propagateEvent( const core::ExecParams* params /* PARAMS FIRST  = core::ExecParams::defaultInstance()*/, Event* );

    /// @}


    /// @name Notifications for graph change listeners
    /// @{

    virtual void notifyAddSlave(core::objectmodel::BaseObject* master, core::objectmodel::BaseObject* slave);
    virtual void notifyRemoveSlave(core::objectmodel::BaseObject* master, core::objectmodel::BaseObject* slave);
    virtual void notifyMoveSlave(core::objectmodel::BaseObject* previousMaster, core::objectmodel::BaseObject* master, core::objectmodel::BaseObject* slave);

    /// @}

    friend std::ostream SOFA_CORE_API & operator << (std::ostream& out, const BaseContext& c );
};

} // namespace objectmodel

} // namespace core

} // namespace sofa

#endif


