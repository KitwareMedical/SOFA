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
//
// C++ Interface: Node
//
// Description:
//
//
// Author: The SOFA team </www.sofa-framework.org>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SOFA_SIMULATION_COMMON_NODE_H
#define SOFA_SIMULATION_COMMON_NODE_H

#include <sofa/core/ExecParams.h>
#include <sofa/core/objectmodel/Context.h>
// moved from GNode (27/04/08)
#include <sofa/core/objectmodel/BaseNode.h>
#include <sofa/core/objectmodel/ConfigurationSetting.h>
#include <sofa/core/BehaviorModel.h>
#include <sofa/core/objectmodel/ContextObject.h>
#include <sofa/core/CollisionModel.h>
#include <sofa/core/visual/VisualModel.h>
#include <sofa/core/visual/VisualManager.h>
#include <sofa/core/visual/Shader.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/Mapping.h>
#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/BaseInteractionForceField.h>
#include <sofa/core/behavior/Mass.h>
#include <sofa/core/behavior/BaseProjectiveConstraintSet.h>
#include <sofa/core/behavior/BaseConstraintSet.h>
#include <sofa/core/topology/Topology.h>
#include <sofa/core/topology/BaseTopologyObject.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/behavior/LinearSolver.h>
#include <sofa/core/behavior/OdeSolver.h>
#include <sofa/core/behavior/ConstraintSolver.h>
#include <sofa/core/behavior/BaseAnimationLoop.h>
#include <sofa/core/visual/VisualLoop.h>
#include <sofa/core/collision/Pipeline.h>
#include <sofa/core/loader/BaseLoader.h>
#include <sofa/core/objectmodel/Event.h>

#include <sofa/simulation/common/common.h>
#include <sofa/simulation/common/MutationListener.h>
#include <sofa/simulation/common/VisitorScheduler.h>
#include <sofa/simulation/common/xml/Element.h>

namespace sofa
{
namespace simulation
{
class Visitor;
}
}
using sofa::simulation::Visitor;
using sofa::simulation::VisitorScheduler;

#include <sofa/helper/system/thread/CTime.h>
#include <string>
#include <stack>

namespace sofa
{

namespace core
{
namespace visual
{
class VisualParams;
} // namespace visual
} // namespace core

namespace simulation
{

/**
   Implements the object (component) management of the core::Context.
   Contains objects in lists and provides accessors.
   The other nodes are not visible (unknown scene graph).

   @author The SOFA team </www.sofa-framework.org>
 */
class SOFA_SIMULATION_COMMON_API Node : public core::objectmodel::BaseNode, public sofa::core::objectmodel::Context
{

public:
    SOFA_ABSTRACT_CLASS2(Node, BaseNode, Context);

    typedef sofa::core::visual::DisplayFlags DisplayFlags;
protected:
    Node(const std::string& name="");

    virtual ~Node();
public:
    /// Create, add, then return the new child of this Node
    virtual Node::SPtr createChild(const std::string& nodeName)=0;

    /// @name High-level interface
    /// @{

    /// Parse the given description to assign values to this object's fields and potentially other parameters
    virtual void parse ( sofa::core::objectmodel::BaseObjectDescription* arg );

    /// Initialize the components
    void init(const core::ExecParams* params);
    bool isInitialized() {return initialized;}
    /// Apply modifications to the components
    void reinit(const core::ExecParams* params);
    /// Do one step forward in time
    void animate(const core::ExecParams* params /* PARAMS FIRST */, double dt);
    /// Draw the objects in an OpenGl context
    void glDraw(core::visual::VisualParams* params);
    /// @}

    /// @name Visitor handling
    /// @{

    /// Execute a recursive action starting from this node.
    /// This method bypasses the actionScheduler of this node if any.
    virtual void doExecuteVisitor(Visitor* action)=0;

    /// Execute a recursive action starting from this node
    void executeVisitor( simulation::Visitor* action);

    /// Execute a recursive action starting from this node
    void execute(simulation::Visitor& action)
    {
        simulation::Visitor* p = &action;
        executeVisitor(p);
    }

    /// Execute a recursive action starting from this node
    void execute(simulation::Visitor* p)
    {
        executeVisitor(p);
    }

    /// Execute a recursive action starting from this node
    template<class Act, class Params>
    void execute(const Params* params)
    {
        Act action(params);
        simulation::Visitor* p = &action;
        executeVisitor(p);
    }

    /// Execute a recursive action starting from this node
    template<class Act>
    void execute(core::visual::VisualParams* vparams)
    {
        Act action(vparams);
        simulation::Visitor* p = &action;
        executeVisitor(p);
    }
    /// @}

    /// @name Component containers
    /// @{
    // methods moved from GNode (27/04/08)

    /// Sequence class to hold a list of objects. Public access is only readonly using an interface similar to std::vector (size/[]/begin/end).
    /// UPDATE: it is now an alias for the Link pointer container
    template < class T, bool strong = false >
    class Sequence : public MultiLink<Node, T, BaseLink::FLAG_DOUBLELINK|(strong ? BaseLink::FLAG_STRONGLINK : BaseLink::FLAG_DUPLICATE)>
    {
    public:
        typedef MultiLink<Node, T, BaseLink::FLAG_DOUBLELINK|(strong ? BaseLink::FLAG_STRONGLINK : BaseLink::FLAG_DUPLICATE)> Inherit;
        typedef T pointed_type;
        typedef typename Inherit::DestPtr value_type;
        //typedef TPtr value_type;
        typedef typename Inherit::const_iterator const_iterator;
        typedef typename Inherit::const_reverse_iterator const_reverse_iterator;
        typedef const_iterator iterator;
        typedef const_reverse_iterator reverse_iterator;

        Sequence(const BaseLink::InitLink<Node>& init)
            : Inherit(init)
        {
        }

        value_type operator[](unsigned int i) const
        {
            return this->get(i);
        }

        /// Swap two values in the list. Uses a const_cast to violate the read-only iterators.
        void swap( iterator a, iterator b )
        {
            value_type& wa = const_cast<value_type&>(*a);
            value_type& wb = const_cast<value_type&>(*b);
            value_type tmp = *a;
            wa = *b;
            wb = tmp;
        }
    };

    /// Class to hold 0-or-1 object. Public access is only readonly using an interface similar to std::vector (size/[]/begin/end), plus an automatic convertion to one pointer.
    /// UPDATE: it is now an alias for the Link pointer container
    template < class T, bool duplicate = true >
    class Single : public SingleLink<Node, T, BaseLink::FLAG_DOUBLELINK|(duplicate ? BaseLink::FLAG_DUPLICATE : BaseLink::FLAG_NONE)>
    {
    public:
        typedef SingleLink<Node, T, BaseLink::FLAG_DOUBLELINK|(duplicate ? BaseLink::FLAG_DUPLICATE : BaseLink::FLAG_NONE)> Inherit;
        typedef T pointed_type;
        typedef typename Inherit::DestPtr value_type;
        //typedef TPtr value_type;
        typedef typename Inherit::const_iterator const_iterator;
        typedef typename Inherit::const_reverse_iterator const_reverse_iterator;
        typedef const_iterator iterator;
        typedef const_reverse_iterator reverse_iterator;

        Single(const BaseLink::InitLink<Node>& init)
            : Inherit(init)
        {
        }

        T* operator->() const
        {
            return this->get();
        }
        T& operator*() const
        {
            return *this->get();
        }
        operator T*() const
        {
            return this->get();
        }
    };

    Sequence<Node,true> child;
    typedef Sequence<Node,true>::iterator ChildIterator;

    Sequence<core::objectmodel::BaseObject,true> object;
    typedef Sequence<core::objectmodel::BaseObject,true>::iterator ObjectIterator;
    typedef Sequence<core::objectmodel::BaseObject,true>::reverse_iterator ObjectReverseIterator;

    Single<core::behavior::BaseAnimationLoop> animationManager;
    Single<core::visual::VisualLoop> visualLoop;

    Sequence<core::BehaviorModel> behaviorModel;
    Sequence<core::BaseMapping> mapping;

    Sequence<core::behavior::OdeSolver> solver;
    Sequence<core::behavior::ConstraintSolver> constraintSolver;
    Sequence<core::behavior::LinearSolver> linearSolver;

    Single<core::topology::Topology> topology;
    Single<core::topology::BaseMeshTopology> meshTopology;
    Sequence<core::topology::BaseTopologyObject> topologyObject;

    Single<core::BaseState> state;
    Single<core::behavior::BaseMechanicalState> mechanicalState;
    Single<core::BaseMapping> mechanicalMapping;
    Single<core::behavior::BaseMass> mass;
    Sequence<core::behavior::BaseForceField> forceField;
    Sequence<core::behavior::BaseInteractionForceField> interactionForceField;
    Sequence<core::behavior::BaseProjectiveConstraintSet> projectiveConstraintSet;
    Sequence<core::behavior::BaseConstraintSet> constraintSet;
    Sequence<core::objectmodel::ContextObject> contextObject;
    Sequence<core::objectmodel::ConfigurationSetting> configurationSetting;

    Sequence<core::visual::Shader> shaders;
    Sequence<core::visual::VisualModel> visualModel;
    Sequence<core::visual::VisualManager> visualManager;

    Sequence<core::CollisionModel> collisionModel;
    Single<core::collision::Pipeline> collisionPipeline;

    Sequence<core::objectmodel::BaseObject> unsorted;

    /// @}


    /// @name Set/get objects
    /// @{

    /// Add an object and return this. Detect the implemented interfaces and add the object to the corresponding lists.
    virtual bool addObject(core::objectmodel::BaseObject::SPtr obj);

    /// Remove an object
    virtual bool removeObject(core::objectmodel::BaseObject::SPtr obj);

    /// Move an object from another node
    virtual void moveObject(core::objectmodel::BaseObject::SPtr obj);

    /// Find an object given its name
    core::objectmodel::BaseObject* getObject(const std::string& name) const;

    void* findLinkDestClass(const core::objectmodel::BaseClass* destType, const std::string& path, const core::objectmodel::BaseLink* link);

#ifdef SOFA_SMP
    /// Get first partition
    Iterative::IterativePartition* getFirstPartition();
#endif

    /// Generic object access, given a set of required tags, possibly searching up or down from the current context
    ///
    /// Note that the template wrapper method should generally be used to have the correct return type,
    virtual void* getObject(const sofa::core::objectmodel::ClassInfo& class_info, const sofa::core::objectmodel::TagSet& tags, SearchDirection dir = SearchUp) const=0;

    /// Generic object access, possibly searching up or down from the current context
    ///
    /// Note that the template wrapper method should generally be used to have the correct return type,
    void* getObject(const sofa::core::objectmodel::ClassInfo& class_info, SearchDirection dir = SearchUp) const
    {
        return getObject(class_info, sofa::core::objectmodel::TagSet(), dir);
    }

    /// Generic object access, given a path from the current context
    ///
    /// Note that the template wrapper method should generally be used to have the correct return type,
    virtual void* getObject(const sofa::core::objectmodel::ClassInfo& class_info, const std::string& path) const=0;

    /// Generic list of objects access, given a set of required tags, possibly searching up or down from the current context
    ///
    /// Note that the template wrapper method should generally be used to have the correct return type,
    virtual void getObjects(const sofa::core::objectmodel::ClassInfo& class_info, GetObjectsCallBack& container, const sofa::core::objectmodel::TagSet& tags, SearchDirection dir = SearchUp) const =0;

    /// Generic list of objects access, possibly searching up or down from the current context
    ///
    /// Note that the template wrapper method should generally be used to have the correct return type,
    void getObjects(const sofa::core::objectmodel::ClassInfo& class_info, GetObjectsCallBack& container, SearchDirection dir = SearchUp) const
    {
        getObjects(class_info, container, sofa::core::objectmodel::TagSet(), dir);
    }




    /// List all objects of this node deriving from a given class
    template<class Object, class Container>
    void getNodeObjects(Container* list)
    {
        this->get<Object, Container>(list, Local);
    }

    /// Return an object of this node deriving from a given class, or NULL if not found.
    /// Note that only the first object is returned.
    template<class Object>
    void getNodeObject(Object*& result)
    {
        result = this->get<Object>(Local);
    }

    template<class Object>
    Object* getNodeObject()
    {
        return this->get<Object>(Local);
    }

    /// List all objects of this node and sub-nodes deriving from a given class
    template<class Object, class Container>
    void getTreeObjects(Container* list)
    {
        this->get<Object, Container>(list, SearchDown);
    }

    /// Return an object of this node and sub-nodes deriving from a given class, or NULL if not found.
    /// Note that only the first object is returned.
    template<class Object>
    void getTreeObject(Object*& result)
    {
        result = this->get<Object>(SearchDown);
    }

    template<class Object>
    Object* getTreeObject()
    {
        return this->get<Object>(SearchDown);
    }

    /// Topology
    virtual core::topology::Topology* getTopology() const;

    /// Mesh Topology (unified interface for both static and dynamic topologies)
    virtual core::topology::BaseMeshTopology* getMeshTopology() const;

    /// Degrees-of-Freedom
    virtual core::BaseState* getState() const;

    /// Mechanical Degrees-of-Freedom
    virtual core::behavior::BaseMechanicalState* getMechanicalState() const;

    /// Shader
    virtual core::visual::Shader* getShader() const;
    virtual core::visual::Shader* getShader(const sofa::core::objectmodel::TagSet& t) const;

    /// @name Solvers and main algorithms
    /// @{

    virtual core::behavior::BaseAnimationLoop* getAnimationLoop() const;
    virtual core::behavior::OdeSolver* getOdeSolver() const;
    virtual core::collision::Pipeline* getCollisionPipeline() const;
    virtual core::visual::VisualLoop* getVisualLoop() const;

    /// @}

    /// Remove odesolvers and mastercontroler
    virtual void removeControllers();

    /// Find a child node given its name
    Node* getChild(const std::string& name) const;

    /// Get a descendant node given its name
    Node* getTreeNode(const std::string& name) const;

    /// Get children nodes
    virtual Children getChildren() const;

    BaseContext* getRootContext() const
    {
        return getRoot()->getContext();
    }

    Node* setDebug(bool);
    bool getDebug() const;
    // debug
    void printComponents();

    const BaseContext* getContext() const;
    BaseContext* getContext();

    /// Update the whole context values, based on parent and local ContextObjects
    virtual void updateContext();

    /// Update the simulation context values(gravity, time...), based on parent and local ContextObjects
    virtual void updateSimulationContext();

    /// Called during initialization to corectly propagate the visual context to the children
    virtual void initVisualContext() {}

    /// Propagate an event
    virtual void propagateEvent(const core::ExecParams* params /* PARAMS FIRST  = sofa::core::ExecParams::defaultInstance()*/, core::objectmodel::Event* event);

    /// Update the visual context values, based on parent and local ContextObjects
    virtual void updateVisualContext();

    Single<VisitorScheduler> actionScheduler;

    // VisitorScheduler can use doExecuteVisitor() method
    friend class VisitorScheduler;

    /// Must be called after each graph modification. Do not call it directly, apply an InitVisitor instead.
    virtual void initialize();

    virtual void bwdInit();

    /// Called after initialization to set the default value of the visual context.
    virtual void setDefaultVisualContextValue();

    template <class RealObject>
    static Node::SPtr create( RealObject* obj, sofa::simulation::xml::Element<sofa::core::objectmodel::BaseNode>*& arg);


    static Node::SPtr create( const std::string& name );

protected:
    bool debug_;
    bool initialized;

    virtual void doAddObject(core::objectmodel::BaseObject::SPtr obj);
    virtual void doRemoveObject(core::objectmodel::BaseObject::SPtr obj);


    std::stack<Visitor*> actionStack;

    virtual void notifyAddChild(Node::SPtr node);
    virtual void notifyRemoveChild(Node::SPtr node);
    virtual void notifyMoveChild(Node::SPtr node, Node* prev);
    virtual void notifyAddObject(core::objectmodel::BaseObject::SPtr obj);
    virtual void notifyRemoveObject(core::objectmodel::BaseObject::SPtr obj);
    virtual void notifyMoveObject(core::objectmodel::BaseObject::SPtr obj, Node* prev);


    BaseContext* _context;

    helper::vector<MutationListener*> listener;


public:

    virtual void notifyAddSlave(core::objectmodel::BaseObject* master, core::objectmodel::BaseObject* slave);
    virtual void notifyRemoveSlave(core::objectmodel::BaseObject* master, core::objectmodel::BaseObject* slave);
    virtual void notifyMoveSlave(core::objectmodel::BaseObject* previousMaster, core::objectmodel::BaseObject* master, core::objectmodel::BaseObject* slave);

    virtual void addListener(MutationListener* obj);
    virtual void removeListener(MutationListener* obj);

    // Added by FF to model component dependencies
    /// Pairs representing component dependencies. First must be initialized before second.
    Data < sofa::helper::vector < std::string > > depend;
    /// Sort the components according to the dependencies expressed in Data depend.
    void sortComponents();


};

}

}

#endif
