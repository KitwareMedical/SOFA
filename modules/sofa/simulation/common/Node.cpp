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
// C++ Implementation: Node
//
// Description:
//
//
// Author: The SOFA team </www.sofa-framework.org>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <sofa/simulation/common/Node.h>
#include <sofa/simulation/common/Simulation.h>
#include <sofa/simulation/common/xml/NodeElement.h>
#include <sofa/simulation/common/PropagateEventVisitor.h>
#include <sofa/simulation/common/UpdateMappingEndEvent.h>
#include <sofa/simulation/common/AnimateVisitor.h>
#include <sofa/simulation/common/DeactivatedNodeVisitor.h>
#include <sofa/simulation/common/InitVisitor.h>
#include <sofa/simulation/common/MechanicalVisitor.h>
#include <sofa/simulation/common/VisualVisitor.h>
#include <sofa/simulation/common/UpdateMappingVisitor.h>

#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/Factory.inl>
#include <sofa/simulation/common/xml/Element.inl>
#include <iostream>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>

//#define DEBUG_VISITOR
//#define DEBUG_LINK

namespace sofa
{

namespace simulation
{
using std::cerr;
using std::endl;
using core::objectmodel::BaseNode;
using core::objectmodel::BaseObject;

Node::Node(const std::string& name)
    : core::objectmodel::BaseNode()
    , sofa::core::objectmodel::Context()
    , child(initLink("child", "Child nodes"))
    , object(initLink("object","All objects attached to this node"))

    , animationManager(initLink("animationLoop","The AnimationLoop attached to this node (only valid for root node)"))
    , visualLoop(initLink("visualLoop", "The VisualLoop attached to this node (only valid for root node)"))

    , behaviorModel(initLink("behaviorModel", "The BehaviorModel attached to this node (only valid for root node)"))
    , mapping(initLink("mapping", "The (non-mechanical) Mapping(s) attached to this node (only valid for root node)"))

    , solver(initLink("odeSolver", "The OdeSolver(s) attached to this node (controlling the mechanical time integration of this branch)"))
    , constraintSolver(initLink("constraintSolver", "The ConstraintSolver(s) attached to this node"))
    , linearSolver(initLink("linearSolver", "The LinearSolver(s) attached to this node"))
    , topology(initLink("topology", "The Topology attached to this node"))
    , meshTopology(initLink("meshTopology", "The MeshTopology / TopologyContainer attached to this node"))
    , topologyObject(initLink("topologyObject", "The topology-related objects attached to this node"))
    , state(initLink("state", "The State attached to this node (storing vectors such as position, velocity)"))
    , mechanicalState(initLink("mechanicalState", "The MechanicalState attached to this node (storing all state vectors)"))
    , mechanicalMapping(initLink("mechanicalMapping", "The MechanicalMapping attached to this node"))
    , mass(initLink("mass", "The Mass attached to this node"))
    , forceField(initLink("forceField", "The (non-interaction) ForceField(s) attached to this node"))
    , interactionForceField(initLink("interactionForceField", "The InteractionForceField(s) attached to this node"))
    , projectiveConstraintSet(initLink("projectiveConstraintSet", "The ProjectiveConstraintSet(s) attached to this node"))
    , constraintSet(initLink("constraintSet", "The ConstraintSet(s) attached to this node"))
    , contextObject(initLink("contextObject", "The ContextObject(s) attached to this node"))
    , configurationSetting(initLink("configurationSetting", "The ConfigurationSetting(s) attached to this node"))

    , shaders(initLink("shaders", "The shaders attached to this node"))
    , visualModel(initLink("visualModel", "The VisualModel(s) attached to this node"))
    , visualManager(initLink("visualManager", "The VisualManager(s) attached to this node"))

    , collisionModel(initLink("collisionModel", "The CollisionModel(s) attached to this node"))
    , collisionPipeline(initLink("collisionPipeline", "The collision Pipeline attached to this node"))

    , unsorted(initLink("unsorted", "The remaining objects attached to this node"))

    , actionScheduler(initLink("visitorScheduler", "The VisitorScheduler attached to this node (deprecated)"))
    , debug_(false)
    , initialized(false)
    , depend(initData(&depend,"depend","Dependencies between the nodes.\nname 1 name 2 name3 name4 means that name1 must be initialized before name2 and name3 before name4"))
{
    _context = this;
    setName(name);
    //f_printLog.setValue(true);
}


Node::~Node()
{
}

void Node::parse( sofa::core::objectmodel::BaseObjectDescription* arg )
{
    Inherit1::parse( arg );
    static const char* oldVisualFlags[] =
    {
        "showAll",
        "showVisual",
        "showBehavior",
        "showCollision",
        "showMapping",
        "showVisualModels",
        "showBehaviorModels",
        "showCollisionModels",
        "showBoundingCollisionModels",
        "showMappings",
        "showMechanicalMappings",
        "showForceFields",
        "showInteractionForceFields",
        "showWireFrame",
        "showNormals",
        NULL
    };
    std::string oldFlags;
    for (unsigned int i=0; oldVisualFlags[i]; ++i)
    {
        const char* str = arg->getAttribute(oldVisualFlags[i], NULL);
        if (str == NULL || !*str) continue;
        bool val;
        if (str[0] == 'T' || str[0] == 't')
            val = true;
        else if (str[0] == 'F' || str[0] == 'f')
            val = false;
        else if ((str[0] >= '0' && str[0] <= '9') || str[0] == '-')
            val = (atoi(str) != 0);
        else continue;
        if (!oldFlags.empty()) oldFlags += ' ';
        if (val) oldFlags += oldVisualFlags[i];
        else { oldFlags += "hide"; oldFlags += oldVisualFlags[i]+4; }
    }
    if (!oldFlags.empty())
    {
        serr << "Deprecated visual flags attributes used. Instead, add the following object within the Node:" << sendl;
        serr << "<VisualStyle displayFlags=\"" << oldFlags << "\" />" << sendl;

        sofa::core::objectmodel::BaseObjectDescription objDesc("displayFlags","VisualStyle");
        objDesc.setAttribute("displayFlags", oldFlags.c_str());
        sofa::core::objectmodel::BaseObject::SPtr obj = sofa::core::ObjectFactory::CreateObject(this, &objDesc);
    }
}

/// Initialize the components of this node and all the nodes which depend on it.
void Node::init(const core::ExecParams* params)
{
    //     cerr<<"Node::init() begin node "<<getName()<<endl;
    execute<simulation::InitVisitor>(params);

    //     cerr<<"Node::init() end node "<<getName()<<endl;
}

/// ReInitialize the components of this node and all the nodes which depend on it.
void Node::reinit(const core::ExecParams* params)
{
    sofa::simulation::DeactivationVisitor deactivate(params, isActive());
    deactivate.execute( this );
}

/// Do one step forward in time
void Node::animate(const core::ExecParams* params /* PARAMS FIRST */, double dt)
{
    simulation::AnimateVisitor vis(params /* PARAMS FIRST */, dt);
    //cerr<<"Node::animate, start execute"<<endl;
    execute(vis);
    //cerr<<"Node::animate, end execute"<<endl;
}

void Node::glDraw(core::visual::VisualParams* vparams)
{
    execute<simulation::VisualUpdateVisitor>(vparams);
    execute<simulation::VisualDrawVisitor>(vparams);
}


/// Add an object. Detect the implemented interfaces and add the object to the corresponding lists.
bool Node::addObject(BaseObject::SPtr obj)
{
    notifyAddObject(obj);
    doAddObject(obj);
    return true;
}

/// Remove an object
bool Node::removeObject(BaseObject::SPtr obj)
{
    notifyRemoveObject(obj);
    doRemoveObject(obj);
    return true;
}

/// Move an object from another node
void Node::moveObject(BaseObject::SPtr obj)
{
    Node* prev = dynamic_cast<Node*>(obj->getContext());
    if (prev==NULL)
    {
        obj->getContext()->removeObject(obj);
        addObject(obj);
    }
    else
    {
        notifyMoveObject(obj,prev);
        prev->doRemoveObject(obj);
        doAddObject(obj);
    }
}



void Node::notifyAddChild(Node::SPtr node)
{
    for (helper::vector<MutationListener*>::const_iterator it = listener.begin(); it != listener.end(); ++it)
        (*it)->addChild(this, node.get());
}


void Node::notifyRemoveChild(Node::SPtr node)
{
    for (helper::vector<MutationListener*>::const_iterator it = listener.begin(); it != listener.end(); ++it)
        (*it)->removeChild(this, node.get());
}


void Node::notifyMoveChild(Node::SPtr node, Node* prev)
{
    for (helper::vector<MutationListener*>::const_iterator it = listener.begin(); it != listener.end(); ++it)
        (*it)->moveChild(prev, this, node.get());
}


void Node::notifyAddObject(core::objectmodel::BaseObject::SPtr obj)
{
    for (helper::vector<MutationListener*>::const_iterator it = listener.begin(); it != listener.end(); ++it)
        (*it)->addObject(this, obj.get());
}

void Node::notifyRemoveObject(core::objectmodel::BaseObject::SPtr obj)
{
    for (helper::vector<MutationListener*>::const_iterator it = listener.begin(); it != listener.end(); ++it)
        (*it)->removeObject(this, obj.get());
}

void Node::notifyMoveObject(core::objectmodel::BaseObject::SPtr obj, Node* prev)
{
    for (helper::vector<MutationListener*>::const_iterator it = listener.begin(); it != listener.end(); ++it)
        (*it)->moveObject(prev, this, obj.get());
}


void Node::notifyAddSlave(core::objectmodel::BaseObject* master, core::objectmodel::BaseObject* slave)
{
    for (helper::vector<MutationListener*>::const_iterator it = listener.begin(); it != listener.end(); ++it)
        (*it)->addSlave(master, slave);
}

void Node::notifyRemoveSlave(core::objectmodel::BaseObject* master, core::objectmodel::BaseObject* slave)
{
    for (helper::vector<MutationListener*>::const_iterator it = listener.begin(); it != listener.end(); ++it)
        (*it)->removeSlave(master, slave);
}

void Node::notifyMoveSlave(core::objectmodel::BaseObject* previousMaster, core::objectmodel::BaseObject* master, core::objectmodel::BaseObject* slave)
{
    for (helper::vector<MutationListener*>::const_iterator it = listener.begin(); it != listener.end(); ++it)
        (*it)->moveSlave(previousMaster, master, slave);
}

void Node::addListener(MutationListener* obj)
{
    // make sure we don't add the same listener twice
    helper::vector< MutationListener* >::iterator it = listener.begin();
    while (it != listener.end() && (*it)!=obj)
        ++it;
    if (it == listener.end())
        listener.push_back(obj);
}

void Node::removeListener(MutationListener* obj)
{
    helper::vector< MutationListener* >::iterator it = listener.begin();
    while (it != listener.end() && (*it)!=obj)
        ++it;
    if (it != listener.end())
        listener.erase(it);
}


/// Find an object given its name
core::objectmodel::BaseObject* Node::getObject(const std::string& name) const
{
    for (ObjectIterator it = object.begin(), itend = object.end(); it != itend; ++it)
        if ((*it)->getName() == name)
            return it->get();
    return NULL;
}

void* Node::findLinkDestClass(const core::objectmodel::BaseClass* destType, const std::string& path, const core::objectmodel::BaseLink* link)
{
    std::string pathStr;
    if (link)
    {
        if (!link->parseString(path,&pathStr))
            return NULL;
    }
    else
    {
        if (!BaseLink::ParseString(path,&pathStr,NULL,this))
            return NULL;
    }
#ifdef DEBUG_LINK
    std::cout << "LINK: Looking for " << destType->className << "<" << destType->templateName << "> " << pathStr << " from Node " << getName() << std::endl;
#endif
    std::size_t ppos = 0;
    std::size_t psize = pathStr.size();
    if (ppos == psize || (ppos == psize-2 && pathStr[ppos] == '[' && pathStr[ppos+1] == ']')) // self-reference
    {
#ifdef DEBUG_LINK
        std::cout << "  self-reference link." << std::endl;
#endif
        if (!link || !link->getOwnerBase()) return destType->dynamicCast(this);
        return destType->dynamicCast(link->getOwnerBase());
    }
    Node* node = this;
    BaseObject* master = NULL;
    bool based = false;
    if (ppos < psize && pathStr[ppos] == '[') // relative index in the list of objects
    {
        if (pathStr[psize-1] != ']')
        {
            serr << "Invalid index-based path \"" << path << "\"" << sendl;
            return NULL;
        }
        int index = atoi(pathStr.c_str()+ppos+1);
#ifdef DEBUG_LINK
        std::cout << "  index-based path to " << index << std::endl;
#endif
        ObjectReverseIterator it = object.rbegin();
        ObjectReverseIterator itend = object.rend();
        if (link && link->getOwnerBase())
        {
            // index from last
            Base* b = link->getOwnerBase();
            while (it != itend && *it != b)
                ++it;
        }
        while (it != itend && index < 0)
        {
            ++it;
            ++index;
        }
        if (it == itend)
            return NULL;
#ifdef DEBUG_LINK
        std::cout << "  found " << it->get()->getTypeName() << " " << it->get()->getName() << "." << std::endl;
#endif
        return destType->dynamicCast(it->get());
    }
    else if (ppos < psize && pathStr[ppos] == '/') // absolute path
    {
#ifdef DEBUG_LINK
        std::cout << "  absolute path" << std::endl;
#endif
        node = dynamic_cast<Node*>(this->getRoot());
        if (!node) return NULL;
        ++ppos;
        based = true;
    }
    while(ppos < psize)
    {
        if ((ppos+1 < psize && pathStr.substr(ppos,2) == "./")
            || pathStr.substr(ppos) == ".")
        {
            // this must be this node
#ifdef DEBUG_LINK
            std::cout << "  to current node" << std::endl;
#endif
            ppos += 2;
            based = true;
        }
        else if ((ppos+2 < psize && pathStr.substr(ppos,3) == "../") // relative
                || pathStr.substr(ppos) == "..")
        {
            ppos += 3;
            if (master)
            {
                master = master->getMaster();
#ifdef DEBUG_LINK
                std::cout << "  to master object " << master->getName() << std::endl;
#endif
            }
            else
            {
                Parents parents = node->getParents();
                if (parents.empty()) return NULL;
                node = dynamic_cast<Node*>(parents[0]); // TODO: explore other parents
                if (!node) return NULL;
#ifdef DEBUG_LINK
                std::cout << "  to parent node " << node->getName() << std::endl;
#endif
            }
            based = true;
        }
        else if (pathStr[ppos] == '/')
        {
            // extra /
#ifdef DEBUG_LINK
            std::cout << "  extra '/'" << std::endl;
#endif
            ppos += 1;
        }
        else
        {
            std::size_t p2pos = pathStr.find('/',ppos);
            if (p2pos == std::string::npos) p2pos = psize;
            std::string name = pathStr.substr(ppos,p2pos-ppos);
            ppos = p2pos+1;
            if (master)
            {
#ifdef DEBUG_LINK
                std::cout << "  to slave object " << name << std::endl;
#endif
                master = master->getSlave(name);
                if (!master) return NULL;
            }
            else
            {
                for (;;)
                {
                    BaseObject* obj = node->getObject(name);
                    Node* child = node->getChild(name);
                    if (child)
                    {
                        node = child;
#ifdef DEBUG_LINK
                        std::cout << "  to child node " << name << std::endl;
#endif
                        break;
                    }
                    else if (obj)
                    {
                        master = obj;
#ifdef DEBUG_LINK
                        std::cout << "  to object " << name << std::endl;
#endif
                        break;
                    }
                    if (based) return NULL;
                    // this can still be found from an ancestor node
                    Parents parents = node->getParents();
                    if (parents.empty()) return NULL;
                    node = dynamic_cast<Node*>(parents[0]); // TODO: explore other parents
                    if (!node) return NULL;
#ifdef DEBUG_LINK
                    std::cout << "  looking in ancestor node " << node->getName() << std::endl;
#endif
                }
            }
            based = true;
        }
    }
    if (master)
    {
#ifdef DEBUG_LINK
        std::cout << "  found " << master->getTypeName() << " " << master->getName() << "." << std::endl;
#endif
        return destType->dynamicCast(master);
    }
    else
    {
        void* r = destType->dynamicCast(node);
        if (r)
        {
#ifdef DEBUG_LINK
            std::cout << "  found node " << node->getName() << "." << std::endl;
#endif
            return r;
        }
        for (ObjectIterator it = node->object.begin(), itend = node->object.end(); it != itend; ++it)
        {
            BaseObject* obj = it->get();
            void *o = destType->dynamicCast(obj);
            if (!o) continue;
#ifdef DEBUG_LINK
            std::cout << "  found " << obj->getTypeName() << " " << obj->getName() << "." << std::endl;
#endif
            if (!r) r = o;
            else return NULL; // several objects are possible, this is an ambiguous path
        }
        if (r) return r;
        // no object found, we look in parent nodes if the searched class is one of the known standard single components (state, topology, ...)
        if (destType->hasParent(sofa::core::BaseState::GetClass()))
            return destType->dynamicCast(node->getState());
        else if (destType->hasParent(core::topology::BaseMeshTopology::GetClass()))
            return destType->dynamicCast(node->getMeshTopology());
        else if (destType->hasParent(core::topology::Topology::GetClass()))
            return destType->dynamicCast(node->getTopology());
        else if (destType->hasParent(core::visual::Shader::GetClass()))
            return destType->dynamicCast(node->getShader());
        else if (destType->hasParent(core::behavior::BaseAnimationLoop::GetClass()))
            return destType->dynamicCast(node->getAnimationLoop());
        else if (destType->hasParent(core::behavior::OdeSolver::GetClass()))
            return destType->dynamicCast(node->getOdeSolver());
        else if (destType->hasParent(core::collision::Pipeline::GetClass()))
            return destType->dynamicCast(node->getCollisionPipeline());
        else if (destType->hasParent(core::visual::VisualLoop::GetClass()))
            return destType->dynamicCast(node->getVisualLoop());

        return NULL;
    }
}

/// Add an object. Detect the implemented interfaces and add the object to the corresponding lists.
void Node::doAddObject(BaseObject::SPtr sobj)
{
    notifyAddObject(sobj);
    //sobj->setContext(this);
    this->setObjectContext(sobj);
    object.add(sobj);
    BaseObject* obj = sobj.get();
    int inserted=0;
    inserted+= animationManager.add(dynamic_cast< core::behavior::BaseAnimationLoop* >(obj));
    inserted+= solver.add(dynamic_cast< core::behavior::OdeSolver* >(obj));
    inserted+= linearSolver.add(dynamic_cast< core::behavior::BaseLinearSolver* >(obj));
    inserted+= constraintSolver.add(dynamic_cast< core::behavior::ConstraintSolver* >(obj));
    inserted+= visualLoop.add(dynamic_cast< core::visual::VisualLoop* >(obj));
    inserted+= state.add(dynamic_cast< core::BaseState* >(obj));
    inserted+= mechanicalState.add(dynamic_cast< core::behavior::BaseMechanicalState* >(obj));
    core::BaseMapping* bmap = dynamic_cast< core::BaseMapping* >(obj);
    bool isMechanicalMapping = false;
    if(bmap)
    {
        isMechanicalMapping = bmap->isMechanical();
        if(isMechanicalMapping)
            inserted += mechanicalMapping.add(bmap);
        else
            inserted += mapping.add(bmap);
    }

    inserted+= mass.add(dynamic_cast< core::behavior::BaseMass* >(obj));
    inserted+= topology.add(dynamic_cast< core::topology::Topology* >(obj));
    inserted+= meshTopology.add(dynamic_cast< core::topology::BaseMeshTopology* >(obj));
    inserted+= topologyObject.add(dynamic_cast< core::topology::BaseTopologyObject* >(obj));
    inserted+= shaders.add(dynamic_cast< sofa::core::visual::Shader* >(obj));

    bool isInteractionForceField = interactionForceField.add(dynamic_cast< core::behavior::BaseInteractionForceField* >(obj));
    inserted+= isInteractionForceField;
    if (!isInteractionForceField)
        forceField.add(dynamic_cast< core::behavior::BaseForceField* >(obj));
    inserted+= projectiveConstraintSet.add(dynamic_cast< core::behavior::BaseProjectiveConstraintSet* >(obj));
    inserted+= constraintSet.add(dynamic_cast< core::behavior::BaseConstraintSet* >(obj));
    inserted+= behaviorModel.add(dynamic_cast< core::BehaviorModel* >(obj));
    inserted+= visualModel.add(dynamic_cast< core::visual::VisualModel* >(obj));
    inserted+= visualManager.add(dynamic_cast< core::visual::VisualManager* >(obj));
    inserted+= collisionModel.add(dynamic_cast< core::CollisionModel* >(obj));
    inserted+= contextObject.add(dynamic_cast< core::objectmodel::ContextObject* >(obj));
    inserted+= configurationSetting.add(dynamic_cast< core::objectmodel::ConfigurationSetting* >(obj));
    inserted+= collisionPipeline.add(dynamic_cast< core::collision::Pipeline* >(obj));
    inserted+= actionScheduler.add(dynamic_cast< VisitorScheduler* >(obj));

    if ( inserted==0 )
    {
        //cerr<<"Node::doAddObject, object "<<obj->getName()<<" is unsorted"<<endl;
        unsorted.add(obj);
    }

}

/// Remove an object
void Node::doRemoveObject(BaseObject::SPtr sobj)
{
    //if (sobj->getContext()==this)
    //{
    //    sobj->setContext(NULL);
    //}
    this->clearObjectContext(sobj);
    object.remove(sobj);
    BaseObject* obj = sobj.get();
    animationManager.remove(dynamic_cast< core::behavior::BaseAnimationLoop* >(obj));
    solver.remove(dynamic_cast< core::behavior::OdeSolver* >(obj));
    linearSolver.remove(dynamic_cast< core::behavior::BaseLinearSolver* >(obj));
    constraintSolver.remove(dynamic_cast< core::behavior::ConstraintSolver* >(obj));
    visualLoop.remove(dynamic_cast< core::visual::VisualLoop* >(obj));
    state.remove(dynamic_cast< core::BaseState* >(obj));
    mechanicalState.remove(dynamic_cast< core::behavior::BaseMechanicalState* >(obj));
    mechanicalMapping.remove(dynamic_cast< core::BaseMapping* >(obj));
    mass.remove(dynamic_cast< core::behavior::BaseMass* >(obj));
    topology.remove(dynamic_cast< core::topology::Topology* >(obj));
    meshTopology.remove(dynamic_cast< core::topology::BaseMeshTopology* >(obj));
    topologyObject.remove(dynamic_cast< core::topology::BaseTopologyObject* >(obj));
    shaders.remove(dynamic_cast<sofa::core::visual::Shader* >(obj));

    forceField.remove(dynamic_cast< core::behavior::BaseForceField* >(obj));
    interactionForceField.remove(dynamic_cast< core::behavior::BaseInteractionForceField* >(obj));
    projectiveConstraintSet.remove(dynamic_cast< core::behavior::BaseProjectiveConstraintSet* >(obj));
    constraintSet.remove(dynamic_cast< core::behavior::BaseConstraintSet* >(obj));
    mapping.remove(dynamic_cast< core::BaseMapping* >(obj));
    behaviorModel.remove(dynamic_cast< core::BehaviorModel* >(obj));
    visualModel.remove(dynamic_cast< core::visual::VisualModel* >(obj));
    visualManager.remove(dynamic_cast< core::visual::VisualManager* >(obj));
    collisionModel.remove(dynamic_cast< core::CollisionModel* >(obj));
    contextObject.remove(dynamic_cast<core::objectmodel::ContextObject* >(obj));
    configurationSetting.remove(dynamic_cast<core::objectmodel::ConfigurationSetting* >(obj));
    collisionPipeline.remove(dynamic_cast< core::collision::Pipeline* >(obj));

    actionScheduler.remove(dynamic_cast< VisitorScheduler* >(obj));

    unsorted.remove(obj);
}


/// Topology
core::topology::Topology* Node::getTopology() const
{
    // return this->topology;
    if (this->topology)
        return this->topology;
    else
        return get<core::topology::Topology>(SearchParents);
}

/// Mesh Topology (unified interface for both static and dynamic topologies)
core::topology::BaseMeshTopology* Node::getMeshTopology() const
{
    if (this->meshTopology)
        return this->meshTopology;
    else
        return get<core::topology::BaseMeshTopology>(SearchParents);
}

/// Degrees-of-Freedom
core::BaseState* Node::getState() const
{
    // return this->state;
    if (this->state)
        return this->state;
    else
        return get<core::BaseState>(SearchParents);
}

/// Mechanical Degrees-of-Freedom
core::behavior::BaseMechanicalState* Node::getMechanicalState() const
{
    // return this->mechanicalModel;
    if (this->mechanicalState)
        return this->mechanicalState;
    else
        return get<core::behavior::BaseMechanicalState>(SearchParents);
}

/// Shader
core::visual::Shader* Node::getShader() const
{
    if (!shaders.empty())
        return *shaders.begin();
    else
        return get<core::visual::Shader>(SearchParents);
}
core::visual::Shader* Node::getShader(const sofa::core::objectmodel::TagSet& t) const
{
    if(t.empty())
        return getShader();
    else // if getShader is Tag filtered
    {
        for(Sequence<core::visual::Shader>::iterator it = shaders.begin(), iend=shaders.end(); it!=iend; it++)
        {
            if ( (*it)->getTags().includes(t) )
                return (*it);
        }
        return get<core::visual::Shader>(t,SearchParents);
    }
}

core::behavior::BaseAnimationLoop* Node::getAnimationLoop() const
{
    if (animationManager)
        return animationManager;
    else
        return get<core::behavior::BaseAnimationLoop>(SearchParents);
}

core::behavior::OdeSolver* Node::getOdeSolver() const
{
    if (!solver.empty())
        return solver[0];
    else
        return get<core::behavior::OdeSolver>(SearchParents);
}

core::collision::Pipeline* Node::getCollisionPipeline() const
{
    if (collisionPipeline)
        return collisionPipeline;
    else
        return get<core::collision::Pipeline>(SearchParents);
}

core::visual::VisualLoop* Node::getVisualLoop() const
{
    if (visualLoop)
        return visualLoop;
    else
        return get<core::visual::VisualLoop>(SearchParents);
}

/// Find a child node given its name
Node* Node::getChild(const std::string& name) const
{
//    cerr<<"Node::getChild, in "<< getName() << ", looking for " << name << endl;
    for (ChildIterator it = child.begin(), itend = child.end(); it != itend; ++it)
    {
//        cerr<<"Node::getChild, see " << (*it)->getName() << endl;
        if ((*it)->getName() == name)
            return it->get();
    }
    return NULL;
}

/// Get a descendant node given its name
Node* Node::getTreeNode(const std::string& name) const
{
    Node* result = NULL;
    result = getChild(name);
    for (ChildIterator it = child.begin(), itend = child.end(); result == NULL && it != itend; ++it)
        result = (*it)->getTreeNode(name);
    return result;
}

/// Get parent node (or NULL if no hierarchy or for root node)
sofa::core::objectmodel::BaseNode::Children Node::getChildren() const
{
    Children list_children;
    list_children.reserve(child.size());
    for (ChildIterator it = child.begin(), itend = child.end(); it != itend; ++it)
        list_children.push_back(it->get());
    return list_children;
}

Node* Node::setDebug(bool b)
{
    debug_=b;
    return this;
}

bool Node::getDebug() const
{
    return debug_;
}




void Node::removeControllers()
{
    removeObject(*animationManager.begin());
    typedef Sequence<core::behavior::OdeSolver> Solvers;
    Solvers solverRemove = solver;
    for ( Solvers::iterator i=solverRemove.begin(), iend=solverRemove.end(); i!=iend; i++ )
        removeObject( *i );
}


core::objectmodel::BaseContext* Node::getContext()
{
    return _context;
}
const core::objectmodel::BaseContext* Node::getContext() const
{
    return _context;
}

// void Node::setContext( core::objectmodel::BaseContext* c )
// {
//     _context=c;
// 	for( ObjectIterator i=object.begin(), iend=object.end(); i!=iend; i++ )
// 		(*i)->setContext(c);
// }


void Node::setDefaultVisualContextValue()
{
    /// @TODO: This method is now broken because getShow*() methods never return -1
    /*
    if (getShowVisualModels() == -1)            setShowVisualModels(true);
    if (getShowBehaviorModels() == -1)          setShowBehaviorModels(false);
    if (getShowCollisionModels() == -1)         setShowCollisionModels(false);
    if (getShowBoundingCollisionModels() == -1) setShowBoundingCollisionModels(false);
    if (getShowMappings() == -1)                setShowMappings(false);
    if (getShowMechanicalMappings() == -1)      setShowMechanicalMappings(false);
    if (getShowForceFields() == -1)             setShowForceFields(false);
    if (getShowInteractionForceFields() == -1)  setShowInteractionForceFields(false);
    if (getShowWireFrame() == -1)               setShowWireFrame(false);
    if (getShowNormals() == -1)                 setShowNormals(false);
    #ifdef SOFA_SMP
    if (showProcessorColor_.getValue() == -1)                 showProcessorColor_.setValue(false);
    #endif
    */
}

void Node::bwdInit()
{
    if (mechanicalMapping && !mechanicalMapping->isMechanical())
    {
        // BUGFIX: the mapping was configured as not mechanical -> remove it from mechanicalMapping and put it in mapping
        core::BaseMapping* bmap = mechanicalMapping.get();
        mapping.add(bmap);
        mechanicalMapping.remove(bmap);
    }
    //printComponents();
}

void Node::initialize()
{
    initialized = true;  // flag telling is the node is initialized
    //cerr<<"Node::initialize()"<<endl;

    initVisualContext();
    sortComponents();
    //     // Put the OdeSolver, if any, in first position. This makes sure that the OdeSolver component is initialized only when all its sibling and children components are already initialized.
    //     /// @todo Putting the solver first means that it will be initialized *before* any sibling or childrens. Is that what we want? -- Jeremie A.
    //     Sequence<BaseObject>::iterator i=object.begin(), iend=object.end();
    //     for ( ; i!=iend && dynamic_cast<core::behavior::OdeSolver*>(*i)==NULL; i++ ) // find the OdeSolver
    //         {}
    //     if ( i!=iend && !object.empty() ) // found
    //     {
    //         // put it first
    //         // BUGFIX 01/12/06 (Jeremie A.): do not modify the order of the other objects
    //         // object.swap( i, object.begin() );
    //         while (i!=object.begin())
    //         {
    //             Sequence<BaseObject>::iterator i2 = i;
    //             --i;
    //             object.swap(i, i2);
    //         }
    //     }

    //
    updateSimulationContext();

    // this is now done by the InitVisitor
    //for (Sequence<Node>::iterator it = child.begin(); it != child.end(); it++) {
    //    (*it)->init();
    //}

}

void Node::updateContext()
{
    updateSimulationContext();
    updateVisualContext();
    if ( debug_ ) std::cerr<<"Node::updateContext, node = "<<getName()<<", updated context = "<< *static_cast<core::objectmodel::Context*>(this) << endl;
}

void Node::updateSimulationContext()
{
    for ( unsigned i=0; i<contextObject.size(); ++i )
    {
        contextObject[i]->init();
        contextObject[i]->apply();
        //       cerr<<"Node::updateContext, modified by node = "<<contextObject[i]->getName()<<endl;
    }
}

void Node::updateVisualContext()
{
    // Apply local modifications to the context
    for ( unsigned i=0; i<contextObject.size(); ++i )
    {
        contextObject[i]->init();
        contextObject[i]->apply();
    }

    if ( debug_ ) std::cerr<<"Node::updateVisualContext, node = "<<getName()<<", updated context = "<< *static_cast<core::objectmodel::Context*>(this) << endl;
}

/// Execute a recursive action starting from this node
void Node::executeVisitor(Visitor* action)
{
    if (!this->isActive()) return;

    if (!action->execParams()->checkValidStorage())
    {
        std::cerr << "IN " << sofa::core::objectmodel::BaseClass::decodeClassName(typeid(*action)) << " at " << this->getPathName() << std::endl;
    }

#ifdef DEBUG_VISITOR
    static int level = 0;
    for (int i=0; i<level; ++i) std::cerr << ' ';
    std::cerr << ">" << sofa::core::objectmodel::BaseClass::decodeClassName(typeid(*action)) << " on " << this->getPathName();
    //     if (MechanicalVisitor* v = dynamic_cast<MechanicalVisitor*>(action))
    if (!action->getInfos().empty())
        std::cerr << "  : " << action->getInfos();
    std::cerr << std::endl;
    ++level;
#endif

    if (actionScheduler)
        actionScheduler->executeVisitor(this,action);
    else
        doExecuteVisitor(action);

#ifdef DEBUG_VISITOR
    --level;
    for (int i=0; i<level; ++i) std::cerr << ' ';
    std::cerr << "<" << sofa::core::objectmodel::BaseClass::decodeClassName(typeid(*action)) << " on " << this->getPathName();
    std::cerr << std::endl;
#endif

}

/// Propagate an event
void Node::propagateEvent(const core::ExecParams* params /* PARAMS FIRST */, core::objectmodel::Event* event)
{
    simulation::PropagateEventVisitor act(params /* PARAMS FIRST */, event);
    this->executeVisitor(&act);
}





void Node::printComponents()
{
    using namespace sofa::core::behavior;
    using core::BaseMapping;
    using core::topology::Topology;
    using core::topology::BaseTopologyObject;
    using core::topology::BaseMeshTopology;
    using core::visual::Shader;
    using core::BehaviorModel;
    using core::visual::VisualModel;
    using core::visual::VisualLoop;
    using core::CollisionModel;
    using core::objectmodel::ContextObject;
    using core::collision::Pipeline;

    serr<<"BaseAnimationLoop: ";
    for ( Single<BaseAnimationLoop>::iterator i=animationManager.begin(), iend=animationManager.end(); i!=iend; i++ )
        serr<<(*i)->getName()<<" ";
    serr<<sendl<<"OdeSolver: ";
    for ( Sequence<OdeSolver>::iterator i=solver.begin(), iend=solver.end(); i!=iend; i++ )
        serr<<(*i)->getName()<<" ";
    serr<<sendl<<"LinearSolver: ";
    for ( Sequence<BaseLinearSolver>::iterator i=linearSolver.begin(), iend=linearSolver.end(); i!=iend; i++ )
        serr<<(*i)->getName()<<" ";
    serr<<sendl<<"ConstraintSolver: ";
    for ( Sequence<ConstraintSolver>::iterator i=constraintSolver.begin(), iend=constraintSolver.end(); i!=iend; i++ )
        serr<<(*i)->getName()<<" ";
    serr<<"VisualLoop: ";
    for ( Single<VisualLoop>::iterator i=visualLoop.begin(), iend=visualLoop.end(); i!=iend; i++ )
        serr<<(*i)->getName()<<" ";
    serr<<sendl<<"InteractionForceField: ";
    for ( Sequence<BaseInteractionForceField>::iterator i=interactionForceField.begin(), iend=interactionForceField.end(); i!=iend; i++ )
        serr<<(*i)->getName()<<" ";
    serr<<sendl<<"ForceField: ";
    for ( Sequence<BaseForceField>::iterator i=forceField.begin(), iend=forceField.end(); i!=iend; i++ )
        serr<<(*i)->getName()<<" ";
    serr<<sendl<<"State: ";
    for ( Single<BaseState>::iterator i=state.begin(), iend=state.end(); i!=iend; i++ )
        serr<<(*i)->getName()<<" ";
    serr<<sendl<<"MechanicalState: ";
    for ( Single<BaseMechanicalState>::iterator i=mechanicalState.begin(), iend=mechanicalState.end(); i!=iend; i++ )
        serr<<(*i)->getName()<<" ";
    serr<<sendl<<"Mechanical Mapping: ";
    for ( Single<BaseMapping>::iterator i=mechanicalMapping.begin(), iend=mechanicalMapping.end(); i!=iend; i++ )
        serr<<(*i)->getName()<<" ";
    serr<<sendl<<"Mapping: ";
    for ( Sequence<BaseMapping>::iterator i=mapping.begin(), iend=mapping.end(); i!=iend; i++ )
        serr<<(*i)->getName()<<" ";
    serr<<sendl<<"Topology: ";
    for ( Single<Topology>::iterator i=topology.begin(), iend=topology.end(); i!=iend; i++ )
        serr<<(*i)->getName()<<" ";
    serr<<sendl<<"MeshTopology: ";
    for ( Single<BaseMeshTopology>::iterator i=meshTopology.begin(), iend=meshTopology.end(); i!=iend; i++ )
        serr<<(*i)->getName()<<" ";
    serr<<sendl<<"Shader: ";
    for ( Sequence<Shader>::iterator i=shaders.begin(), iend=shaders.end(); i!=iend; i++ )
        serr<<(*i)->getName()<<" ";
    serr<<sendl<<"ProjectiveConstraintSet: ";
    for ( Sequence<BaseProjectiveConstraintSet>::iterator i=projectiveConstraintSet.begin(), iend=projectiveConstraintSet.end(); i!=iend; i++ )
        serr<<(*i)->getName()<<" ";
    serr<<sendl<<"ConstraintSet: ";
    for ( Sequence<BaseConstraintSet>::iterator i=constraintSet.begin(), iend=constraintSet.end(); i!=iend; i++ )
        serr<<(*i)->getName()<<" ";
    serr<<sendl<<"BehaviorModel: ";
    for ( Sequence<BehaviorModel>::iterator i=behaviorModel.begin(), iend=behaviorModel.end(); i!=iend; i++ )
        serr<<(*i)->getName()<<" ";
    serr<<sendl<<"VisualModel: ";
    for ( Sequence<VisualModel>::iterator i=visualModel.begin(), iend=visualModel.end(); i!=iend; i++ )
        serr<<(*i)->getName()<<" ";
    serr<<sendl<<"CollisionModel: ";
    for ( Sequence<CollisionModel>::iterator i=collisionModel.begin(), iend=collisionModel.end(); i!=iend; i++ )
        serr<<(*i)->getName()<<" ";
    serr<<sendl<<"ContextObject: ";
    for ( Sequence<ContextObject>::iterator i=contextObject.begin(), iend=contextObject.end(); i!=iend; i++ )
        serr<<(*i)->getName()<<" ";
    serr<<sendl<<"Pipeline: ";
    for ( Single<Pipeline>::iterator i=collisionPipeline.begin(), iend=collisionPipeline.end(); i!=iend; i++ )
        serr<<(*i)->getName()<<" ";
    serr<<sendl<<"VisitorScheduler: ";
    for ( Single<VisitorScheduler>::iterator i=actionScheduler.begin(), iend=actionScheduler.end(); i!=iend; i++ )
        serr<<(*i)->getName()<<" ";
    serr<<sendl;
}

/** @name Dependency graph
This graph reflects the dependencies between the components. It is used internally to ensure that the initialization order is comform to the dependencies.
*/
/// @{
// Vertices
struct component_t
{
    typedef boost::vertex_property_tag kind;
};
typedef boost::property<component_t, BaseObject::SPtr> VertexProperty;

// Graph
typedef ::boost::adjacency_list < ::boost::vecS, ::boost::vecS, ::boost::bidirectionalS, VertexProperty > DependencyGraph;

void Node::sortComponents()
{
    if (depend.getValue().empty())
        return;
    typedef DependencyGraph::vertex_descriptor Vertex;
    DependencyGraph dependencyGraph;
    // map vertex->component
    boost::property_map< DependencyGraph, component_t >::type  component_from_vertex = boost::get( component_t(), dependencyGraph );
    // map component->vertex
    std::map< BaseObject::SPtr, Vertex > vertex_from_component;

    // build the graph
    for (int i = object.size() - 1; i >= 0; i--) // in the reverse order for a final order more similar to the current one
    {
        Vertex v = add_vertex( dependencyGraph );
        component_from_vertex[v] = object[i];
        vertex_from_component[object[i]] = v;
    }
    assert( depend.getValue().size()%2 == 0 ); // must contain only pairs
    for ( unsigned i=0; i<depend.getValue().size(); i+=2 )
    {
        BaseObject* o1 = getObject( depend.getValue()[i] );
        BaseObject* o2 = getObject( depend.getValue()[i+1] );
        if ( o1==NULL ) cerr<<"Node::sortComponent, could not find object called "<<depend.getValue()[i]<<endl;
        else if ( o2==NULL ) cerr<<"Node::sortComponent, could not find object called "<<depend.getValue()[i+1]<<endl;
        else
        {
            boost::add_edge( vertex_from_component[o1], vertex_from_component[o2], dependencyGraph );
            //cerr<<"Node::sortComponents, added edge "<<o1->getName()<<" -> "<<o2->getName()<<endl;
        }
    }

    // sort the components according to the dependencies
    typedef std::vector< Vertex > container;
    container c;
    boost::topological_sort(dependencyGraph, std::back_inserter(c));

    // remove all the components
    for ( container::reverse_iterator ii=c.rbegin(); ii!=c.rend(); ++ii)
    {
        removeObject(component_from_vertex[*ii]);
    }

    // put the components in the right order
    //cerr << "Node::sortComponents, New component order: ";
    for ( container::reverse_iterator ii=c.rbegin(); ii!=c.rend(); ++ii)
    {
        addObject(component_from_vertex[*ii]);
        //cerr << component_from_vertex[*ii]->getName() << " ";
    }
    //cerr << endl;

}

#ifdef SOFA_SMP
Iterative::IterativePartition* Node::getFirstPartition()
{
    if(is_partition())
        return partition_;
    for (sofa::simulation::Node::ChildIterator it= child.begin(); it != child.end(); ++it)
    {
        sofa::simulation::Node *g=static_cast<sofa::simulation::Node *>(*it);
        if(g)
        {
            Iterative::IterativePartition* p= g->getFirstPartition();
            if(p)
                return p;
        }
    }
    return NULL;
}
#endif


template <class RealObject>
Node::SPtr Node::create( RealObject*, sofa::simulation::xml::Element<sofa::core::objectmodel::BaseNode>*& arg)
{
//    Node::SPtr obj=getSimulation()->createNewGraph(arg->getName());
    Node::SPtr obj=getSimulation()->createNewNode(arg->getName());
    obj->parse(arg);
    return obj;
}

Node::SPtr Node::create( const std::string& name )
{
    return getSimulation()->createNewNode(name);
}


SOFA_DECL_CLASS(Node)
//create method of Node called if the user wants the default node. The object created will depend on the simulation currently in use.
helper::Creator<xml::NodeElement::Factory, Node> NodeClass("default");

}

}
