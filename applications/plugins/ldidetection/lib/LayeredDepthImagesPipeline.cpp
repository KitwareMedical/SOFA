//
// C++ Interface: LayeredDepthImagesPipeline
//
// Description: Collision Pipeline performing collision detection and response
//
//
// Author: Francois Faure, Sebastien Barbier, Jeremie Allard, Florent Falipou
//
// Licence: QPL, See LICENCE.txt file that comes with this distribution
// Copyright: See COPYING file that comes with this distribution
//
//
// #include <sofa/component/collision/LayeredDepthImagesPipeline.h>
#include "LayeredDepthImagesPipeline.h"

#include <sofa/core/CollisionModel.h>
#include <sofa/simulation/common/Node.h>
#include <sofa/core/ObjectFactory.h>

#include <sofa/core/collision/ContactManager.h>

#include <sofa/component/collision/DefaultContactManager.h>
#include <sofa/component/collision/CubeModel.h>
#include <sofa/component/collision/RayModel.h>
#include <sofa/component/collision/TriangleModel.h>
#include <sofa/component/visualmodel/VisualModelImpl.h>
#include <sofa/simulation/common/Visitor.h>

#include <sofa/simulation/common/Simulation.h>

#include <sofa/core/collision/CollisionGroupManager.h>
#include <sofa/core/collision/BroadPhaseDetection.h>

#include "LDIDetection.h"

//#include <GL/gl.h>

#define VERBOSE(a) if (f_printLog.getValue()) a;



namespace sofa
{

namespace component
{

namespace collision
{

using namespace core;
using namespace core::objectmodel;
using namespace core::collision;
using namespace sofa::defaulttype;
typedef sofa::core::collision::ContactManager::ContactVector ContactVector;

sofa::core::ObjectFactory::ClassEntry* classTriangleModel;


SOFA_DECL_CLASS(LayeredDepthImagesPipeline);

int LayeredDepthImagesPipelineClass = core::RegisterObject("The LDI modeling pipeline")
        .add< LayeredDepthImagesPipeline >()
        .addAlias("LDIPipeline")
        .addLicense("QPL")
        .addAuthor("Francois Faure, Sebastien Barbier, Jeremie Allard, Florent Falipou");
;

LayeredDepthImagesPipeline::LayeredDepthImagesPipeline():
    contactResponse     (initData(&contactResponse,        true,
                                  "contactResponse",      "Use Pressure Based Force Field to do the collision response"))
  , bSelfCollision      (initData(&bSelfCollision,         false,
                                  "selfcollision",        "Detect self collision"))
  , Kpressure           (initData(&Kpressure,              1.0f,
                                  "pressure",             "Constant for Collision Response"))
  , Kselfpressure       (initData(&Kselfpressure,          1.0f,
                                  "pressureself",         "Constant for Self Collision Response"))
  , Kviscosity          (initData(&Kviscosity,             0.0f,
                                  "viscosity",            "Constant for fluid friction"))
  , resolution          (initData(&resolution,             (unsigned int)64,
                                  "resolution",           "Max Resolution of the depth texture"))
  , resolutionPixel     (initData(&resolutionPixel,        0.0f,
                                  "resolutionPixel",      "Resolution of the Pixels in millimeter:\nonly used if strictly positive"))
  , subdivision         (initData(&subdivision,            Vec<3,unsigned int>(),
                                  "subdivision",          "Subdivision of the volume, to provides multiple constraints by contact region"))
  , minSubdivision      (initData(&minSubdivision,         Vec<3,unsigned int>(),
                                  "minSubdivision",       "Minimal Subdivision of the volume in number of pixel used"))
  , bAdaptiveResolution(initData(&bAdaptiveResolution,     0.0,
                                 "adaptiveResolution",    "Use an Adaptive resolution for the depth peeling:\nonly used if strictly positive"))
  , depthBB             (initData(&depthBB,               (unsigned int)4,
                                  "depthBB",              "Depth of Hierarchical BoundingBoxes used for Broad detection"))
  , showBBox            (initData(&showBBox,              false,
                                  "showBBox",             "Draw the rendering Bounding Boxes"))
  , showInteractionPairs(initData(&showInteractionPairs,  false,
                                  "showInteractionPairs", "Draw the interactions found between pixels"))
  , showDepthPixels     (initData(&showDepthPixels,       false,
                                  "showDepthPixels",      "Draw the region corresponding to a pixel found in interaction"))
  , bGPUCollisionVolume (initData(&bGPUCollisionVolume,   false,
                                  "GPUCollisionVolume",   "Use GPU to compite collision volumes when available (G80+)"))
  , bGlobalRendering    (initData(&bGlobalRendering,      false,
                                  "globalRendering",      "Several models rendered in one step"))
  , bUseGlobalBBox      (initData(&bUseGlobalBBox,        false,
                                  "useGlobalBBox",        "Use the global Bounding Box"))
  , bGlobalBBox         (initData(&bGlobalBBox,
                                  "globalBBox",           "Global Bounding Box where the collisions will be computed"))
  , isInitialized (false)
{
    sofa::core::ObjectFactory::AddAlias("TriangleMeshModel", "TriangleModelInRegularGrid", true, &classTriangleModel);
    sofa::core::ObjectFactory::AddAlias("TriangleSetModel" , "TriangleModelInRegularGrid", true, &classTriangleModel);
    sofa::core::ObjectFactory::AddAlias("TriangleMesh"     , "TriangleModelInRegularGrid", true, &classTriangleModel);
    sofa::core::ObjectFactory::AddAlias("TriangleSet"      , "TriangleModelInRegularGrid", true, &classTriangleModel);
    sofa::core::ObjectFactory::AddAlias("TriangleModel"    , "TriangleModelInRegularGrid", true, &classTriangleModel);
    sofa::core::ObjectFactory::AddAlias("Triangle"         , "TriangleModelInRegularGrid", true, &classTriangleModel);
}

typedef simulation::Visitor::ctime_t ctime_t;

void LayeredDepthImagesPipeline::updateDataVisibility()
{


}

void LayeredDepthImagesPipeline::init()
{
    this->getContext()->get( LDI);
    if( ! LDI )
    {
        serr << "LayeredDepthImagesPipeline::init(). Can't find the LDI Detection component in the scene graph." << sendl;
        exit(0);
    }
    sofa::simulation::PipelineImpl::init();
    updateDataVisibility();
}
void LayeredDepthImagesPipeline::reinit()
{
    //Initialization of the LDI detection by Glew init and compilation of the shader for the depth peeling operations

    if (bAdaptiveResolution.getValue() > 0) LDI->setAdaptiveLength(bAdaptiveResolution.getValue());
    else                                    LDI->setAdaptiveLength(0);

    LDI->setResolution(resolution.getValue());
    LDI->setResolutionPixel(resolutionPixel.getValue());
    //Create the different vectors, compile the shader...
    LDI->setK(Kpressure.getValue(),Kselfpressure.getValue());
    LDI->setSelfCollision(bSelfCollision.getValue());
    LDI->setViscosity(Kviscosity.getValue());
    LDI->setSubdivision(subdivision.getValue(), minSubdivision.getValue());
    LDI->setDepthBB(depthBB.getValue());
    LDI->setVerbose(f_printLog.getValue());
    LDI->setGPUCollisionVolume(bGPUCollisionVolume.getValue());
    if (bUseGlobalBBox.getValue()) LDI->setBBox(bGlobalBBox.getValue());

    updateDataVisibility();
    LDI->initialize();
}

helper::set< std::string > LayeredDepthImagesPipeline::getResponseList() const
{

    helper::set< std::string > listResponse;
    listResponse.insert("LDI");
    listResponse.insert("LDIConstraint");

    return listResponse;

}

//***********************************************************************************************************************************
void LayeredDepthImagesPipeline::doCollisionReset()
{
    core::objectmodel::BaseContext* scene = getContext();
//    simulation::Node* node = dynamic_cast<simulation::Node*>(scene);
//    ctime_t t0 = 0;
    const std::string category = "collision";
    VERBOSE(std::cout << "Reset collisions"<<std::endl);
    // clear all contacts
    if (contactManager!=NULL)
    {
        //	    const sofa::helper::vector<Contact*>& contacts = contactManager->getContacts();
        const ContactVector& contacts = contactManager->getContacts();
        for (ContactVector::const_iterator it = contacts.begin(); it!=contacts.end(); it++)
        {
            (*it)->removeResponse();
        }
    }
    // clear all collision groups
    if (groupManager!=NULL)
    {
        groupManager->clearGroups(scene);
    }
}

//***********************************************************************************************************************************
void LayeredDepthImagesPipeline::doCollisionDetection(const sofa::helper::vector<core::CollisionModel*>& collisionModels)
{
#ifdef SOFA_DUMP_VISITOR_INFO
    sofa::simulation::Visitor::printNode("reInit");
#endif
    // std::cerr<<"LayeredDepthImagesPipeline::doCollisionDetection"<<std::endl;
    core::objectmodel::BaseContext* scene = getContext();
    simulation::Node* node = dynamic_cast<simulation::Node*>(scene);
//    ctime_t t0 = 0;
    const std::string category = "collision";

    //***********************************************************************************************************************************
    //LDI Detection
    //Initialization of the LDI detection by Glew init and compilation of the shader for the depth peeling operations

    if (!isInitialized )
    {
        reinit();isInitialized = true;
    }

    LDI->clear();
    //Bounding Box for a pair of collision model
    sofa::helper::vector  < sofa::helper::vector< TriangleModel*> > list_interaction;
    sofa::helper::vector  < BoundingBox >     list_BBox;

#ifdef SOFA_DUMP_VISITOR_INFO
    sofa::simulation::Visitor::printCloseNode("reInit");
#endif

#ifdef SOFA_DUMP_VISITOR_INFO
    sofa::simulation::Visitor::printNode("rayPick");
#endif
    VERBOSE(std::cout << "Find the Colliding Zone\n");
    sofa::helper::vector<CollisionModel*> vectBoundingVolume;
    bool RayPick = false;
    const bool continuous = intersectionMethod->useContinuous();
    const double dt       = getContext()->getDt();

    sofa::helper::vector<CollisionModel*>::const_iterator it;;
    sofa::helper::vector<CollisionModel*>::const_iterator itEnd = collisionModels.end();
    int nActive = 0;
//    bool desactivateMapping=true;

    for (it = collisionModels.begin(); it != itEnd; it++)
    {
        if (RayModel* rayPick = dynamic_cast< RayModel* >(*it))
        {
            RayPick=((CollisionModel*)rayPick)->getSize();
//            desactivateMapping = ((CollisionModel*)rayPick)->getSize() == 0;
            break;
        }
    }
    //Ray Pick doesn't work if depthBB == 0
#ifdef SOFA_DUMP_VISITOR_INFO
    sofa::simulation::Visitor::printCloseNode("rayPick");
#endif

#ifdef SOFA_DUMP_VISITOR_INFO
    sofa::simulation::Visitor::printNode("ComputeBoundingTree");
#endif
    const int RayPickPossible = RayPick && depthBB.getValue()==0;
    for (it = collisionModels.begin(); it != itEnd; it++)
    {

        if (!(*it)->isActive() /* || !(*it)->isSimulated()*/) continue;
        if (continuous)
            (*it)->computeContinuousBoundingTree(dt, 4*RayPickPossible + depthBB.getValue());
        else
            (*it)->computeBoundingTree(RayPickPossible + depthBB.getValue());
        vectBoundingVolume.push_back ((*it)->getFirst());
        ++nActive;
    }
#ifdef SOFA_DUMP_VISITOR_INFO
    sofa::simulation::Visitor::printCloseNode("ComputeBoundingTree");
#endif
    VERBOSE(std::cout << "Computed "<<nActive<<" BBoxs"<<std::endl);
#ifdef SOFA_DUMP_VISITOR_INFO
    sofa::simulation::Visitor::printNode("LDIBroadPhase");
#endif

    LDI->BroadPhase( collisionModels, list_interaction, list_BBox, node, this);
#ifdef SOFA_DUMP_VISITOR_INFO
    sofa::simulation::Visitor::printCloseNode("LDIBroadPhase");
#endif
    //In list_interaction we have the list of the pairs of elements (or single element, depending if the self collision has been activated)...
    //... that can enter in collision
    //list_BBox contains the BoundingBox where the collision could occur.
    //*******************************************************************************************************************************

#ifdef SOFA_DUMP_VISITOR_INFO
    sofa::simulation::Visitor::printNode("DepthPeeling");
#endif
    //*******************************************************************************************************************************
    //Do the Depth Peeling for each pair of colliding elements
    if (!bGlobalRendering.getValue())
    {
        for (unsigned int i=0;i<list_interaction.size();i++)
        {
            LDI->NarrowPhase( list_BBox[i], list_interaction[i], node, this );
        }
    }
    else
    {
        sofa::helper::vector< TriangleModel* > all_collisionModels = LDI->getAllCollisionModels();
        if (list_BBox.size() != 0)
        {
            BoundingBox globalBB=list_BBox[0];
            sofa::helper::vector< bool > BB_used(list_BBox.size(),false);

            for (unsigned int i=1;i<list_BBox.size();++i)
            {

                LDI->unionAABB(list_BBox[i],globalBB, globalBB);
            }
            LDI->NarrowPhase(globalBB, all_collisionModels, node, this);
        }
    }
    LDI->EndNarrowPhase( node, this );
    VERBOSE(std::cout << "\n\n\n\n");

#ifdef SOFA_DUMP_VISITOR_INFO
    sofa::simulation::Visitor::printCloseNode("DepthPeeling");
#endif



    //***********************************************************************************************************************************
    //Support for Ray Interaction

    if (broadPhaseDetection==NULL || !RayPick || nActive<2) return; // can't go further
#ifdef SOFA_DUMP_VISITOR_INFO
    sofa::simulation::Visitor::printNode("RayPickHack");
#endif
    VERBOSE(std::cout << "BroadPhaseDetection "<<broadPhaseDetection->getName()<<std::endl);
    broadPhaseDetection->beginBroadPhase();
    broadPhaseDetection->addCollisionModels(vectBoundingVolume);  // detection is done there
    broadPhaseDetection->endBroadPhase();
#ifdef SOFA_DUMP_VISITOR_INFO
    sofa::simulation::Visitor::printCloseNode("RayPickHack");
#endif

    // then we start the narrow phase
    if (narrowPhaseDetection==NULL) return; // can't go further
#ifdef SOFA_DUMP_VISITOR_INFO
    sofa::simulation::Visitor::printNode("NarrowPhase");
#endif
    VERBOSE(std::cout << "NarrowPhaseDetection "<<narrowPhaseDetection->getName()<<std::endl);
    narrowPhaseDetection->beginNarrowPhase();

    sofa::helper::vector<std::pair<CollisionModel*, CollisionModel*> >& vectCMPair = broadPhaseDetection->getCollisionModelPairs();
    sofa::helper::vector<std::pair<CollisionModel*, CollisionModel*> > mouseInteraction;

    for (unsigned int i=0;i<vectCMPair.size();i++)
    {
        if (vectCMPair[i].first->getContext()->getName() =="mouse" ||
                vectCMPair[i].second->getContext()->getName() == "mouse")
            mouseInteraction.push_back(vectCMPair[i]);
    }
    VERBOSE(std::cout << mouseInteraction.size()<<" colliding model pairs"<<std::endl);
    narrowPhaseDetection->addCollisionPairs(mouseInteraction);
    narrowPhaseDetection->endNarrowPhase();
#ifdef SOFA_DUMP_VISITOR_INFO
    sofa::simulation::Visitor::printCloseNode("NarrowPhase");
#endif
}

void LayeredDepthImagesPipeline::doCollisionResponse()
{
    core::objectmodel::BaseContext* scene = getContext();
//    simulation::Node* node = dynamic_cast<simulation::Node*>(scene);
//    ctime_t t0 = 0;
    const std::string category = "collision";

    // then we start the creation of contacts
    if (contactManager==NULL) return; // can't go further
    //if (LDI->getDetectionOutputs().size() == 0) return;
    VERBOSE(std::cout << "Create Contacts "<<contactManager->getName()<<std::endl);

    VERBOSE(std::cout << "object interactions : "   <<LDI->getDetectionOutputs().size()<<std::endl);
    NarrowPhaseDetection::DetectionOutputMap &list_collision = LDI->getDetectionOutputs();
    const NarrowPhaseDetection::DetectionOutputMap &list_mouse = narrowPhaseDetection->getDetectionOutputs();

    if (list_mouse.size() == 1) list_collision.insert( *list_mouse.begin() );

    contactManager->createContacts(list_collision);

    const ContactVector& contacts = contactManager->getContacts();


    // finally we start the creation of collisionGroup

    // First we remove all contacts with static objects and directly add them
    ContactVector notStaticContacts;

    for (ContactVector::const_iterator it = contacts.begin(); it!=contacts.end(); it++)
    {
        Contact* c = it->get();
        if (!c->getCollisionModels().first->isSimulated())
        {
            c->createResponse(c->getCollisionModels().second->getContext());
        }
        else if (!c->getCollisionModels().second->isSimulated())
        {
            c->createResponse(c->getCollisionModels().first->getContext());
        }
        else if (c->getCollisionModels().first->getContext() == c->getCollisionModels().second->getContext()) // self-collision
        {
            c->createResponse(c->getCollisionModels().first->getContext());
        }
        else
            notStaticContacts.push_back(c);
    }


    if (groupManager==NULL)
    {
        VERBOSE(std::cout << "Linking all contacts to Scene"<<std::endl);
        for (ContactVector::const_iterator it = notStaticContacts.begin(); it!=notStaticContacts.end(); it++)
        {
            (*it)->createResponse(scene);
        }
    }
    else
    {
        VERBOSE(std::cout << "Create Groups "<<groupManager->getName()<<std::endl);
        groupManager->createGroups(scene, notStaticContacts);

    }

    VERBOSE(std::cout << "Response Done "<<std::endl);



}






} // namespace collision

} // namespace component

} // namespace sofa

