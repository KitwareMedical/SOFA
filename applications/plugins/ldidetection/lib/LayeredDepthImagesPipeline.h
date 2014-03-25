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
#ifndef SOFA_COMPONENT_COLLISION_LDIPIPELINE_H
#define SOFA_COMPONENT_COLLISION_LDIPIPELINE_H

/* #include <sofa/component/collision/LDIDetection.h> */
/* #include <sofa/component/collision/DepthPeeling.h> */

#include "LDIDetection.h"
#include "DepthPeeling.h"
#include <sofa/core/collision/Pipeline.h>
#include <sofa/simulation/common/PipelineImpl.h>
//#include <sofa/core/VisualModel.h>
#include <sofa/component/collision/TriangleModel.h>

//#include "headers.h"



namespace sofa
{

  namespace component
  {

    namespace collision
    {
      using namespace sofa::defaulttype;
      typedef std::pair< Vector3, Vector3 > BoundingBox;

      class SOFA_LDIDETECTION_API LayeredDepthImagesPipeline : public sofa::simulation::PipelineImpl
	{

      public:
          SOFA_CLASS(LayeredDepthImagesPipeline,sofa::simulation::PipelineImpl);
 	Data<bool>           contactResponse; 
	Data<bool>           bSelfCollision;
	Data<bool>           pressureResponse;
	Data<float>          Kpressure;
	Data<float>          Kselfpressure;	
	Data<float>          Kviscosity;
	Data<unsigned int>   resolution;
	Data<float>   resolutionPixel;
	Data<Vec<3,unsigned int> >  subdivision;
	Data<Vec<3,unsigned int> >  minSubdivision;
	Data<double>         bAdaptiveResolution;
	Data<unsigned int>   depthBB;
	Data<bool>           showBBox;	
	Data<bool>           showInteractionPairs;
	Data<bool>           showDepthPixels;	
        Data<bool>           bGPUCollisionVolume;
	Data<bool>           bGlobalRendering;
	Data<bool>           bUseGlobalBBox;
	Data<Vec<6,float> >  bGlobalBBox;
        LayeredDepthImagesPipeline();
        ~LayeredDepthImagesPipeline(){};

  // -- Pipeline interface 
  /// get the set of response available with the current collision pipeline
   helper::set< std::string > getResponseList() const;
	// -- VisualModel interface
	void draw() { LDI->draw(showBBox.getValue(), showInteractionPairs.getValue(), showDepthPixels.getValue());}; 
	void initTextures() { }
	void update() { }

	void init();
	void reinit();
	void updateDataVisibility();
      protected:

	// -- Pipeline interface

	/// Remove collision response from last step
	virtual void doCollisionReset();
	/// Detect new collisions. Note that this step must not modify the simulation graph
	virtual void doCollisionDetection(const sofa::helper::vector<core::CollisionModel*>& collisionModels);
	/// Add collision response in the simulation graph
	virtual void doCollisionResponse();
	

      private:
	LDIDetection* LDI;
	bool isInitialized;
      };

    } // namespace collision

  } // namespace component

} // namespace sofa

#endif
