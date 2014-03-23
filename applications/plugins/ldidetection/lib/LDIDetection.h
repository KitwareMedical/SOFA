//
// C++ Interface: LDIDetection
//
// Description: Implements the two phases of the collision detection, BroadPhase and NarrowPhase.
//
//
// Author: Francois Faure, Sebastien Barbier, Jeremie Allard, Florent Falipou
//
// Licence: QPL, See LICENCE.txt file that comes with this distribution
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SOFA_COMPONENT_COLLISION_LDIDETECTION_H
#define SOFA_COMPONENT_COLLISION_LDIDETECTION_H

//#include <sofa/component/collision/DepthPeeling.h>
#include "DepthPeeling.h"
#include "DepthPeelingUtility.h"


//#include <sofa/core/VisualModel.h>
#include <sofa/component/collision/TriangleModel.h>

#include <sofa/core/collision/DetectionOutput.h>
#include <sofa/core/collision/CollisionAlgorithm.h>

#include <sofa/core/collision/ContactManager.h>
#include <sofa/core/CollisionElement.h>

#include <sofa/component/collision/RayModel.h>

#include <sofa/simulation/common/Node.h>

#include <sofa/defaulttype/Vec.h>
namespace sofa
{

  namespace component
  {

    namespace collision
    {
      using namespace sofa::defaulttype;
      using namespace core::collision;


      /***************************************************************************************/

      struct LDIDetectionOutput
      {
        typedef int64_t ContactId;
        /// Pair of colliding elements.
        std::pair<core::CollisionElementIterator, core::CollisionElementIterator> elem;

        /// Unique id of the contact for the given pair of collision models.
        ContactId id;
        /// Contact points on the surface of each model. They are expressed in the local coordinate system of the model if any is defined..
        Vec<2, float> point[2];
        /// Normal of the contact, pointing outward from the first model
        int normal;
        /// Volume described by this output
        float value;

        float subdivisionFactor;
        double *volume;

        unsigned int index;

        /// If using a continuous collision detection, estimated of time of contact.
        float dS;
        float dx;
        float dy;


        float K;
        float viscosity;
        Vec3d dVelocity;

        LDIDetectionOutput()
        : elem (),id ( 0 ),normal ( 0 ), subdivisionFactor(1.0f), dS ( 0 ),K ( 0.0f ), dVelocity ( Vec3d() )
        {
        }
      };

      /***************************************************************************************/

      class LDIDetectionOutputVector : public DetectionOutputVector
      {
        public:

          ~LDIDetectionOutputVector()
          {
          }

          unsigned int size() const
          {
            return Vector.size();
          }

          void clear()
          {
            Vector.clear();
          }

          sofa::helper::vector< LDIDetectionOutput > Vector;

      };



      typedef std::pair< Vector3, Vector3 > BoundingBox;

      class LDIDetection: public CollisionAlgorithm
      {
        public :
        SOFA_CLASS(LDIDetection,CollisionAlgorithm);
          typedef NarrowPhaseDetection::DetectionOutputMap    DetectionOutputMap;

          LDIDetection ( bool verbose=false, unsigned int resolution_=32 )
            :   bGPUCollisionVolume ( false ),  Kpressure ( 200 ), selfCollision ( false ), adaptiveLength ( 0.0 ), frame ( 0 ), text ( verbose ), resolution ( resolution_ ), depth_BB ( 4 ), initialized ( false ), first_triangle ( NULL ), current_triangle ( NULL )
          {
						
          };

          ~LDIDetection()
          {
            if ( first_triangle ) delete first_triangle;
            if ( current_triangle ) delete current_triangle;
          };

          virtual void changeInstance ( Instance /*inst*/ )
          {
            // change the collision pipeline instance.
            // TODO
            // We should store all the data in maps depending on the instance and switch the context here...
          }


          void init();

          ///Allocation of memory for the textures, layers... Done only once for a kind of collision detection
          void initialize();

          ///Clear the vector containing the detection output vectors
          void clear();

          ///Optimization step: calculate hierarchical bounding box to determine intersection.
          ///Only the intersecting bounding box will be kept for further collision detection
          ///  * collisionModels contains the initial list of all collisionModels: only those with a triangle model description will be used
          ///  * the other lists will contain the result from the search
          ///    --> list_collision contains a list of all the collisions we need to study. They are stored as index correponding to a triangleModel of list_triangleModel
          ///    --> list_BBox contains the intersection BB where the collision will be searched. It is stored in the same order as list_collision
          void BroadPhase ( const sofa::helper::vector< core::CollisionModel* >               &collisionModels,
                            sofa::helper::vector  < sofa::helper::vector< TriangleModel*> >   &list_collision,
                            sofa::helper::vector  < BoundingBox >                             &list_BBox,
                            simulation::Node* timelog = NULL, core::objectmodel::BaseObject* timeobj = NULL );

          ///We make at the same time a depth peeling and an exploration of the texture resulting of the depth peeling.
          ///  * interval is the intervals of coordinates in the world coordinates, where we will look after a collision
          ///  * BB is the bounding box where the depth peeling will be performed
          ///  * list_index contains the index of the collision model we want to study.
          ///  * timelog is an optionnal pointer used only to log computation times used by each step
          void NarrowPhase ( const BoundingBox &BB,
                             const sofa::helper::vector< TriangleModel*> &list_index,
                             simulation::Node* timelog = NULL, core::objectmodel::BaseObject* timeobj = NULL );

          void EndNarrowPhase ( simulation::Node* timelog = NULL, core::objectmodel::BaseObject* timeobj = NULL );

          bool retrieveInformationFromPixel( const unsigned int &layerDepth, const unsigned int &index, const bool &globalRendering,
                                             bool &front, unsigned int &indexModel, 
                                             float &alpha, float &beta) const ;

          void draw ( bool BBox, bool interactions, bool regions );

          void draw();

          //Once modified, the LDI Detection must be initialized
          //Interface to set the parameters of the detection
          void                setResolution ( const unsigned int r ) {resolution=r; initialized=false;}
          void                setResolutionPixel ( const float r )  {resolutionPixel=r; initialized=false;}
          void                setVerbose ( bool verbose )           { text=verbose; }
          void                setOrientation ( int orient )         { current_orientation = orient%3; }
          void                setK ( float kp, float kselfp )       { Kpressure = kp;Kselfpressure = kselfp; }
          void                setViscosity ( float v )              { Kviscosity=v; }
          void                setSelfCollision ( bool b )           { selfCollision = b; }
          void                setDepthBB ( unsigned int d )         { depth_BB = d; }
          void                setGPUCollisionVolume ( bool val )    { bGPUCollisionVolume = val; }
          void                setAdaptiveLength ( double l )        { adaptiveLength = l;}
          void                setSubdivision ( const Vec<3,unsigned int> &s,const Vec<3,unsigned int> &minS )   
          {            
            for (unsigned int i=0;i<3;++i)
              if (s[i]==0) subDivisionBB[i]=1;
              else         subDivisionBB[i]=s[i];
            minSubDivisionBB=minS;
          }
          void                setBBox ( Vec6f val )                 { globalBBox.first = Vector3 ( val[0],val[1],val[2] ); globalBBox.second = Vector3 ( val[3],val[4],val[5] ); }

          unsigned int        getResolution() const              {return resolution;}
          DetectionOutputMap &getDetectionOutputs()             {return detectionOutputs; }


          sofa::helper::vector< TriangleModel* > getAllCollisionModels() {return all_collisionModels;};
          const sofa::helper::vector< sofa::helper::vector< float> >& getDeepthLayers() { return depth_layer;}

          void crossingAABB ( const BoundingBox &b0, const BoundingBox &b1, bool &isCrossing, BoundingBox &resultingBBox ) const;
          void unionAABB ( const BoundingBox &b0, const BoundingBox &b1, BoundingBox &resultingBBox ) const;
					
         
          
         

          
        protected :

          bool validForNarrowPhase( const TriangleModel *m1, const TriangleModel *m2) const;

#ifdef SOFA_GPU_CUDA
          void NarrowPhaseCUDA ( const BoundingBox &BB,const  sofa::helper::vector< TriangleModel*> &models, const sofa::helper::vector<unsigned int > &index_models, bool autoCollision,
                                 LDIDetectionOutputVector *outputs, simulation::Node* timelog, core::objectmodel::BaseObject* timeobj );
#endif
          void NarrowPhaseCPU ( const BoundingBox &BB,const  sofa::helper::vector< TriangleModel*> &models, const sofa::helper::vector<unsigned int > &index_models, bool autoCollision,
                                LDIDetectionOutputVector *outputs, simulation::Node* timelog, core::objectmodel::BaseObject* timeobj );

          //Working on the AABBs
          inline void constructBoundingHierarchy ( TriangleModel *m, sofa::helper::vector< sofa::helper::vector< BoundingBox > > &boundingHierarchy ) const;
          void        findIntersectionBoundingBox ( const BoundingBox  &m0BB, const BoundingBox &m1BB,
                  bool &isCrossing,BoundingBox &result ) const;
          inline BoundingBox computeCorrectVisualizationBB ( const BoundingBox &bb,  const sofa::helper::vector< unsigned int > index_models ) const;

          inline bool inCollision ( const bool same_object,
                                    const Vec3d &normal0,const Vec3d &normal1, int &SC_ray_cast ) const ;


          //Elements need to perform the depth peeling

          unsigned char current_orientation;                   //index between 0 and 2 corresponding to one of the three orientation

          sofa::helper::vector< sofa::helper::vector< float> > depth_layer; //layers of the depth peeling


          //Elements used for the collision Detection
          DetectionOutputMap detectionOutputs;
          sofa::helper::vector< TriangleModel* > all_collisionModels; //list of the triangle models present in the scene
          sofa::helper::vector< BoundingBox    > all_collisionBB;     //list of the triangle models present in the scene
          sofa::helper::vector  < sofa::helper::vector< TriangleModel*> >  all_collisionPair; //list of the pairs of colliding elements
          sofa::helper::vector< BoundingBox >    visualization_BBox;
          std::set< TriangleModel* > index_used;                        //index of the list "all_collisionModels" corresponding to triangle models used in the scene

          bool bGPUCollisionVolume;

          float Kpressure;     //constant used for the computation of the resulting pressure to perform the collision response
          float Kselfpressure;
          float Kviscosity;
          float attenuation; //damping coeff
          //global boolean: if self collision must be operated. Should be better if this boolean was present in each collision model.
          bool selfCollision;
          double adaptiveLength;
          DepthPeelingUtility depthPeelingProcess;
          Vec<3,unsigned int> subDivisionBB; 
          Vec<3,unsigned int> minSubDivisionBB; 


          struct NarrowPhaseTest
          {
            core::CollisionModel* model1;
            core::CollisionModel* model2;
            int result_id[3];
            float K;
            float Kviscosity;
          };

          sofa::helper::vector<NarrowPhaseTest> narrowPhaseTests;


        private:
          //parameters
          int frame;
          bool text;
          unsigned int resolution;
          float resolutionPixel;

          unsigned int depth_BB;
          BoundingBox globalBBox;
          bool initialized;
          unsigned int id_collision;
          Triangle *first_triangle;
          Triangle *current_triangle;

          
      };
    }
  }
}

#endif
