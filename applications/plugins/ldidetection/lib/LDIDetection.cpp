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
// #include <sofa/component/collision/LDIDetection.h>
#include "LDIDetection.h"
#include <sofa/helper/system/gl.h>
#include <sofa/helper/gl/template.h>
#include <sofa/helper/gl/BasicShapes.h>
#include <math.h>
#include <sofa/core/CollisionModel.h>
#include <sofa/simulation/common/Node.h>
#include <sofa/core/ObjectFactory.h>


#include "DepthPeelingUtility.h"

#include <sofa/component/collision/CubeModel.h>
#include <sofa/component/visualmodel/VisualModelImpl.h>

#ifdef SOFA_DUMP_VISITOR_INFO
#include <sofa/simulation/common/Visitor.h>
#endif
#include <sofa/simulation/common/Simulation.h>

namespace sofa
{

  namespace component
  {

    namespace collision
    {
#define EPSILON 0.25
#define X  (current_orientation+1)%3
#define Y  (current_orientation+2)%3
#define exact_volume 0
      using namespace core;
      using namespace core::objectmodel;
      using namespace core::collision;
      using namespace sofa::defaulttype;

      //***********************************************************************************************************************************
      // Factory declaration
      SOFA_DECL_CLASS ( LDIDetection );

      int LDIDetectionClass = core::RegisterObject ( "The LDI colision detection" )
                              .add< LDIDetection >()
                              .addLicense ( "QPL" )
                              .addAuthor ( "Francois Faure, Sebastien Barbier, Jeremie Allard, Florent Falipou" );


      void LDIDetection::init()
      {
       
      }

      //***********************************************************************************************************************************
      // Memory allocation
      void LDIDetection::initialize()
      {
        if ( text ) std::cout << "Initialize\n";
        depthPeelingProcess.init ( resolution, resolutionPixel, bGPUCollisionVolume );
        depth_layer.clear();

        current_orientation = 0;

        first_triangle   = new Triangle();
        current_triangle = new Triangle();

        initialized = true;
      }

      //***********************************************************************************************************************************
      ///Broad Phase :
      ///We search pairs of collision models (made of triangle) whose BB could be in intersection
      ///We construct:
      ///             * a list of pair of elements (or single element for self collision) that can collide (list_Collision)
      ///             * a list of bounding box (list_BBox) corresponding to the parameters of the orthographic camera that will be used to perform the depth peeling
      //*******************************************************************************************************************************

      void LDIDetection::BroadPhase ( const sofa::helper::vector< core::CollisionModel* > &collisionModels,
                                      sofa::helper::vector  < sofa::helper::vector< TriangleModel*> >   &list_collision,
                                      sofa::helper::vector  < BoundingBox >                             &list_BBox,
                                      simulation::Node* /*timelog*/, core::objectmodel::BaseObject* /*timeobj*/ )
      {



#ifdef SOFA_DUMP_VISITOR_INFO
        sofa::simulation::Visitor::printNode("BroadPhase");
#endif

       

//        simulation::Visitor::ctime_t t0 = 0;
        if ( text ) std::cout << "****************************************************************************\nInteractions studied : \n";


        sofa::helper::vector< std::pair< TriangleModel *, sofa::helper::vector< BoundingBox > > > list_triangleModel;
        //indicate if the list has been entirely created
        bool list_triangleModel_initialized = false;

        //depthPeelingProcess.setNbModels(collisionModels.size());

        bool isCrossing;
        BoundingBox intersectionBoundingBox;

        //Hierarchical Bounding Box for the two models in collision
        sofa::helper::vector< BoundingBox > modelBB[2];


        std::map< unsigned int, unsigned int> corresponding_index;
        std::set< TriangleModel* >::iterator item_found;

        TriangleModel *mesh;
        //For each collision models
        for ( unsigned int col_model=0; col_model < collisionModels.size(); col_model++ )
        {

          //First, we find one triangle model:
          if ( ! ( mesh = dynamic_cast< TriangleModel *> ( collisionModels[col_model] ) ) ) continue;
          if ( !mesh->isActive() ) continue;
          if ( frame==0 ) index_used.insert ( mesh );
          //if the complete list of all the triangle models has been created, we just have to read it. No need to compute again the Bounding Box
          if ( !list_triangleModel_initialized )
          {
            //Get the Hierarchical Bounding Box for the current TriangleModel
            sofa::helper::vector< BoundingBox > hierarchicalBB;
            //std::cout << "constructBoundingHierarchy..."<<std::endl;

            //if ( model->getPrevious() ) model = model->getPrevious(); // skip the first layer of bounding cubes
            dynamic_cast<CubeModel*>(mesh->getPrevious())->getBoundingTree( hierarchicalBB);

            //constructBoundingHierarchy ( mesh, hierarchicalBB );
            /*
            std::cout << "hierarchicalBB =";
            for (unsigned i=0;i<hierarchicalBB.size();++i)
            {
                std::cout << " [";
                for (unsigned j=0;j<hierarchicalBB[i].size();++j)
              std::cout << " <" << hierarchicalBB[i][j].first<<">-<" << hierarchicalBB[i][j].second<<">";
                std::cout << "]";
            }
            std::cout << std::endl;
            */
            if ( hierarchicalBB.size() == 0 ) {continue;}
            //Store the new triangle model in the list with its own hierarchical Bounding Box
            corresponding_index[col_model] = list_triangleModel.size();
            list_triangleModel.push_back ( std::make_pair ( mesh, hierarchicalBB ) );
            //List of all the triangle models present in the current scene
            all_collisionBB.push_back ( hierarchicalBB[0] );
            modelBB[0] = hierarchicalBB;
            //If the self collision is activated, and the mesh is not static, we add it in the list of the collision to perform
            if ( selfCollision && mesh->isSimulated() && mesh->canCollideWith ( mesh ) )
            {
              item_found = index_used.find ( mesh );
              if ( item_found == index_used.end() )
              {
                index_used.insert ( mesh );
                sofa::helper::vector< TriangleModel*> selfcollision_vec ( 1,mesh );
                list_collision.push_back ( selfcollision_vec );
                list_BBox.push_back ( hierarchicalBB[0] );
              }
            }
          }
          else
          {
            //second element: the hierarchical bounding box
            modelBB[0] = list_triangleModel[corresponding_index[col_model]].second;
          }

          //Then we try to compare its bounding box, with the bounding box of another triangle collision model
          for ( unsigned int col_modelCompare = col_model+1; col_modelCompare< collisionModels.size(); col_modelCompare++ )
          {

            TriangleModel *mesh2;
            if ( ! ( mesh2 = dynamic_cast< TriangleModel *> ( collisionModels[col_modelCompare] ) ) ) continue;

            if ( !mesh2->isActive() ) continue;


            if ( !list_triangleModel_initialized )
            {

              //We only work with triangle model
              //Get the Hierarchical Bounding Box for the first TriangleModel
              sofa::helper::vector< BoundingBox > hierarchicalBB;

						//if ( model->getPrevious() ) model = model->getPrevious(); // skip the first layer of bounding cubes
							dynamic_cast<CubeModel*>(mesh2->getPrevious())->getBoundingTree( hierarchicalBB);
							//constructBoundingHierarchy ( mesh2, hierarchicalBB );

              if ( hierarchicalBB.size() == 0 ) {continue;}
              corresponding_index[col_modelCompare] = list_triangleModel.size();
              list_triangleModel.push_back ( std::make_pair ( mesh2, hierarchicalBB ) );
              all_collisionBB.push_back ( hierarchicalBB[0] );

              modelBB[1] = hierarchicalBB;

              if ( selfCollision && mesh2->isSimulated() && mesh2->canCollideWith ( mesh2 ) )
              {
                item_found = index_used.find ( mesh2 );
                if ( item_found == index_used.end() )
                {
                  index_used.insert ( mesh2 );

                  sofa::helper::vector< TriangleModel*> selfcollision_vec ( 1,mesh2 );
                  list_collision.push_back ( selfcollision_vec );

                  list_BBox.push_back ( hierarchicalBB[0] );

                }
              }
            }
            else
            {
              //second element: the hierarchical bounding box
              modelBB[1] = list_triangleModel[corresponding_index[col_modelCompare]].second;
            }

            if (!validForNarrowPhase(mesh, mesh2)) continue;
            

            if ( text ) std::cout << "Collision between : {"
                                  << mesh->getName()
                                  << " -- "
                                  << mesh2->getName()  << "}";

            //using the hierarchical bounding box of the two triangle models, we find an intersection.
            //if one exists, then isCrossing is true, and :
            //    * intersectionBoundingBox contains relevant information on the intersecting AABB
            findIntersectionBoundingBox ( modelBB[0][0], modelBB[1][0], isCrossing, intersectionBoundingBox );

            if ( isCrossing )
            {
//               if ( intersectionBoundingBox.second[0] - intersectionBoundingBox.first[0]<0.00001 ||
//                    intersectionBoundingBox.second[1] - intersectionBoundingBox.first[1]<0.00001 ||
//                    intersectionBoundingBox.second[2] - intersectionBoundingBox.first[2]<0.00001 ) continue;

              // Avoid to detect the collisions between a cutting model and its corresponding collision model.
              if( mesh->getTopology() == mesh2->getTopology()) continue;

              //We update the list of collision to study
              //we memorize the pair of triangle model that could collide, and the AABB where it could happen
              sofa::helper::vector< TriangleModel*> collision_pair;
              collision_pair.push_back ( mesh );
              collision_pair.push_back ( mesh2 );

              index_used.insert ( mesh );
              index_used.insert ( mesh2 );

              list_collision.push_back ( collision_pair );
              list_BBox.push_back ( intersectionBoundingBox );

              if ( text ) std::cout
                << " will be calculated between "
                << list_triangleModel[corresponding_index[col_model]].first->getName()
                << " -- "
                << list_triangleModel[corresponding_index[col_modelCompare]].first->getName()
                << " =>"
                <<  " [" << intersectionBoundingBox.first << "] [" << intersectionBoundingBox.second << "]\n";
            }
            else
            {
              if ( text ) std::cout << " Not Found\n";
            }
          }
          //Once all the collision models visited, we know that the list_triangleModel contains all the information needed
          //no need to recompute hierarchical boundingboxes
          list_triangleModel_initialized = true;
        }

        std::set< TriangleModel* >::const_iterator item;

        for ( item = index_used.begin(); item != index_used.end(); item++ )
        {
          all_collisionModels.push_back ( ( *item ) );
        }

        depthPeelingProcess.loadModels ( all_collisionModels );
        //visualization_BBox = list_BBox;
        visualization_BBox.clear();
        all_collisionPair = list_collision;

#ifdef SOFA_DUMP_VISITOR_INFO
        sofa::simulation::Visitor::printCloseNode("BroadPhase");
#endif
      }     

   
      bool LDIDetection::validForNarrowPhase( const TriangleModel *mesh1, const TriangleModel *mesh2) const
      {
        
        //If the two collision models are not simulated, there is no need to send them to the narrow phase
        if (!mesh1->isSimulated() && !mesh2->isSimulated()) return false;
        //Special case of Pushing & Cutting pair of triangle models
        if ( (mesh1->hasTag(Tag ( "Cutting" )) && mesh2->hasTag(Tag ( "Pushing" )) ) ||
             (mesh1->hasTag(Tag ( "Pushing" )) && mesh2->hasTag(Tag ( "Cutting" )) ) ) return false;
        
        const sofa::core::objectmodel::TagSet &tags = mesh1->getTags();
        
        //If one of the mesh has no tags, we validate them
        if (mesh1->getTags().empty() || mesh2->getTags().empty()) return true;
        
        //If both of them have a set of tags, we need to find at least one compatible

        for (sofa::core::objectmodel::TagSet::const_iterator it=tags.begin(); it!=tags.end(); ++it)
          {
            if (mesh2->hasTag(*it)) return true;
          }

        return false;
      }


      //************************************************************************************************************************************
      /// Narrow Phase : given a bounding box of visualization and a list of indexes corresponding to the triangle models stored in the VBOs
      ///                fill the map containing the set of collision outputs.
      void LDIDetection::NarrowPhase ( const BoundingBox &BB_visu,
                                       const  sofa::helper::vector< TriangleModel*> &models,
                                       simulation::Node* timelog, core::objectmodel::BaseObject* timeobj )
      {
#ifdef SOFA_DUMP_VISITOR_INFO
        sofa::simulation::Visitor::printNode("NarrowPhase");
#endif
        if ( text )
        {
          std::cout << "Using DepthPeeling Method \n";
          std::cout << "-------------------------------------------------------------------------\n";
          std::cout << "Collision studied in [" << BB_visu.first << "] [" << BB_visu.second << "] : ";
          for ( unsigned int i=0;i<models.size();i++ )
          {
            std::cout << models[i]->getName() << // " " << (*it_list_index)<<
            "\t";
          }
          std::cout << "\n";
        }
        visualization_BBox.push_back ( BB_visu );

        sofa::helper::vector<unsigned int > index_models;
        if ( models.size() <= 2 ) index_models.resize ( 2,0 );
        else index_models.resize ( models.size(),0 );

        for ( unsigned int i=0,found=0;i<models.size() && found <index_models.size();++i )
        {
          for ( unsigned int j=0;j<all_collisionModels.size();++j )
          {
            if ( all_collisionModels[j] == models[i] ) {index_models[found] = j;found++;}
          }

        }

        bool globalRendering = models.size() > 2;

        //if only one model is given as input, it means that a self collision detection if wanted

        bool autoCollision = models.size() == 1;

        //Create the Vector that will be used by the contact manager to create the collision response
        LDIDetectionOutputVector *outputs=NULL;

        if ( !globalRendering )
        {
          DetectionOutputVector * &o=detectionOutputs[std::make_pair ( models[0],models[autoCollision ? 0 : 1] ) ];
          if ( o == NULL )
            o = new LDIDetectionOutputVector();
          outputs = dynamic_cast< LDIDetectionOutputVector *> ( o );
        }
        else
        {
          for ( unsigned int i=0;i<all_collisionPair.size();++i )
          {
            DetectionOutputVector * &o=detectionOutputs[std::make_pair ( all_collisionPair[i][0],all_collisionPair[i][1] ) ];
            if ( o == NULL )
              o = new LDIDetectionOutputVector();
          }
        }


#ifdef SOFA_GPU_CUDA
        if ( !globalRendering && bGPUCollisionVolume )
        {
          NarrowPhaseCUDA ( BB_visu,models, index_models, autoCollision,  outputs, timelog, timeobj );
        }
        else
#endif
        {
          NarrowPhaseCPU ( BB_visu,models, index_models, autoCollision, outputs,  timelog, timeobj );
        }

#ifdef SOFA_DUMP_VISITOR_INFO
        sofa::simulation::Visitor::printCloseNode("NarrowPhase");
#endif
      }

#ifdef SOFA_GPU_CUDA
      void LDIDetection::NarrowPhaseCUDA ( const BoundingBox &BB_visu,const sofa::helper::vector< TriangleModel*> &models, const sofa::helper::vector<unsigned int > &index_models,bool autoCollision,
                                           LDIDetectionOutputVector *outputs,
                                           simulation::Node* timelog, core::objectmodel::BaseObject* timeobj )
      {
        simulation::Node::ctime_t t0 = 0;
        NarrowPhaseTest t;
        BoundingBox BB = BB_visu;
        t.model1 = models[0];
        t.model2 = ( models.size() ==1 ) ? models[0] : models[1];
        float K;
        if ( autoCollision ) K = Kselfpressure;
        else                K = Kpressure;
        t.K = K;
        t.Kviscosity = Kviscosity;
        for ( int orientation=0; orientation<3; ++orientation )
        {
          current_orientation = orientation;
          BB = computeCorrectVisualizationBB ( BB_visu, index_models );
          depthPeelingProcess.doCollisionVolume ( models, t.result_id[orientation], BB, orientation, timelog, timeobj );
          if ( t.result_id[orientation] < -1 )
          {
            int ncolls = -1-t.result_id[orientation];
            const float* collisions = depthPeelingProcess.readCollisionVolume ( // t.result_id[orientation]
                                      );
            if ( collisions!=NULL )
            {
              outputs->Vector.resize ( outputs->size() +ncolls );
              LDIDetectionOutput *detection = &* ( outputs->Vector.end()-ncolls );
              TriangleModel *first_model, *current_model;
              //Values to save redondant computations
              const float depth_scale_factor = depthPeelingProcess.getDepthScaleFactor();
              Vector3 size_pixel = ( BB.second - BB.first ) / ( ( float ) resolution );
              const double dx = size_pixel[ ( current_orientation+1 ) %3];
              const double dy = size_pixel[ ( current_orientation+2 ) %3];
              const float dS = ( float ) ( dx*dy );
              for ( int i=0;i<ncolls;i++ )
              {
                int index_texture = i*8;
                float first_depth = collisions[index_texture+0];
                Vec2f first_bary ( collisions[index_texture+1], collisions[index_texture+2] );
                int first_index_triangle   = sofa::helper::rnear ( MAXELEMENT * collisions[index_texture+3] );
                int first_index_model; // = (int) ( floor(depthPeelingProcess.getNbModels()*depth_layer[layer_depth][index_texture]+0.5) -1);
                if ( first_bary[0] > 0.5f )
                {
                  first_index_model = 1;
                  first_bary[0] -= 0.625f;
                }
                else
                {
                  first_index_model = 0;
                  first_bary[0] -= 0.125f;
                }
                first_model = models[first_index_model];
                bool first_front = false;
                if ( first_bary[1] > 0.5f )
                {
                  first_front = true;
                  first_bary[1] -= 0.625f;
                }
                else
                {
                  first_front = false;
                  first_bary[1] -= 0.125f;
                }
                first_bary *= 4.0f;
                index_texture += 4;
                float current_depth = collisions[index_texture+0];
                Vec2f bary ( collisions[index_texture+1], collisions[index_texture+2] );
                int current_index_triangle = sofa::helper::rnear ( MAXELEMENT * collisions[index_texture+3] );
                int index_model; // = (int) ( floor(depthPeelingProcess.getNbModels()*depth_layer[layer_depth][index_texture]+0.5) -1);
                if ( bary[0] > 0.5f )
                {
                  index_model = 1;
                  bary[0] -= 0.625f;
                }
                else
                {
                  index_model = 0;
                  bary[0] -= 0.125f;
                }
                current_model = models[index_model];
                bool front = false;
                if ( bary[1] > 0.5f )
                {
                  front = true;
                  bary[1] -= 0.625f;
                }
                else
                {
                  front = false;
                  bary[1] -= 0.125f;
                }
                bary *= 4.0f;

                {
                  // COLLISION
                  if ( ( int ) first_index_triangle   >= first_model->getSize() )   {std::cout << "ERROR in triangle1 indices!!!! : \n";break;}
                  if ( ( int ) current_index_triangle >= current_model->getSize() ) {std::cout << "ERROR in triangle2 indices!!!! : \n";break;}


                  bool pair_inverted = ( first_model != models[0] );
                  detection->id = id_collision++;
                  detection->dx= ( float ) dx; detection->dy= ( float ) dy;
                  if ( !pair_inverted )
                  {
                    //Same coordinates along in the plane normal to the direction of the camera
                    detection->point[0][0] = first_bary[0];detection->point[0][1] = first_bary[1];
                    detection->point[1][0] =       bary[0];detection->point[1][1] =       bary[1];
                    detection->elem = std::pair< core::CollisionElementIterator, core::CollisionElementIterator > ( CollisionElementIterator ( first_model, first_index_triangle ), CollisionElementIterator ( current_model, current_index_triangle ) );
                    detection->dS = dS;
                  }
                  else
                  {
                    //Same coordinates along in the plane normal to the direction of the camera
                    detection->point[0][0] =       bary[0];detection->point[0][1] =       bary[1];
                    detection->point[1][0] = first_bary[0];detection->point[1][1] = first_bary[1];
                    detection->elem = std::pair< core::CollisionElementIterator, core::CollisionElementIterator > ( CollisionElementIterator ( current_model, current_index_triangle ), CollisionElementIterator ( first_model, first_index_triangle ) );
                    detection->dS = -dS;
                  }
                  float K ;
                  if ( autoCollision ) K = Kselfpressure;
                  else                K = Kpressure;

                  detection->normal      = current_orientation;
                  detection->K = K;
                  detection->viscosity = Kviscosity;
                  detection->value = fabs ( ( current_depth - first_depth ) * depth_scale_factor * dS );
                  detection->volume = NULL;
                  ++detection;
                  //if (detection->value < 0) std::cerr << "ERROR: negative volume between colliding layers" << std::endl;
                }
              }

            }
          }
        }
        narrowPhaseTests.push_back ( t );
      }
#endif

      void LDIDetection::NarrowPhaseCPU ( const BoundingBox &BB_visu,const  sofa::helper::vector< TriangleModel*> &models, const sofa::helper::vector<unsigned int > &index_models,bool autoCollision,
                                          LDIDetectionOutputVector *outputs,
                                          simulation::Node* /*timelog*/, core::objectmodel::BaseObject* /*timeobj*/ )
      {
//        simulation::Visitor::ctime_t t0 = 0;
        BoundingBox BB = BB_visu;
        bool globalRendering = models.size() > 2;
        double noCollidingVolume[3]={0,0,0};
        //Keep in memory the last triangles from the at the same coordinate (x,y) but at a different depth
        std::deque< std::pair< int, int > > list_triangle_collision; // list of two cells containing < numero of the layer, the index of the collision model >
        current_orientation--;
        //we will proceed at 3 renders max in the 3 different axes of the AABB used.

        Vec<3,float> unitaryBB;
        Vec<3,float> volumeIdBB;
        for ( int counter=0;counter<3;counter++ )
        {
          current_orientation = ( current_orientation+1 ) %3;
          //given an orientation, we compute the current X and Y

          //****************************************************************************************************************
          unsigned int resolution_used;

          BB = computeCorrectVisualizationBB ( BB_visu, index_models );

          resolution_used= depthPeelingProcess.doDepthPeeling ( models, depth_layer,BB, current_orientation, !globalRendering, adaptiveLength );



          //            if (text)depthPeelingProcess.printDepthPeeling(depth_layer);

//           const float depth_scale_factor = depthPeelingProcess.getDepthScaleFactor();

          //Values to save redondant computationss
          Vector3 size_pixel = ( BB.second - BB.first ) / ( ( float ) resolution_used );
          size_pixel[current_orientation]=BB.second[current_orientation] - BB.first[current_orientation];
          const double dx = size_pixel[X];
          const double dy = size_pixel[Y];
          const float dS = ( float ) ( dx*dy );

          Vec<3,unsigned int> s;

          for (unsigned int i=0;i<3;++i)
            {
                s[i]=resolution_used/subDivisionBB[i];
                if (minSubDivisionBB[i] != 0 && s[i] < minSubDivisionBB[i]) s[i]=minSubDivisionBB[i];
            }

          unitaryBB = size_pixel.linearProduct(s);

          Vector3 centerBB=(BB.second+BB.first)*0.5;
          for (unsigned int i=0;i<3;++i)
            {
              //Nb of independant Volume
              int nbIdPossible = (resolution_used/(float)s[i]);
              if (resolution_used/(float)s[i] != 0) nbIdPossible+=2;
              volumeIdBB[i] = centerBB[i]-(unitaryBB[i]*nbIdPossible*0.5);
            }


          if ( text )     std::cout << "\tUsing dS=" << dS   << "\n";

          //int SC_ray_cast;
          TriangleModel *first_model, *current_model;
          //int index_triangle[2];

          sofa::helper::vector< std::pair< unsigned int,unsigned  int > > &pixel_loc = depthPeelingProcess.pixel_location;


          if ( pixel_loc.size() <= 1 ) continue; //if only one layer, or less, there can't occur collision
          //Exploration of the interval

          for ( unsigned int current_index=pixel_loc[1].first; current_index<pixel_loc[1].second;current_index++ )
          {
            //SC_ray_cast = 0;
            const int index_texture = current_index*4;

            int first_index_model = -1; first_model = NULL;
            Vec2f first_bary;
            unsigned int first_layer_depth = 0;
            bool first_front = false;
            int inside_counter = 0;
            int inside_model=-1;
            int layerPassed=0;
            for ( unsigned int layer_depth=0; layer_depth<pixel_loc.size();layer_depth++ )
            {

              unsigned int index_model=0;
              bool front=false;
              Vec2f bary;

              if (!retrieveInformationFromPixel(layer_depth, index_texture, globalRendering,
                                                front, index_model, bary[0], bary[1]))
                break;

              layerPassed ++;
              if (globalRendering)
                {
                  index_model = ( int ) ( depthPeelingProcess.getNbModels() *index_model ) -1;
                  current_model = depthPeelingProcess.getModel ( ( unsigned int ) index_model );
                }
              else current_model =  models[index_model];

//              bool volumetricDetection=false;

              if ( inside_counter ==0 ) inside_model=-1;

              if ( front )
              {
                if ( inside_counter==0 ) inside_model=first_index_model;
                ++inside_counter;
              }
              else
              {
                --inside_counter;
                bool collision;
                // else
                if ( autoCollision )
                  {
                    collision = first_front && inside_counter > 0;
                    //                     volumetricDetection=collision;
                  }
                else
                {
                  collision = first_front && index_model != (unsigned int)first_index_model;
//                  volumetricDetection = collision;
                  if ( collision && globalRendering )
                  {
                    collision = detectionOutputs.find ( std::make_pair ( first_model, current_model ) ) != detectionOutputs.end() ||
                                detectionOutputs.find ( std::make_pair ( current_model, first_model ) ) != detectionOutputs.end();
                  }

                  if (layerPassed > 1 && inside_counter > 0 && inside_model >=0 && first_front &&
                      index_model == (unsigned int)first_index_model && (unsigned int)inside_model != index_model) //object complety inside another: need to store the volume to perform the good response
                    {
//                       noCollidingVolume[current_orientation] += fabs((depth_layer[layer_depth][index_texture] - depth_layer[first_layer_depth][index_texture]) * depth_scale_factor * dS);
//                      volumetricDetection = true;
                    }

                }

                if ( layer_depth == 0 )
                  {
                    collision = false;
                    continue;
                  }

                if ( collision )
                {
                  // COLLISION

                  unsigned int first_index_triangle   = ( unsigned int ) ceil ( MAXELEMENT * depth_layer[first_layer_depth][index_texture+3] );
                  unsigned int current_index_triangle = ( unsigned int ) ceil ( MAXELEMENT * depth_layer[      layer_depth][index_texture+3] );

                  if ( ( int ) first_index_triangle   >= first_model->getSize() )   {std::cout << "ERROR in triangle1 indices!!!! : \n";break;}
                  if ( ( int ) current_index_triangle >= current_model->getSize() ) {std::cout << "ERROR in triangle2 indices!!!! : \n";break;}

                  {

//                    bool pair_inverted = ( first_model != models[0] );

                    if ( globalRendering )
                    {
                      if ( detectionOutputs.find ( std::make_pair ( first_model, current_model ) ) != detectionOutputs.end() )
                      {
                        DetectionOutputVector * &o=detectionOutputs[std::make_pair ( first_model, current_model ) ];
                        outputs = dynamic_cast< LDIDetectionOutputVector *> ( o );
//                        pair_inverted = false;
                      }
                      else if ( detectionOutputs.find ( std::make_pair ( current_model, first_model ) ) != detectionOutputs.end() )
                      {
                        DetectionOutputVector * &o=detectionOutputs[std::make_pair ( current_model, first_model ) ];
                        outputs = dynamic_cast< LDIDetectionOutputVector *> ( o );
//                        pair_inverted = true;
                      }
                    }
                    if ( outputs != NULL )
                    {
                      Triangle current_t ( current_model, current_index_triangle );
                      Triangle first_t ( first_model, first_index_triangle );
                      Vector3 P[2];
                      P[0] =    current_t.p1()* ( 1-bary[0]-bary[1] ) + current_t.p2()*bary[0] + current_t.p3()*bary[1];
                      P[1] =    first_t.p1()* ( 1-first_bary[0]-first_bary[1] ) + first_t.p2()*first_bary[0] + first_t.p3()*first_bary[1];



                      Vec<3,int> indexInBB;
                      indexInBB[current_orientation] = (P[0][current_orientation]-volumeIdBB[current_orientation])/unitaryBB[current_orientation];
                      indexInBB[X]                   = (P[0][X]-volumeIdBB[X])/unitaryBB[X];
                      indexInBB[Y]                   = (P[0][Y]-volumeIdBB[Y])/unitaryBB[Y];

                      float initDepth = P[0][current_orientation];
                      float endDepth  = P[1][current_orientation];

                      int subdiv=indexInBB[current_orientation];

                      float PenetrationReference=endDepth-initDepth;
                      while (initDepth < endDepth)
                        {
                          float penetrationUsed = 0.0f;
                          float endPortion = volumeIdBB[current_orientation]+(subdiv+1)*unitaryBB[current_orientation];
                          if (endPortion > endDepth) penetrationUsed = endDepth - initDepth;
                          else                       penetrationUsed = endPortion - initDepth;




                          outputs->Vector.resize ( outputs->size() +1 );
                          LDIDetectionOutput *detection = & ( outputs->Vector.back() );
                          detection->value = penetrationUsed;
                          detection->subdivisionFactor = penetrationUsed / PenetrationReference;
                          detection->index = subDivisionBB[0]*subDivisionBB[1]*indexInBB[0] + subDivisionBB[0]*indexInBB[1] + indexInBB[2];

                          detection->id = id_collision++;
                          detection->dx= ( float ) dx; detection->dy= ( float ) dy;


                          detection->point[0][0] = first_bary[0];detection->point[0][1] = first_bary[1];
                          detection->point[1][0] =       bary[0];detection->point[1][1] =       bary[1];    
                          detection->elem = std::pair< core::CollisionElementIterator, core::CollisionElementIterator > ( CollisionElementIterator ( first_model, first_index_triangle ), CollisionElementIterator ( current_model, current_index_triangle ) );
                          detection->dS = dS;

                          float K ;
                          if ( autoCollision ) K = Kselfpressure;
                          else                 K = Kpressure;
                          detection->normal    = current_orientation;
                          detection->K = K;
                          detection->viscosity = Kviscosity;
                          detection->volume=NULL;
//                           if (detection->value < 0) std::cerr << "ERROR: negative volume between colliding layers" << std::endl;

                          subdiv++;
                          initDepth=endPortion;
                        }
                    }
                  }
                }
              }
              first_index_model = index_model; first_model = current_model;
              first_bary = bary;
              first_layer_depth = layer_depth;
              first_front = front;
            }
          }
          if ( exact_volume && outputs->size() != 0 && noCollidingVolume[current_orientation] != 0 )
          {
            std::cout << ( char ) ( 'X' - 0 + current_orientation ) << " : noCollidingVolume = " << noCollidingVolume[current_orientation] << " size : " << outputs->size() << "\n";
            outputs->Vector[outputs->size()-1].volume=&noCollidingVolume[current_orientation];
          }
        }

//        if (!cuttingModels.empty()) voxelDetection->removeVoxels();
      }

      void LDIDetection::EndNarrowPhase ( simulation::Node* /*timelog*/, core::objectmodel::BaseObject* /*timeobj*/ )
      {
//        simulation::Visitor::ctime_t t0 = 0;


        if ( !narrowPhaseTests.empty() )
        {

#ifdef SOFA_GPU_CUDA
          if ( bGPUCollisionVolume )
          {
//    depthPeelingProcess.beginReadCollisionVolume();
            for ( sofa::helper::vector<NarrowPhaseTest>::iterator t=narrowPhaseTests.begin(); t!=narrowPhaseTests.end(); ++t )
            {
              const float* results[3];
              float volume = 0.0f;
              for ( int orientation=0; orientation<3; ++orientation )
              {
                if ( t->result_id[orientation]<0 ) continue;
                results[orientation] = depthPeelingProcess.readCollisionVolume ( // t->result_id[orientation]
                                       );
                if ( results[orientation]==NULL ) continue;
                float vol = results[orientation][0];
                if ( vol > 0 )
                {
                  std::cout << "Results["<<orientation<<"]=";
                  std::cout << " " << vol << std::endl;
                  int i0 = 1;
                  int nbv = static_cast<TriangleModel*> ( t->model1 )->getMechanicalState()->getSize();
                  for ( int v=0;v<nbv;v++ )
                    if ( results[orientation][i0+v] != 0.0f ) std::cout << " "<<v<<"="<<results[orientation][i0+v];
                  std::cout << std::endl;
                  i0 += nbv;
                  nbv = static_cast<TriangleModel*> ( t->model2 )->getMechanicalState()->getSize();
                  for ( int v=0;v<nbv;v++ )
                    if ( results[orientation][i0+v] != 0.0f ) std::cout << " "<<v<<"="<<results[orientation][i0+v];
                  std::cout << std::endl;
                  if ( vol > volume ) volume = vol;
                }
              }

              volume = 0;

              if ( volume != 0 )
              {
                //Create the Vector that will be used by the contact manager to create the collision response
                LDIDetectionOutputVector *outputs=NULL;
                {
                  DetectionOutputVector * &o=detectionOutputs[std::make_pair ( t->model1,t->model2 ) ];
                  if ( o == NULL )
                    o = new LDIDetectionOutputVector();
                  outputs = dynamic_cast< LDIDetectionOutputVector *> ( o );
                }
                outputs->Vector.resize ( outputs->size() +1 );
                LDIDetectionOutput *detection = &* ( outputs->Vector.end()-1 );
                detection->id = 0;
                detection->elem = std::pair< core::CollisionElementIterator, core::CollisionElementIterator > ( CollisionElementIterator ( t->model1, 0 ), CollisionElementIterator ( t->model2, 0 ) );
                detection->dS = 0.0f;
                detection->normal      = -1;
                detection->K = t->K;
                detection->viscosity = t->Kviscosity;
                detection->value = volume;
              }
            }
            depthPeelingProcess.endReadCollisionVolume();
          }
#endif
          narrowPhaseTests.clear();
        }
        ++frame;

        {
          //clean the output vector

          for ( unsigned int i=0;i<all_collisionPair.size();++i )
          {
            TriangleModel *pairCollision[2];
            if ( all_collisionPair[i].size() == 1 )
            {
              pairCollision[0] = pairCollision[1] = all_collisionPair[i][0];
            }
            else
            {
              pairCollision[0] = all_collisionPair[i][0];
              pairCollision[1] = all_collisionPair[i][1];
            }
            if ( detectionOutputs.find ( std::make_pair ( pairCollision[0],pairCollision[1] ) ) != detectionOutputs.end() )
            {
              DetectionOutputVector * &o=detectionOutputs[std::make_pair ( pairCollision[0],pairCollision[1] ) ];

              if ( !o && o->size() ==0 )
              {
                detectionOutputs.erase ( std::make_pair ( pairCollision[0],pairCollision[1] ) );
              }
            }
          }
        }
      }


      //***********************************************************************************************************************************
      void LDIDetection::draw ( bool BBox, bool interactions, bool regions )
      {
				//voxelDetection->drawParsedHexas();
				//voxelDetection->drawRemovedHexas();
				//voxelDetection->drawCollisionTriangles();
				//voxelDetection->clearDebugVectors();

        /*/ Display the position of the removed Hexas
        sofa::helper::vector<Vector3>& detectedVoxels = voxelDetection->removedVoxelsCenters;
        glColor4f ( 1,0,1,1 );
        glPointSize ( 10 );
        glBegin ( GL_POINTS );
        for( sofa::helper::vector<Vector3>::iterator it = detectedVoxels.begin(); it != detectedVoxels.end(); it++)
        {
            helper::gl::glVertexT ( *it);
        }
        glEnd();
        //*/

        if ( !BBox && !interactions && ! regions ) return;
        GLdouble mvMatrix[16];
        glGetDoublev ( GL_MODELVIEW_MATRIX, mvMatrix );

        // Calculate viewpoint
//        double eye[3];
//        eye[0] = ( ( -mvMatrix[12]*mvMatrix[0] ) + ( -mvMatrix[13]*mvMatrix[1] ) +
//                   ( -mvMatrix[14]*mvMatrix[2] ) );
//        eye[1] = ( ( -mvMatrix[12]*mvMatrix[4] ) + ( -mvMatrix[13]*mvMatrix[5] ) +
//                   ( -mvMatrix[14]*mvMatrix[6] ) );
//        eye[2] = ( ( -mvMatrix[12]*mvMatrix[8] ) + ( -mvMatrix[13]*mvMatrix[9] ) +
//                   ( -mvMatrix[14]*mvMatrix[10] ) );

        glDisable ( GL_LIGHTING );

        DetectionOutputMap::const_iterator it;

        if ( BBox )
        {
          glBegin ( GL_LINES );
          for ( unsigned int i=0;i<visualization_BBox.size();i++ )
          {
            glColor4f ( 1,1,1,1 );
            const Vec3d p[8] =
            {
              visualization_BBox[i].first,
              Vec3d ( visualization_BBox[i].second[0],visualization_BBox[i].first[1],visualization_BBox[i].first[2] ),
              Vec3d ( visualization_BBox[i].second[0],visualization_BBox[i].second[1],visualization_BBox[i].first[2] ),
              Vec3d ( visualization_BBox[i].first[0],visualization_BBox[i].second[1],visualization_BBox[i].first[2] ),
              Vec3d ( visualization_BBox[i].first[0],visualization_BBox[i].first[1],visualization_BBox[i].second[2] ),
              Vec3d ( visualization_BBox[i].second[0],visualization_BBox[i].first[1],visualization_BBox[i].second[2] ),
              visualization_BBox[i].second,
              Vec3d ( visualization_BBox[i].first[0],visualization_BBox[i].second[1],visualization_BBox[i].second[2] )
            };
            helper::gl::glVertexT ( p[0] );
            helper::gl::glVertexT ( p[1] );

            helper::gl::glVertexT ( p[1] );
            helper::gl::glVertexT ( p[2] );

            helper::gl::glVertexT ( p[2] );
            helper::gl::glVertexT ( p[3] );

            helper::gl::glVertexT ( p[3] );
            helper::gl::glVertexT ( p[0] );

            helper::gl::glVertexT ( p[4] );
            helper::gl::glVertexT ( p[5] );

            helper::gl::glVertexT ( p[5] );
            helper::gl::glVertexT ( p[6] );

            helper::gl::glVertexT ( p[6] );
            helper::gl::glVertexT ( p[7] );

            helper::gl::glVertexT ( p[7] );
            helper::gl::glVertexT ( p[4] );

            helper::gl::glVertexT ( p[0] );
            helper::gl::glVertexT ( p[4] );

            helper::gl::glVertexT ( p[1] );
            helper::gl::glVertexT ( p[5] );

            helper::gl::glVertexT ( p[2] );
            helper::gl::glVertexT ( p[6] );

            helper::gl::glVertexT ( p[3] );
            helper::gl::glVertexT ( p[7] );

          }
          glEnd();
        }


        if ( interactions )
        {
          // Draw intersection pairs
          glBegin ( GL_LINES );
          for ( it=detectionOutputs.begin(); it!= detectionOutputs.end(); it++ )
          {
            if ( ! ( dynamic_cast<TriangleModel*> ( ( *it ).first.first ) && dynamic_cast<TriangleModel*> ( ( *it ).first.second ) ) ) continue; //handle mouse interaction : not drawn here
            LDIDetectionOutputVector* detected = dynamic_cast<LDIDetectionOutputVector*> ( ( *it ).second );
            for ( unsigned int i=0;i<detected->size();i++ )
            {
              LDIDetectionOutput* o=&detected->Vector[i];

              Triangle t[2] = { static_cast<Triangle> ( o->elem.first ), static_cast<Triangle> ( o->elem.second ) };

              Vec3d p[2]    = { t[0].p1() * ( 1-o->point[0][0]-o->point[0][1] ) + t[0].p2() *o->point[0][0] + t[0].p3() *o->point[0][1],
                                t[1].p1() * ( 1-o->point[1][0]-o->point[1][1] ) + t[1].p2() *o->point[1][0] + t[1].p3() *o->point[1][1]
                              };

              glColor4f ( ( float ) ( o->normal==0 ), ( float ) ( o->normal==1 ), ( float ) ( o->normal==2 ),0.2f );
              helper::gl::glVertexT ( p[0] );
              helper::gl::glVertexT ( p[1] );
            }
          }
          glEnd();
        }


        if ( regions )
        {
          //Sort the quads to draw
          std::multimap< double, std::pair< Vec3d, LDIDetectionOutput*> > sorted_quad;
          for ( it=detectionOutputs.begin(); it!= detectionOutputs.end(); it++ )
          {
            if ( ! ( dynamic_cast<TriangleModel*> ( ( *it ).first.first ) && dynamic_cast<TriangleModel*> ( ( *it ).first.second ) ) ) continue; //handle mouse interaction : not drawn here
            LDIDetectionOutputVector* detected = dynamic_cast<LDIDetectionOutputVector*> ( ( *it ).second );
            for ( unsigned int i=0;i<detected->size();i++ )
            {
              LDIDetectionOutput* o=&detected->Vector[i];

              Triangle t[2] = { static_cast<Triangle> ( o->elem.first ), static_cast<Triangle> ( o->elem.second ) };
              Vec3d p[2]    = { t[0].p1() * ( 1-o->point[0][0]-o->point[0][1] ) + t[0].p2() *o->point[0][0] + t[0].p3() *o->point[0][1],
                                t[1].p1() * ( 1-o->point[1][0]-o->point[1][1] ) + t[1].p2() *o->point[1][0] + t[1].p3() *o->point[1][1]
                              };
              sorted_quad.insert ( std::make_pair ( p[0][o->normal], std::make_pair ( p[0], o ) ) );
              sorted_quad.insert ( std::make_pair ( p[1][o->normal], std::make_pair ( p[1], o ) ) );
            }
          }


          glDisable ( GL_DEPTH_TEST );
          glEnable ( GL_BLEND );
          glBlendFunc ( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
          glBegin ( GL_QUADS );
          // Draw intersection pairs

          for ( std::multimap<double, std::pair< Vec3d, LDIDetectionOutput*> >::const_iterator it=sorted_quad.begin(); it != sorted_quad.end();it++ )
          {
            int normal = ( *it ).second.second->normal;
            float quad_x = ( *it ).second.second->dx/2.0f;
            float quad_y = ( *it ).second.second->dy/2.0f;
            Vec3d p= ( *it ).second.first;
            glColor4f ( normal==0,normal==1,normal==2,0.25 );

            Vec3d axe=Vec3d();
            axe[ ( normal+1 ) %3] = -quad_x;axe[ ( normal+2 ) %3] = -quad_y;helper::gl::glVertexT ( p+axe );
            axe[ ( normal+1 ) %3] = -quad_x;axe[ ( normal+2 ) %3] =  quad_y;helper::gl::glVertexT ( p+axe );
            axe[ ( normal+1 ) %3] =  quad_x;axe[ ( normal+2 ) %3] =  quad_y;helper::gl::glVertexT ( p+axe );
            axe[ ( normal+1 ) %3] =  quad_x;axe[ ( normal+2 ) %3] = -quad_y;helper::gl::glVertexT ( p+axe );
          }

          glEnd();
          glDisable ( GL_BLEND );
          glEnable ( GL_DEPTH_TEST );
        }


      }

      bool LDIDetection::retrieveInformationFromPixel( const unsigned int &layerDepth, const unsigned int &index, const bool &globalRendering,
                                                       bool &front, unsigned int &indexModel,
                                                       float &alpha, float &beta) const
      {
        Vec2f bary ( depth_layer[layerDepth][index+1], depth_layer[layerDepth][index+2] );
        if ( bary[0] == 0 ) return false;
        if ( bary[0] > 0.5f )
          {
            if ( !globalRendering ) indexModel = 1;
            bary[0] -= 0.625f;
          }
        else
          {
            if ( !globalRendering ) indexModel = 0;
            bary[0] -= 0.125f;
          }

        if ( globalRendering ) indexModel = depth_layer[layerDepth][index];


        front = false;
        if ( bary[1] > 0.5f )
          {
            front = true;
            bary[1] -= 0.625f;
          }
        else
          {
            front = false;
            bary[1] -= 0.125f;
          }
        bary *= 4.0f;
        alpha=bary[0]; beta=bary[1];
        return true;
      }
      //***********************************************************************************************************************************
      //Detection from the information given by the normals
      bool LDIDetection::inCollision ( const bool same_object, const Vec3d &normal0,const Vec3d &normal1, int &SC_ray_cast ) const
      {

        const float sign[2] = { ( float )-normal0[current_orientation],
                                ( float )-normal1[current_orientation]
                              };

        if ( same_object )
        {
          //self collision : before creating a self collision response, we need to find two entering sides of the triangles (normals oriented in an opposite direction)
          if ( SC_ray_cast == 0 )
          {
            if ( sign[0]< 0 && sign[1] < 0 ) SC_ray_cast ++;
            return false;
          }
          SC_ray_cast = 0;
        }
        return ( sign[0]< 0 && sign[1] > 0 );
      }




      //***********************************************************************************************************************************
      //Given a triangle model, we construct its hierarchical bounding boxes:
      //   the first element of the vector will be the bounding boxes at the given level
      //   the second and last element will be the whole bounding box of the triangle model
      void LDIDetection::constructBoundingHierarchy ( TriangleModel *m, sofa::helper::vector< sofa::helper::vector< BoundingBox > > &boundingHierarchy ) const
      {

        CollisionModel *model;
        sofa::helper::vector< BoundingBox > bounding;

        //Calculate the hierarchical bounding boxes of the triangle model
        //m->computeBoundingTree(depth_BB);
        model = m->getPrevious();
        if ( model->getPrevious() ) model = model->getPrevious(); // skip the first layer of bounding cubes


        CubeModel *last_model=NULL;
        int layer = 0;
        //We explorate the hierarchy of collision models, and calculate at each depth the bounding boxes
        //We want to keep only the lowest and highest level of the hierarchy
        while ( model != NULL )
        {
          if ( CubeModel *BB = dynamic_cast< CubeModel *> ( model ) )
          {
            last_model = BB;
            if ( BB->getNumberCells() != 0 && ( boundingHierarchy.empty() || model->getPrevious() ==NULL ) ) // first or last
            {
              BB->getBoundingTree ( bounding );
              boundingHierarchy.push_back ( bounding );
            }
            ++layer;
          }
          model = model->getPrevious();
        }
        if ( last_model == NULL ) return;

        //last_model->getBoundingTree( bounding );
        //boundingHierarchy.push_back( bounding );
      }


      //***********************************************************************************************************************************
      //Given two bounding boxes, we determine if they intersect, and if isCrossing is true, result contains the resulting AABB
      void LDIDetection::findIntersectionBoundingBox ( const BoundingBox &m0BB, const BoundingBox  &m1BB,
              bool &isCrossing, BoundingBox &result ) const
      {
        //************************************************************************************************
        //First test with the whole bounding box:
        crossingAABB ( m0BB, m1BB , isCrossing, result );
//         if (isCrossing)
//           {
//             if (m0BB.first[0]  < result.first [0]) result.first[0]  = m0BB.first[0];
//             if (m0BB.first[1]  < result.first [1]) result.first[1]  = m0BB.first[1];
//             if (m0BB.first[2]  < result.first [2]) result.first[2]  = m0BB.first[2];
//             if (m0BB.second[0] > result.second[0]) result.second[0] = m0BB.second[0];
//             if (m0BB.second[1] > result.second[1]) result.second[1] = m0BB.second[1];
//             if (m0BB.second[2] > result.second[2]) result.second[2] = m0BB.second[2];
//           }
       }


      //***********************************************************************************************************************************
      //Make the intersection of the Bounding Box. isCrossing gives the information about the presence of this intersection
      //resultingBBox is relevant if only isCrossing is true.
      void LDIDetection::crossingAABB ( const  BoundingBox &b0, const BoundingBox &b1, bool &isCrossing, BoundingBox &resultingBBox ) const
      {
        //If it exists an intersection between the two bounding box, we calcule the resulting BB
        if ( b0.first[0]<b1.second[0] && b0.first[1]<b1.second[1] && b0.first[2]<b1.second[2] &&
             b1.first[0]<b0.second[0] && b1.first[1]<b0.second[1] && b1.first[2]<b0.second[2] )
        {
          resultingBBox.first [0] = ( b0.first[0]  < b1.first[0] ) ?b1.first[0]:b0.first[0];
          resultingBBox.first [1] = ( b0.first[1]  < b1.first[1] ) ?b1.first[1]:b0.first[1];
          resultingBBox.first [2] = ( b0.first[2]  < b1.first[2] ) ?b1.first[2]:b0.first[2];
          resultingBBox.second[0] = ( b0.second[0] > b1.second[0] ) ?b1.second[0]:b0.second[0];
          resultingBBox.second[1] = ( b0.second[1] > b1.second[1] ) ?b1.second[1]:b0.second[1];
          resultingBBox.second[2] = ( b0.second[2] > b1.second[2] ) ?b1.second[2]:b0.second[2];
          isCrossing = ( resultingBBox.first[0]<resultingBBox.second[0] &&
                         resultingBBox.first[1]<resultingBBox.second[1] &&
                         resultingBBox.first[2]<resultingBBox.second[2] );
        }
        else
        {
          isCrossing = false;
        }
      }
      //***********************************************************************************************************************************
      //Make the intersection of the Bounding Box. isCrossing gives the information about the presence of this intersection
      //resultingBBox is relevant if only isCrossing is true.
      void LDIDetection::unionAABB ( const  BoundingBox &b0, const BoundingBox &b1, BoundingBox &resultingBBox ) const
      {
        //If it exists an intersection between the two bounding box, we calcule the resulting BB

        resultingBBox.first [0] = ( b0.first[0]  < b1.first[0] ) ?b0.first[0]:b1.first[0];
        resultingBBox.first [1] = ( b0.first[1]  < b1.first[1] ) ?b0.first[1]:b1.first[1];
        resultingBBox.first [2] = ( b0.first[2]  < b1.first[2] ) ?b0.first[2]:b1.first[2];
        resultingBBox.second[0] = ( b0.second[0] > b1.second[0] ) ?b0.second[0]:b1.second[0];
        resultingBBox.second[1] = ( b0.second[1] > b1.second[1] ) ?b0.second[1]:b1.second[1];
        resultingBBox.second[2] = ( b0.second[2] > b1.second[2] ) ?b0.second[2]:b1.second[2];

      }

			

			

      BoundingBox LDIDetection::computeCorrectVisualizationBB ( const BoundingBox &bb, const sofa::helper::vector< unsigned int > //index_models
                                                              ) const
      {
        BoundingBox correctBB=bb;
        const Vec3d s = bb.second - bb.first;
//         if (exact_volume)
//         {
//           correctBB.first[current_orientation]  = std::min(all_collisionBB[index_models[0]].first[current_orientation],
//               all_collisionBB[index_models[1]].first[current_orientation])-0.1*s[current_orientation];

//           correctBB.second[current_orientation] = std::max(all_collisionBB[index_models[0]].second[current_orientation],
//               all_collisionBB[index_models[1]].second[current_orientation])+0.1*s[current_orientation];
//         }
        if ( s[X] == 0 || s[Y] == 0 ) return correctBB; //the pixel will not be set as squared: degenerated case
        if ( s[X] > s[Y] )
        {
          correctBB.second[Y] = correctBB.first[Y] + s[X];
        }
        else
        {
          correctBB.second[X] = correctBB.first[X] + s[Y];
        }
        return correctBB;
      }

      void LDIDetection::clear()
      {
        DetectionOutputMap::iterator it;

        for ( it = detectionOutputs.begin(); it != detectionOutputs.end();it++ )
        {
          if ( dynamic_cast< RayModel* > ( ( *it ).first.first )  || dynamic_cast< RayModel* > ( ( *it ).first.second ) ) continue;
//          delete ( *it ).second;
          ( *it ).second->release();
        }

        detectionOutputs.clear();
        id_collision = 0;
        index_used.clear();
        all_collisionModels.clear();
        all_collisionBB.clear();
        depthPeelingProcess.nextFrame();
      }

      void LDIDetection::draw()
      {

      }


    }
  }
}
