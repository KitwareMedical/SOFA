//
// C++ Interface: DepthPeelingUtility
//
// Description: Implements the Depth Peeling algorithm 
//
//
// Author: Francois Faure, Sebastien Barbier, Jeremie Allard, Florent Falipou
//
// Licence: QPL, See LICENCE.txt file that comes with this distribution
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef _DEPTHPEELINGUTILITY_H_
#define _DEPTHPEELINGUTILITY_H_

#include "DepthPeeling.h"
#include <sofa/helper/vector.h>
#include <sofa/component/collision/TriangleModel.h>
#include <sofa/component/topology/RegularGridTopology.h>
#include <sofa/simulation/common/Node.h>
#ifdef SOFA_GPU_CUDA
#include <sofa/gpu/cuda/CudaTypes.h>
//#include <cudpp.h>
#endif

namespace sofa
{

  namespace component
  {

    namespace collision
    {

//number max of triangles : 2^20
      
#define MAXELEMENT 2097152.0f
//#define MAXELEMENT 1048576.0f
#define MAXLAYER   30


//#define MAXELEMENT 4096.0f


      using namespace sofa::component::collision;
      
      
      typedef std::pair< Vector3, Vector3 > BoundingBox;
      //*******************************************************************************************
      //Class to easily perform a depth peeling
      class SOFA_LDIDETECTION_API DepthPeelingUtility
      {
	public:

    typedef core::behavior::MechanicalState<Vec3Types>::VecCoord VecCoord;

	  DepthPeelingUtility();
	  ~DepthPeelingUtility()
	  {
	    if (!initialized)return;
	    if (PBO.size() > 0) glDeleteBuffersARB(2, &PBO[0]);
	    //if (VBO.size() >0)  glDeleteBuffersARB(VBO.size(), &(VBO[0]));
	    
	    
	    glDeleteTextures(3, shadowMap);	    
	    glDeleteTextures(2, idsMap);	    
	    glDeleteFramebuffersEXT(2, FBO);
	  }
	  
	  
	  /// Initialization: Init GLEW, and reserve memory space for the different textures
	  void init(unsigned int _resolution, float _resolutionPixel, bool initCUDA);	  
	  
	  void setCUDA(bool b){useCUDA=b;}
	  void         setNbModels(int n){nbModels = n;};
	  unsigned int getNbModels() const{return nbModels;}
	  TriangleModel* getModel(unsigned int index);
	  
	  TriangleModel* getModel(float index);
	  
	  /// Load one model in a VBO
	  void loadModels( const sofa::helper::vector<TriangleModel *> &model);
	  
	  /// Load models in a VBO
	  void loadModel( TriangleModel * &model);
	  
	  /// Do a complete Depth Peeling of ONE object (load in VBO, and get the layers)
	  unsigned int doDepthPeeling( TriangleModel *  _models, 
                  	       sofa::helper::vector< sofa::helper::vector< float> > &layer, 
			  const unsigned char orientation=0, bool usePairShader = true, double adaptiveLength=0.0);	  
	  
	  /// Do a Depth Peeling of a scene containing objects that MUST have been already loaded! The Bounding Box will be automatically calculated  
	  unsigned int doDepthPeeling( const sofa::helper::vector< TriangleModel *>  &_models, 
			       sofa::helper::vector< sofa::helper::vector< float> > &layer, 
	  const unsigned char orientation=0, bool usePairShader = false, double adaptiveLength=0.0);	  
	      
	  /// Do a Depth Peeling of a scene containing objects that MUST have been already loaded given the region of interest
	  unsigned int doDepthPeeling( const sofa::helper::vector< TriangleModel *>  &_models, 
			        sofa::helper::vector< sofa::helper::vector< float> > &layer,
	                       const BoundingBox &bb,
	                       const unsigned char orientation=0, bool usePairShader = false, double adaptiveLength=0.0 );
	  
	  
	  /// Generate a texture 3d given the layers generated by a depth peeling
	  void doVoxelize(const sofa::helper::vector< sofa::helper::vector< float> > &layer,
			  sofa::helper::vector< sofa::helper::vector< sofa::helper::vector< unsigned char > > > &texture3d) ;
	 
	  /// Generate a texture 3d given a list of TriangleModels
	  void doVoxelize(const sofa::helper::vector< TriangleModel *>  &_models, 
			  sofa::helper::vector< sofa::helper::vector< sofa::helper::vector< unsigned char > > > &texture3d,
                          const unsigned char orientation=0 )  ;
	  
	  /// Generate a texture 3d given ONE TriangleModel
	  void doVoxelize( TriangleModel *  _models, 
			  sofa::helper::vector< sofa::helper::vector< sofa::helper::vector< unsigned char > > > &texture3d,
                          const unsigned char orientation=0 )  ;
	  
	  /// Compute the bounding box of a triangle model
	 static BoundingBox constructBoundingBox(  TriangleModel *m) ;
	 		 
	 //DEBUG Functions
	 /// ASCII display of the layers
	 void printDepthPeeling(const sofa::helper::vector< sofa::helper::vector< float> > &layer);
	 
	 /// ASCII display of the texture 3d generated
	 static void printVoxelization(const sofa::helper::vector< sofa::helper::vector< sofa::helper::vector< unsigned char > > > &texture3d);
	 
	 //for each level of the depth peeling, store the index of the first pixel found, and the last pixel : when they are equal, depth peeling is over for this level
	 sofa::helper::vector< std::pair< unsigned int,unsigned  int > > pixel_location;
	 
      float getDepthScaleFactor() const { return (float)(sceneMaxBBox[current_orientation] - sceneMinBBox[current_orientation]); }	  

          void doCollisionVolume( const sofa::helper::vector< TriangleModel *>  &_models, 
                                  int &result_id,
                                  const BoundingBox &bb,
                                  const unsigned char orientation=0, simulation::Node* timelog = NULL, core::objectmodel::BaseObject* timeobj = NULL );

#ifdef SOFA_GPU_CUDA
          void beginReadCollisionVolume();
          const float* readCollisionVolume();
          void endReadCollisionVolume();  

#endif
	protected:
	   	  	  
	  void setResolutionTexture( const BoundingBox &bb);
	      
	     // -- Function only used to retrieve information about the coordinates and number of triangle of the current Triangle Model: it is only used for CPU->GPU transfert( creation of a VBO )
	  void fillArrays( TriangleModel *m);
 
	  /// Allocation of GPU memory
          void initFBO();

#ifdef SOFA_GPU_CUDA
	  void initFBOVolume();	  
#endif
	  /// Render the objects specified
	  void renderScene (  const  sofa::helper::vector< TriangleModel *>  &_models, bool usePairShader = false ) const;
	  
	  bool computeInterval( const unsigned int index_layer,
				const  sofa::helper::vector< sofa::helper::vector< float> > &layer);
	  
	  
	  void initDepthPeeling(const BoundingBox &visualizationBB);
	  void initRenderView( const BoundingBox &boundingBox );
	  void endRenderView();
	  void setDepthBufferTest() const;
	  void resetDepthBufferTest() const;

	  inline void createVBO( TriangleModel *m, int counter_model, unsigned int indexVBO);
	  
	  /// Information about the orientation of the triangle : is it a front face?
	  inline bool frontFace( Triangle &t) const;
	  	  
	  void fillTexture3d(unsigned int x, unsigned int y, 
			     const sofa::helper::vector< std::pair < std::pair< int, unsigned int> , bool > > &ray,
	                     sofa::helper::vector< sofa::helper::vector< sofa::helper::vector< unsigned char > > > &texture3d) const ;
 
	  float getDepth(const float i, const float alpha, const float beta);
	  float getDepth(const float z);
	  bool isGridTopology(TriangleModel* m) const;
	      
	  bool initialized;
	  DepthPeelingShader depthPeeling_shader;
	  DepthPeelingShaderRenderPair depthPeeling_shader_pair;
 
          PeelingShader* shader;
	  
    std::map<  TriangleModel *, unsigned int> modelsId;
    struct ModelInfo
    {
        TriangleModel* model;
        GLuint vbo; ///< id of the VBO of this model
        int lastUpdateTS; ///< last time this model was updated
        int lastNbTriangles; ///< number of triangles at last update
        float* uploadBuffer;
        std::string grid; ///< if not empty: name of the 2x2x2 grid-mapping
        component::topology::GridTopology* gridTopology;
        core::behavior::MechanicalState< Vec3Types >* gridState;
        ModelInfo()
        : model(NULL), vbo(0), lastUpdateTS(-1), lastNbTriangles(0), uploadBuffer(NULL), gridTopology(NULL), gridState(NULL)
        {
        }
    };
    
    sofa::helper::vector<ModelInfo> modelsInfo;
    
    int currentTS; ///< Current time-stamp
public:
    void nextFrame() { ++currentTS; }
protected:
	  //std::map<  TriangleModel *, unsigned int> modelsIdVBO;
	  std::map<  std::string , unsigned int>    modelsGridVBO;	
	  //std::map< TriangleModel*, unsigned int >  staticObject;
	  //std::map< TriangleModel*, unsigned int >  gridObject;	  
	  
	  //std::map< TriangleModel*, unsigned int >  staticObjectIdVBO;
	  //std::map< TriangleModel*, unsigned int >  gridObjectIdVBO;	  
	  
	  unsigned int nbModels;
	  unsigned int nbMaxTriangles;
	  float VOI[3]; // size of the volume of interest
	  //Parameters
	  unsigned int resolution;	 
	  float        resolutionPixel;	 
	  unsigned int resolutionMax;	  
	  unsigned int number_layer ;               //current number of layers present computated
	
	  
	  
	  //OpenGL elements
	  float *init_shadowmap; //used to clear the shadow map
	  float *init_idsmap;    //used to clear the identity map
	  sofa::helper::vector< float >  array_identity;
          GLuint idVBOIdentity;
	  sofa::helper::vector< float >  array_coord   ;
	  GLuint                idsMap[2];          //id of the texture that will contain information about the identity of the triangles
	  GLuint                shadowMap[3];       //id of the texture that will contain information about the depth of the last drawn depth peel
	  GLuint                FBO[2];             //id of the FBOs
	  sofa::helper::vector< GLuint > PBO   ;             //id of the PBO
	  //sofa::helper::vector< GLuint > VBO;                //ids of all the VBO containing the collision models	 
	  GLdouble mvMatrix[16];  //model view matrix for adaptive resolution
	  
          // Geometry-shader based collision volume computation
          GLuint                idLTexture;                         ///<id of the texture that will contain LDI layers
          GLuint                VBOpixels;                          ///<VBO with one dummy position per pixel in LDI
          GLuint                FBOVolume;                          ///<id of the FBO used to render volume contributions
          sofa::helper::vector<GLuint>  idsVolumeResultTexture;     ///<id of the textures that will contain collision volumes
          sofa::helper::vector<GLuint>  idsVolumeResultPBO;         ///<id of the PBOs that will contain collision volumes
          sofa::helper::vector< GLuint > IBO;                       ///<ids of all the VBO containing the triangle indices of collision models
          sofa::helper::vector< GLuint > idsIndexBuffer;            ///<ids of all the buffer textures containing the triangle indices of collision models
          GLuint                geometryQuery;                      ///<query used to know how many points have been generated by the geometry shaders
          //GLuint                geometryBuffer;                     ///<buffer object for geometry transform feedback
          int                   prevVolumeResultSize;               ///<previous size of the results computed in the volume textures
          int                   volumeResultSize;                   ///<size of the results computed in the volume textures
          int volumeResultCurrent;                                  ///<index of the currently mapped volume result
          const float* volumeResultMapping;                         ///<current mapping of volume result textures
          enum { VOLUME_TEXTURE_RESOLUTION = 1024, VOLUME_TEXTURE_SIZE=VOLUME_TEXTURE_RESOLUTION*VOLUME_TEXTURE_RESOLUTION };

#ifdef SOFA_GPU_CUDA
          // CUDA based collision volume computation
          // Note that idLTexture, FBOVolume and IBO are reused from above
          //GLuint                idLTexture;                         ///<id of the texture that will contain LDI layers
          //GLuint                FBOVolume;                          ///<id of the FBO used to render volume contributions
          GLuint                  idLPBO;                           ///<id of the PBO used to copy layers from OpenGL to CUDA
          //sofa::helper::vector< GLuint > IBO;                       ///<ids of all the VBO containing the triangle indices of collision models
          sofa::gpu::cuda::CudaVector<unsigned int> scan_in, scan_out;
          sofa::gpu::cuda::CudaVector<unsigned int> partial_sums;
          sofa::gpu::cuda::CudaVector<float> output_collisions;
          //CUDPPScanConfig scan_config;
#endif

          GLuint occlusionQuery[2]; ///<id of occlusion query, or 0 if not used/available

	  unsigned char ping_pong;
	  unsigned char current_orientation;       //orientation for the depth peeling
	  
	  bool useCUDA;

	  Vector3 sceneMinBBox;
          Vector3 sceneMaxBBox;
      };

	  struct DepthCmp {
		 bool operator()( const float s1, const float s2 ) const {
			return s1 <= s2;
		}
  };
            
    } //collision
  } //component
} //sofa


#endif
