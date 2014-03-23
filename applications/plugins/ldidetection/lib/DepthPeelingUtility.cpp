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

#include "DepthPeelingUtility.h"
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/component/collision/CubeModel.h>
#include <sofa/component/collision/TriangleModel.h>
#include <sofa/core/VecId.h>
#include <sofa/simulation/common/Node.h>
#include <sofa/component/topology/RegularGridTopology.h>
#include <sofa/component/topology/MeshTopology.h>
#include <sofa/simulation/common/Simulation.h>
#include <sofa/helper/accessor.h>
#ifdef SOFA_GPU_CUDA
#include <sofa/gpu/cuda/mycuda.h>
#endif

#include <string.h>

#include <algorithm>
using std::cerr;
using std::endl;

namespace sofa
{

  namespace component
  {

    namespace collision
    {

      // Use a single FBO for ping-pong instead of 2
#define USE_SINGLE_FBO

      // Use occlusion queries
#define USE_OCCLUSION_QUERY



      //#define CHECK_GL_ERRORS

#ifdef CHECK_GL_ERRORS

#define GLCHECK()							\
      {									\
	GLenum err;							\
	while ((err = glGetError()) != GL_NO_ERROR)			\
	  {								\
	    std::cerr << __FILE__ << '(' << __LINE__ << "): ERROR GL"<<err<<": "; \
	    switch(err)							\
	      {								\
	      case GL_INVALID_ENUM: std::cerr << "GL_INVALID_ENUM"; break; \
	      case GL_INVALID_VALUE: std::cerr << "GL_INVALID_VALUE"; break; \
	      case GL_INVALID_OPERATION: std::cerr << "GL_INVALID_OPERATION"; break; \
	      case GL_STACK_OVERFLOW: std::cerr << "GL_STACK_OVERFLOW"; break; \
	      case GL_STACK_UNDERFLOW: std::cerr << "GL_STACK_UNDERFLOW"; break; \
	      case GL_OUT_OF_MEMORY: std::cerr << "GL_OUT_OF_MEMORY"; break; \
	      default: std::cerr << "unknown error 0x"<<std::hex<<err<<std::dec; break; \
	      }								\
	    std::cerr << std::endl;					\
	  }								\
      }

#else
#define GLCHECK() {}
#endif

      DepthPeelingUtility::DepthPeelingUtility(): initialized(false), nbModels(0), nbMaxTriangles(0), resolution(64), number_layer(0), idVBOIdentity(0)
      {
	resolutionPixel = resolution;
	useCUDA=false;
      }

      //***********************************************************************************************************************************
      void DepthPeelingUtility::init(unsigned int _resolution, float _resolutionPixel, bool bCUDA)
      {

	useCUDA = bCUDA;
	resolutionPixel = _resolutionPixel;
	if (initialized && _resolution == resolutionMax) return; //no need to recreate identical things
	resolutionMax = _resolution;
	resolution = resolutionMax;
#ifdef SOFA_HAVE_GLEW

	GLenum err = glewInit();
	if (GLEW_OK != err)
	  {
	    /* Problem: glewInit failed, something is seriously wrong. */
	    fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
	  }
#endif

	init_shadowmap = new float[resolutionMax*resolutionMax];
	memset(init_shadowmap, 1, resolutionMax*resolutionMax*sizeof(float));


	depthPeeling_shader.load(); // = DepthPeelingShader(true);
	depthPeeling_shader_pair.load();


	shader = &depthPeeling_shader_pair;

	//VBO.reserve(200);

	//Creation of the PBOs
	if (PBO.size() != 0)
	  glDeleteBuffersARB(PBO.size(), &PBO[0]);

	PBO.resize(2);
	glGenBuffersARB(2, &PBO[0]);
	for (unsigned int i=0;i<2;i++)
	  {
	    glBindBufferARB(GL_PIXEL_PACK_BUFFER_ARB, PBO[i]);
	    glBufferDataARB(GL_PIXEL_PACK_BUFFER_ARB,  resolutionMax*resolutionMax*4*sizeof(float),NULL, GL_STREAM_READ_ARB);
	    glBindBufferARB(GL_PIXEL_PACK_BUFFER_ARB,0);
	  }

	glGenTextures( 3, shadowMap );
	glGenTextures( 2, idsMap);
	glGenFramebuffersEXT(2, FBO);
	initFBO();
	
#ifdef SOFA_GPU_CUDA
	if (useCUDA)	  initFBOVolume();
#endif

#ifdef SOFA_HAVE_GLEW
	if (GLEW_ARB_color_buffer_float)
	  {
	    glClampColorARB(GL_CLAMP_VERTEX_COLOR_ARB, GL_FALSE);
	    glClampColorARB(GL_CLAMP_FRAGMENT_COLOR_ARB, GL_FALSE);
	    glClampColorARB(GL_CLAMP_READ_COLOR_ARB, GL_FALSE);
	  }
#endif

#ifdef USE_OCCLUSION_QUERY
#ifdef SOFA_HAVE_GLEW
	if (GLEW_ARB_occlusion_query)
	  {
	    std::cout << "DepthPeeling: using occlusion query" << std::endl;
	    if (!occlusionQuery[0])
	      glGenQueriesARB(2,occlusionQuery);
	  }
	else
#endif
#endif
	  occlusionQuery[0] = 0;
	initialized = true;
	GLCHECK();
      }



      //***********************************************************************************************************************************
      void DepthPeelingUtility::setResolutionTexture( const BoundingBox &bb)
      {
	if (resolutionPixel <= 0) return; //fixed resolution wanted by the used
	SReal distance = bb.second[current_orientation]-bb.first[current_orientation];
	resolution = (unsigned int) (distance/(resolutionPixel/1000.0f));
	resolution += resolution%2;

	if (resolution > resolutionMax)  resolution = resolutionMax;
      }



      //***********************************************************************************************************************************
      void DepthPeelingUtility::loadModels( const sofa::helper::vector<TriangleModel *> &_models)
      {

	//const bool text = false;
	//one VBO for the coordinates of the points, and another to keep the information of the depth, index of the collision model and numero of triangle

	if (nbModels < _models.size()) nbModels = _models.size();

	std::vector<int> modelsToUpdate;
	bool updateVBOIdentity = false;
	//unsigned int nbVBO = 0;
	std::set< std::string> names;
	for (unsigned int i=0;i<_models.size();i++)
	  {
	    assert(_models[i]->getSize() < MAXELEMENT);

	    if (!modelsId.count(_models[i]))
	      {
		modelsId[_models[i]] = modelsInfo.size();
		ModelInfo info;
		info.model = _models[i];
		if (isGridTopology(_models[i]))
		  {
		    if (TriangleModel *meshmodel = dynamic_cast< TriangleModel* >(_models[i]) )
		      {
			if (static_cast< sofa::component::topology::MeshTopology *>(meshmodel->getTopology())->getFilename() != "")
			  info.grid =static_cast< sofa::component::topology::MeshTopology *>(meshmodel->getTopology())->getFilename();
		      }
		  }
		modelsInfo.push_back(info);
		if (_models[i]->getSize() > (int)nbMaxTriangles)
		  {
		    nbMaxTriangles = _models[i]->getSize();
		    updateVBOIdentity = true;
		  }
	      }

	    int id = modelsId[_models[i]];
	    ModelInfo& info = modelsInfo[id];
	    if (info.vbo == 0 || (info.model->isMoving() && info.lastUpdateTS != currentTS))
	      modelsToUpdate.push_back(id);
	  }
	if (updateVBOIdentity)
	  {
	    if (!idVBOIdentity)
	      glGenBuffers(1,&idVBOIdentity);
	    if (array_identity.size() < 4*3*nbMaxTriangles)
	      {
		//each triangle has 3 vertices and we store on RGBA
		const unsigned int initial_size = array_identity.size()/12;
		array_identity.resize(4*3*nbMaxTriangles);
		const float b0 = 0.125f;
		const float b1 = 0.125f+0.25f;
		for (unsigned int i=initial_size;i<nbMaxTriangles;i++)
		  {
		    const int index = 3*i;
		    //Load the barycentric coordinates for each triangle
		    const float idTriangle = i/MAXELEMENT;
		    array_identity[(index)*4+0] = 0.0f;
		    array_identity[(index)*4+1] = b0;
		    array_identity[(index)*4+2] = b0;
		    array_identity[(index)*4+3] = idTriangle;

		    array_identity[(index+1)*4+0] = 0.0f;
		    array_identity[(index+1)*4+1] = b1;
		    array_identity[(index+1)*4+2] = b0;
		    array_identity[(index+1)*4+3] = idTriangle;

		    array_identity[(index+2)*4+0] = 0.0f;
		    array_identity[(index+2)*4+1] = b0;
		    array_identity[(index+2)*4+2] = b1;
		    array_identity[(index+2)*4+3] = idTriangle;
		  }
	      }
	    std::cout << "Uploading identity array with "<<nbMaxTriangles<<" triangles"<<std::endl;
	    glBindBufferARB(GL_ARRAY_BUFFER_ARB, idVBOIdentity);
	    glBufferDataARB(GL_ARRAY_BUFFER_ARB, nbMaxTriangles*3*4*sizeof(float), &(array_identity[0]), GL_STATIC_DRAW_ARB);
	    glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
	  }
	if (modelsToUpdate.empty()) return;
	for (unsigned int i=0;i<modelsToUpdate.size();++i)
	  {
	    int id = modelsToUpdate[i];
	    ModelInfo& info = modelsInfo[id];
	    if (!info.vbo)
	      glGenBuffers(1,&info.vbo);
	    int nbt = info.model->getSize();
	    glBindBufferARB(GL_ARRAY_BUFFER_ARB, info.vbo);
	    if (nbt != info.lastNbTriangles)
	      {
		std::cout << "(Re)Creating VBO for "<<info.model->getName()<<" with "<<nbt<<" triangles"<<std::endl;
#ifdef SOFA_GPU_CUDA
		if (useCUDA)
		  {
		    if (info.lastNbTriangles)
		      sofa::gpu::cuda::mycudaGLUnregisterBufferObject(info.vbo);

		    if (info.uploadBuffer)
		      sofa::gpu::cuda::mycudaFreeHost(info.uploadBuffer);
		    void* ptr = NULL;
		    sofa::gpu::cuda::mycudaMallocHost(&ptr, nbt*3*3*sizeof(float));
		    info.uploadBuffer = static_cast<float*>(ptr);
		    //glBufferDataARB(GL_ARRAY_BUFFER_ARB, nbt*3*3*sizeof(float), NULL, (info.model->isMoving()?GL_DYNAMIC_DRAW_ARB:GL_STATIC_DRAW_ARB));
		  }
		else
		  glBufferDataARB(GL_ARRAY_BUFFER_ARB, nbt*3*3*sizeof(float), NULL, (info.model->isMoving()?GL_DYNAMIC_DRAW_ARB:GL_STATIC_DRAW_ARB));
#else
		glBufferDataARB(GL_ARRAY_BUFFER_ARB, nbt*3*3*sizeof(float), NULL, (info.model->isMoving()?GL_DYNAMIC_DRAW_ARB:GL_STATIC_DRAW_ARB));
#endif
	      }
	    else
	      {
		//std::cout << "Updating VBO for "<<info.model->getName()<<" with "<<nbt<<" triangles"<<std::endl;
	      }

	    {
	      float* data=NULL;
#ifdef SOFA_GPU_CUDA
	      if (useCUDA)
		data = info.uploadBuffer;
	      else
		data = static_cast<float*>(glMapBufferARB(GL_ARRAY_BUFFER_ARB, GL_WRITE_ONLY_ARB));
#else
	      data = static_cast<float*>(glMapBufferARB(GL_ARRAY_BUFFER_ARB, GL_WRITE_ONLY_ARB));
#endif
	      const sofa::core::behavior::MechanicalState< Vec3Types >* mstate = info.model->getMechanicalState();
	      const Vec3Types::VecCoord& x = *mstate->getX();
	      for (int i=0;i<nbt;i++)
		{
		  Triangle t(info.model, i);
		  int v;
		  v = t.p1Index();
		  *(data++) = (float)x[v][0];
		  *(data++) = (float)x[v][1];
		  *(data++) = (float)x[v][2];
		  v = t.p2Index();
		  *(data++) = (float)x[v][0];
		  *(data++) = (float)x[v][1];
		  *(data++) = (float)x[v][2];
		  v = t.p3Index();
		  *(data++) = (float)x[v][0];
		  *(data++) = (float)x[v][1];
		  *(data++) = (float)x[v][2];
		}
	    }
#ifdef SOFA_GPU_CUDA
	    if (useCUDA)
	      {
		if (nbt != info.lastNbTriangles)
		  {
		    glBufferDataARB(GL_ARRAY_BUFFER_ARB, nbt*3*3*sizeof(float), info.uploadBuffer, (info.model->isMoving()?GL_DYNAMIC_DRAW_ARB:GL_STATIC_DRAW_ARB));
		    sofa::gpu::cuda::mycudaGLRegisterBufferObject(info.vbo);
		  }
		else
		  {
		    void* devicePtr = NULL;
		    sofa::gpu::cuda::mycudaGLMapBufferObject(&devicePtr, info.vbo);
		    sofa::gpu::cuda::mycudaMemcpyHostToDevice(devicePtr, info.uploadBuffer, nbt*3*3*sizeof(float));
		    sofa::gpu::cuda::mycudaGLUnmapBufferObject(info.vbo);
		  }
	      }
	    else
	      glUnmapBufferARB(GL_ARRAY_BUFFER_ARB);
#else
	    glUnmapBufferARB(GL_ARRAY_BUFFER_ARB);
#endif
            glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
	    info.lastNbTriangles = nbt;

	    info.lastUpdateTS = currentTS;
	  }


      GLCHECK();
    }

    //***********************************************************************************************************************************
    TriangleModel* DepthPeelingUtility::getModel(unsigned int index)
    {
      if (index < modelsInfo.size())
	return modelsInfo[index].model;
      else
	return NULL;
    }

    //***********************************************************************************************************************************
    TriangleModel* DepthPeelingUtility::getModel(float index)
    {
      int i = (int)(index *  modelsInfo.size())-1;
      if (i< 0 || i >= (int) modelsInfo.size())
	{
	  std::cout << " PROBLEM !!! " << i <<"\n";
	  return NULL;
	}
      //TODO Verify if i is a good correction!!
      //      return modelsInfo[index].model;
      return modelsInfo[i].model;

    }

    //***********************************************************************************************************************************
    void DepthPeelingUtility::loadModel( TriangleModel * &_model)
    {
      sofa::helper::vector<TriangleModel *> vecModel; vecModel.push_back(_model);
      loadModels(vecModel);
    }

    //***********************************************************************************************************************************
    unsigned int DepthPeelingUtility::doDepthPeeling( TriangleModel *  _models,
						      sofa::helper::vector< sofa::helper::vector< float> >  &layer,
						      const unsigned char orientation, bool usePairShader, double adaptiveLength )
    {
      const BoundingBox bb = constructBoundingBox(_models);
      const sofa::helper::vector< TriangleModel *> vecModel(1,_models);
      return doDepthPeeling(vecModel, layer, bb, orientation, usePairShader, adaptiveLength);
    }


    //***********************************************************************************************************************************
    unsigned int DepthPeelingUtility::doDepthPeeling( const sofa::helper::vector< TriangleModel *>  &_models,
						      sofa::helper::vector< sofa::helper::vector< float> >  &layer,
						      const unsigned char orientation, bool usePairShader, double adaptiveLength)
    {
      BoundingBox bb = constructBoundingBox(_models[0]);
      for (unsigned int i=1;i<_models.size();i++)
	{
	  BoundingBox temp = constructBoundingBox(_models[i]);
	  if (temp.first[0] < bb.first[0]) bb.first[0] = temp.first[0];
	  if (temp.first[1] < bb.first[1]) bb.first[1] = temp.first[1];
	  if (temp.first[2] < bb.first[2]) bb.first[2] = temp.first[2];

	  if (temp.second[0] > bb.second[0]) bb.second[0] = temp.second[0];
	  if (temp.second[1] > bb.second[1]) bb.second[1] = temp.second[1];
	  if (temp.second[2] > bb.second[2]) bb.second[2] = temp.second[2];
	}
      return doDepthPeeling(_models, layer, bb, orientation, usePairShader,adaptiveLength);

    };


    //***********************************************************************************************************************************


    unsigned int DepthPeelingUtility::doDepthPeeling( const sofa::helper::vector< TriangleModel *>  &_models,
						      sofa::helper::vector< sofa::helper::vector< float> >  &layer,
						      const BoundingBox &bb,
						      const unsigned char orientation,
						      bool usePairShader, double adaptiveLength)
    {
      resolution = resolutionMax;
// 	{
// 	  initFBO();
// 	}
      DepthPeelingShader *shader_In_Use = static_cast< DepthPeelingShader * >(shader);

      if (!initialized) init(resolutionMax, resolutionPixel, useCUDA);

      if (_models.size() == 0) return 0;

      loadModels(_models);


      setResolutionTexture(bb);

      if (adaptiveLength != 0.0)
	{
	  //compute the resolution used
	  Vector3 positionRender = (bb.second+bb.first)*0.5;

	  glGetDoublev(GL_MODELVIEW_MATRIX, mvMatrix);

	  // Calculate viewpoint
	  Vector3 eye;
	  eye[0] = ((-mvMatrix[12]*mvMatrix[0])+(-mvMatrix[13]*mvMatrix[1])+
		    (-mvMatrix[14]*mvMatrix[2]));
	  eye[1] = ((-mvMatrix[12]*mvMatrix[4])+(-mvMatrix[13]*mvMatrix[5])+
		    (-mvMatrix[14]*mvMatrix[6]));
	  eye[2] = ((-mvMatrix[12]*mvMatrix[8])+(-mvMatrix[13]*mvMatrix[9])+
		    (-mvMatrix[14]*mvMatrix[10]));



	  const double distance = (positionRender-eye).norm();
	  const int valueMax = (int)(log((float)resolutionMax)/log(2.f));
	  const int value = (int)(distance*valueMax/adaptiveLength);

	  if (resolution/pow(2.f,value) < 8) resolution = 8;
	  else                               resolution = (int)(resolution/pow(2.f,value));

	}

      VOI[0]=(float)(bb.second[(orientation+1)%3] - bb.first[(orientation+1)%3]);
      VOI[1]=(float)(bb.second[(orientation+2)%3] - bb.first[(orientation+2)%3]);
      VOI[2]=(float)(bb.second[(orientation)%3]   - bb.first[(orientation)%3])*1.2f;
      current_orientation = orientation;
      assert(orientation < 3);


      bool init=true;

      GLint viewport[4];
      //Save the current viewport
      glGetIntegerv(GL_VIEWPORT, viewport);
      initDepthPeeling(bb);

      if (usePairShader) shader = &depthPeeling_shader_pair;
      else shader = &depthPeeling_shader;


      glUseProgramObjectARB(*(shader->getShader()));
      glEnable(GL_VERTEX_PROGRAM_TWO_SIDE);
      // 	glDisable(GL_BLEND);
      // 	glDisable(GL_ALPHA);
      if (occlusionQuery[0])
	{

	  // Depth Peeling algorithm using occlusion queries and asynchronous transfers:
	  // 1- L = 0
	  // 2- LOOP:
	  //   3- Render layer L on fbo[L&1]
	  //   4- IF L>=1 and layer L-1 is not empty: Start transfer of layer L-1 to PBO[(L-1)&1]
	  //   5- IF L>=2: Read layer L-2 from PBO[(L-2)&1]
	  //   6- IF L>=1 and layer L-1 is empty: BREAK
	  //   7- L=L+1
	  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, FBO[0]);
	  glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);
	  glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
#ifdef USE_SINGLE_FBO
	  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_TEXTURE_2D, shadowMap[0], 0);
#endif
	  unsigned int prev_layer_size = layer.size();
	  while (true)
	    {
	      //
	      //   3- Render layer L on texture L&1
	      //
	      //std::cout << "Render layer "<<number_layer<<std::endl;


	      glUniform1fARB(shader_In_Use->getLocationResolution(), resolution/(float)resolutionMax);

	      glEnable(GL_TEXTURE_2D);
	      //If NO depth peeling has ever been performed yet, we initialize the depth map.
	      if (init)	  { glBindTexture( GL_TEXTURE_2D, shadowMap[2]);           }
	      else 	  { glBindTexture( GL_TEXTURE_2D, shadowMap[1-(number_layer&1)]); }

	      glUniform1iARB(shader_In_Use->getLocation(),0);

	      glBeginQueryARB(GL_SAMPLES_PASSED_ARB, occlusionQuery[(number_layer&1)]);

	      renderScene (_models, usePairShader);	// render scene culling first layer with old depth buffer

	      glEndQueryARB(GL_SAMPLES_PASSED_ARB);
	      glDisable(GL_TEXTURE_2D);
	      glBindTexture ( GL_TEXTURE_2D, 0 );
	      //
	      //   4- IF L>=1 and layer L-1 is not empty: Start transfer of layer L-1 to PBO[(L-1)&1]
	      //

	      bool prev_layer_empty = false;
	      GLint npix;
	      if (number_layer >= 1)
		{
		  npix = -1;
		  glGetQueryObjectivARB(occlusionQuery[(number_layer-1)&1],GL_QUERY_RESULT_ARB,&npix);
		  //std::cout << "Layer "<<number_layer-1<<" contains "<<npix<<" pixels"<<std::endl;
		  prev_layer_empty = (npix==0);
		}
	      // flip the active FBO
#ifdef USE_SINGLE_FBO
	      glReadBuffer(((number_layer+1)&1) ? GL_COLOR_ATTACHMENT1_EXT:GL_COLOR_ATTACHMENT0_EXT);
	      glDrawBuffer(((number_layer+1)&1) ? GL_COLOR_ATTACHMENT1_EXT:GL_COLOR_ATTACHMENT0_EXT);
	      glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_TEXTURE_2D, shadowMap[((number_layer+1)&1)], 0);
#else
	      glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, FBO[((number_layer+1)&1)]);
	      glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);
	      glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
#endif
	      if (number_layer >= 1 && !prev_layer_empty)
		{
		  //std::cout << "Start transfer of layer "<<number_layer-1<<" to PBO"<<std::endl;
		  //glBindBufferARB(GL_PIXEL_PACK_BUFFER_ARB, PBO[ ((number_layer-1)&1) ]);
		  glReadPixels(0,0,
			       resolution,resolution,
			       GL_RGBA, GL_FLOAT, 0);
		}

	      // flip the active PBO
	      glBindBufferARB(GL_PIXEL_PACK_BUFFER_ARB, PBO[(number_layer&1)]);

	      //
	      //   5- IF L>=2: Read layer L-2 from PBO[(L-2)&1]
	      //
	      if (number_layer >= 2)
		{
		  //std::cout << "Read layer "<<number_layer-2<<" from PBO" << " : " << resolution <<std::endl;
		  if (layer.size() == number_layer-2)
		    layer.push_back(sofa::helper::vector< float >(4*resolutionMax*resolutionMax));
		  const void *mem = glMapBufferARB(GL_PIXEL_PACK_BUFFER_ARB, GL_READ_ONLY_ARB);
		  memcpy(&(layer[number_layer-2][0]), mem, resolution*resolution*4*sizeof(float));
		  glUnmapBufferARB(GL_PIXEL_PACK_BUFFER_ARB);
		  pixel_location.push_back(std::make_pair(0,resolution*resolution));
		}

	      //
	      //   6- IF L>=1 and layer L-1 is empty: BREAK
	      //

	      if (number_layer>=1 && prev_layer_empty)
		{
		  //std::cout << "BREAK"<<std::endl;
		  break;
		}
	      else if (number_layer > MAXLAYER)
		{
		  break;
		}

	      //
	      //   7- L=L+1
	      //

	      ++number_layer;
	      init = false;
	    }
	  //std::cout << "DepthPeelingUtility::doDepthPeeling: generated "<<pixel_location.size()<<" layers."<<std::endl;
	  number_layer -= 2; // the last 2 layers are empty
	  glBindBufferARB(GL_PIXEL_PACK_BUFFER_ARB,0);
	  if (prev_layer_size != layer.size())
	    std::cout << "DepthPeeling: max layers = "<<layer.size()<<std::endl;
	}
      else
	{
#ifdef USE_SINGLE_FBO
	  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, FBO[0]);
#endif
	  while (true)
	    {
	      //std::cout << "DP layer "<<number_layer<<std::endl;
	      if (layer.size() == number_layer)
		layer.push_back(sofa::helper::vector< float >(4*resolutionMax*resolutionMax));


	      //Need to perform a Depth Peeling
	      ping_pong = 1 - ping_pong;
#ifdef USE_SINGLE_FBO
	      glReadBuffer(ping_pong ? GL_COLOR_ATTACHMENT1_EXT:GL_COLOR_ATTACHMENT0_EXT);
	      glDrawBuffer(ping_pong ? GL_COLOR_ATTACHMENT1_EXT:GL_COLOR_ATTACHMENT0_EXT);
	      glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_TEXTURE_2D, shadowMap[ping_pong], 0);
#else
	      glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, FBO[ping_pong]);
	      glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);
	      glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
#endif

	      glUniform1fARB(shader_In_Use->getLocationResolution(), resolution/(float)resolutionMax);

	      glEnable(GL_TEXTURE_2D);
	      //If NO depth peeling has ever been performed yet, we initialize the depth map.
	      if (init)	  { glBindTexture( GL_TEXTURE_2D, shadowMap[2]);           }
	      else 	  { glBindTexture( GL_TEXTURE_2D, shadowMap[1-ping_pong]); }

	      glUniform1iARB(shader_In_Use->getLocation(),0);

	      renderScene (_models, usePairShader);	// render scene culling first layer with old depth buffer

	      glDisable(GL_TEXTURE_2D);
	      glBindTexture ( GL_TEXTURE_2D, 0 );

	      glBindBufferARB(GL_PIXEL_PACK_BUFFER_ARB, PBO[ ping_pong ]);

	      glReadPixels(0,0,
			   resolution,resolution,
			   GL_RGBA, GL_FLOAT, 0);
	      void *mem = glMapBufferARB(GL_PIXEL_PACK_BUFFER_ARB, GL_READ_ONLY_ARB);
	      memcpy(&(layer[number_layer][0]), mem, resolution*resolution*4*sizeof(float));
	      glUnmapBufferARB(GL_PIXEL_PACK_BUFFER_ARB);
	      glBindBufferARB(GL_PIXEL_PACK_BUFFER_ARB,0);

	      // manually look for rendered pixels
	      ++number_layer;

	      //	    glReadBuffer(GL_NONE);

	      bool is_active = computeInterval( number_layer-1, layer);
	      if (!is_active) break;
	      init = false;
	    }
	}
      glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
      glReadBuffer(GL_BACK);
      glDrawBuffer(GL_BACK);
      glUseProgramObjectARB(0);
      glDisable(GL_VERTEX_PROGRAM_TWO_SIDE);

      endRenderView();
      glViewport(viewport[0],viewport[1],viewport[2],viewport[3]);
      
      return resolution;
    }


    struct layer_sort_fn
    {
      bool operator()(const Vec<4,float>& a, const Vec<4,float>& b) { return a[0] < b[0]; }
    };

    //***********************************************************************************************************************************
    void DepthPeelingUtility::doVoxelize(const sofa::helper::vector< TriangleModel *>  &_models,
					 sofa::helper::vector< sofa::helper::vector< sofa::helper::vector< unsigned char > > > &texture3d,
					 const unsigned char orientation )
    {
      sofa::helper::vector< sofa::helper::vector< float> > layer;
      doDepthPeeling(_models, layer, orientation, true);
      doVoxelize(layer, texture3d);
    }


    //***********************************************************************************************************************************
    void DepthPeelingUtility::doVoxelize( TriangleModel *  _models,
					  sofa::helper::vector< sofa::helper::vector< sofa::helper::vector< unsigned char > > > &texture3d,
					  const unsigned char orientation )
    {
      sofa::helper::vector< sofa::helper::vector< float> > layer;
      doDepthPeeling(_models, layer, orientation, true);
      doVoxelize(layer, texture3d);
    }

    //***********************************************************************************************************************************
    //Configure the FBO
    void DepthPeelingUtility::initFBO()
    {
      //Initialization of the first FBO

      glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, FBO[0]);
      //Texture of identification :
      glBindTexture(GL_TEXTURE_2D, idsMap[0]);
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, resolutionMax, resolutionMax,
		   0, GL_RGBA, GL_FLOAT, NULL);
      glTexParameteri(GL_TEXTURE_2D,
		      GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D,
		      GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D,
		      GL_TEXTURE_WRAP_S, GL_CLAMP);
      glTexParameteri(GL_TEXTURE_2D,
		      GL_TEXTURE_WRAP_T, GL_CLAMP);
      //Attach this texture to the Color info
      glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, idsMap[0], 0);


      //Shadow Map for depth test
	{
	  glBindTexture(GL_TEXTURE_2D, shadowMap[0]);
	  glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT32, resolutionMax, resolutionMax,
		       0, GL_DEPTH_COMPONENT, GL_FLOAT, init_shadowmap);
	  glTexParameteri(GL_TEXTURE_2D,
			  GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	  glTexParameteri(GL_TEXTURE_2D,
			  GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	  glTexParameteri(GL_TEXTURE_2D,
			  GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	  glTexParameteri(GL_TEXTURE_2D,
			  GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE_ARB, GL_COMPARE_R_TO_TEXTURE_ARB);
	  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC_ARB, GL_LEQUAL);
	  //Attach this texture to the Depth info
	  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_TEXTURE_2D, shadowMap[0], 0);
	}

      {
	//Check FBO
	GLenum status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
	if (status == GL_FRAMEBUFFER_COMPLETE_EXT) ;//std::cout << "FBO OK"<<std::endl;
	else
	  {
	    std::cerr << "ERROR: FBO status "<<status<<std::endl;
	  }
      }
      //glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
      //End

      //******************************************************************************************************

      //Initialization of the second FBO
#ifndef USE_SINGLE_FBO
      glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, FBO[1]);
#endif
      //Texture of identification :
      glBindTexture(GL_TEXTURE_2D, idsMap[1]);
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, resolutionMax, resolutionMax,
		   0, GL_RGBA, GL_FLOAT, NULL);
      glTexParameteri(GL_TEXTURE_2D,
		      GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D,
		      GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D,
		      GL_TEXTURE_WRAP_S, GL_CLAMP);
      glTexParameteri(GL_TEXTURE_2D,
		      GL_TEXTURE_WRAP_T, GL_CLAMP);
      //Attach this texture to the Color info
#ifdef USE_SINGLE_FBO
      glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT1_EXT, GL_TEXTURE_2D, idsMap[1], 0);
#else
      glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, idsMap[1], 0);
#endif

	{
	  //Shadow Map for depth test
	  glBindTexture(GL_TEXTURE_2D, shadowMap[1]);
	  glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT32, resolutionMax, resolutionMax,
		       0, GL_DEPTH_COMPONENT, GL_FLOAT, init_shadowmap);
	  glTexParameteri(GL_TEXTURE_2D,
			  GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	  glTexParameteri(GL_TEXTURE_2D,
			  GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	  glTexParameteri(GL_TEXTURE_2D,
			  GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	  glTexParameteri(GL_TEXTURE_2D,
			  GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE_ARB, GL_COMPARE_R_TO_TEXTURE_ARB);
	  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC_ARB, GL_LEQUAL);
	  //Attach this texture to the Depth info
	  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_TEXTURE_2D, shadowMap[1], 0);
	}

      //Check FBO
      {
	GLenum status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
	if (status == GL_FRAMEBUFFER_COMPLETE_EXT) ;//std::cout << "FBO OK"<<std::endl;
	else
	  {
	    std::cerr << "ERROR: FBO status "<<status<<std::endl;
	  }
      }

      glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);

	{
	  //Shadow Map Initialization for depth test
	  glBindTexture(GL_TEXTURE_2D, shadowMap[2]);
	  glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT32, resolutionMax, resolutionMax,
		       0, GL_DEPTH_COMPONENT, GL_FLOAT, init_shadowmap);
	  glTexParameteri(GL_TEXTURE_2D,
			  GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	  glTexParameteri(GL_TEXTURE_2D,
			  GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	  glTexParameteri(GL_TEXTURE_2D,
			  GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	  glTexParameteri(GL_TEXTURE_2D,
			  GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE_ARB, GL_COMPARE_R_TO_TEXTURE_ARB);
	  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC_ARB, GL_LEQUAL);

	  //glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_TEXTURE_2D, shadowMap[2], 0);
	  //Attach this texture to the Depth info
	  //End
	}
      //delete [] init_shadowmap;
      //delete [] init_stencilmap;
    }

    //***********************************************************************************************************************************
    bool DepthPeelingUtility::isGridTopology(TriangleModel* m) const
    {
      using namespace core::objectmodel;
      sofa::component::topology::RegularGridTopology* grid;
      m->getContext()->get(grid, core::objectmodel::BaseContext::SearchUp);
      if (grid != NULL)
	{
	  if ( grid->getNx() == 2 && grid->getNy() == 2 && grid->getNz() == 2)
	    {
	      return true;
	    }
	}
      return false;
    }

    //***********************************************************************************************************************************
    void DepthPeelingUtility::renderScene (  const  sofa::helper::vector< TriangleModel *>  &_models, bool usePairShader ) const
    {

      glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);


      glBindBufferARB( GL_ARRAY_BUFFER_ARB, idVBOIdentity);
      glColorPointer ( 4, GL_FLOAT, 0, 0);
      glEnableClientState( GL_COLOR_ARRAY );

      // 	glEnableClientState( GL_VERTEX_ARRAY ); //should be ENABLED, but already done and kept on in the viewer

      for (unsigned int i=0;i<_models.size();i++)
	{
	  std::map< TriangleModel*,unsigned int>::const_iterator iter = modelsId.find(_models[i]);
	  if (iter == modelsId.end()) continue;

	  //get the index of the VBO corresponding to a model
	  unsigned int index_model=iter->second;
	  const ModelInfo& info = modelsInfo[index_model];

	  float IdVBO;
	  if (usePairShader)
	    IdVBO = i*1.0f;
	  else
	    IdVBO = (index_model+1)/((float)nbModels);


	  //SPECIAL CASE of a regular grid 2x2x2 object. We need to pass the coordinates of the 8 points.
	  //And find the good VBO (as several "bunny" object are stored in a unique VBO for example
	  Vec3f bb_grid[8];
	  float RegGrid = 1.0f;
	  if (info.gridState != NULL && info.gridTopology != NULL)
	    {
	      RegGrid = 0.0f;
	      const sofa::helper::vector< Vector3 > &dof = *info.gridState->getX();

	      bb_grid[0] = dof[ info.gridTopology->point( 0,0,0) ];
	      bb_grid[1] = dof[ info.gridTopology->point( info.gridTopology->getNx()-1,0,0)];
	      bb_grid[2] = dof[ info.gridTopology->point( 0,info.gridTopology->getNy()-1,0) ];
	      bb_grid[3] = dof[ info.gridTopology->point( info.gridTopology->getNx()-1,info.gridTopology->getNy()-1,0) ];
	      bb_grid[4] = dof[ info.gridTopology->point( 0,0,info.gridTopology->getNz()-1) ];
	      bb_grid[5] = dof[ info.gridTopology->point( info.gridTopology->getNx()-1,0,info.gridTopology->getNz()-1) ];
	      bb_grid[6] = dof[ info.gridTopology->point( 0,info.gridTopology->getNy()-1, info.gridTopology->getNz()-1) ];
	      bb_grid[7] = dof[ info.gridTopology->point( info.gridTopology->getNx()-1,info.gridTopology->getNy()-1, info.gridTopology->getNz()-1) ];
	      glUniform3fvARB(shader->getLocationBB(),8,(GLfloat*) &bb_grid[0]);
	    }

	  glUniform1fARB(shader->getLocationRegGrid(),RegGrid);

	  glUniform1fARB(shader->getLocationIdModels(),(GLfloat) IdVBO);
	  //Render the VBO

	  //glBindBufferARB( GL_ARRAY_BUFFER_ARB, VBO[2*indexVBO+IDENTITY]);
	  //glColorPointer ( 4, GL_FLOAT, 0, 0);

	  glBindBufferARB( GL_ARRAY_BUFFER_ARB, info.vbo);
	  glVertexPointer( 3, GL_FLOAT, 0, 0);

	  glDrawArrays(GL_TRIANGLES, 0, 3*_models[i]->getSize());
	}

      glDisableClientState( GL_COLOR_ARRAY );

      // 	glDisableClientState( GL_VERTEX_ARRAY );  //should be DISABLED, but a brutal initialization in the viewer was done before


      glBindBufferARB( GL_ARRAY_BUFFER_ARB, 0);

      glFlush();
    }


    //***********************************************************************************************************************************
    //Compute the interval where the information in the texture is present. If the interval is empty, it means that the whole cell is empty.
    bool DepthPeelingUtility::computeInterval( const unsigned int index_layer,
					       const sofa::helper::vector< sofa::helper::vector< float> >  &layer )
    {

      unsigned int init, end;
      unsigned int initial= 0;
      unsigned int final  = resolution;
      unsigned int i;

      if (index_layer == 0)
	{
	  init = 0;
	  end=resolution*resolution;
	}
      else
	{
	  init = pixel_location[index_layer-1].first;
	  end  = pixel_location[index_layer-1].second;
	}

      for (i=init;i<end;i++)
	{

	  if (layer[index_layer][(i*4)] != 0)
	    {
	      initial=i;
	      break;
	    }
	}
      if (i==end) return false;

      for (i=end -1;i>=init;i--)
	{
	  if (layer[index_layer][(i*4)] !=0)
	    {
	      final=i;
	      break;
	    }
	}
      if (initial == final)
	{
	  return false;
	}
      else
	{
	  final++;
	  pixel_location.push_back(std::make_pair( initial,final ));
	  return true;
	}
    }


    //***********************************************************************************************************************************
    void DepthPeelingUtility::initDepthPeeling(const BoundingBox &visualizationBB)
    {
      number_layer = 0;
      pixel_location.clear();

      ping_pong = 1;

      //Create the new viewport
      glViewport(0,0, resolution, resolution);


      initRenderView( visualizationBB );
    }


    //***********************************************************************************************************************************
    void DepthPeelingUtility::initRenderView( const BoundingBox &bb )
    {
      BoundingBox boundingBox = bb;

      //Compute the World Bounding Box
      sceneMinBBox = boundingBox.first;
      sceneMaxBBox = boundingBox.second;
      Vector3 size_BB = sceneMaxBBox - sceneMinBBox;
      sceneMinBBox[current_orientation] -= 0.1*size_BB[current_orientation];
      sceneMaxBBox[current_orientation]+= 0.1*size_BB[current_orientation];
      size_BB = sceneMaxBBox - sceneMinBBox;


      setDepthBufferTest();

      glMatrixMode   ( GL_PROJECTION );                   // select the projection matrix
      glPushMatrix   ();                                  // store the projection matrix
      glLoadIdentity ();                                  // reset the projection matrix

      // set up an ortho screen
      glOrtho( -size_BB[(current_orientation+1)%3],0,
	       -size_BB[(current_orientation+2)%3],0,
	       0,size_BB[current_orientation]);

      glMatrixMode   ( GL_MODELVIEW );                    // select the modelview matrix
      glPushMatrix   ();                                  // store the modelview matrix
      glLoadIdentity ();                                  // reset the modelview matrix
      //Set the camera

      Vec<3, int> orientation_depth_peeling = Vec<3,int>();
      Vec<3, int> orientation2_depth_peeling = Vec<3,int>();
      orientation_depth_peeling[current_orientation] = -1;

      orientation2_depth_peeling[(current_orientation+2)%3] = 1;

      gluLookAt(sceneMaxBBox.x(),
		sceneMaxBBox.y(),
		sceneMaxBBox.z(),

		sceneMaxBBox.x() +orientation_depth_peeling[0],
		sceneMaxBBox.y() +orientation_depth_peeling[1],
		sceneMaxBBox.z() +orientation_depth_peeling[2],

		+orientation2_depth_peeling[0],
		+orientation2_depth_peeling[1],
		+orientation2_depth_peeling[2]);

    }



    //***********************************************************************************************************************************
    void DepthPeelingUtility::endRenderView()
    {
      //End: back to initial state
      glMatrixMode   ( GL_MODELVIEW );                    // select the modelview matrix
      glPopMatrix   ();                                  // store the modelview matrix

      glMatrixMode   ( GL_PROJECTION );                    // select the modelview matrix
      glPopMatrix   ();                                  // store the modelview matrix
      resetDepthBufferTest();
    }


    //***********************************************************************************************************************************
    void DepthPeelingUtility::setDepthBufferTest() const
    {
      glClearColor   ( 0.0, 0.0, 0.0, 0.0 );
      glClearStencil(0);
      glClearDepth(1.0f);
      glEnable       ( GL_DEPTH_TEST );
      glDepthFunc    ( GL_LEQUAL );

      glHint ( GL_POLYGON_SMOOTH_HINT,         GL_NICEST );
      glHint ( GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST );
    }


    //***********************************************************************************************************************************
    void DepthPeelingUtility::resetDepthBufferTest() const
    {
      glDepthFunc(GL_LESS);
    }


    //***********************************************************************************************************************************
    BoundingBox DepthPeelingUtility::constructBoundingBox(  TriangleModel *m)
    {
      core::CollisionModel *model;
      using component::collision::CubeModel;
      //Calculate the hierarchical bounding boxes of the triangle model
      m->computeBoundingTree(6);
      model = m->getPrevious();

      CubeModel *last_model=NULL;
      while (model != NULL)
	{
	  if (CubeModel *BB = dynamic_cast< CubeModel *> (model) )
	    last_model = BB;

	  model = model->getPrevious();
	}

      sofa::helper::vector< BoundingBox > bounding;
      last_model->getBoundingTree( bounding );
      return bounding[0];
    }


    void DepthPeelingUtility::fillArrays( TriangleModel *m)
    {
      using namespace sofa::core::objectmodel;
      const int size = m->getSize();

      //compute the bounding box of the mechanical object
      sofa::helper::vector< sofa::core::behavior::MechanicalState< Vec3Types >* > mecha;


      Vector3 boundingbox[2] = {Vector3(0,0,0),Vector3(0,0,0)};
      Vector3 sizeBB(1,1,1);
      //looking for a regular grid 2x2x2

	  sofa::component::topology::RegularGridTopology* grid;
	  m->getContext()->get(grid, core::objectmodel::BaseContext::SearchUp);
	  if (isGridTopology(m))
	    {
#ifndef SOFA_FLOAT
            {
              core::behavior::MechanicalState< Vec3dTypes >* mstate;
	      grid->getContext()->get(mstate);

	      if (mstate != NULL)
		{
 		  const sofa::helper::vector< Vec3d > &dof = *mstate->getX();
      

 		  boundingbox[0] = dof[ grid->point(0,0,0) ];
 		  boundingbox[1] = dof[ grid->point( grid->getNx()-1,grid->getNy()-1, grid->getNz()-1) ];
 		  sizeBB = boundingbox[1] - boundingbox[0];
		}
            }
#endif
#ifndef SOFA_DOUBLE
            {
	      core::behavior::MechanicalState< Vec3fTypes >* mstate;
	      grid->getContext()->get(mstate);

	      if (mstate != NULL)
		{
 		  const sofa::helper::vector< Vec3f > &dof = *mstate->getX();
      
 		  boundingbox[0] = dof[ grid->point(0,0,0) ];
 		  boundingbox[1] = dof[ grid->point( grid->getNx()-1,grid->getNy()-1, grid->getNz()-1) ];
 		  sizeBB = boundingbox[1] - boundingbox[0];
		}
            }
#endif

	    }
      sizeBB[0] = 1.0/sizeBB[0];sizeBB[1] = 1.0/sizeBB[1];sizeBB[2] = 1.0/sizeBB[2];
      for ( int i=0;i<3*size; i+=3)
	{
	  //For each triangle of the model, we store the coordinates of the vertices and information about each of them
	  Triangle t(m, i/3);
	  //Point 1
	  array_coord[i*3  ]     = (float) ((t.p1()[0]-boundingbox[0][0])*sizeBB[0]);
	  array_coord[i*3+1]     = (float) ((t.p1()[1]-boundingbox[0][1])*sizeBB[1]);
	  array_coord[i*3+2]     = (float) ((t.p1()[2]-boundingbox[0][2])*sizeBB[2]);
	  //Point 2
	  array_coord[(i+1)*3+0] = (float) ((t.p2()[0]-boundingbox[0][0])*sizeBB[0]);
	  array_coord[(i+1)*3+1] = (float) ((t.p2()[1]-boundingbox[0][1])*sizeBB[1]);
	  array_coord[(i+1)*3+2] = (float) ((t.p2()[2]-boundingbox[0][2])*sizeBB[2]);
	  //Point 3
	  array_coord[(i+2)*3+0] = (float) ((t.p3()[0]-boundingbox[0][0])*sizeBB[0]);
	  array_coord[(i+2)*3+1] = (float) ((t.p3()[1]-boundingbox[0][1])*sizeBB[1]);
	  array_coord[(i+2)*3+2] = (float) ((t.p3()[2]-boundingbox[0][2])*sizeBB[2]);
	}
    }

    //***********************************************************************************************************************************
    void DepthPeelingUtility::printDepthPeeling(const  sofa::helper::vector< sofa::helper::vector< float> >  &layer)
    {
      if ((void *) shader == &depthPeeling_shader_pair)
	{
	  for (unsigned int level=0;level < layer.size();level++)
	    {
	      for (int y=resolution-1;y>=0;y--)
		{
		  std::cout << "\t";
		  for (unsigned int x=0;x<resolution;x++)
		    {
		      if (layer[level][((y*resolution + x)<<2)+1]==0)
			std::cout << '_';
		      else
			{

			  int index_model = (layer[level][((y*resolution + x)<<2)+1]>0.5?1:0);
			  bool front = (layer[level][((y*resolution + x)<<2)+2]>0.5);
			  std::cout << (char)((front ? 'A':'a')+index_model);
			}
		    }
		  std::cout << "\n";
		}
	      std::cout << "\n\n";
	    }
	}
      else
	{
	  for (unsigned int level=0;level < layer.size();level++)
	    {
	      for (int y=resolution-1;y>=0;y--)
		{
		  std::cout << "\t";
		  for (unsigned int x=0;x<resolution;x++)
		    {
		      int index_model = (int)(nbModels*layer[level][((y*resolution + x)<<2)]) ;

		      if (index_model<=0) std::cout << "_";
		      else                        std::cout << index_model;

		    }
		  std::cout << "\n";
		}
	      std::cout << "\n\n";
	    }
	}
    }


    //***********************************************************************************************************************************
    void DepthPeelingUtility::printVoxelization(const sofa::helper::vector< sofa::helper::vector< sofa::helper::vector< unsigned char > > > &texture3d)
    {
      for (unsigned int z=0;z<texture3d.size();z++)
	{
	  for (unsigned int x=0;x<texture3d[0].size();x++)
	    {
	      for (int y=texture3d[0][0].size()-1;y>=0;y--)
		{
		  if ((int)(texture3d[z][x][y]) == 0) std::cout << " ";
		  else std::cout << (int)(texture3d[z][x][y]);
		}
	      std::cout << "\n";
	    }
	  std::cout << "###############################################################\n\n";
	}
    }



    //***********************************************************************************************************************************
    // Determine the depth of a given pixel from the gpu"s texture
    float DepthPeelingUtility::getDepth(const float i, const float alpha, const float beta)
    {
      TriangleModel *m = getModel((unsigned int)(i-1));
      Triangle t(m,(int)((i - (int)(i))*1.1f*m->getSize()+0.5f));
      return (float)((1-alpha-beta)*t.p1()[current_orientation] + alpha*t.p2()[current_orientation] + beta*t.p3()[current_orientation]);
    }

    //***********************************************************************************************************************************
    // Determine the depth of a given pixel from the gpu"s texture
    float DepthPeelingUtility::getDepth(const float z)
    {
      return (float)(sceneMinBBox[current_orientation]+(sceneMaxBBox[current_orientation]-sceneMinBBox[current_orientation])*z);
    }


    //***********************************************************************************************************************************
    void DepthPeelingUtility::doVoxelize(const sofa::helper::vector< sofa::helper::vector< float> > &layer,
					 sofa::helper::vector< sofa::helper::vector< sofa::helper::vector< unsigned char > > > &texture3d)
    {
      bool init = true;
      float min = -1;
      float max = -1;

      for (unsigned int level=0;level<layer.size()-1;level++)
	{
	  for (int y=resolution-1;y>=0;y--)
	    {
	      for (unsigned int x=0;x<resolution;x++)
		{
		  if (layer[level][((y*resolution + x)<<2)] != 0)
		    {
		      float depth =  getDepth(nbModels*layer[level][((y*resolution + x)<<2)],layer[level][((y*resolution + x)<<2)+1],layer[level][((y*resolution + x)<<2)+2]);
		      if (depth > max) max = depth;
		      if (init || depth < min) {min = depth;init = false;}
		    }
		}
	    }
	}
      const float dx=VOI[0]/((float) resolution);
      const float dy=VOI[1]/((float) resolution);

      unsigned int X_texture, Y_texture;
      if (VOI[0] > VOI[1])
	{
	  X_texture = resolution;
	  Y_texture = (unsigned int) ceil(X_texture*(VOI[1]/VOI[0]));
	}
      else
	{
	  Y_texture = resolution;
	  X_texture = (unsigned int) ceil(Y_texture*(VOI[0]/VOI[1]));
	}

      const float step = (dx>dy?dx:dy);
      unsigned int depth_texture = (unsigned int ) ((max-min)/step)+1;

      //Initialization of the texture3d
      texture3d.resize(depth_texture);
      for (unsigned int z=0;z<depth_texture;z++)
	{
	  texture3d[z].resize(X_texture);
	  for (unsigned int x=0;x<X_texture;x++)
	    {
	      texture3d[z][x].resize(Y_texture);
	    }
	}

      //Filling the texture3d
      for (int y=resolution-1;y>=0;y--)
	{
	  for (unsigned int x=0;x<resolution;x++)
	    {
	      sofa::helper::vector< std::pair < std::pair< int, unsigned int> , bool > > ray;
	      for (unsigned int level=0;level < layer.size()-1;level++)
		{
		  int index_model =(int)(nbModels*layer[level][((y*resolution + x)<<2)])-1 ;
		  if (index_model<0) break;

		  unsigned int depth = (unsigned int)( (getDepth(nbModels*layer[level][((y*resolution + x)<<2)],layer[level][((y*resolution + x)<<2)+1],layer[level][((y*resolution + x)<<2)+2])- min)/step);

		  //for (std::map< TriangleModel *, unsigned int>::const_iterator it=models.begin(); it != models.end(); it++)
		  {
		    //if ((int)(*it).second == index_model)
		    if ((unsigned)index_model < modelsInfo.size())
		      {
			TriangleModel *m = modelsInfo[index_model].model;
			unsigned int index_triangle = (unsigned int)( (m->getSize()*1.1f)*(nbModels*layer[level][(y*resolution + x)<<2] - (int)(nbModels*layer[level][(y*resolution + x)<<2]))+0.5f);

			Triangle t(m,index_triangle);

			ray.push_back(std::make_pair( std::make_pair(index_model,depth) , frontFace(t) ) );
			break;
		      }
		  }
		}
	      fillTexture3d(x,y,ray,texture3d);

	    }
	}
    }


    //***********************************************************************************************************************************
    bool DepthPeelingUtility::frontFace(Triangle &t) const
    {
      Vector3 normal_triangle=t.n();
      if (normal_triangle[current_orientation] > 0) return true;
      else return false;
    }



    //***********************************************************************************************************************************
    void DepthPeelingUtility::fillTexture3d(unsigned int _x, unsigned int _y, const sofa::helper::vector< std::pair < std::pair< int, unsigned int> , bool > > &ray,
					    sofa::helper::vector< sofa::helper::vector< sofa::helper::vector< unsigned char > > > &texture3d) const
    {
      //Conversion from LDI coordinates to texture3d coordinates
      float factor[2]= {texture3d[0].size()/((float)resolution), texture3d[0][0].size()/((float)resolution)};

      unsigned int x = (unsigned int)(_x*factor[0]);
      unsigned int y = (unsigned int)(_y*factor[1]);
      for (unsigned int i=0;i<ray.size();i++)
	{
	  //Depth of the first element
	  unsigned int depth_final = ray[i].first.second;
	  //FRONT
	  if (ray[i].second == 1)
	    {
	      for (unsigned int j=i+1;j<ray.size();j++)
		{
		  //if same element
		  if (ray[j].first.first == ray[i].first.first)
		    {
		      if (ray[j].second == 0) depth_final = ray[j].first.second;
		      break;
		    }
		}
	    }
	  //  	  std::cout << texture3d.size() << " " << texture3d[0].size() << " " << texture3d[0][0].size()
	  //  	            <<  " Filling : " << ray[i].first.second << " " << depth_final << " \t\t" << i << "\n";
	  for (int z=ray[i].first.second; z> (int)depth_final;z--)
	    {
	      texture3d[z][x][y] = ray[i].first.first +1;
	    }
	}
    }

#ifdef SOFA_GPU_CUDA
    //***********************************************************************************************************************************      
    void DepthPeelingUtility::initFBOVolume()
    {
      glGenTextures(1, &idLTexture);
      glBindTexture(GL_TEXTURE_2D_ARRAY_EXT, idLTexture);
      glTexImage3D(GL_TEXTURE_2D_ARRAY_EXT, 0, GL_RGBA32F_ARB, resolutionMax, resolutionMax, MAXLAYER,
		   0, GL_RGBA, GL_FLOAT, NULL);
      glTexParameteri(GL_TEXTURE_2D_ARRAY_EXT,
		      GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D_ARRAY_EXT,
		      GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D_ARRAY_EXT,
		      GL_TEXTURE_WRAP_S, GL_CLAMP);
      glTexParameteri(GL_TEXTURE_2D_ARRAY_EXT,
		      GL_TEXTURE_WRAP_T, GL_CLAMP);
      glBindTexture(GL_TEXTURE_2D_ARRAY_EXT, 0);

      glGenFramebuffersEXT(1, &FBOVolume);

      glGenQueriesARB(1, &geometryQuery);
      //glGenBuffers(1,&geometryBuffer);
      //glBindBufferARB(GL_TRANSFORM_FEEDBACK_BUFFER_NV, geometryBuffer);
      //glBufferDataARB(GL_TRANSFORM_FEEDBACK_BUFFER_NV,  1024*1024, NULL, GL_STREAM_READ_ARB);
      //glBindBufferARB(GL_TRANSFORM_FEEDBACK_BUFFER_NV,0);

      glGenBuffers(1,&VBOpixels);
      glBindBufferARB(GL_ARRAY_BUFFER_ARB, VBOpixels);
      glBufferDataARB(GL_ARRAY_BUFFER_ARB,  resolutionMax*resolutionMax*2*sizeof(short), NULL, GL_STATIC_DRAW_ARB);
      glBindBufferARB(GL_ARRAY_BUFFER_ARB,0);

      volumeResultSize = 0;
      prevVolumeResultSize = 0;
      volumeResultCurrent = -1;
      volumeResultMapping = NULL;

      glGenBuffers(1,&idLPBO);
      glBindBufferARB(GL_PIXEL_PACK_BUFFER_ARB, idLPBO);
      glBufferDataARB(GL_PIXEL_PACK_BUFFER_ARB, resolutionMax*resolutionMax*MAXLAYER*4*sizeof(float), NULL, GL_DYNAMIC_DRAW_ARB);
      glBindBufferARB(GL_PIXEL_PACK_BUFFER_ARB, 0);
      scan_in.resize(resolutionMax*resolutionMax);
      //scan_out.resize(resolution*resolution);
      scan_out.resize(2*((resolutionMax*resolutionMax+BSIZE-1)/BSIZE));
      partial_sums.resize(((resolutionMax*resolutionMax+BSIZE-1)/BSIZE));
      output_collisions.resize(16*1024*2*4);
      /*scan_config.direction = CUDPP_SCAN_FORWARD;
	scan_config.exclusivity = CUDPP_SCAN_EXCLUSIVE;
	scan_config.op = CUDPP_ADD;
	scan_config.datatype = CUDPP_UINT;
	scan_config.maxNumElements = resolutionMax*resolutionMax;
	scan_config.maxNumRows = 1;
	scan_config.rowPitch = 0;
	cudppInitializeScan(&scan_config);*/
    }

    extern "C"
    {
      void CollisionVolume_count(const void* layers, void* counts, void* count_per_bloc, const int nlayers, const int npixels, const int bsize, bool self);
      void CollisionVolume_write(const void* layers, const void* counts, const void* bloc_input_pos, const void* bloc_output_pos, void* collisions, const int nlayers, const int npixels, const int bsize, const int nblocs, bool self);
    }

    void DepthPeelingUtility::doCollisionVolume( const sofa::helper::vector< TriangleModel *>  &_models,
						 int &result_id,
						 const BoundingBox &bb,
						 const unsigned char orientation,
						 simulation::Node* timelog, core::objectmodel::BaseObject* timeobj)
    {
      if (!initialized) init(resolutionMax, resolutionPixel, true);
      if (_models.size() == 0) return;
      simulation::Node::ctime_t t0 = 0;
      if (timelog) t0 = timelog->startTime();
      //if (models.size() == 0)
      {
        loadModels(_models);
        if (timelog) t0 = timelog->endTime(t0, "collision/NarrowPhase/GPU/LoadModels", timeobj);
      }

      resolution = resolutionMax;

      GLCHECK();

      VOI[0]=(float)(bb.second[(orientation+1)%3] - bb.first[(orientation+1)%3]);
      VOI[1]=(float)(bb.second[(orientation+2)%3] - bb.first[(orientation+2)%3]);
      VOI[2]=(float)(bb.second[(orientation)%3]   - bb.first[(orientation)%3])*1.2f;
      current_orientation = orientation;
      assert(orientation < 3);
      bool init=true;

      GLint viewport[4];
      //Save the current viewport
      glGetIntegerv(GL_VIEWPORT, viewport);
      initDepthPeeling(bb);

      shader = &depthPeeling_shader_pair;

      //////////////////////////// PROBLEM !!

      glUseProgramObjectARB(*(shader->getShader()));
      glEnable(GL_VERTEX_PROGRAM_TWO_SIDE);

      if (!occlusionQuery[0])
	{
	  std::cerr << "ERROR: DepthPeelingUtility::doCollisionVolume requires occlusion queries support." << std::endl;
	  result_id = -1;
	  return;
	}
      else
	{
	  // Depth Peeling algorithm using occlusion queries and asynchronous transfers:
	  // 1- L = 0
	  // 2- LOOP:
	  //   3- Render layer L on fbo[L&1]
	  //   6- IF L>=1 and layer L-1 is empty: BREAK
	  //   7- L=L+1
	  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, FBO[0]);
#ifdef USE_SINGLE_FBO
	  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_TEXTURE_2D, shadowMap[0], 0);
#endif
	  glReadBuffer(GL_COLOR_ATTACHMENT2_EXT);
	  glDrawBuffer(GL_COLOR_ATTACHMENT2_EXT);
	  // set the active PBO

	  if (useCUDA)
	    glBindBufferARB(GL_PIXEL_PACK_BUFFER_ARB, idLPBO);

	  while (number_layer < MAXLAYER)
	    {
	      //
	      //   3- Render layer L on texture L&1
	      //
	      //std::cout << "Render layer "<<number_layer<<std::endl;
	      glFramebufferTextureLayerEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT2_EXT, idLTexture, 0, number_layer);

	      glEnable(GL_TEXTURE_2D);
	      //If NO depth peeling has ever been performed yet, we initialize the depth map.
	      if (init)	  { glBindTexture( GL_TEXTURE_2D, shadowMap[2]);           }
	      else 	  { glBindTexture( GL_TEXTURE_2D, shadowMap[1-(number_layer&1)]); }

	      glUniform1iARB(static_cast< DepthPeelingShader*>(shader)->getLocation(),0);
	      glUniform1fARB(shader->getLocationResolution(), 1.0f); //resolution/(float)resolutionMax);

	      glBeginQueryARB(GL_SAMPLES_PASSED_ARB, occlusionQuery[(number_layer&1)]);

	      renderScene (_models, true);	// render scene culling first layer with old depth buffer

	      glEndQueryARB(GL_SAMPLES_PASSED_ARB);
	      glDisable(GL_TEXTURE_2D);
	      glBindTexture ( GL_TEXTURE_2D, 0 );

	      glReadPixels(0,0,
			   resolution,resolution,
			   GL_RGBA, GL_FLOAT, ((char*)NULL)+(number_layer*resolution*resolution*4*sizeof(float)));

	      //
	      //   4- IF L>=1 and layer L-1 is not empty: Start transfer of layer L-1 to PBO[(L-1)&1]
	      //

	      bool prev_layer_empty = false;
	      GLint npix;
	      if (number_layer >= 1)
		{
		  npix = -1;
		  glGetQueryObjectivARB(occlusionQuery[(number_layer-1)&1],GL_QUERY_RESULT_ARB,&npix);
		  //std::cout << "Layer "<<number_layer-1<<" contains "<<npix<<" pixels"<<std::endl;
		  prev_layer_empty = (npix==0);
		}
	      // flip the active FBO
#ifdef USE_SINGLE_FBO
	      glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_TEXTURE_2D, shadowMap[((number_layer+1)&1)], 0);
#else
	      glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, FBO[((number_layer+1)&1)]);
#endif

	      //
	      //   6- IF L>=1 and layer L-1 is empty: BREAK
	      //

	      if (number_layer>=1 && prev_layer_empty)
		{
		  --number_layer;
		  //std::cout << "BREAK"<<std::endl;
		  break;
		}

	      ++number_layer;
	      init = false;
	    }
	  glBindBufferARB(GL_PIXEL_PACK_BUFFER_ARB, 0);
	}

      glUseProgramObjectARB(0);
      glDisable(GL_VERTEX_PROGRAM_TWO_SIDE);

	glBindBufferARB(GL_PIXEL_PACK_BUFFER_ARB, 0);

      endRenderView();

      GLCHECK();
      // Now process the rendered layers

      if (timelog) t0 = timelog->endTime(t0, "collision/NarrowPhase/GPU/LDI", timeobj);

      glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
      glReadBuffer(GL_BACK);
      glDrawBuffer(GL_BACK);
      glViewport(viewport[0],viewport[1],viewport[2],viewport[3]);
      GLCHECK();

      if (timelog) t0 = timelog->endTime(t0, "collision/NarrowPhase/GPU", timeobj);

      result_id = -1;
      sofa::gpu::cuda::mycudaGLRegisterBufferObject(idLPBO);
       void* devLPtr = NULL;
       sofa::gpu::cuda::mycudaGLMapBufferObject(&devLPtr, idLPBO);

      //if (timelog) t0 = timelog->endTime(t0, "collision/NarrowPhase/GPU/CUDA/Begin", timeobj);

      CollisionVolume_count(devLPtr, scan_in.deviceWrite(), partial_sums.deviceWrite(), number_layer, resolution*resolution, BSIZE, _models.size()==1);
      const unsigned int* count_sums = partial_sums.hostRead();



      //if (timelog) t0 = timelog->endTime(t0, "collision/NarrowPhase/GPU/CUDA/Count", timeobj);
      int first_bloc = -1;
      int last_bloc = -2;
      int nblocs = 0;
      int ncolls = 0;
      for (unsigned int i=0;i<partial_sums.size();i++)
	{
	  if (count_sums[i] != 0)
	    {
	      if (first_bloc == -1) first_bloc = i;
	      last_bloc = i;
	      ++nblocs;
	      //vol += vol_sums[2*i+1];
	      ncolls += count_sums[i];
	    }
	}
      //if (nblocs) std::cout << "ncolls="<<ncolls<<" nblocs="<<nblocs<<" "<<first_bloc<<"-"<<last_bloc<<std::endl;
      if (nblocs > 0)
	{

#ifdef CHECK_SCAN
	  const unsigned int* p_scan_in = scan_in.hostRead();
	  int nerr = 0;
	  for (unsigned int b = 0; b < partial_sums.size();++b)
	    {
	      unsigned int ctotal = partial_sums[b];
	      unsigned int prev = 0;
	      for (unsigned int t = 0; t < BSIZE; ++t)
		{
		  unsigned int c = p_scan_in[b*BSIZE+t];
		  if (t==0 && c != 0)
		    { if (!(nerr++)) std::cerr << "doCollisionVolume("<<_models[0]->getName()<<"-"<<(_models.size()==1?std::string("SELF"):_models[1]->getName())<<"):\n";  std::cerr << "SCAN ERROR: bloc "<<b<<"/"<<partial_sums.size()<<" thread "<<t<<" count="<<c<<std::endl; }
		  if (c > ctotal)
		    { if (!(nerr++)) std::cerr << "doCollisionVolume("<<_models[0]->getName()<<"-"<<(_models.size()==1?std::string("SELF"):_models[1]->getName())<<"):\n";  std::cerr << "SCAN ERROR: bloc "<<b<<"/"<<partial_sums.size()<<" thread "<<t<<" count="<<c<<" while bloc total="<<ctotal<<std::endl; }
		  if (c < prev)
		    { if (!(nerr++)) std::cerr << "doCollisionVolume("<<_models[0]->getName()<<"-"<<(_models.size()==1?std::string("SELF"):_models[1]->getName())<<"):\n";  std::cerr << "SCAN ERROR: bloc "<<b<<"/"<<partial_sums.size()<<" thread "<<t<<" count="<<c<<" while thread "<<t-1<<" count="<<prev<<std::endl; }
		  prev = c;
		}
	    }
#endif
	  scan_out.fastResize(nblocs*2);
          unsigned int* bloc_input_pos = scan_out.hostWriteAt(0);
          unsigned int* bloc_output_pos = scan_out.hostWriteAt(nblocs);
	  int pos = 0;
	  int bloc = 0;
	  for (unsigned int i=0;i<partial_sums.size();i++)
	    {
	      if (count_sums[i] != 0)
		{
		  bloc_input_pos[bloc] = i*BSIZE;
		  bloc_output_pos[bloc] = pos;
		  ++bloc;
		  pos += count_sums[i];
		}
	    }
	  output_collisions.fastResize(ncolls*2*4);
	  //        CollisionVolume_write(((const char*)devLPtr)+(first_bloc*BSIZE*4*sizeof(float)), scan_out.deviceRead(first_bloc*BSIZE), output_collisions.deviceWrite(), number_layer, resolution*resolution, BSIZE, (last_bloc-first_bloc+1), _models.size()==1);
	  CollisionVolume_write(devLPtr, scan_in.deviceRead(), scan_out.deviceRead(0), scan_out.deviceRead(nblocs), output_collisions.deviceWrite(), number_layer, resolution*resolution, BSIZE, nblocs, _models.size()==1);
	  output_collisions.hostRead();
	  //if (timelog) t0 = timelog->endTime(t0, "collision/NarrowPhase/GPU/CUDA/Write", timeobj);
	  result_id = -1-ncolls;
	}
      else
        result_id = -1;
      
      sofa::gpu::cuda::mycudaGLUnmapBufferObject(idLPBO);      
      sofa::gpu::cuda::mycudaGLUnregisterBufferObject(idLPBO);
      if (timelog) t0 = timelog->endTime(t0, "collision/NarrowPhase/GPU/CUDA", timeobj);

    }

    void DepthPeelingUtility::beginReadCollisionVolume()
    {
    }

    const float* DepthPeelingUtility::readCollisionVolume()
    {
      return output_collisions.hostRead();
    }

    void DepthPeelingUtility::endReadCollisionVolume()
    {
    }

#endif
  }
}
}
