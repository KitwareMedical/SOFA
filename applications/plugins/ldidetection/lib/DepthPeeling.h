/*************************************************************************
* Shaders.h                                                              *
* Auteurs        : Christian Boucheny, Sebastien Barbier                 *
* Version        : 2.0                                                   *
* Date           : 06.06.07                                              *
* Specifications : Load a shader on the GPU                              *
**************************************************************************/
#ifndef _DEPTHPEELING_H_
#define _DEPTHPEELING_H_

#include "initldidetection.h"
#ifdef SOFA_HAVE_GLEW 
#include <GL/glew.h>
#else
#include <sofa/helper/system/gl.h>
#include <sofa/helper/system/glu.h>
#endif
#include <vector>
#include <fstream>

namespace sofa
{

namespace component
{

namespace collision
{

enum VBO_tag {POSITION,IDENTITY};


//*******************************************************************************************
//GLSL basic shader
class SOFA_LDIDETECTION_API GLSL_Shader
{
public:
    GLSL_Shader();
    virtual ~GLSL_Shader(){}

    virtual bool init(std::string vshader_name, std::string gshader_name, std::string fshader_name);

    GLhandleARB *getShader(){ return program_shader; }
    
    bool isReady() { return ready; }

    GLint getGeometryInputType() { return geometry_input_type; }
    void  setGeometryInputType(GLint v) { geometry_input_type = v; }

    GLint getGeometryOutputType() { return geometry_output_type; }
    void  setGeometryOutputType(GLint v) { geometry_output_type = v; }

    GLint getGeometryVerticesOut() { return geometry_vertices_out; }
    void  setGeometryVerticesOut(GLint v) { geometry_vertices_out = v; }
    
    virtual int getLocationBB() const {return 0; }
    virtual int getLocationRegGrid() const {return 0; }
    virtual int getLocationIdModels()const{return 0;};

protected:

    GLhandleARB *compileShader(GLint target, const std::string& source);

    GLhandleARB *createProgram();

    bool attachShader(GLhandleARB *program, GLhandleARB * shader);

    bool linkProgram(GLhandleARB *program);

    std::string getStringFromFile( const char* filename);

    GLhandleARB *program_shader;
    GLhandleARB *vertex_shader;
    GLhandleARB *geometry_shader;
    GLhandleARB *fragment_shader;
    
    bool ready;
    
    GLint geometry_input_type;
    GLint geometry_output_type;
    GLint geometry_vertices_out;

};

      //*******************************************************************************************
      //Specialized Shader for Peeling purpose
      class SOFA_LDIDETECTION_API PeelingShader: public GLSL_Shader
{
    public:
        PeelingShader();
        ~PeelingShader(){}
        
        virtual bool init(std::string vshader_name, std::string gshader_name, std::string fshader_name);

        int getLocationBB() const {return location_BB; }
        int getLocationRegGrid() const {return location_RegGrid; }
        int getLocationIdModels()const{return location_IdModels;};
        int getLocationResolution() const {return location_resolution;}

    protected:
	
        int location_BB;
        int location_RegGrid;
        int location_IdModels;
        int location_resolution;
};

      //*******************************************************************************************
      //Specialized Shader for DepthPeeling purpose
      class SOFA_LDIDETECTION_API DepthPeelingShader: public PeelingShader
      {
      public:
          DepthPeelingShader();
	~DepthPeelingShader(){}
        
        virtual bool load();
        virtual bool init(std::string vshader_name, std::string gshader_name, std::string fshader_name);

	int getLocation() const {return location_depthMap; }

      protected:
	
	int location_depthMap;
      };
    
    //*******************************************************************************************
    //Specialized Shader for DepthPeeling purpose, when at most 2 models are rendered
    class SOFA_LDIDETECTION_API DepthPeelingShaderRenderPair: public DepthPeelingShader
    {
    public:
        virtual bool load();
    };
    
    
    
} //collision
} //component
} //sofa


#endif
