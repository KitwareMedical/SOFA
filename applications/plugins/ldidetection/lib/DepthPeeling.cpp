/*************************************************************************
* Shaders.h                                                              *
* Auteurs        : Christian Boucheny, Sebastien Barbier                 *
* Version        : 2.0                                                   *
* Date           : 06.06.07                                              *
* Specifications : Load a shader on the GPU                              *
**************************************************************************/
// #include <sofa/component/collision/DepthPeeling.h>
#include "DepthPeeling.h"
#include "DepthPeelingUtility.h"

//#include "headers.h"
#include <fstream>
#include <iostream>
#include <sstream>

/**************************************
              SHADERS
***************************************/

namespace sofa
{

namespace component
{

namespace collision
{

GLhandleARB* GLSL_Shader::compileShader(GLint target, const std::string& source)
{
    const char* stype = "";
    if (target == GL_VERTEX_SHADER_ARB) stype = "vertex";
    else if (target == GL_FRAGMENT_SHADER_ARB) stype = "fragment";
#ifdef GL_GEOMETRY_SHADER_EXT
    else if (target == GL_GEOMETRY_SHADER_EXT) stype = "geometry";
#endif

    GLhandleARB shader = glCreateShaderObjectARB(target);

    const char* src = source.c_str();
    
    glShaderSourceARB(shader, 1, &src, NULL);
    
    glCompileShaderARB(shader);

    GLint compiled = 0, length = 0, laux = 0;
    glGetObjectParameterivARB(shader, GL_OBJECT_COMPILE_STATUS_ARB, &compiled);
    if (!compiled) std::cerr << "ERROR: Compilation of "<<stype<<" shader failed:\n"<<src<<std::endl;
    //    else std::cout << "Compilation of "<<stype<<" shader OK" << std::endl;
    glGetObjectParameterivARB(shader, GL_OBJECT_INFO_LOG_LENGTH_ARB, &length);
    if (length)
    {
        GLcharARB *logString = (GLcharARB *)malloc(length * sizeof(GLcharARB));
        glGetInfoLogARB(shader, length, &laux, logString);
        std::cerr << logString << std::endl;
        free(logString);
    }
    if (compiled)
        return new GLhandleARB(shader);
    else
    {
        return NULL;
    }
}

GLhandleARB *GLSL_Shader::createProgram()
{
    return new GLhandleARB(glCreateProgramObjectARB());
}

bool GLSL_Shader::attachShader(GLhandleARB *program, GLhandleARB * shader)
{
    if (program == NULL || shader == NULL) return false;
    glAttachObjectARB(*program, *shader);
    return true;
}

bool GLSL_Shader::linkProgram(GLhandleARB *program)
{
    if (program == NULL) return false;
    glLinkProgramARB(*program);

    GLint compiled = 0, length = 0, laux = 0;
    glGetObjectParameterivARB(*program, GL_OBJECT_LINK_STATUS_ARB, &compiled);
    if (!compiled) std::cerr << "ERROR: Link of program shader failed:\n"<<std::endl;
    //    else std::cout << "Link of program shader OK" << std::endl;
    glGetObjectParameterivARB(*program, GL_OBJECT_INFO_LOG_LENGTH_ARB, &length);
    if (length)
    {
        GLcharARB *logString = (GLcharARB *)malloc(length * sizeof(GLcharARB));
        glGetInfoLogARB(*program, length, &laux, logString);
        std::cerr << logString << std::endl;
        free(logString);
    }
    return (bool) compiled;
}

GLSL_Shader::GLSL_Shader()
: program_shader(NULL)
, vertex_shader(NULL)
, geometry_shader(NULL)
, fragment_shader(NULL)
, geometry_input_type(-1)
, geometry_output_type(-1)
, geometry_vertices_out(-1)
{
}

bool GLSL_Shader::init(std::string vshader_name, std::string gshader_name, std::string fshader_name)
{
    program_shader = createProgram();

    bool result = true;

    vertex_shader = compileShader(GL_VERTEX_SHADER_ARB, vshader_name);
    result &= attachShader(program_shader, vertex_shader);

    if (!gshader_name.empty())
    {
#ifdef GL_GEOMETRY_SHADER_EXT
        geometry_shader = compileShader(GL_GEOMETRY_SHADER_EXT, gshader_name);
        if (geometry_input_type != -1) glProgramParameteriEXT(*program_shader, GL_GEOMETRY_INPUT_TYPE_EXT, geometry_input_type );
        if (geometry_output_type != -1) glProgramParameteriEXT(*program_shader, GL_GEOMETRY_OUTPUT_TYPE_EXT, geometry_output_type );
        if (geometry_vertices_out != -1) glProgramParameteriEXT(*program_shader, GL_GEOMETRY_VERTICES_OUT_EXT, geometry_vertices_out );
#endif
        result &= attachShader(program_shader, geometry_shader);
    }

    fragment_shader = compileShader(GL_FRAGMENT_SHADER_ARB, fshader_name);
    result &= attachShader(program_shader, fragment_shader);

    if (result)
        result &= linkProgram(program_shader);

    ready = result;

    return result;
}

//**********************************************************************************
//Constructors

PeelingShader::PeelingShader()
: location_BB(0)
, location_RegGrid(0)
, location_resolution(0)
{
}

bool PeelingShader::init(std::string vshader_name, std::string gshader_name, std::string fshader_name)
{
    bool r = GLSL_Shader::init(vshader_name, gshader_name, fshader_name);
    if (r)
    {
        glUseProgramObjectARB(*program_shader);
        
        location_BB       = glGetUniformLocationARB(*program_shader, "BB");
        location_RegGrid    = glGetUniformLocationARB(*program_shader, "RegGrid");
	location_IdModels = glGetUniformLocationARB(*program_shader, "IdModels");
	location_resolution = glGetUniformLocationARB(*program_shader, "resolution");
        
        glUniform1fARB (location_resolution,1.0f);
        
        glUseProgramObjectARB(0);
    }
    return r;
}

DepthPeelingShader::DepthPeelingShader()
: location_depthMap(0)
{
}

bool DepthPeelingShader::init(std::string vshader_name, std::string gshader_name, std::string fshader_name)
{
    bool r = PeelingShader::init(vshader_name, gshader_name, fshader_name);
    if (r)
    {
        glUseProgramObjectARB(*program_shader);
        
        location_depthMap = glGetUniformLocationARB(*program_shader, "depthMap");
        
        glUniform1iARB (location_depthMap,0);
        
        glUseProgramObjectARB(0);
    }
    return r;
}

bool DepthPeelingShader::load()
{
    return init(
        std::string(
          "varying	vec3	   pt;\n\
	  varying vec4  id;\n\
	  uniform vec3 BB[8];\n\
	  uniform float RegGrid;\n\
	  uniform float IdModels;\n\
	  uniform float resolution;\n\
	  void	main ()\n\
	  {\n\
	  float a = gl_Vertex.x;\n\
	  float b = gl_Vertex.y;\n\
	  float c = gl_Vertex.z;\n\
	  vec4 bary;\n\
	  bary.xyz = (1.-a)*(1.-b)*(1.-c)*BB[0] + a*(1.-b)*(1.-c)*BB[1] + (1.-a)*b*(1.-c)*BB[2] + a*b*(1.-c)*BB[3] + (1.-a)*(1.-b)*c*BB[4] + a*(1.-b)*c*BB[5] + (1.-a)*b*c*BB[6] + a*b*c*BB[7];\n\
	  bary             = gl_ModelViewProjectionMatrix*vec4(bary.xyz,1.);\n\
	  vec4 pos = ftransform();\n\
	  pos = RegGrid*pos+(1. - RegGrid)*bary;\n\
	  pt = (pos.xyz/pos.w + vec3(1.0))*0.5;\n\
	  pt = pt*vec3(resolution, resolution, 1.)+ vec3(0., 0., -0.00001);\n\
          //pt              = (pos.xyz / pos.w + vec3 ( 1.0 )) * 0.5;\n\
	  id = vec4(IdModels,gl_Color.g,gl_Color.b,gl_Color.a);\n\
          gl_FrontColor = vec4(0.0,0.0,0.5,0.0);\n\
          gl_BackColor = vec4(0.0,0.0,0.0,0.0);\n\
	  gl_Position     = pos;\n\
        }"),
        std::string(),
        std::string(
"#extension GL_ARB_draw_buffers: enable\n\
varying	vec3 pt;\n\
varying vec4 id;\n\
uniform sampler2DShadow	depthMap;\n\
void main (void){\n\
//const vec3 eps   = vec3 ( 0, 0, -0.00001 );\n\
vec4       depth = shadow2D  ( depthMap, pt);\n\
if ( depth.x > 0.0 ) {		// reject previously visible fragment\n\
discard;\n\
}\n\
gl_FragData[0] = id+gl_Color;\n\
        }"));
}
    
    bool DepthPeelingShaderRenderPair::load()
    {
        return init(
                    std::string(
                                "varying	vec3	   pt;\n\
                                varying vec4  id;\n\
                                uniform vec3 BB[8];\n\
                                uniform float RegGrid;\n\
                                uniform float IdModels;\n\
	                        uniform float resolution;\n\
                                void	main ()\n\
                                {\n\
                                  vec4 pos;\n\
                                  if (RegGrid == 0.0) { \n\
                                    float a = gl_Vertex.x;\n\
                                    float b = gl_Vertex.y;\n\
                                    float c = gl_Vertex.z;\n\
                                    vec4 bary;\n\
                                    bary.xyz = (1.-a)*(1.-b)*(1.-c)*BB[0] + a*(1.-b)*(1.-c)*BB[1] + (1.-a)*b*(1.-c)*BB[2] + a*b*(1.-c)*BB[3] + (1.-a)*(1.-b)*c*BB[4] + a*(1.-b)*c*BB[5] + (1.-a)*b*c*BB[6] + a*b*c*BB[7];\n\
                                    pos             = gl_ModelViewProjectionMatrix*vec4(bary.xyz,1.);\n\
                                  } else { \n\
                                    pos = ftransform();\n\
                                  } \n\
                                  pt = (pos.xyz/pos.w + vec3(1.0))*0.5;\n\
	                          pt = pt*vec3(resolution, resolution, 1.)+ vec3(0., 0., -0.00001);\n\
	                           //pt              = (pos.xyz /* / pos.w */ ) * 0.5 + vec3 ( 0.5, 0.5, 0.5 -0.00001 );\n\
                                  id = vec4(pt.z,gl_Color.g,gl_Color.b,gl_Color.a);\n\
                                  gl_FrontColor = vec4(0.0,IdModels*0.5,0.5,0.0);\n\
                                  gl_BackColor = vec4(0.0,IdModels*0.5,0.0,0.0);\n\
                                  gl_Position     = pos;\n\
                                }"),
                    std::string(),
                    std::string(
                                "#extension GL_ARB_draw_buffers: enable\n\
                                varying	vec3 pt;\n\
                                varying vec4 id;\n\
                                uniform sampler2DShadow	depthMap;\n\
                                void main (void){\n\
                                vec4       depth = shadow2D  ( depthMap, pt  );\n\
                                if ( depth.x > 0.0 ) {		// reject previously visible fragment\n\
                                discard;\n\
                                }\n\
                                gl_FragData[0] = id+gl_Color;\n\
                                }"));
    }




template<class T>
static std::string toString(const T& t)
{
    std::ostringstream o;
    o << t;
    return o.str();
}

} //collision
} //component
} //sofa

//gl_FragData[0] = vec4(identity.rg, gl_FragCoord.z,1);
