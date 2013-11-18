/*******************************************************************************
* CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
* version 0.1                                                                  *
* Copyright (C) 2009-2012, IGG Team, LSIIT, University of Strasbourg           *
*                                                                              *
* This library is free software; you can redistribute it and/or modify it      *
* under the terms of the GNU Lesser General Public License as published by the *
* Free Software Foundation; either version 2.1 of the License, or (at your     *
* option) any later version.                                                   *
*                                                                              *
* This library is distributed in the hope that it will be useful, but WITHOUT  *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
* for more details.                                                            *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this library; if not, write to the Free Software Foundation,      *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
*                                                                              *
* Web site: http://cgogn.unistra.fr/                                           *
* Contact information: cgogn@unistra.fr                                        *
*                                                                              *
*******************************************************************************/

#include "Utils/Shaders/shaderSimpleColor.h"

namespace CGoGN
{

namespace Utils
{

#include "shaderSimpleColor.vert"
#include "shaderSimpleColor.frag"

//std::string ShaderSimpleColor::vertexShaderText =
//		"ATTRIBUTE vec3 VertexPosition, VertexNormal;\n"
//		"uniform mat4 ModelViewProjectionMatrix;\n"
////		"INVARIANT_POS;\n"
//		"void main ()\n"
//		"{\n"
//		"	gl_Position = ModelViewProjectionMatrix * vec4 (VertexPosition, 1.0);\n"
//		"}";
//
//
//std::string ShaderSimpleColor::fragmentShaderText =
//		"PRECISON;\n"
//		"uniform vec4 color;\n"
//		"FRAG_OUT_DEF;\n"
//		"void main()\n"
//		"{\n"
//		"	gl_FragColor = color;\n"
//		"}";


ShaderSimpleColor::ShaderSimpleColor(bool black_is_transparent)
{
	m_nameVS = "ShaderSimpleColor_vs";
	m_nameFS = "ShaderSimpleColor_fs";
	m_nameGS = "ShaderSimpleColor_gs";

	// chose GL defines (2 or 3)
	// and compile shaders
	std::string glxvert(*GLSLShader::DEFINES_GL);
	glxvert.append(vertexShaderText);

	std::string glxfrag(*GLSLShader::DEFINES_GL);
	if (black_is_transparent)
		glxfrag.append("#define BLACK_TRANSPARENCY 1\n");
	glxfrag.append(fragmentShaderText);

	loadShadersFromMemory(glxvert.c_str(), glxfrag.c_str());

	*m_unif_color = glGetUniformLocation(this->program_handler(),"color");

	//Default values
	Geom::Vec4f color(0.1f, 0.9f, 0.1f, 0.0f);
	setColor(color);
}

void ShaderSimpleColor::setColor(const Geom::Vec4f& color)
{
	m_color = color;
	bind();
	glUniform4fv(*m_unif_color, 1, color.data());
	unbind();
}

unsigned int ShaderSimpleColor::setAttributePosition(VBO* vbo)
{
	m_vboPos = vbo;
	bind();
	unsigned int id = bindVA_VBO("VertexPosition", vbo);
	unbind();
	return id;
}

void ShaderSimpleColor::restoreUniformsAttribs()
{
	*m_unif_color = glGetUniformLocation(this->program_handler(), "color");
	bind();
	glUniform4fv(*m_unif_color, 1, m_color.data());
	bindVA_VBO("VertexPosition", m_vboPos);
	unbind();
}

} // namespace Utils

} // namespace CGoGN
