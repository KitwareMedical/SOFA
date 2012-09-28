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
#ifndef SOFA_GPU_OPENCL_OPENCLPROGRAM_H
#define SOFA_GPU_OPENCL_OPENCLPROGRAM_H

#include "myopencl.h"
#include "OpenCLCommon.h"

#include <string>
#include <cstdio>
#include <stdlib.h>
#include <iostream>
#include <map>

namespace sofa
{

namespace gpu
{

namespace opencl
{

class OpenCLProgram
{
    cl_program _program;
    std::string _filename;
    std::string  _source,_inputs,_operation;
    std::map<std::string,std::string> _types;

public:
    OpenCLProgram();
    OpenCLProgram(std::string filename, std::string srcPrefix = stringBSIZE);
    OpenCLProgram(std::string *source);
    OpenCLProgram(std::string *source,std::map<std::string,std::string> *types);
    OpenCLProgram(std::string filename, std::string srcPrefix,std::map<std::string,std::string> *types);
    OpenCLProgram(std::string *source,std::string *operation);
    OpenCLProgram(std::string *source,std::string *operation, std::map<std::string,std::string> *types);

    ~OpenCLProgram();


    void setSource(std::string  s) {_source=s;}

    void setSourceFile(std::string filename, std::string srcPrefix = stringBSIZE);

    void setInputs(std::string s) {_inputs=s;}
    void setTypes(std::map<std::string,std::string> types) {_types=types;}


    void addMacro(std::string name,std::string method);

    void addMacros(std::string* sources,std::string option);


    void createProgram();

    std::string createTypes(std::map<std::string,std::string> *types);

    void createProgram(std::string * s);

    void createProgram(std::string *source,std::map<std::string,std::string> *types);

    void createProgram(std::string *source,std::string *operation);

    void createProgram(std::string *source,std::string *operation, std::map<std::string,std::string> *types);

    void* program() {return _program;}

    void buildProgram();

    void buildProgram(char* flags);

    /// Load a source kernel file, appending its content into dest
    static bool loadSource(const std::string& file_source, std::string* dest);

    std::string buildLog(int device);

    std::string source();

    std::string sourceLog();

};

}

}

}

#endif // SOFA_GPU_OPENCL_OPENCLPROGRAM_H
