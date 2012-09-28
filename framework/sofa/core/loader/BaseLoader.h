/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *
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
*                              SOFA :: Framework                              *
*                                                                             *
* Authors: The SOFA Team (see Authors.txt)                                    *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_CORE_LOADER_BASELOADER_H
#define SOFA_CORE_LOADER_BASELOADER_H

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/objectmodel/DataFileName.h>
#include <sofa/core/objectmodel/BaseObjectDescription.h>

//#include <stdlib.h>
#include <string>
#include <cstring>
//#include <iostream>
//#include <sofa/helper/vector.h>

namespace sofa
{

namespace core
{

namespace loader
{

bool SOFA_CORE_API canLoad(const char* filename);

class BaseLoader : public virtual objectmodel::BaseObject
{
public:
    SOFA_ABSTRACT_CLASS(BaseLoader, objectmodel::BaseObject);
protected:
    ///Constructor
    BaseLoader(): m_filename(initData(&m_filename,"filename","Filename of the object")) {}

    ///Destructor
    virtual ~BaseLoader() { }
public:
    // virtual bool load(const char *filename) = 0;

    virtual bool load() = 0;

    virtual void parse(sofa::core::objectmodel::BaseObjectDescription *arg)
    {
        objectmodel::BaseObject::parse(arg);
        if (canLoad())
            load(/*m_filename.getFullPath().c_str()*/);
        else
            sout << "Doing nothing" << sendl;
    }


    virtual bool canLoad()
    {
        std::string cmd;

        // -- Check filename field:
        if(m_filename.getValue() == "")
        {
            serr << "Error: MeshLoader: No file name given." << sendl;
            return false;
        }


        // -- Check if file exist:
        const char* filename = m_filename.getFullPath().c_str();
        std::string sfilename (filename);

        if (!sofa::helper::system::DataRepository.findFile(sfilename))
        {
            serr << "Error: MeshLoader: File '" << m_filename << "' not found. " << sendl;
            return false;
        }

        std::ifstream file(filename);

        // -- Check if file is readable:
        if (!file.good())
        {
            serr << "Error: MeshLoader: Cannot read file '" << m_filename << "'." << sendl;
            return false;
        }

        // -- Step 2.2: Check first line.
        file >> cmd;
        if (cmd.empty())
        {
            serr << "Error: MeshLoader: Cannot read first line in file '" << m_filename << "'." << sendl;
            file.close();
            return false;
        }

        file.close();
        return true;
    };


    void setFilename(std::string f)
    {
        m_filename.setValue(f);
    }

    const std::string &getFilename()
    {
        return m_filename.getValue();
    }


    static void skipToEOL(FILE* f)
    {
        int ch;
        while ((ch = fgetc(f)) != EOF && ch != '\n') ;
    }


    static bool readLine(char* buf, int size, FILE* f)
    {
        buf[0] = '\0';
        if (fgets(buf, size, f) == NULL)
            return false;
        if ((int)strlen(buf)==size-1 && buf[size-1] != '\n')
            skipToEOL(f);
        return true;
    }

    /*
    template<class T>
    static bool canCreate(T*& obj, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg)
    {

      return BaseMapping::canCreate(obj, context, arg);
    }
    */
    sofa::core::objectmodel::DataFileName m_filename;

};



} // namespace loader

} // namespace core

} // namespace sofa

#endif
