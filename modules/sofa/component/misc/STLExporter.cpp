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

#include "STLExporter.h"

#include <sstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>

#include <sofa/core/ObjectFactory.h>

#include <sofa/core/objectmodel/Event.h>
#include <sofa/simulation/common/AnimateBeginEvent.h>
#include <sofa/simulation/common/AnimateEndEvent.h>
#include <sofa/core/objectmodel/KeypressedEvent.h>
#include <sofa/core/objectmodel/KeyreleasedEvent.h>


namespace sofa
{

namespace component
{

namespace misc
{

SOFA_DECL_CLASS(STLExporter)

int STLExporterClass = core::RegisterObject("Read State vectors from file at each timestep")
        .add< STLExporter >();

STLExporter::STLExporter()
    : stepCounter(0)
    , stlFilename( initData(&stlFilename, "filename", "output STL file name"))
    , m_fileFormat( initData(&m_fileFormat, (bool)true, "binaryformat", "if true, save in binary format, otherwise in ascii"))
    , m_position( initData(&m_position, "position", "points coordinates"))
    , m_triangle( initData(&m_triangle, "triangle", "triangles indices"))
    , m_quad( initData(&m_quad, "quad", "quads indices"))
    , exportEveryNbSteps( initData(&exportEveryNbSteps, (unsigned int)0, "exportEveryNumberOfSteps", "export file only at specified number of steps (0=disable)"))
    , exportAtBegin( initData(&exportAtBegin, (bool)false, "exportAtBegin", "export file at the initialization"))
    , exportAtEnd( initData(&exportAtEnd, (bool)false, "exportAtEnd", "export file when the simulation is finished"))
{
}

STLExporter::~STLExporter()
{
    if (outfile)
        delete outfile;
}

void STLExporter::init()
{
    sofa::core::objectmodel::BaseContext* context = this->getContext();
    
    context->get(topology, sofa::core::objectmodel::BaseContext::Local);
    context->get(mstate, sofa::core::objectmodel::BaseContext::Local);
    context->get(vmodel, sofa::core::objectmodel::BaseContext::Local);
    
    // Test if the position has not been modified
    if(!m_position.isSet())
    {
        sofa::core::objectmodel::BaseData* pos = NULL;
        sofa::core::objectmodel::BaseData* tri = NULL;
        sofa::core::objectmodel::BaseData* qua = NULL;
        if(vmodel)
        {
            pos = vmodel->findField("position");
            tri = vmodel->findField("triangles");
            qua = vmodel->findField("quads");
        }
        else if(topology)
        {
            pos = topology->findField("position");
            tri = topology->findField("triangles");
            qua = topology->findField("quads");
        }
        else
        {
            serr << "STLExporter : error, STLExporter needs VisualModel or Topology" << sendl;
            return;
        }
        
        m_position.setParent(pos);
        m_triangle.setParent(tri);
        m_quad.setParent(qua);
    }

    // Activate the listening to the event in order to be able to export file at the nth-step
    if(exportEveryNbSteps.getValue() != 0)
        this->f_listening.setValue(true);
    
    nbFiles = 0;

}

void STLExporter::writeSTL()
{
    std::string filename = stlFilename.getFullPath();
    if (maxStep)
    {
        std::ostringstream oss;
        oss.width(5);
        oss.fill('0');
        oss << nbFiles;
        filename += oss.str();
    }
    filename += ".stl";

    outfile = new std::ofstream(filename.c_str());
    if( !outfile->is_open() )
    {
        serr << "Error creating file " << filename << sendl;
        delete outfile;
        outfile = NULL;
        return;
    }
    
    helper::ReadAccessor< Data< vector< sofa::component::topology::Triangle > > > triangleIndices = m_triangle;  
    helper::ReadAccessor< Data< vector< sofa::component::topology::Quad > > > quadIndices = m_quad;  
    helper::ReadAccessor<Data<defaulttype::Vec3Types::VecCoord> > positionIndices = m_position;
    
    vector< sofa::component::topology::Triangle > vecTri;
        
    if(positionIndices.empty())
    {
        serr << "STLExporter::writeSTL : error, no positions in topology" << sendl;
        return;
    }
    if(!triangleIndices.empty())
    {
        for(unsigned int i=0;i<triangleIndices.size();i++)
        {
            vecTri.push_back(triangleIndices[i]);
        }
    }
    else if(!quadIndices.empty())
    {
        sofa::component::topology::Triangle tri;
        for(unsigned int i=0;i<quadIndices.size();i++)
        {
            for(int j=0;j<3;j++)
            {
                tri[j] = quadIndices[i][j];
            }
            vecTri.push_back(tri);
            tri[0] = quadIndices[i][0];
            tri[1] = quadIndices[i][2];
            tri[2] = quadIndices[i][3];
            vecTri.push_back(tri);
        }
    }
    else
    {
        serr << "STLExporter::writeSTLBinary : error, neither triangles nor quads" << sendl;
        return;
    }
    
    /* Get number of facets */
    const int nbt = vecTri.size();
    
    // Sets the floatfield format flag for the str stream to fixed
    std::cout.precision(6);
    
    /* solid */
    *outfile << "solid Exported from Sofa" << std::endl;
    
    
    for(int i=0;i<nbt;i++)
    {
        /* normal */
        *outfile << "facet normal 0 0 0" << std::endl;
        *outfile << "outer loop" << std::endl;
        for (int j=0;j<3;j++)
        {
            /* vertices */
            *outfile << "vertex " << std::fixed << positionIndices[ vecTri[i][j] ] << std::endl;
        }
        *outfile << "endloop" << std::endl;
        *outfile << "endfacet" << std::endl;
    }
    
    /* endsolid */
    *outfile << "endsolid Exported from Sofa" << std::endl;
    
    outfile->close();
    std::cout << filename << " written" << std::endl;
    nbFiles++;
}

void STLExporter::writeSTLBinary()
{
    std::string filename = stlFilename.getFullPath();
    if (maxStep)
    {
        std::ostringstream oss;
        oss.width(5);
        oss.fill('0');
        oss << nbFiles;
        filename += oss.str();
    }
    filename += ".stl";
    
    outfile = new std::ofstream(filename.c_str(), std::ios::out | std::ios::binary);
    if( !outfile->is_open() )
    {
        serr << "Error creating file " << filename << sendl;
        delete outfile;
        outfile = NULL;
        return;
    }
    
    helper::ReadAccessor< Data< vector< sofa::component::topology::Triangle > > > triangleIndices = m_triangle;  
    helper::ReadAccessor< Data< vector< sofa::component::topology::Quad > > > quadIndices = m_quad;  
    helper::ReadAccessor<Data<defaulttype::Vec3Types::VecCoord> > positionIndices = m_position;
    
    vector< sofa::component::topology::Triangle > vecTri;
    
    if(positionIndices.empty())
    {
        serr << "STLExporter::writeSTLBinary : error, no positions in topology" << sendl;
        return;
    }
    if(!triangleIndices.empty())
    {
        for(unsigned int i=0;i<triangleIndices.size();i++)
        {
            vecTri.push_back(triangleIndices[i]);
        }
    }
    else if(!quadIndices.empty())
    {
        sofa::component::topology::Triangle tri;
        for(unsigned int i=0;i<quadIndices.size();i++)
        {
            for(int j=0;j<3;j++)
            {
                tri[j] = quadIndices[i][j];
            }
            vecTri.push_back(tri);
            tri[0] = quadIndices[i][0];
            tri[1] = quadIndices[i][2];
            tri[2] = quadIndices[i][3];
            vecTri.push_back(tri);
        }
    }
    else
    {
        serr << "STLExporter::writeSTLBinary : error, neither triangles nor quads" << sendl;
        return;
    }
    
    // Sets the floatfield format flag for the str stream to fixed
    std::cout.precision(6);
    
    /* Creating header file */
    char* buffer = new char [80];
    // Cleaning buffer
    for(int i=0;i<80;i++)
    {
        buffer[i]='\0';
    }
    strcpy(buffer, "Exported from Sofa");
    outfile->write(buffer,80);
        
    /* Number of facets */
    const unsigned int nbt = vecTri.size();
    outfile->write((char*)&nbt,4);
    
    // Parsing facets
    for(unsigned long i=0;i<nbt;i++)
    {
        /* normal */
        float nul = 0.; // normals are set to 0
        outfile->write((char*)&nul, 4);
        outfile->write((char*)&nul, 4);
        outfile->write((char*)&nul, 4);
        for (int j=0;j<3;j++)
        {
            /* vertices */
            float iOne = (float)positionIndices[ vecTri[i][j] ][0];
            float iTwo = (float)positionIndices[ vecTri[i][j] ][1];
            float iThree = (float)positionIndices[ vecTri[i][j] ][2];
            outfile->write( (char*)&iOne, 4);
            outfile->write( (char*)&iTwo, 4);
            outfile->write( (char*)&iThree, 4);
        }
        
        /* Attribute byte count */
        // attribute count is currently not used, it's garbage
        unsigned int zero = 0;
        outfile->write((char*)&zero, 2);
    }
    
    outfile->close();
    std::cout << filename << " written" << std::endl;
    nbFiles++;
}

void STLExporter::handleEvent(sofa::core::objectmodel::Event *event)
{
    if (sofa::core::objectmodel::KeypressedEvent* ev = dynamic_cast<sofa::core::objectmodel::KeypressedEvent*>(event))
    {
        std::cout << "key pressed " << std::endl;
        switch(ev->getKey())
        {

        case 'E':
        case 'e':
            if(m_fileFormat.getValue())
                writeSTLBinary();
            else
                writeSTL();
            break;

        case 'F':
        case 'f':
            break;
        }
    }


    if ( /*simulation::AnimateEndEvent* ev =*/  dynamic_cast<simulation::AnimateEndEvent*>(event))
    {
        maxStep = exportEveryNbSteps.getValue();
        if (maxStep == 0) return;

        stepCounter++;
        if(stepCounter % maxStep == 0)
        {
            if(m_fileFormat.getValue())
                writeSTLBinary();
            else
                writeSTL();
        }
    }
}

void STLExporter::cleanup()
{
    if (exportAtEnd.getValue())
        (m_fileFormat.getValue()) ? writeSTLBinary() : writeSTL();

}

void STLExporter::bwdInit()
{
    if (exportAtBegin.getValue())
        (m_fileFormat.getValue()) ? writeSTLBinary() : writeSTL();
}

} // namespace misc

} // namespace component

} // namespace sofa

