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
#include <sofa/core/ObjectFactory.h>
#include <sofa/component/loader/MeshXspLoader.h>
#include <sofa/core/visual/VisualParams.h>
#include <iostream>

namespace sofa
{

namespace component
{

namespace loader
{

using namespace sofa::defaulttype;

SOFA_DECL_CLASS(MeshXspLoader)

int MeshXspLoaderClass = core::RegisterObject("Specific mesh loader for Xsp file format.")
        .add< MeshXspLoader >()
        ;

MeshXspLoader::MeshXspLoader() : MeshLoader()
    , gravity(initData(&gravity,"gravity","Gravity coordinates loaded in this mesh."))
    , viscosity(initData(&viscosity,"viscosity","viscosity values loaded in this mesh."))
{
    gravity.setPersistent(false);
    viscosity.setPersistent(false);
}


bool MeshXspLoader::load()
{

    sout << "Loading Xsp file: " << m_filename << sendl;

    std::string cmd;
    bool fileRead = false;

    // -- Loading file
    const char* filename = m_filename.getFullPath().c_str();
    std::ifstream file(filename);

    if (!file.good())
    {
        serr << "Error: MeshXspLoader: Cannot read file '" << m_filename << "'." << sendl;
        return false;
    }


    // -- Check first line.
    file >> cmd;

    // -- Reading file version
    if (cmd == "Xsp")
    {
        float version = 0.0f;
        file >> version;

        if (version == 3.0)
            fileRead = readXsp(file, false);
        else if (version == 4.0)
            fileRead = readXsp(file, true);
    }
    else
    {
        serr << "Error: MeshXspLoader: File '" << m_filename << "' finally appears not to be a Xsp file." << sendl;
        file.close();
        return false;

    }

    file.close();
    return fileRead;
}



bool MeshXspLoader::readXsp (std::ifstream &file, bool vector_spring)
{
    sout << "Reading Xsp file: " << vector_spring << sendl;


    std::string cmd;
    int npoints = 0;
    int nlines = 0;

    int totalNumMasses;
    int totalNumSprings;

    file >> cmd;

    // then find out number of masses and springs
    if (cmd == "numm")
    {
        file >> totalNumMasses;
        npoints=totalNumMasses;
    }

    if (cmd=="nums")
    {
        file >> totalNumSprings;
        nlines=totalNumSprings;
    }


    helper::vector<sofa::defaulttype::Vector3>& my_positions = *(positions.beginEdit());
    helper::vector<Edge >& my_edges = *(edges.beginEdit());

    helper::vector <Vector3>& my_gravity = *(gravity.beginEdit());
    helper::vector <double>& my_viscosity = *(viscosity.beginEdit());

    while (!file.eof())
    {
        file  >> cmd;
        if (cmd=="mass")
        {
            int index;
            char location;
            double px,py,pz,vx,vy,vz,mass=0.0,elastic=0.0;
            bool fixed=false;
            file >> index >> location >> px >> py >> pz >> vx >> vy >> vz >> mass >> elastic;
            if (mass < 0)
            {
                // fixed point initialization
                mass = -mass;
                fixed = true;
            }
            my_positions.push_back(Vector3(px, py, pz));
        }
        else if (cmd=="lspg")	// linear springs connector
        {
            int	index;
            Edge m;
            double ks=0.0,kd=0.0,initpos=-1;
            double restx=0.0,resty=0.0,restz=0.0;
            if (vector_spring)
                file >> index >> m[0] >> m[1] >> ks >> kd >> initpos >> restx >> resty >> restz;
            else
                file >> index >> m[0] >> m[1] >> ks >> kd >> initpos;
            --m[0];
            --m[1];

            addEdge(&my_edges, m);
        }
        else if (cmd == "grav")
        {
            double gx,gy,gz;
            file >> gx >> gy >> gz;
            my_gravity.push_back(Vector3(gx, gy, gz));
        }
        else if (cmd == "visc")
        {
            double visc;
            file >> visc;
            my_viscosity.push_back (visc);
        }
        else if (cmd == "step")
        {
        }
        else if (cmd == "frce")
        {
        }
        else if (cmd[0] == '#')	// it's a comment
        {
        }
        else		// it's an unknown keyword
        {
            std::cerr << "Unknown MassSpring keyword:" << cmd;
            positions.endEdit();
            edges.endEdit();
            gravity.endEdit();
            viscosity.endEdit();
            return false;
        }
    }

    positions.endEdit();
    edges.endEdit();
    gravity.endEdit();
    viscosity.endEdit();

    return true;
}

} // namespace loader

} // namespace component

} // namespace sofa

