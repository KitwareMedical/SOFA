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
#include <sofa/component/loader/MeshGmshLoader.h>
#include <sofa/core/visual/VisualParams.h>
#include <iostream>

namespace sofa
{

namespace component
{

namespace loader
{

using namespace sofa::defaulttype;

SOFA_DECL_CLASS(MeshGmshLoader)

int MeshGmshLoaderClass = core::RegisterObject("Specific mesh loader for Gmsh file format.")
        .add< MeshGmshLoader >()
        ;

bool MeshGmshLoader::load()
{
    sout << "Loading Gmsh file: " << m_filename << sendl;

    std::string cmd;
    bool fileRead = false;
    unsigned int gmshFormat = 0;

    // -- Loading file
    const char* filename = m_filename.getFullPath().c_str();
    std::ifstream file(filename);

    if (!file.good())
    {
        serr << "Error: MeshGmshLoader: Cannot read file '" << m_filename << "'." << sendl;
        return false;
    }


    // -- Looking for Gmsh version of this file.
    std::getline(file, cmd); //Version
    std::istringstream versionReader(cmd);
    std::string version;
    versionReader >> version;
    if (version == "$MeshFormat") // Reading gmsh 2.0 file
    {
        gmshFormat = 2;
//      sout << "Gmsh format 2.0" << sendl;
        std::string line;
        std::getline(file, line); // we don't care about this line (2 0 8)
        std::getline(file, cmd); // end Version
        std::istringstream endMeshReader(cmd);
        std::string endMesh;
        endMeshReader >> endMesh;

        if (endMesh != std::string("$EndMeshFormat") ) // it should end with $EndMeshFormat
        {
            serr << "Closing File" << sendl;
            file.close();
            return false;
        }
        else
        {
            std::getline(file, cmd); // First Command
        }
    }
    else
    {
        //sout << "Gmsh format 1.0" << sendl;
        gmshFormat = 1;
    }


    std::istringstream nodeReader(cmd);
    std::string node;
    nodeReader >> node;
    // -- Reading file
    if (node == "$NOD" || node == "$Nodes") // Gmsh format
    {
        fileRead = readGmsh(file, gmshFormat);
    }
    else //if it enter this "else", it means there is a problem before in the factory or in canLoad()
    {
        serr << "Error: MeshGmshLoader: File '" << m_filename << "' finally appears not to be a Gmsh file." << sendl;
        file.close();
        return false;
    }

    return fileRead;
}



bool MeshGmshLoader::readGmsh(std::ifstream &file, const unsigned int gmshFormat)
{
    sout << "Reading Gmsh file: " << gmshFormat << sendl;

    std::string cmd;
    std::string line;

    unsigned int npoints = 0;
    unsigned int nelems = 0;

    unsigned int nlines = 0;
    unsigned int ntris = 0;
    unsigned int nquads = 0;
    unsigned int ntetrahedra = 0;
    unsigned int ncubes = 0;

    // --- Loading Vertices ---
    file >> npoints; //nb points

    helper::vector<sofa::defaulttype::Vector3>& my_positions = *(positions.beginEdit());

    std::vector<unsigned int> pmap; // map for reordering vertices possibly not well sorted
    for (unsigned int i=0; i<npoints; ++i)
    {
        unsigned int index = i;
        double x,y,z;
        file >> index >> x >> y >> z;

        my_positions.push_back(Vector3(x, y, z));

        if (pmap.size() <= index)
            pmap.resize(index+1);

        pmap[index] = i; // In case of hole or switch
        //sout << "pmap[" << index << "] = " << pmap[index] << sendl;
    }
    positions.endEdit();

    file >> cmd;
    if (cmd != "$ENDNOD" && cmd != "$EndNodes")
    {
        serr << "Error: MeshGmshLoader: '$ENDNOD' or '$EndNodes' expected, found '" << cmd << "'" << sendl;
        file.close();
        return false;
    }


    // --- Loading Elements ---
    file >> cmd;
    if (cmd != "$ELM" && cmd != "$Elements")
    {
        serr << "Error: MeshGmshLoader: '$ELM' or '$Elements' expected, found '" << cmd << "'" << sendl;
        file.close();
        return false;
    }

    file >> nelems; //Loading number of Element

    helper::vector<Edge>& my_edges = *(edges.beginEdit());
    helper::vector<Triangle>& my_triangles = *(triangles.beginEdit());
    helper::vector<Quad>& my_quads = *(quads.beginEdit());
    helper::vector<Tetrahedron>& my_tetrahedra = *(tetrahedra.beginEdit());
    helper::vector<Hexahedron>& my_hexahedra = *(hexahedra.beginEdit());


    for (unsigned int i=0; i<nelems; ++i) // for each elem
    {
        int index, etype, rphys, relem, nnodes, ntags, tag;

        if (gmshFormat==1)
        {
            // version 1.0 format is
            // elm-number elm-type reg-phys reg-elem number-of-nodes <node-number-list ...>
            file >> index >> etype >> rphys >> relem >> nnodes;
        }
        else if (gmshFormat == 2)
        {
            // version 2.0 format is
            // elm-number elm-type number-of-tags < tag > ... node-number-list
            file >> index >> etype >> ntags;

            for (int t=0; t<ntags; t++)
            {
                file >> tag;
            }


            switch (etype)
            {
            case 15: //point
                nnodes = 1;
                break;
            case 1: // Line
                nnodes = 2;
                break;
            case 2: // Triangle
                nnodes = 3;
                break;
            case 3: // Quad
                nnodes = 4;
                break;
            case 4: // Tetra
                nnodes = 4;
                break;
            case 5: // Hexa
                nnodes = 8;
                break;
            default:
                serr << "Error: MeshGmshLoader: Elements of type 1, 2, 3, 4, 5, or 6 expected. Element of type " << etype << " found." << sendl;
                nnodes = 0;
            }
        }
        //store real index of node and not line index

        helper::vector <unsigned int> nodes;
        nodes.resize (nnodes);

        for (int n=0; n<nnodes; ++n)
        {
            int t = 0;
            file >> t;
            nodes[n] = (((unsigned int)t)<pmap.size())?pmap[t]:0;
            //sout << "nodes[" << n << "] = " << nodes[n] << sendl;
        }
        Hexahedron hexa;
        switch (etype)
        {
        case 1: // Line
            addEdge(&my_edges, Edge(nodes[0], nodes[1]));
            ++nlines;
            break;
        case 2: // Triangle
            addTriangle(&my_triangles, Triangle(nodes[0], nodes[1], nodes[2]));
            ++ntris;
            break;
        case 3: // Quad
            addQuad(&my_quads, Quad(nodes[0], nodes[1], nodes[2], nodes[3]));
            ++nquads;
            break;
        case 4: // Tetra
            addTetrahedron(&my_tetrahedra, Tetrahedron(nodes[0], nodes[1], nodes[2], nodes[3]));
            ++ntetrahedra;
            break;
        case 5: // Hexa

            for (unsigned int n=0; n<8; n++)
                hexa[n] = nodes[n];
            addHexahedron(&my_hexahedra,hexa);
            ++ncubes;
            break;

        default:
            //if the type is not handled, skip rest of the line
            std::string tmp;
            std::getline(file, tmp);
        }
    }

    edges.endEdit();
    triangles.endEdit();
    quads.endEdit();
    tetrahedra.endEdit();
    hexahedra.endEdit();

    file >> cmd;
    if (cmd != "$ENDELM" && cmd!="$EndElements")
    {
        serr << "Error: MeshGmshLoader: '$ENDELM' or '$EndElements' expected, found '" << cmd << "'" << sendl;
        file.close();
        return false;
    }

    // 	sout << "Loading topology complete:";
    // 	if (npoints>0) sout << ' ' << npoints << " points";
    // 	if (nlines>0)  sout << ' ' << nlines  << " lines";
    // 	if (ntris>0)   sout << ' ' << ntris   << " triangles";
    // 	if (nquads>0)  sout << ' ' << nquads  << " quads";
    // 	if (ntetrahedra>0) sout << ' ' << ntetrahedra << " tetrahedra";
    // 	if (ncubes>0)  sout << ' ' << ncubes  << " cubes";
    // 	sout << sendl;

    file.close();
    return true;
}


} // namespace loader

} // namespace component

} // namespace sofa

