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
#include <sofa/helper/io/MeshOBJ.h>
#include <sofa/helper/system/FileRepository.h>
#include <sofa/helper/system/SetDirectory.h>

namespace sofa
{

namespace helper
{

namespace io
{

using namespace sofa::defaulttype;
using namespace sofa::core::loader;

SOFA_DECL_CLASS(MeshOBJ)

Creator<Mesh::FactoryMesh,MeshOBJ> MeshOBJClass("obj");

void MeshOBJ::init (std::string filename)
{
    if (!sofa::helper::system::DataRepository.findFile(filename))
    {
        std::cerr << "File " << filename << " not found " << std::endl;
        return;
    }
    std::ifstream file(filename.c_str());
    readOBJ(file, filename);
    file.close();
}

void MeshOBJ::readOBJ (std::ifstream &file, const std::string &filename)
{
    vector< vector<int> > vertNormTexIndices;
    vector<int>vIndices, nIndices, tIndices;
    int vtn[3];
    Vec3d result;
    Vec3d texCoord;
    Vec3d normal;
    int nbf = facets.size();

    std::string line;

    PrimitiveGroup curGroup;

    while( std::getline(file,line) )
    {
        if (line.empty()) continue;
        std::istringstream values(line);
        std::string token;
        values >> token;
        if (token == "#")
        {
            /* comment */
        }
        else if (token == "v")
        {
            /* vertex */
            values >> result[0] >> result[1] >> result[2];
            vertices.push_back(result);
        }
        else if (token == "vn")
        {
            /* normal */
            values >> result[0] >> result[1] >> result[2];
            normals.push_back(result);
        }
        else if (token == "vt")
        {
            /* texcoord */
            values >> result[0] >> result[1];
            texCoords.push_back(result);
        }
        else if (token == "mtllib")
        {
            while (!values.eof())
            {
                std::string materialLibaryName;
                values >> materialLibaryName;
                std::string mtlfile = sofa::helper::system::SetDirectory::GetRelativeFromFile(materialLibaryName.c_str(), filename.c_str());
                readMTL(mtlfile.c_str());
            }
        }
        else if (token == "usemtl" || token == "g")
        {
            // end of current group
            curGroup.nbp = nbf - curGroup.p0;
            if (curGroup.nbp > 0) groups.push_back(curGroup);
            curGroup.p0 = nbf;
            curGroup.nbp = 0;
            if (token == "usemtl")
            {
                curGroup.materialId = -1;
                values >> curGroup.materialName;
                vector<Material>::iterator it = materials.begin();
                vector<Material>::iterator itEnd = materials.end();
                for (; it != itEnd; it++)
                {
                    if (it->name == curGroup.materialName)
                    {
                        // std::cout << "Using material "<<it->name<<std::endl;
                        (*it).activated = true;
                        if (!material.activated)
                            material = *it;
                        curGroup.materialId = it - materials.begin();
                        break;
                    }
                }
            }
            else if (token == "g")
            {
                curGroup.groupName.clear();
                while (!values.eof())
                {
                    std::string g;
                    values >> g;
                    if (!curGroup.groupName.empty())
                        curGroup.groupName += " ";
                    curGroup.groupName += g;
                }
            }

        }
        else if (token == "l" || token == "f")
        {
            /* face */
            vIndices.clear();
            nIndices.clear();
            tIndices.clear();
            vertNormTexIndices.clear();

            while (!values.eof())
            {
                std::string face;
                values >> face;
                if (face.empty()) continue;
                for (int j = 0; j < 3; j++)
                {
                    vtn[j] = 0;
                    std::string::size_type pos = face.find('/');
                    std::string tmp = face.substr(0, pos);
                    if (pos == std::string::npos)
                        face = "";
                    else
                    {
                        face = face.substr(pos + 1);
                    }

                    if (!tmp.empty())
                        vtn[j] = atoi(tmp.c_str()) - 1; // -1 because the numerotation begins at 1 and a vector begins at 0
                }

                vIndices.push_back(vtn[0]);
                nIndices.push_back(vtn[1]);
                tIndices.push_back(vtn[2]);
            }

            vertNormTexIndices.push_back (vIndices);
            vertNormTexIndices.push_back (nIndices);
            vertNormTexIndices.push_back (tIndices);
            facets.push_back(vertNormTexIndices);
            ++nbf;
        }
        else
        {
            // std::cerr << "readObj : Unknown token for line " << line << std::endl;
        }
    }

    // end of current group
    curGroup.nbp = nbf - curGroup.p0;
    if (curGroup.nbp > 0) groups.push_back(curGroup);

    // announce the model statistics
    // std::cout << " Vertices: " << vertices.size() << std::endl;
    // std::cout << " Normals: " << normals.size() << std::endl;
    // std::cout << " Texcoords: " << texCoords.size() << std::endl;
    // std::cout << " Triangles: " << facets.size() << std::endl;

    if (vertices.size()>0)
    {
        // compute bbox
        Vector3 minBB = vertices[0];
        Vector3 maxBB = vertices[0];
        for (unsigned int i=1; i<vertices.size(); i++)
        {
            Vector3 p = vertices[i];
            for (int c=0; c<3; c++)
            {
                if (minBB[c] > p[c])
                    minBB[c] = p[c];
                if (maxBB[c] < p[c])
                    maxBB[c] = p[c];
            }
        }

        // 		std::cout << "BBox: <"<<minBB[0]<<','<<minBB[1]<<','<<minBB[2]<<">-<"<<maxBB[0]<<','<<maxBB[1]<<','<<maxBB[2]<<">\n";
    }

}

// -----------------------------------------------------
// readMTL: read a wavefront material library file
//
//    model - properly initialized GLMmodel structure
//    name  - name of the material library
// -----------------------------------------------------
void MeshOBJ::readMTL(const char* filename)
{
    FILE* file;
    char buf[128];
    file = fopen(filename, "r");
    Material *mat = NULL;
    if (!file);//std::cerr << "readMTL() failed: can't open material file " << filename << std::endl;
    else
    {
        /* now, read in the data */
        while (fscanf(file, "%s", buf) != EOF)
        {

            switch (buf[0])
            {
            case '#':
                /* comment */
                /* eat up rest of line */
                if ( fgets(buf, sizeof(buf), file) == NULL)
                {
                    if (feof (file) )
                        std::cerr << "Error: MeshOBJ: fgets function has encountered end of file." << std::endl;
                    else
                        std::cerr << "Error: MeshOBJ: fgets function has encountered an error." << std::endl;
                }
                break;
            case 'n':
                /* newmtl */
                if (mat != NULL)
                {
                    materials.push_back(*mat);
                    delete mat;
                    mat = NULL;
                }
                mat = new Material();
                if ( fgets(buf, sizeof(buf), file) == NULL)
                {
                    if (feof (file) )
                        std::cerr << "Error: MeshOBJ: fgets function has encountered end of file." << std::endl;
                    else
                        std::cerr << "Error: MeshOBJ: fgets function has encountered an error." << std::endl;
                }
                sscanf(buf, "%s %s", buf, buf);
                mat->name = buf;
                break;
            case 'N':
                switch (buf[1])
                {
                case 'i':
                {
                    float optical_density;
                    if (fscanf(file, "%f", &optical_density) == EOF)
                        std::cerr << "Error: MeshOBJ: fscanf has encountered an error" << std::endl;

                    break;
                }
                case 's':
                    if( fscanf(file, "%f", &mat->shininess) == EOF)
                        std::cerr << "Error: MeshOBJ: fscanf has encountered an error" << std::endl;
                    // wavefront shininess is from [0, 1000], so scale for OpenGL
                    //mat->shininess /= 1000.0;
                    //mat->shininess *= 128.0;
                    mat->useShininess = true;
                    break;
                default:
                    /* eat up rest of line */
                    if ( fgets(buf, sizeof(buf), file) == NULL)
                    {
                        if (feof (file) )
                            std::cerr << "Error: MeshOBJ: fgets function has encountered end of file." << std::endl;
                        else
                            std::cerr << "Error: MeshOBJ: fgets function has encountered an error." << std::endl;
                    }
                    break;
                }
                break;
            case 'K':
                switch (buf[1])
                {
                case 'd':
                    if( fscanf(file, "%f %f %f", &mat->diffuse[0], &mat->diffuse[1], &mat->diffuse[2]) == EOF)
                        std::cerr << "Error: MeshOBJ: fscanf has encountered an error" << std::endl;
                    mat->useDiffuse = true;
                    /*std::cout << mat->name << " diffuse = "<<mat->diffuse[0]<<' '<<mat->diffuse[1]<<'*/ /*'<<mat->diffuse[2]<<std::endl;*/
                    break;
                case 's':
                    if( fscanf(file, "%f %f %f", &mat->specular[0], &mat->specular[1], &mat->specular[2]) == EOF)
                        std::cerr << "Error: MeshOBJ: fscanf has encountered an error" << std::endl;
                    mat->useSpecular = true;
                    /*std::cout << mat->name << " specular = "<<mat->specular[0]<<' '<<mat->specular[1]<<'*/ /*'<<mat->specular[2]<<std::endl;*/
                    break;
                case 'a':
                    if( fscanf(file, "%f %f %f", &mat->ambient[0], &mat->ambient[1], &mat->ambient[2]) == EOF)
                        std::cerr << "Error: MeshOBJ: fscanf has encountered an error" << std::endl;
                    mat->useAmbient = true;
                    /*std::cout << mat->name << " ambient = "<<mat->ambient[0]<<' '<<mat->ambient[1]<<'*/ /*'<<mat->ambient[2]<<std::endl;*/
                    break;
                default:
                    /* eat up rest of line */
                    if (fgets(buf, sizeof(buf), file) == NULL)
                    {
                        if (feof (file) )
                            std::cerr << "Error: MeshOBJ: fgets function has encountered end of file." << std::endl;
                        else
                            std::cerr << "Error: MeshOBJ: fgets function has encountered an error." << std::endl;
                    }

                    break;
                }
                break;
            case 'd':
            case 'T':
                // transparency value
                if( fscanf(file, "%f", &mat->diffuse[3]) == EOF)
                    std::cerr << "Error: MeshOBJ: fscanf has encountered an error" << std::endl;
                break;

            case 'm':
            {
                //texture map
                char charFilename[128] = {0};
                if (fgets(charFilename, sizeof(charFilename), file)==NULL)
                {
                    std::cerr << "Error: MeshOBJ: fgets has encountered an error" << std::endl;
                }
                else
                {
                    mat->useTexture = true;

                    //store the filename of the texture map in the material

                    std::string stringFilename(charFilename);
                    //delete carriage return from the string assuming the next property of the .mtl file is at the next line
                    stringFilename.erase(stringFilename.end()-1, stringFilename.end());
                    stringFilename.erase(stringFilename.begin(), stringFilename.begin()+1);
                    mat->textureFilename = stringFilename;
                }
            }
            break;
            case 'b':
            {
                //bump mapping texture map
                char charFilename[128] = {0};
                if (fgets(charFilename, sizeof(charFilename), file)==NULL)
                {
                    std::cerr << "Error: MeshOBJ: fgets has encountered an error" << std::endl;
                }
                else
                {
                    mat->useBumpMapping = true;

                    //store the filename of the texture map in the material

                    std::string stringFilename(charFilename);
                    //delete carriage return from the string assuming the next property of the .mtl file is at the next line
                    stringFilename.erase(stringFilename.end()-1, stringFilename.end());
                    stringFilename.erase(stringFilename.begin(), stringFilename.begin()+1);
                    mat->bumpTextureFilename = stringFilename;
                }
            }
            break;
            default:
                /* eat up rest of line */
                if ( fgets(buf, sizeof(buf), file) == NULL)
                {
                    if (feof (file) )
                        std::cerr << "Error: MeshOBJ: fgets function has encountered end of file." << std::endl;
                    else
                        std::cerr << "Error: MeshOBJ: fgets function has encountered an error." << std::endl;
                }
                break;
            }

        }
        fclose(file);
    }

    if (mat != NULL)
    {
        materials.push_back(*mat);
        delete mat;
        mat = NULL;
    }
}

} // namespace io

} // namespace helper

} // namespace sofa

