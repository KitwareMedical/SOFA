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
#include <sofa/component/collision/RemovePrimitivePerformer.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/topology/TopologicalMapping.h>
#include <sofa/helper/Factory.inl>
#include <sofa/helper/system/glut.h>

#if 0
#ifdef SOFA_DEV
#include <sofa/component/collision/TetrahedronModel.h>
#endif // SOFA_DEV
#endif
#include <sofa/simulation/common/Simulation.h>

namespace sofa
{

namespace component
{

namespace collision
{

using namespace sofa::core::topology;

template <class DataTypes>
RemovePrimitivePerformer<DataTypes>::RemovePrimitivePerformer(BaseMouseInteractor *i)
    :TInteractionPerformer<DataTypes>(i)
    ,firstClick (0)
    ,surfaceOnVolume(false)
    ,volumeOnSurface(false)
    ,topo_curr(NULL)
{}


/// Functions called in framework of the mouse Interactor
//***************************************************************************************************************

template <class DataTypes>
void RemovePrimitivePerformer<DataTypes>::start()
{
    if (topologicalOperation != 0)
    {
        if (firstClick)
            firstClick = false;
        else
            firstClick = true;
    }
}


template <class DataTypes>
void RemovePrimitivePerformer<DataTypes>::execute()
{
    // - STEP 1: Get body picked and Mstate associated
    picked=this->interactor->getBodyPicked();
    if (!picked.body) return;

    mstateCollision = dynamic_cast< core::behavior::MechanicalState<DataTypes>*    >(picked.body->getContext()->getMechanicalState());
    if (!mstateCollision)
    {
        std::cerr << "incompatible MState during Mouse Interaction " << std::endl;
        return;
    }

    // - STEP 1: Checking type of operation
    if (topologicalOperation == 0) // normal case, remove directly one element
    {
        core::CollisionElementIterator collisionElement( picked.body, picked.indexCollisionElement);
        core::CollisionModel* model = collisionElement.getCollisionModel();

        sofa::core::topology::TopologyModifier* topologyModifier;
        picked.body->getContext()->get(topologyModifier);

        // Handle Removing of topological element (from any type of topology)
        if(topologyModifier)
            topologyChangeManager.removeItemsFromCollisionModel(model, (int)picked.indexCollisionElement);

        picked.body=NULL;
        this->interactor->setBodyPicked(picked);
    }
    else // second case remove a zone of element
    {
        if (firstClick) // first click detected => creation of the zone
        {
            if (!createElementList())
                return;
        }
        else // second clic removing zone stored in selectedElem
        {
            if (selectedElem.empty())
                return;

            core::CollisionElementIterator collisionElement( picked.body, picked.indexCollisionElement);

            sofa::core::topology::TopologyModifier* topologyModifier;
            picked.body->getContext()->get(topologyModifier);

            // Problem of type takeng by functions called: Converting selectedElem <unsigned int> in <int>
            helper::vector<int> ElemList_int;
            ElemList_int.resize(selectedElem.size());
            for (unsigned int i = 0; i<selectedElem.size(); ++i)
                ElemList_int[i] = selectedElem[i];

            // Creating model of collision
            core::CollisionModel::SPtr model;
            if (surfaceOnVolume) // In the case of deleting a volume from a surface an volumique collision model is needed (only tetra available for the moment)
            {
#if 0
#ifdef SOFA_DEV
                model = sofa::core::objectmodel::New<TetrahedronModel>();
                //model->setContext(topo_curr->getContext());
                topo_curr->getContext()->addObject(model);
#endif // SOFA_DEV
#endif
            }
            else // other cases, collision model from pick is taken
            {
                model = collisionElement.getCollisionModel();
            }

            // Handle Removing of topological element (from any type of topology)
            if(topologyModifier) topologyChangeManager.removeItemsFromCollisionModel(model.get(),ElemList_int );
            picked.body=NULL;
            this->interactor->setBodyPicked(picked);

            if (surfaceOnVolume) // In the case of deleting a volume from a surface an volumique collision model is needed (only tetra available for the moment)
            {
#if 0
#ifdef SOFA_DEV
                topo_curr->getContext()->removeObject(model);
#endif // SOFA_DEV
#endif
            }
            selectedElem.clear();
        }
    }
}


template <class DataTypes>
void RemovePrimitivePerformer<DataTypes>::end()
{
    std::cout << "RemovePrimitivePerformer::end()" << std::endl;
    //	firstClick = true;
}



//***************************************************************************************************************
// Internal functions

// ** Creating a list of elements concerned by the removal operation **
template <class DataTypes>
bool RemovePrimitivePerformer<DataTypes>::createElementList()
{
    // - STEP 1: Looking for current topology type
    topo_curr = picked.body->getContext()->getMeshTopology();
    if (topo_curr->getNbHexahedra())
        topoType = HEXAHEDRON;
    else if (topo_curr->getNbTetrahedra())
        topoType = TETRAHEDRON;
    else if (topo_curr->getNbQuads())
        topoType = QUAD;
    else if (topo_curr->getNbTriangles())
        topoType = TRIANGLE;
    else
    {
        std::cerr << "Error: No topology has been found." << std::endl;
        return false;
    }

    // Initialization of first element
    selectedElem.clear();
    selectedElem.resize (1);
    selectedElem[0] = picked.indexCollisionElement;

    // - STEP 2: Looking for type of zone to remove
    if (!volumicMesh) // Surfacique case
    {
        volumeOnSurface = false;
        sofa::core::topology::TopologyObjectType topoTypeTmp = topoType;

        // - STEP 3: Looking for tricky case
        if (topoType == TETRAHEDRON || topoType == HEXAHEDRON) // special case: removing a surface volume on the mesh (tetra only for the moment)
        {
            // looking for mapping VolumeToSurface
            simulation::Node *node_curr = dynamic_cast<simulation::Node*>(topo_curr->getContext());
            std::vector< core::objectmodel::BaseObject * > listObject;
            node_curr->get<core::objectmodel::BaseObject>(&listObject, core::objectmodel::BaseContext::SearchRoot);

            for(unsigned int i=0; i<listObject.size(); ++i) // loop on all components to find mapping
            {
                sofa::core::topology::TopologicalMapping *topoMap = dynamic_cast<sofa::core::topology::TopologicalMapping *>(listObject[i]);
                if (topoMap)
                {
                    // Mapping found: 1- looking for volume, 2- looking for surface element on border, 3- looking for correspondant ID element in surfacique mesh
                    const BaseMeshTopology::TrianglesInTetrahedron& tetraTri = topo_curr->getTrianglesInTetrahedron(selectedElem[0]);

                    int volTmp = -1;
                    std::map<unsigned int, unsigned int> MappingMap = topoMap->getGlob2LocMap();
                    std::map<unsigned int, unsigned int>::iterator it;

                    for (unsigned int j = 0; j<4; ++j)
                    {
                        it = MappingMap.find (tetraTri[j]);
                        if ( it != MappingMap.end())
                        {
                            volTmp = (*it).second;
                            break;
                        }
                    }

                    if (volTmp == -1)
                    {
                        std::cerr << "Error while looking for corresponding element on surfacique mesh." << std::endl;
                        return false;
                    }

                    // Surfacique element has been found, computation will be done on surfacique mesh => switch temporary all variables to surface
                    selectedElem[0] = (unsigned int)volTmp;
                    volumeOnSurface = true;
                    topo_curr = topoMap->getTo();
                    topoType = TRIANGLE;
                }
            }

            if (!volumeOnSurface)
            {
                std::cerr << "Error: Trying to remove a volume at the surface of the mesh without using mapping volume to surface mesh. This case is not handle." << std::endl;
                return false;
            }
        }


        // - STEP 4: Loop on getNeighboorElements and getElementInZone until no more nighboor are in zones
        // Initialization
        bool end = false;
        VecIds tmp = getNeighboorElements (selectedElem);
        VecIds tmp2;

        while (!end) // Creating region of interest
        {
            tmp2 = getElementInZone (tmp);
            tmp.clear();

            if (tmp2.empty())
                end = true;

            for (unsigned int t = 0; t<tmp2.size(); ++t)
                selectedElem.push_back (tmp2[t]);

            tmp = getNeighboorElements (tmp2);
            tmp2.clear ();
        }


        // - STEP 5: Postprocessing: zone using surface element has been found, extract volumes behind (small error on boundary regarding barycentric points)
        if (volumeOnSurface)
        {
            // Get dofs on surface
            for (unsigned int i = 0; i<selectedElem.size(); ++i)
            {
                helper::vector<unsigned int> elem;

                switch ( topoType ) // Get surfacique elements as array of vertices
                {
                case QUAD:
                {
                    const BaseMeshTopology::Quad& quad = topo_curr->getQuad(selectedElem[i]);
                    elem.resize(4);
                    for (unsigned int j = 0; j<4; ++j)
                        elem[j] = quad[j];
                    break;
                }
                case TRIANGLE:
                {
                    const BaseMeshTopology::Triangle& tri = topo_curr->getTriangle(selectedElem[i]);
                    elem.resize(3);
                    for (unsigned int j = 0; j<3; ++j)
                        elem[j] = tri[j];
                    break;
                }
                default:
                    break;
                }

                // Pattern fill vector without redundancy
                for (unsigned int j = 0; j<elem.size(); ++j)
                {
                    bool dofFind = false;
                    unsigned int Selem = elem[j];

                    for (unsigned int k = 0; k<tmp2.size(); ++k)
                        if (tmp2[j] == Selem)
                        {
                            dofFind = true;
                            break;
                        }

                    if (!dofFind)
                        tmp2.push_back(Selem);
                }
            }

            // Switching variables to initial topology (topotype, topology) clear list of surfacique elements selected
            topo_curr = picked.body->getMeshTopology();
            topoType = topoTypeTmp;
            selectedElem.clear();

            // Get Volumique elements from list of vertices in tmp2
            for (unsigned int i = 0; i<tmp2.size(); ++i)
            {
                helper::vector<unsigned int> elem;

                switch ( topoType )
                {
                case HEXAHEDRON:
                {
                    const BaseMeshTopology::HexahedraAroundVertex& hexaV = topo_curr->getHexahedraAroundVertex(tmp2[i]);
                    for (unsigned int j = 0; j<hexaV.size(); ++j)
                        elem.push_back(hexaV[j]);

                    break;
                }
                case TETRAHEDRON:
                {
                    const BaseMeshTopology::TetrahedraAroundVertex& tetraV = topo_curr->getTetrahedraAroundVertex(tmp2[i]);
                    for (unsigned int j = 0; j<tetraV.size(); ++j)
                        elem.push_back(tetraV[j]);

                    break;
                }
                default:
                    break;
                }

                // Pattern fill vector without redundancy
                for (unsigned int j = 0; j<elem.size(); ++j)
                {
                    bool Vfind = false;
                    unsigned int VelemID = elem[j];

                    for (unsigned int k = 0; k<selectedElem.size(); ++k) // Check if not already insert
                        if (selectedElem[k] == VelemID)
                        {
                            Vfind = true;
                            break;
                        }

                    if (!Vfind)
                        selectedElem.push_back (VelemID);
                }
            }
        }

    }
    else // - STEP 2: Volumique case
    {
        surfaceOnVolume = false;

        // - STEP 3: Looking for tricky case
        if (topoType == TRIANGLE || topoType == QUAD) // Special case: removing a volumique zone on the mesh while starting at the surface
        {
            // looking for mapping VolumeToSurface
            simulation::Node *node_curr = dynamic_cast<simulation::Node*>(topo_curr->getContext());
            std::vector< core::objectmodel::BaseObject * > listObject;
            node_curr->get<core::objectmodel::BaseObject>(&listObject, core::objectmodel::BaseContext::Local);

            for(unsigned int i=0; i<listObject.size(); ++i) // loop on all components to find mapping (only tetra for the moment)
            {
                sofa::core::topology::TopologicalMapping *topoMap = dynamic_cast<sofa::core::topology::TopologicalMapping *>(listObject[i]);
                if (topoMap)
                {
                    // Mapping found: 1- get surface element ID in volumique topology, 2- get volume element ID behind surface element, 3- switching all variables to volumique case
                    unsigned int volTmp = (topoMap->getLoc2GlobVec()).getValue()[selectedElem[0]];
                    topo_curr = topoMap->getFrom();
                    selectedElem[0] = topo_curr->getTetrahedraAroundTriangle(volTmp)[0];
                    surfaceOnVolume = true;
                    topoType = TETRAHEDRON;
                }
            }

            if (!surfaceOnVolume)
            {
                std::cerr << "Error: Trying to remove a volume using a surfacique mesh without mapping to volumique mesh." << std::endl;
                return false;
            }
        }

        // - STEP 4: Loop on getNeighboorElements and getElementInZone until no more nighboor are in zones
        // Initialization
        bool end = false;
        VecIds tmp = getNeighboorElements (selectedElem);
        VecIds tmp2;

        while (!end) // Creating region of interest
        {
            tmp2 = getElementInZone (tmp);

            tmp.clear();

            if (tmp2.empty())
                end = true;

            for (unsigned int t = 0; t<tmp2.size(); ++t)
                selectedElem.push_back (tmp2[t]);

            tmp = getNeighboorElements (tmp2);
            tmp2.clear ();
        }

    }

    return true;
}



// ** Return a vector of elements directly neighboor of a given list of elements **
template <class DataTypes>
sofa::helper::vector <unsigned int> RemovePrimitivePerformer<DataTypes>::getNeighboorElements(VecIds& elementsToTest)
{
    VecIds vertexList;
    VecIds neighboorList;


    // - STEP 1: get list of element vertices
    for (unsigned int i = 0; i<elementsToTest.size(); ++i)
    {
        helper::vector<unsigned int> elem;

        switch ( topoType ) // Get element as array of vertices
        {
        case HEXAHEDRON:
        {
            const BaseMeshTopology::Hexa& hexa = topo_curr->getHexahedron(elementsToTest[i]);
            elem.resize(8);
            for (unsigned int j = 0; j<8; ++j)
                elem[j] = hexa[j];
            break;
        }
        case TETRAHEDRON:
        {
            const BaseMeshTopology::Tetra& tetra = topo_curr->getTetrahedron(elementsToTest[i]);
            elem.resize(4);
            for (unsigned int j = 0; j<4; ++j)
                elem[j] = tetra[j];
            break;
        }
        case QUAD:
        {
            const BaseMeshTopology::Quad& quad = topo_curr->getQuad(elementsToTest[i]);
            elem.resize(4);
            for (unsigned int j = 0; j<4; ++j)
                elem[j] = quad[j];
            break;
        }
        case TRIANGLE:
        {
            const BaseMeshTopology::Triangle& tri = topo_curr->getTriangle(elementsToTest[i]);
            elem.resize(3);
            for (unsigned int j = 0; j<3; ++j)
                elem[j] = tri[j];
            break;
        }
        default:
            break;
        }

        // Pattern fill vector without redundancy
        for (unsigned int j = 0; j<elem.size(); ++j) // Insert vertices for each element
        {
            bool Vfind = false;
            unsigned int VelemID = elem[j];

            for (unsigned int k = 0; k<vertexList.size(); ++k) // Check if not already insert
                if (vertexList[k] == VelemID)
                {
                    Vfind = true;
                    break;
                }

            if (!Vfind)
                vertexList.push_back (VelemID);
        }
    }

    // - STEP 2: get list of element around vertices previously obtained
    for (unsigned int i = 0; i<vertexList.size(); ++i)
    {
        VecIds elemAroundV;

        switch ( topoType ) // Get elements around vertices as array of ID
        {
        case HEXAHEDRON:
        {
            elemAroundV = topo_curr->getHexahedraAroundVertex (vertexList[i]);
            break;
        }
        case TETRAHEDRON:
        {
            elemAroundV = topo_curr->getTetrahedraAroundVertex (vertexList[i]);
            break;
        }
        case QUAD:
        {
            elemAroundV = topo_curr->getQuadsAroundVertex (vertexList[i]);
            break;
        }
        case TRIANGLE:
        {
            elemAroundV = topo_curr->getTrianglesAroundVertex (vertexList[i]);
            break;
        }
        default:
            break;
        }

        // Pattern fill vector without redundancy + checking not insert in input selectedElem
        for (unsigned int j = 0; j<elemAroundV.size(); ++j)  // Insert each element as new neighboor
        {
            bool Efind = false;
            unsigned int elemID = elemAroundV[j];

            for (unsigned int k = 0; k<neighboorList.size(); ++k) // Check if not already insert
                if (neighboorList[k] == elemID)
                {
                    Efind = true;
                    break;
                }

            if (!Efind)
                for (unsigned int k = 0; k<selectedElem.size(); ++k) // Check if not in selected list
                    if (selectedElem[k] == elemID)
                    {
                        Efind = true;
                        break;
                    }

            if (!Efind)
                neighboorList.push_back (elemID);
        }
    }

    return neighboorList;
}


// ** Function testing if elements are in the range of a given zone **
template <class DataTypes>
sofa::helper::vector <unsigned int> RemovePrimitivePerformer<DataTypes>::getElementInZone(VecIds& elementsToTest)
{
    // - STEP 0: Compute appropriate scale from BB:  selectorScale = 100 => zone = all mesh
    Vec<3, SReal> sceneMinBBox, sceneMaxBBox;
    core::objectmodel::BaseNode* root = dynamic_cast<core::objectmodel::BaseNode*>(mstateCollision->getContext());
    if (root) root = root->getRoot();
    if (root) { sceneMinBBox = root->f_bbox.getValue().minBBox(); sceneMaxBBox = root->f_bbox.getValue().maxBBox(); }
    else      { sceneMinBBox = mstateCollision->getContext()->f_bbox.getValue().minBBox(); sceneMaxBBox = mstateCollision->getContext()->f_bbox.getValue().maxBBox(); }
    Real BB_size = (sceneMaxBBox - sceneMinBBox).norm();
    if (BB_size == 0)
    {
        std::cerr << "Error while computing Boundingbox size, size return null." << std::endl;
        BB_size = 1; // not to crash program
    }
    Real zone_size = (BB_size*selectorScale)/200;
    Real dist;
    Coord center = picked.point;

    // - STEP 2: Compute baryCoord of elements in list:
    const VecCoord& X = *mstateCollision->getX();

    VecCoord baryCoord;
    baryCoord.resize (elementsToTest.size());

    for (unsigned int i = 0; i<elementsToTest.size(); ++i)
    {
        unsigned int N = 1;

        switch ( topoType ) // get element as array of vertices and sum the coordinates
        {
        case HEXAHEDRON:
        {
            const BaseMeshTopology::Hexa& hexa = topo_curr->getHexahedron(elementsToTest[i]);
            baryCoord[i] = X[hexa[0]] + X[hexa[1]] + X[hexa[2]] + X[hexa[3]] +
                    X[hexa[4]] + X[hexa[5]] + X[hexa[6]] + X[hexa[7]];
            N = 8;

            break;
        }
        case TETRAHEDRON:
        {
            const BaseMeshTopology::Tetra& tetra = topo_curr->getTetrahedron(elementsToTest[i]);
            baryCoord[i] = X[tetra[0]] + X[tetra[1]] + X[tetra[2]] + X[tetra[3]];
            N = 4;

            break;
        }
        case QUAD:
        {
            const BaseMeshTopology::Quad& quad = topo_curr->getQuad(elementsToTest[i]);
            baryCoord[i] = X[quad[0]] + X[quad[1]] + X[quad[2]] + X[quad[3]];
            N = 4;

            break;
        }
        case TRIANGLE:
        {
            const BaseMeshTopology::Triangle& tri = topo_curr->getTriangle(elementsToTest[i]);
            baryCoord[i] = X[tri[0]] + X[tri[1]] + X[tri[2]];
            N = 3;

            break;
        }
        default:
            break;
        }

        for (unsigned int j = 0; j<center.size(); ++j) // divided each coordinate by N (number of vertices)
            baryCoord[i][j] = baryCoord[i][j]/N;

    }


    VecIds elemInside;
    // - STEP 3: Test if barycentric points are inside the zone
    for (unsigned int i = 0; i<elementsToTest.size(); ++i)
    {
        //compute distance from barycenter to center zone
        dist = (baryCoord[i] - center).norm();

        if (dist < zone_size)
            elemInside.push_back (elementsToTest[i]);
    }

    return elemInside;
}



//***************************************************************************************************************

template <class DataTypes>
void RemovePrimitivePerformer<DataTypes>::draw(const core::visual::VisualParams* )
{
    if (picked.body == NULL) return;

    if (mstateCollision == NULL) return;


    const VecCoord& X = *mstateCollision->getX();
    //core::topology::BaseMeshTopology* topo = picked.body->getMeshTopology();

    glDisable(GL_LIGHTING);
    glColor3f(0.3,0.8,0.3);


    if (topoType == QUAD || topoType == HEXAHEDRON)
        glBegin (GL_QUADS);
    else
        glBegin (GL_TRIANGLES);


    for (unsigned int i=0; i<selectedElem.size(); ++i)
    {
        helper::vector<unsigned int> elem;

        switch ( topoType )
        {
        case HEXAHEDRON:
        {
            const BaseMeshTopology::Hexa& hexa = topo_curr->getHexahedron(selectedElem[i]);
            Coord coordP[8];

            for (unsigned int j = 0; j<8; j++)
                coordP[j] = X[hexa[j]];

            for (unsigned int j = 0; j<8; ++j)
            {
                glVertex3d(coordP[j][0], coordP[j][1], coordP[j][2]);
                glVertex3d(coordP[(j+1)%4][0], coordP[(j+1)%4][1], coordP[(j+1)%4][2]);
                glVertex3d(coordP[(j+2)%4][0], coordP[(j+2)%4][1], coordP[(j+2)%4][2]);
                glVertex3d(coordP[(j+3)%4][0], coordP[(j+3)%4][1], coordP[(j+3)%4][2]);
            }
            break;
        }
        case TETRAHEDRON:
        {
            const BaseMeshTopology::Tetra& tetra = topo_curr->getTetrahedron(selectedElem[i]);
            Coord coordP[4];

            for (unsigned int j = 0; j<4; j++)
                coordP[j] = X[tetra[j]];

            for (unsigned int j = 0; j<4; ++j)
            {
                glVertex3d(coordP[j][0], coordP[j][1], coordP[j][2]);
                glVertex3d(coordP[(j+1)%4][0], coordP[(j+1)%4][1], coordP[(j+1)%4][2]);
                glVertex3d(coordP[(j+2)%4][0], coordP[(j+2)%4][1], coordP[(j+2)%4][2]);
            }
            break;
        }
        case QUAD:
        {
            const BaseMeshTopology::Quad& quad = topo_curr->getQuad(selectedElem[i]);

            for (unsigned int j = 0; j<4; j++)
            {
                Coord coordP = X[quad[j]];
                glVertex3d(coordP[0], coordP[1], coordP[2]);
            }
            break;
        }
        case TRIANGLE:
        {
            const BaseMeshTopology::Triangle& tri = topo_curr->getTriangle(selectedElem[i]);

            for (unsigned int j = 0; j<3; j++)
            {
                Coord coordP = X[tri[j]];
                glVertex3d(coordP[0] * 1.001, coordP[1] * 1.001, coordP[2] * 1.001);
            }
            for (unsigned int j = 0; j<3; j++)
            {
                Coord coordP = X[tri[j]];
                glVertex3d(coordP[0] * 0.999, coordP[1] * 0.999, coordP[2] * 0.999);
            }

            break;
        }
        default:
            break;
        }



    }
    glEnd();
}


}
}
}

