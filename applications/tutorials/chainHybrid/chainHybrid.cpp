/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU General Public License as published by the Free  *
* Software Foundation; either version 2 of the License, or (at your option)   *
* any later version.                                                          *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    *
* more details.                                                               *
*                                                                             *
* You should have received a copy of the GNU General Public License along     *
* with this program; if not, write to the Free Software Foundation, Inc., 51  *
* Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.                   *
*******************************************************************************
*                            SOFA :: Applications                             *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/

#include "../objectCreator/ObjectCreator.h"

#include <sofa/helper/ArgumentParser.h>
#include <sofa/simulation/tree/TreeSimulation.h>
#include <sofa/simulation/common/Node.h>

#include <sofa/gui/GUIManager.h>
#include <sofa/gui/Main.h>
#include <sofa/helper/system/FileRepository.h>


#include <sofa/component/loader/MeshGmshLoader.h>
#include <sofa/component/topology/MeshTopology.h>
#include <sofa/component/topology/RegularGridTopology.h>


//Using double by default, if you have SOFA_FLOAT in use in you sofa-default.cfg, then it will be FLOAT.
#include <sofa/component/typedef/Sofa_typedef.h>
// ---------------------------------------------------------------------
// ---
// ---------------------------------------------------------------------
using namespace sofa::simulation;
using namespace sofa::component::container;
using namespace sofa::component::loader;
using namespace sofa::component::topology;

Node *createChainHybrid(Node *root)
{
    const std::string visualModel="mesh/torus.obj";
    const std::string collisionModel="mesh/torus_for_collision.obj";

    std::vector<std::string> modelTypes;
    modelTypes.push_back("Triangle");
    modelTypes.push_back("Line");
    modelTypes.push_back("Point");

    //Elements of the scene
    //------------------------------------
    Node::SPtr  chain = root->createChild("Chain");

    //************************************
    //Torus Fixed
    {
        Node::SPtr  torusFixed = sofa::SimpleObjectCreator::CreateObstacle(chain,"mesh/torus_for_collision.obj", "mesh/torus.obj", "gray");
    }
    //************************************
    //Torus FEM
    {
        Node::SPtr  torusFEM = sofa::SimpleObjectCreator::CreateEulerSolverNode(chain,"FEM");

        MeshGmshLoader::SPtr  loaderFEM = New<MeshGmshLoader>();
        loaderFEM->setFilename(sofa::helper::system::DataRepository.getFile("mesh/torus_low_res.msh"));
        loaderFEM->load();
        torusFEM->addObject(loaderFEM);

        MeshTopology::SPtr meshTorusFEM = sofa::core::objectmodel::New<MeshTopology>();
        meshTorusFEM->setSrc("",loaderFEM.get());
        torusFEM->addObject(meshTorusFEM);

        const Deriv3 translation(2.5,0,0);
        const Deriv3 rotation(90,0,0);

        MechanicalObject3::SPtr dofFEM = sofa::core::objectmodel::New<MechanicalObject3>(); dofFEM->setName("FEM Object");
        dofFEM->setTranslation(translation[0],translation[1],translation[2]);
        dofFEM->setRotation(rotation[0],rotation[1],rotation[2]);
        torusFEM->addObject(dofFEM);

        UniformMass3::SPtr uniMassFEM = sofa::core::objectmodel::New<UniformMass3>();
        uniMassFEM->setTotalMass(5); //the whole object will have 5 as given mass
        torusFEM->addObject(uniMassFEM);

        TetrahedronFEMForceField3::SPtr tetraFEMFF = sofa::core::objectmodel::New<TetrahedronFEMForceField3>();
        tetraFEMFF->setName("FEM");
        tetraFEMFF->setComputeGlobalMatrix(false);
        tetraFEMFF->setMethod("large");
        tetraFEMFF->setPoissonRatio(0.3);
        tetraFEMFF->setYoungModulus(1000);
        torusFEM->addObject(tetraFEMFF);

        //Node VISUAL
        Node::SPtr  FEMVisualNode = sofa::SimpleObjectCreator::CreateVisualNodeVec3(torusFEM, dofFEM.get(),visualModel, "red", translation, rotation);

        //Node COLLISION
        Node::SPtr  FEMCollisionNode = sofa::SimpleObjectCreator::CreateCollisionNodeVec3(torusFEM, dofFEM.get(),collisionModel,modelTypes, translation, rotation );
    }
    //************************************
    //Torus Spring
    {
        Node::SPtr  torusSpring = sofa::SimpleObjectCreator::CreateEulerSolverNode(chain,"Spring");

        MeshGmshLoader::SPtr  loaderSpring = New<MeshGmshLoader>();
        loaderSpring->setFilename(sofa::helper::system::DataRepository.getFile("mesh/torus_low_res.msh"));
        loaderSpring->load();

        torusSpring->addObject(loaderSpring);
        loaderSpring->init();

        MeshTopology::SPtr meshTorusSpring = sofa::core::objectmodel::New<MeshTopology>();
        meshTorusSpring->setSrc("",loaderSpring.get());
        torusSpring->addObject(meshTorusSpring);

        const Deriv3 translation(5,0,0);
        const Deriv3 rotation(0,0,0);

        MechanicalObject3::SPtr dofSpring = sofa::core::objectmodel::New<MechanicalObject3>(); dofSpring->setName("Spring Object");

        dofSpring->setTranslation(translation[0],translation[1],translation[2]);
        dofSpring->setRotation(rotation[0],rotation[1],rotation[2]);

        torusSpring->addObject(dofSpring);

        UniformMass3::SPtr uniMassSpring = sofa::core::objectmodel::New<UniformMass3>();
        uniMassSpring->setTotalMass(5); //the whole object will have 5 as given mass
        torusSpring->addObject(uniMassSpring);

        MeshSpringForceField3::SPtr springFF = sofa::core::objectmodel::New<MeshSpringForceField3>();
        springFF->setName("Springs");
        springFF->setStiffness(400);
        springFF->setDamping(0);
        torusSpring->addObject(springFF);


        //Node VISUAL
        Node::SPtr  SpringVisualNode = sofa::SimpleObjectCreator::CreateVisualNodeVec3(torusSpring, dofSpring.get(), visualModel,"green", translation, rotation);

        //Node COLLISION
        Node::SPtr  SpringCollisionNode = sofa::SimpleObjectCreator::CreateCollisionNodeVec3(torusSpring, dofSpring.get(), collisionModel,modelTypes,translation, rotation);
    }
    //************************************
    //Torus FFD
    {
        Node::SPtr  torusFFD = sofa::SimpleObjectCreator::CreateEulerSolverNode(chain,"FFD");

        const Deriv3 translation(7.5,0,0);
        const Deriv3 rotation(90,0,0);

        MechanicalObject3::SPtr dofFFD = sofa::core::objectmodel::New<MechanicalObject3>(); dofFFD->setName("FFD Object");
        dofFFD->setTranslation(translation[0],translation[1],translation[2]);
        dofFFD->setRotation(rotation[0],rotation[1],rotation[2]);
        torusFFD->addObject(dofFFD);

        UniformMass3::SPtr uniMassFFD = sofa::core::objectmodel::New<UniformMass3>();
        uniMassFFD->setTotalMass(5); //the whole object will have 5 as given mass
        torusFFD->addObject(uniMassFFD);

        RegularGridTopology::SPtr gridTopo = sofa::core::objectmodel::New<RegularGridTopology>(6,2,5); //dimension of the grid
        gridTopo->setPos(
            -2.5,2.5,  //Xmin, Xmax
            -0.5,0.5,  //Ymin, Ymax
            -2,2       //Zmin, Zmax
        );
        torusFFD->addObject(gridTopo);

        RegularGridSpringForceField3::SPtr FFDFF = sofa::core::objectmodel::New<RegularGridSpringForceField3>();
        FFDFF->setName("Springs FFD");
        FFDFF->setStiffness(200);
        FFDFF->setDamping(0);
        torusFFD->addObject(FFDFF);

        //Node VISUAL
        Node::SPtr  FFDVisualNode = sofa::SimpleObjectCreator::CreateVisualNodeVec3(torusFFD, dofFFD.get(), visualModel,"yellow", translation);

        //Node COLLISION
        Node::SPtr  FFDCollisionNode = sofa::SimpleObjectCreator::CreateCollisionNodeVec3(torusFFD ,dofFFD.get(),collisionModel,modelTypes, translation);
    }

    //************************************
    //Torus Rigid
    {
        Node::SPtr  torusRigid = sofa::SimpleObjectCreator::CreateEulerSolverNode(chain,"Rigid");

        const Deriv3 translation(10,0,0);
        const Deriv3 rotation(0,0,0);

        MechanicalObjectRigid3::SPtr dofRigid = sofa::core::objectmodel::New<MechanicalObjectRigid3>(); dofRigid->setName("Rigid Object");
        dofRigid->setTranslation(translation[0],translation[1],translation[2]);
        dofRigid->setRotation(rotation[0],rotation[1],rotation[2]);
        torusRigid->addObject(dofRigid);

        UniformMassRigid3::SPtr uniMassRigid = sofa::core::objectmodel::New<UniformMassRigid3>();
        uniMassRigid->setTotalMass(1); //the whole object will have 5 as given mass
        torusRigid->addObject(uniMassRigid);

        //Node VISUAL
        Node::SPtr  RigidVisualNode = sofa::SimpleObjectCreator::CreateVisualNodeRigid(torusRigid, dofRigid.get(), visualModel,"gray");

        //Node COLLISION
        Node::SPtr  RigidCollisionNode = sofa::SimpleObjectCreator::CreateCollisionNodeRigid(torusRigid, dofRigid.get(),collisionModel,modelTypes);
    }
    return root;
}



int main(int argc, char** argv)
{
    glutInit(&argc,argv);

    std::vector<std::string> files;

    sofa::helper::parse("This is a SOFA application. Here are the command line arguments")
    (argc,argv);


    sofa::simulation::setSimulation(new sofa::simulation::tree::TreeSimulation());

    sofa::gui::initMain();
    sofa::gui::GUIManager::Init(argv[0]);

    // The graph root node
    Node::SPtr root = sofa::SimpleObjectCreator::CreateRootWithCollisionPipeline();
    root->setGravity( Coord3(0,0,-10) );

    //Add the objects
    createChainHybrid(root.get());

    root->setAnimate(false);

    sofa::simulation::getSimulation()->init(root.get());


    //=======================================
    // Run the main loop
    sofa::gui::GUIManager::MainLoop(root);

    return 0;
}
