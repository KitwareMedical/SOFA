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
#include <sofa/component/contextobject/Gravity.h>
#include <sofa/component/contextobject/CoordinateSystem.h>
#include <sofa/component/odesolver/EulerSolver.h>
#include <sofa/component/typedef/Sofa_typedef.h>
#include <sofa/component/visualmodel/OglModel.h>
#include <sofa/component/visualmodel/VisualStyle.h>

#include <sofa/core/objectmodel/Context.h>
#include <sofa/core/VecId.h>

#include <sofa/gui/GUIManager.h>
#include <sofa/gui/Main.h>

#include <sofa/helper/ArgumentParser.h>
#include <sofa/helper/system/FileRepository.h>
#include <sofa/helper/system/glut.h>

#include <sofa/simulation/tree/GNode.h>
#include <sofa/simulation/tree/TreeSimulation.h>

#include <iostream>
#include <fstream>

using namespace sofa::simulation::tree;
using sofa::simulation::Node;
using sofa::component::odesolver::EulerSolver;
using sofa::component::topology::MeshTopology;
using sofa::component::visualmodel::OglModel;
using sofa::core::objectmodel::Data;
using sofa::helper::ReadAccessor;
using sofa::helper::WriteAccessor;
using sofa::core::VecId;
// ---------------------------------------------------------------------
// ---
// ---------------------------------------------------------------------
int main(int argc, char** argv)
{
    glutInit(&argc,argv);
    sofa::helper::parse("This is a SOFA application.")
    (argc,argv);

    sofa::gui::initMain();
    sofa::gui::GUIManager::Init(argv[0]);

    // The graph root node : gravity already exists in a GNode by default
    sofa::simulation::setSimulation(new sofa::simulation::tree::TreeSimulation());
    sofa::simulation::Node::SPtr groot = sofa::simulation::getSimulation()->createNewGraph("root");
    groot->setGravity( Coord3(0,-10,0) );

    // One solver for all the graph
    EulerSolver::SPtr solver = sofa::core::objectmodel::New<EulerSolver>();
    solver->setName("solver");
    solver->f_printLog.setValue(false);
    groot->addObject(solver);

    // Tetrahedron degrees of freedom
    MechanicalObject3::SPtr DOF = sofa::core::objectmodel::New<MechanicalObject3>();
    groot->addObject(DOF);
    DOF->resize(4);
    DOF->setName("DOF");
    //get write access to the position vector of mechanical object DOF
    WriteAccessor<Data<VecCoord3> > x = *DOF->write(VecId::position());

    x[0] = Coord3(0,10,0);
    x[1] = Coord3(10,0,0);
    x[2] = Coord3(-10*0.5,0,10*0.866);
    x[3] = Coord3(-10*0.5,0,-10*0.866);

    // Tetrahedron uniform mass
    UniformMass3::SPtr mass = sofa::core::objectmodel::New<UniformMass3>();
    groot->addObject(mass);
    mass->setMass(2);
    mass->setName("mass");

    // Tetrahedron topology
    MeshTopology::SPtr topology = sofa::core::objectmodel::New<MeshTopology>();
    topology->setName("mesh topology");
    groot->addObject( topology );
    topology->addTetra(0,1,2,3);

    // Tetrahedron constraints
    FixedConstraint3::SPtr constraints = sofa::core::objectmodel::New<FixedConstraint3>();
    constraints->setName("constraints");
    groot->addObject(constraints);
    constraints->addConstraint(0);

    // Tetrahedron force field
    TetrahedronFEMForceField3::SPtr fem = sofa::core::objectmodel::New<TetrahedronFEMForceField3>();
    fem->setName("FEM");
    groot->addObject(fem);
    fem->setMethod("polar");
    fem->setUpdateStiffnessMatrix(true);
    fem->setYoungModulus(6);

    // Tetrahedron skin
    Node::SPtr skin = groot.get()->createChild("skin");
    // The visual model
    OglModel::SPtr visual = sofa::core::objectmodel::New<OglModel>();
    visual->setName( "visual" );
    visual->load(sofa::helper::system::DataRepository.getFile("mesh/liver-smooth.obj"), "", "");
    visual->setColor("red");
    visual->applyScale(0.7, 0.7, 0.7);
    visual->applyTranslation(1.2, 0.8, 0);
    skin->addObject(visual);

    // The mapping between the tetrahedron (DOF) and the liver (visual)
    BarycentricMapping3_to_Ext3::SPtr mapping = sofa::core::objectmodel::New<BarycentricMapping3_to_Ext3>();
    mapping->setModels(DOF.get(), visual.get());
    mapping->setName( "mapping" );
    skin->addObject(mapping);

    // Display Flags
    sofa::component::visualmodel::VisualStyle::SPtr style = sofa::core::objectmodel::New<sofa::component::visualmodel::VisualStyle>();
    groot->addObject(style);
    sofa::core::visual::DisplayFlags& flags = *style->displayFlags.beginEdit();
    flags.setShowNormals(false);
    flags.setShowInteractionForceFields(false);
    flags.setShowMechanicalMappings(false);
    flags.setShowCollisionModels(false);
    flags.setShowBoundingCollisionModels(false);
    flags.setShowMappings(false);
    flags.setShowForceFields(true);
    flags.setShowWireFrame(true);
    flags.setShowVisualModels(true);
    flags.setShowBehaviorModels(true);
    style->displayFlags.endEdit();

    // Init the scene
    sofa::simulation::tree::getSimulation()->init(groot.get());
    groot->setAnimate(false);


    //=======================================
    // Run the main loop
    sofa::gui::GUIManager::MainLoop(groot);

    return 0;
}
