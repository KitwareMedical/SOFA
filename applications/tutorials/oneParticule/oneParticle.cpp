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
#include <sofa/helper/ArgumentParser.h>
#include <sofa/simulation/tree/TreeSimulation.h>
#include <sofa/simulation/common/Node.h>
#include <sofa/component/contextobject/Gravity.h>
#include <sofa/component/contextobject/CoordinateSystem.h>
#include <sofa/component/odesolver/EulerSolver.h>
#include <sofa/component/visualmodel/VisualStyle.h>
#include <sofa/core/objectmodel/Context.h>
#include <sofa/core/VecId.h>
#include <sofa/gui/GUIManager.h>
#include <sofa/gui/Main.h>

#include <sofa/helper/system/glut.h>
#include <sofa/helper/accessor.h>



using namespace sofa::simulation::tree;
using sofa::component::odesolver::EulerSolver;
using sofa::core::objectmodel::Data;
using sofa::helper::ReadAccessor;
using sofa::helper::WriteAccessor;
using sofa::core::VecId;

//Using double by default, if you have SOFA_FLOAT in use in you sofa-default.cfg, then it will be FLOAT.
#include <sofa/component/typedef/Sofa_typedef.h>
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

    // The graph root node
    sofa::simulation::setSimulation(new sofa::simulation::tree::TreeSimulation());
    sofa::simulation::Node::SPtr groot = sofa::simulation::getSimulation()->createNewGraph("root");
    groot->setGravity( Coord3(0,-10,0) );

    // One solver for all the graph
    EulerSolver::SPtr solver = sofa::core::objectmodel::New<EulerSolver>();
    solver->setName("solver");
    solver->f_printLog.setValue(false);
    groot->addObject(solver);

    // One node to define the particle
    sofa::simulation::Node::SPtr particule_node = groot.get()->createChild("particle_node");
    // The particule, i.e, its degrees of freedom : a point with a velocity
    MechanicalObject3::SPtr particle = sofa::core::objectmodel::New<MechanicalObject3>();
    particle->setName("particle");
    particule_node->addObject(particle);
    particle->resize(1);
    // get write access the particle positions vector
    WriteAccessor< Data<MechanicalObject3::VecCoord> > positions = *particle->write( VecId::position() );
    positions[0] = Coord3(0,0,0);
    // get write access the particle velocities vector
    WriteAccessor< Data<MechanicalObject3::VecDeriv> > velocities = *particle->write( VecId::velocity() );
    velocities[0] = Deriv3(0,0,0);

    // Its properties, i.e, a simple mass node
    UniformMass3::SPtr mass = sofa::core::objectmodel::New<UniformMass3>();
    mass->setName("mass");
    particule_node->addObject(mass);
    mass->setMass( 1 );

    // Display Flags
    sofa::component::visualmodel::VisualStyle::SPtr style = sofa::core::objectmodel::New<sofa::component::visualmodel::VisualStyle>();
    groot->addObject(style);
    sofa::core::visual::DisplayFlags& flags = *style->displayFlags.beginEdit();
    flags.setShowBehaviorModels(true);
    style->displayFlags.endEdit();

    sofa::simulation::tree::getSimulation()->init(groot.get());
    groot->setAnimate(false);

    //=======================================
    // Run the main loop
    sofa::gui::GUIManager::MainLoop(groot);

    return 0;
}
