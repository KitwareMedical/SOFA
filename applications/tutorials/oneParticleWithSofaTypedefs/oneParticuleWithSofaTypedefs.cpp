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

// launch sofaTypedefs.exe to generate sofa.h
#include <sofa/sofa.h>
#include <sofa/helper/ArgumentParser.h>
#include <sofa/simulation/tree/TreeSimulation.h>
#include <sofa/simulation/tree/GNode.h>
#include <sofa/gui/GUIManager.h>
#include <sofa/gui/Main.h>

using namespace sofa::simulation::tree;
using namespace sofa::component::odesolver;
using namespace sofa::component::container;
using namespace sofa::component::mass;
using namespace sofa::component::visualmodel;
using namespace sofa::component::collision;
using namespace sofa::component::mapping;
using namespace sofa::component::constraintset;
using namespace sofa::component::projectiveconstraintset;
using namespace sofa::component::controller;
using namespace sofa::component::forcefield;
using namespace sofa::component::fem::material;
using namespace sofa::component::fem::forcefield;
using namespace sofa::component::animationloop;
using namespace sofa::component::linearsolver;
using namespace sofa::component::topology;
using namespace sofa::component::interactionforcefield;
using namespace sofa::component::engine;
using namespace sofa::component::behaviormodel::eulerianfluid;
using namespace sofa::component::misc;
using namespace sofa::component::configurationsetting;
using namespace sofa::component::loader;
using namespace sofa::core::objectmodel;
using namespace sofa::core::visual;
using namespace sofa::helper;

/*using sofa::component::odesolver::EulerSolver;
using sofa::component::container::MechanicalObject_Vec3f;
using sofa::component::mass::UniformMass_Vec3f;
using sofa::core::objectmodel::New;
using sofa::core::visual::DisplayFlags;
using sofa::helper::WriteAccessor;*/

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
    GNode::SPtr groot = New<GNode>();
    groot->setName( "root" );
    groot->setGravity( Vec3f(0,-10,0) );

    // One solver for all the graph
    EulerSolver::SPtr solver = New<EulerSolver>();
    solver->setName("solver");
    solver->f_printLog.setValue(false);
    groot->addObject(solver);

    // One node to define the particle
    GNode::SPtr particule_node = New<GNode>("particle_node", groot.get());

    // The particule, i.e, its degrees of freedom : a point with a velocity
    MechanicalObject_Vec3f::SPtr particle = New<MechanicalObject_Vec3f>();
    particle->setName("particle");
    particule_node->addObject(particle);
    particle->resize(1);
    // get write access the particle positions vector
    WriteAccessor< Data<MechanicalObject_Vec3f::VecCoord> > positions = *particle->write( VecId::position() );
    positions[0] = MechanicalObject_Vec3f::Coord(0,0,0);
    // get write access the particle velocities vector
    WriteAccessor< Data<MechanicalObject_Vec3f::VecDeriv> > velocities = *particle->write( VecId::velocity() );
    velocities[0] = MechanicalObject_Vec3f::Deriv(0,0,0);

    // Its properties, i.e, a simple mass node
    UniformMass_Vec3f::SPtr mass = New<UniformMass_Vec3f>();
    mass->setName("mass");
    particule_node->addObject(mass);
    mass->setMass( 1 );

    // Display Flags
    VisualStyle::SPtr style = New<VisualStyle>();
    groot->addObject(style);
    DisplayFlags& flags = *style->displayFlags.beginEdit();
    flags.setShowBehaviorModels(true);
    style->displayFlags.endEdit();

    sofa::simulation::tree::getSimulation()->init(groot.get());
    groot->setAnimate(false);

    //=======================================
    // Run the main loop
    sofa::gui::GUIManager::MainLoop(groot);

    return 0;
}
