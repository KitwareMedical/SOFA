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
*                               SOFA :: Plugins                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef PLUGINS_PIM_SCULPTBODYPERFORMER_H
#define PLUGINS_PIM_SCULPTBODYPERFORMER_H

#include <sofa/component/topology/TetrahedronSetTopologyContainer.h>
#include <sofa/component/misc/MeshTetraStuffing.h>
#include <sofa/component/collision/TriangleModel.h>
#include <sofa/component/forcefield/TetrahedralCorotationalFEMForceField.h>
#include <sofa/component/odesolver/EulerImplicitSolver.h>
#include <sofa/component/linearsolver/CGLinearSolver.h>
#include <sofa/component/mass/UniformMass.h>
#include <sofa/component/mass/DiagonalMass.h>
#include <sofa/component/constraint/FixedConstraint.h>
#include <plugins/CGALPlugin/MeshGenerationFromPolyhedron.h>
#include <sofa/gui/PickHandler.h>
#include <set>
#include <sofa/simulation/tree/GNode.h>
//#include "EventManager.h"
#include "ComputeMeshIntersection.h"
#include <sofa/component/forcefield/StiffSpringForceField.h>
#include <sofa/component/visualmodel/OglModel.h>

namespace plugins
{

namespace pim
{

using sofa::component::linearsolver::GraphScatteredMatrix;
using sofa::component::linearsolver::GraphScatteredVector;
using namespace sofa::component::collision;
using namespace sofa;
using namespace sofa::component;
using namespace sofa::gui;

class SculptBodyPerformerConfiguration
{
public:
    SculptBodyPerformerConfiguration() {}
    void setForce(double f) {force=f;}
    void setScale(double s) {scale=s;}
    void setMass(double m) {mass=m;}
    void setStiffness(double s) {stiffness=s;}
    void setDamping(double d) {damping=d;}
    void setCheckedFix(bool b) {checkedFix = b;}
    void setCheckedInflate(bool b) {checkedInflate = b;}
    void setCheckedDeflate(bool b) {checkedDeflate = b;}

protected:
    SReal force;
    SReal scale;
    SReal mass;
    SReal stiffness;
    SReal damping;
    bool checkedFix, checkedInflate, checkedDeflate;
};

template <class DataTypes>
class SculptBodyPerformer: public TInteractionPerformer<DataTypes>, public SculptBodyPerformerConfiguration
{
    typedef typename DataTypes::Coord                                         Coord;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef sofa::component::topology::PointSubset SetIndex;
    typedef sofa::component::collision::RayModel MouseCollisionModel;

public:
    SculptBodyPerformer(BaseMouseInteractor *i);
    ~SculptBodyPerformer() {};

    void start();
    void execute();
    void end();
    void draw();
    void animate();

protected:

    void computeNeighborhood();
    void inflate();
    unsigned int computeNearestVertexFromPickedPoint();
    Coord computeNormal(const unsigned int& index);
    void updatePosition(const unsigned int& index, const Coord& normal, const double& W);
    void createTetraVolume(BaseMeshTopology* mesh);
    void stopSimulation();
    void saveCurrentStateAndMesh();
    void ComputeIntersectonLayer();
    void CreateNewMechanicalState(MechanicalState<DataTypes>* fatmstate, BaseMeshTopology* fatMesh);
    void UpdateNewMechanicalState(MechanicalState<DataTypes>* fatmstate, BaseMeshTopology* fatMesh);

    sofa::simulation::tree::GNode* createCollisionNode(sofa::simulation::tree::GNode* parentNode);
    void createVisualNode(sofa::simulation::tree::GNode* parentNode);
    sofa::simulation::tree::GNode* createSubsetMapping(sofa::simulation::tree::GNode* parentNode, MechanicalState<DataTypes>* subsetState,
            BaseMeshTopology* subsetTopology, double m, double ymodulus, double damping);
    void resetForce(core::objectmodel::BaseNode* node,  int index);

    std::set<unsigned int> vertexNeighborhood, vertexInInfluenceZone, alreadyCheckedVertex, drawTriangles;
    BodyPicked picked;
    core::behavior::MechanicalState<DataTypes>* mstateCollision;
//        EventManager<DataTypes> eventManager;
    MechanicalState<DataTypes>* createContactSurfaceMapping(sofa::simulation::tree::GNode* parentNode);
    simulation::Node* root, *visualRoot;
    sofa::helper::vector<sofa::core::behavior::OdeSolver*> solvers;
    sofa::helper::vector<sofa::component::linearsolver::CGLinearSolver<sofa::component::linearsolver::GraphScatteredMatrix,sofa::component::linearsolver::GraphScatteredVector>*> linear;
    core::objectmodel::BaseContext* bc;
    ComputeMeshIntersection<DataTypes>* cmi;
    bool startSculpt;
};



#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_COLLISION_SCULPTBODYPERFORMER_CPP)
extern template class SOFA_COMPONENT_MISC_API SculptBodyPerformer<defaulttype::Vec3Types>;
#endif


}
}

#endif
