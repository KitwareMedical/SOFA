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
#ifndef SOFA_COMPONENT_MISC_READTOPOLOGY_H
#define SOFA_COMPONENT_MISC_READTOPOLOGY_H

#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/simulation/common/AnimateBeginEvent.h>
#include <sofa/simulation/common/AnimateEndEvent.h>
#include <sofa/simulation/common/Visitor.h>
#include <sofa/component/component.h>
#include <sofa/core/objectmodel/DataFileName.h>

#ifdef SOFA_HAVE_ZLIB
#include <zlib.h>
#endif

#include <fstream>

namespace sofa
{

namespace component
{

namespace misc
{

/** Read Topology containers informations from file at each timestep
*/
class SOFA_LOADER_API ReadTopology: public core::objectmodel::BaseObject
{
public:
    SOFA_CLASS(ReadTopology,core::objectmodel::BaseObject);

    sofa::core::objectmodel::DataFileName f_filename;
    Data < double > f_interval;
    Data < double > f_shift;
    Data < bool > f_loop;

protected:
    core::topology::BaseMeshTopology* m_topology;
    std::ifstream* infile;
#ifdef SOFA_HAVE_ZLIB
    gzFile gzfile;
#endif
    double nextTime;
    double lastTime;
    double loopTime;

    ReadTopology();

    virtual ~ReadTopology();
public:
    virtual void init();

    virtual void reset();

    void setTime(double time);

    virtual void handleEvent(sofa::core::objectmodel::Event* event);

    void processReadTopology();
    void processReadTopology(double time);

    /// Read the next values in the file corresponding to the last timestep before the given time
    bool readNext(double time, std::vector<std::string>& lines);

    /// Pre-construction check method called by ObjectFactory.
    /// Check that DataTypes matches the MechanicalState.
    template<class T>
    static bool canCreate(T* obj, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg)
    {
        if (dynamic_cast<core::topology::BaseMeshTopology*>(context->getMeshTopology()) == NULL)
            return false;
        return BaseObject::canCreate(obj, context, arg);
    }



    // Tab of 2D elements composition
    Data< helper::vector< helper::fixed_array <unsigned int,2> > > edges;
    Data< helper::vector< helper::fixed_array <unsigned int,3> > > triangles;
    Data< helper::vector< helper::fixed_array <unsigned int,4> > > quads;

    // Tab of 3D elements composition
    Data< helper::vector< helper::fixed_array<unsigned int,4> > > tetrahedra;
    Data< helper::vector< helper::fixed_array<unsigned int,8> > > hexahedra;


};


///Create ReadTopology component in the graph each time needed
class SOFA_LOADER_API ReadTopologyCreator: public Visitor
{
public:
    ReadTopologyCreator(const core::ExecParams* params);
    ReadTopologyCreator(const std::string &n, bool _createInMapping, const core::ExecParams* params, bool i=true, int c=0 );
    virtual Result processNodeTopDown( simulation::Node*  );

    void setSceneName(std::string &n) { sceneName = n;}
    void setCounter(int c) {counterReadTopology = c;};
    void setCreateInMapping(bool b) {createInMapping=b;}
    virtual const char* getClassName() const { return "ReadTopologyCreator"; }
protected:
    std::string sceneName;
    std::string extension;
    bool createInMapping;
    bool init;
    int counterReadTopology; //avoid to have two same files if two Topologies are present with the same name

    void addReadTopology(core::topology::BaseMeshTopology* topology, simulation::Node* gnode);

};

class SOFA_LOADER_API ReadTopologyActivator: public Visitor
{
public:
    ReadTopologyActivator(const core::ExecParams* params, bool active)
        :Visitor(params), state(active) {}
    virtual Result processNodeTopDown( simulation::Node*  );

    bool getTopology() const {return state;};
    void setTopology(bool active) {state=active;};
    virtual const char* getClassName() const { return "ReadTopologyActivator"; }
protected:
    void changeTopologyReader(sofa::component::misc::ReadTopology *rt);

    bool state;
};

class SOFA_LOADER_API ReadTopologyModifier: public simulation::Visitor
{
public:
    ReadTopologyModifier(const core::ExecParams* params /* PARAMS FIRST */, double _time)
        :Visitor(params), time(_time) {}

    virtual Result processNodeTopDown( simulation::Node*  );

    double getTime() const { return time; }
    void setTime(double _time) { time=_time; }
    virtual const char* getClassName() const { return "ReadTopologyModifier"; }
protected:
    void changeTimeReader(sofa::component::misc::ReadTopology *rt) { rt->processReadTopology(time); }

    double time;
};

} // namespace misc

} // namespace component

} // namespace sofa

#endif
