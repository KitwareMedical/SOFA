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
#ifndef SOFA_COMPONENT_MISC_WRITESTATE_H
#define SOFA_COMPONENT_MISC_WRITESTATE_H

#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/BaseMechanicalState.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/simulation/common/AnimateBeginEvent.h>
#include <sofa/simulation/common/AnimateEndEvent.h>
#include <sofa/defaulttype/DataTypeInfo.h>
#include <sofa/simulation/common/Visitor.h>
#include <sofa/component/component.h>

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

/** Write State vectors to file at a given set of time instants
 * A period can be etablished at the last time instant
 * The DoFs to print can be chosen using DOFsX and DOFsV
 * Stop to write the state if the kinematic energy reach a given threshold (stopAt)
 * The energy will be measured at each period determined by keperiod
*/
class SOFA_EXPORTER_API WriteState: public core::objectmodel::BaseObject
{
public:
    SOFA_CLASS(WriteState,core::objectmodel::BaseObject);

    sofa::core::objectmodel::DataFileName f_filename;
    Data < bool > f_writeX;
    Data < bool > f_writeX0;
    Data < bool > f_writeV;
    Data < bool > f_writeF;
    Data < double > f_interval;
    Data < helper::vector<double> > f_time;
    Data < double > f_period;
    Data < helper::vector<unsigned int> > f_DOFsX;
    Data < helper::vector<unsigned int> > f_DOFsV;
    Data < double > f_stopAt;
    Data < double > f_keperiod;

protected:
    core::behavior::BaseMechanicalState* mmodel;
    std::ofstream* outfile;
#ifdef SOFA_HAVE_ZLIB
    gzFile gzfile;
#endif
    unsigned int nextTime;
    double lastTime;
    bool kineticEnergyThresholdReached;
    double timeToTestEnergyIncrease;
    double savedKineticEnergy;


    WriteState();

    virtual ~WriteState();
public:
    virtual void init();

    virtual void reset();

    virtual void handleEvent(sofa::core::objectmodel::Event* event);


    /// Pre-construction check method called by ObjectFactory.
    /// Check that DataTypes matches the MechanicalState.
    template<class T>
    static bool canCreate(T*& obj, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg)
    {
        if (dynamic_cast<core::behavior::BaseMechanicalState*>(context->getMechanicalState()) == NULL)
            return false;
        return BaseObject::canCreate(obj, context, arg);
    }

};

///Create WriteState component in the graph each time needed
class SOFA_EXPORTER_API WriteStateCreator: public Visitor
{
public:
    WriteStateCreator(const core::ExecParams* params);
    WriteStateCreator(const core::ExecParams* params, const std::string &n, bool _recordX, bool _recordV, bool _recordF, bool _createInMapping, int c=0);
    virtual Result processNodeTopDown( simulation::Node*  );

    void setSceneName(std::string &n) { sceneName = n; }
    void setRecordX(bool b) {recordX=b;}
    void setRecordV(bool b) {recordV=b;}
    void setRecordF(bool b) {recordF=b;}
    void setCreateInMapping(bool b) { createInMapping=b; }
    void setCounter(int c) { counterWriteState = c; }
    virtual const char* getClassName() const { return "WriteStateCreator"; }
protected:
    std::string sceneName;
    std::string extension;
    bool recordX,recordV,recordF;
    bool createInMapping;

    int counterWriteState; //avoid to have two same files if two mechanical objects has the same name

    void addWriteState(sofa::core::behavior::BaseMechanicalState*ms, simulation::Node* gnode);

};

class SOFA_EXPORTER_API WriteStateActivator: public simulation::Visitor
{
public:
    WriteStateActivator( const core::ExecParams* params /* PARAMS FIRST */, bool active) : Visitor(params), state(active) {}
    virtual Result processNodeTopDown( simulation::Node*  );

    bool getState() const { return state; }
    void setState(bool active) { state=active; }
    virtual const char* getClassName() const { return "WriteStateActivator"; }
protected:
    void changeStateWriter(sofa::component::misc::WriteState *ws);

    bool state;
};

} // namespace misc

} // namespace component

} // namespace sofa

#endif
