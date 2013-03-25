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

#ifndef SOFA_COMPONENT_MISC_EVALPOINTSDISTANCE_H
#define SOFA_COMPONENT_MISC_EVALPOINTSDISTANCE_H

#include <sofa/core/behavior/MechanicalState.h>

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/core/objectmodel/DataFileName.h>

#include <sofa/helper/vector.h>

#include <fstream>

namespace sofa
{

namespace component
{

namespace misc
{

using namespace core::behavior;

/** Compute the distance between point/node positions in two objects
*/
template<class TDataTypes>
class EvalPointsDistance: public virtual sofa::core::objectmodel::BaseObject
{

public:
    SOFA_CLASS(SOFA_TEMPLATE(EvalPointsDistance,TDataTypes), sofa::core::objectmodel::BaseObject);

    typedef TDataTypes DataTypes;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::Real Real;
    typedef typename defaulttype::Vector3 Vec3;

    /// Rendering of lines between associated points (activation)
    Data < bool > f_draw;
    /// Output file name
    sofa::core::objectmodel::DataFileName f_filename;
    /// Period between outputs
    Data < double > f_period;
    /// Computed distances
    Data < helper::vector<Real> > dist;
    /// Computed distances (mean, min, max, standard deviation)
    Data < double > distMean, distMin, distMax, distDev;
    /// Relative computed distances (mean, min, max, standard deviation)
    Data < double > rdistMean, rdistMin, rdistMax, rdistDev;

    typedef SingleLink<EvalPointsDistance<DataTypes>,core::behavior::MechanicalState<DataTypes>,BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> LinkMState;
    /// First model mechanical state
    LinkMState mstate1;
    /// Second model mechanical state
    LinkMState mstate2;

protected:
    /** Default constructor
    */
    EvalPointsDistance();
    virtual ~EvalPointsDistance();
public:
    /// Init the computation of the distances
    virtual void init();
    /// Update
    virtual void reinit();
    /// Reset the computation of the distances
    virtual void reset() { reinit(); };

    /** Distance computation */

    /// Get the nodes/points coordinates of the two objects and compute the distances
    virtual SReal eval();
    /// Compute the distances between the two objects
    virtual SReal doEval(const VecCoord& x1, const VecCoord& x2, const VecCoord& x0);


    virtual void handleEvent(sofa::core::objectmodel::Event* event);
    virtual void draw(const core::visual::VisualParams* vparams);
    virtual void doDraw(const VecCoord& x1, const VecCoord& x2);

    /// Retrieve the associated MechanicalState (First model)
    core::behavior::MechanicalState<DataTypes>* getMState1() { return mstate1.get(); }
    core::behavior::BaseMechanicalState* getMechModel1() { return mstate1.get(); }

    /// Retrieve the associated MechanicalState (Second model)
    core::behavior::MechanicalState<DataTypes>* getMState2() { return mstate2.get(); }
    core::behavior::BaseMechanicalState* getMechModel2() { return mstate2; }


    /// Pre-construction check method called by ObjectFactory.
    /// Check that DataTypes matches the MechanicalState.
    template<class T>
    static bool canCreate(T*& obj, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg)
    {
        std::string object1 = arg->getAttribute("object1","@./");
        std::string object2 = arg->getAttribute("object2","@./");
        if (!LinkMState::CheckPath(object1, context))
            return false;
        if (!LinkMState::CheckPath(object2, context))
            return false;

        return core::objectmodel::BaseObject::canCreate(obj, context, arg);
    }

    virtual std::string getTemplateName() const
    {
        return templateName(this);
    }

    static std::string templateName(const EvalPointsDistance<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }

    void setPathToMS1(const std::string &o) { mstate1.setPath(o); }
    void setPathToMS2(const std::string &o) { mstate2.setPath(o); }
    void setPeriod(const double& _dt)      { f_period.setValue(_dt); }

protected:
    /// output file
    std::ofstream* outfile;
    /// time value for the distance computations
    double lastTime;

    sofa::defaulttype::BoundingBox box1;
    sofa::defaulttype::BoundingBox box2;
};

} // namespace misc

} // namespace component

} // namespace sofa

#endif
