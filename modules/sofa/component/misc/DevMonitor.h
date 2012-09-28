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
#ifndef SOFA_COMPONENT_MISC_DEVMONITOR_H
#define SOFA_COMPONENT_MISC_DEVMONITOR_H

#include <sofa/helper/vector.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/simulation/common/AnimateEndEvent.h>
#include <sofa/core/DevBaseMonitor.h>

namespace sofa
{

namespace component
{

namespace misc
{

template <class TDataTypes>
class DevMonitor: public virtual core::DevBaseMonitor
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(DevMonitor,TDataTypes), core::DevBaseMonitor);

    typedef TDataTypes DataTypes;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Real Real;
    //typedef double Real;
    typedef typename std::pair< Coord,Real > TData;


    Data < double > f_period;
    Data< sofa::helper::vector< unsigned int > > f_indices;

    DevMonitor():
        f_period( initData(&f_period, 1.0, "period", "period between outputs"))
        , f_indices( initData(&f_indices,"indices","Indices of the points which will be monitored") )
        , lastTime(0)
    {
    }

    sofa::helper::vector<TData> getData()
    {
        sofa::helper::vector<TData> copy;
        copy = data;
        data.clear();
        return copy;
    }

    void handleEvent(sofa::core::objectmodel::Event* event)
    {
        if (dynamic_cast<simulation::AnimateEndEvent*>(event))
        {
            timestamp = getContext()->getTime();
            // write the state using a period
            if (timestamp+getContext()->getDt()/2 >= (lastTime + f_period.getValue()))
            {
                eval();
                lastTime += f_period.getValue();
            }
        }
    }

    virtual std::string getTemplateName() const
    {
        return templateName(this);
    }

    static std::string templateName(const DevMonitor<TDataTypes>* = NULL)
    {
        return TDataTypes::Name();
    }
protected:
    double lastTime;
    double timestamp;
    sofa::helper::vector<TData> data;
};

} // namespace misc

} // namespace component

} // namespace sofa

#endif
