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
/*
 * WiimoteDriver.h
 *
 *  Created on: 8 sept. 2009
 *      Author: froy
 */

#ifndef SOFAVRPNCLIENT_WIIMOTEDRIVER_H_
#define SOFAVRPNCLIENT_WIIMOTEDRIVER_H_

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/DataEngine.h>

#include <VRPNDevice.h>

#include <vrpn/vrpn_Analog.h>

namespace sofavrpn
{

namespace client
{

template<class DataTypes>
class WiimoteDriver : public virtual sofa::core::objectmodel::BaseObject, public virtual sofa::core::DataEngine
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(WiimoteDriver, DataTypes), sofa::core::objectmodel::BaseObject);

    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Point;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::VecCoord VecCoord;

    static const unsigned int WIIMOTE_NUMBER_OF_IR_DOTS;

    //input
    sofa::core::objectmodel::Data<sofa::helper::vector<Real> > f_channels;

    //output
    sofa::core::objectmodel::Data<VecCoord> f_dots;

    //Parameters
    sofa::core::objectmodel::Data<bool> p_viewDots;

    WiimoteDriver();
    virtual ~WiimoteDriver();

//	void init();
//	void reinit();
    void update();
    void draw();

private:

};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFAVRPNCLIENT_WIIMOTEDRIVER_CPP_)
#ifndef SOFA_FLOAT
extern template class SOFA_SOFAVRPNCLIENT_API WiimoteDriver<defaulttype::Vec3dTypes>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
extern template class SOFA_SOFAVRPNCLIENT_API WiimoteDriver<defaulttype::Vec3fTypes>;
#endif //SOFA_DOUBLE
#endif

}

}

#endif /* SOFAVRPNCLIENT_WIIMOTEDRIVER_H_ */
