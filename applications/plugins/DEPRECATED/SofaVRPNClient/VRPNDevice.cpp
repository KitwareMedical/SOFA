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
 * VRPNDevice.cpp
 *
 *  Created on: 8 sept. 2009
 *      Author: froy
 */

#include "VRPNDevice.h"
#include <sofa/simulation/common/AnimateBeginEvent.h>
#include <sofa/simulation/common/AnimateEndEvent.h>
namespace sofavrpn
{

namespace client
{

VRPNDevice::VRPNDevice()
    : deviceName(initData(&deviceName, std::string("Dummy"), "deviceName", "Name of this device"))
    , serverName(initData(&serverName, std::string("127.0.0.1"), "serverName", "VRPN server name"))
    , serverPort(initData(&serverPort, std::string("3883"), "serverPort", "VRPN server port"))
{
    // TODO Auto-generated constructor stub

}

VRPNDevice::~VRPNDevice()
{
    // TODO Auto-generated destructor stub
}

void VRPNDevice::init()
{
    this->reinit();
}

void VRPNDevice::reinit()
{
    deviceURL = deviceName.getValue() + std::string("@") + serverName.getValue() + std::string(":") + serverPort.getValue();

    bool connected = connect();
    if (!connected)
        std::cout << getName() << " : Not Connected" << std::endl;

}

bool VRPNDevice::connect()
{
    std::cout << "Opening: " << deviceURL << "." << std::endl;

    return connectToServer();
}

void VRPNDevice::handleEvent(sofa::core::objectmodel::Event *event)
{
    if (dynamic_cast<sofa::simulation::AnimateBeginEvent *>(event))
    {
        update();
    }
}

}

}
