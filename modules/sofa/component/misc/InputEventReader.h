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
#ifndef SOFA_COMPONENT_MISC_INPUTEVENTREADER_H
#define SOFA_COMPONENT_MISC_INPUTEVENTREADER_H

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/simulation/common/AnimateBeginEvent.h>
#include <sofa/core/objectmodel/MouseEvent.h>
#include <sofa/core/objectmodel/DataFileName.h>

#ifdef __linux__
#include <linux/input.h>
#include <poll.h>
#endif

namespace sofa
{

namespace component
{

namespace misc
{

#ifndef __linux__
struct input_event {};
#endif

/**
 * @brief InputEventReader Class
 *
 * Reads mouse Linux events from file /dev/input/eventX and propagate them as SOFA MouseEvents.
 */
class InputEventReader : public core::objectmodel::BaseObject
{
public:
    SOFA_CLASS(InputEventReader,core::objectmodel::BaseObject);
protected:
    /**
     * @brief Default Constructor.
     */
    InputEventReader();

    /**
     * @brief Default Destructor.
     */
    ~InputEventReader();
public:
    /**
     * @brief SceneGraph callback initialization method.
     */
    void init();

    /**
     * @brief handle an event.
     *
     * At every simulation step transforms the mouse Linux events in SOFA mouse events and propagates them
     */
    virtual void handleEvent(core::objectmodel::Event *event);

private:

    sofa::core::objectmodel::DataFileName filename; ///< file in which the events are read.
    Data<bool> inverseSense; ///< inverse the sense of the mouvement
    Data<bool> p_printEvent;
    Data<char> p_key1, p_key2;
    Data<bool> p_writeEvents;
    sofa::core::objectmodel::DataFileName p_outputFilename;
    std::ifstream* inFile;
    std::ofstream* outFile;

//	Data<double> timeout;
    int fd; ///< desciptor to open and read the file.

    int deplX, deplY; ///< mouse relative deplacements.

    enum PedalState { LEFT_PEDAL, RIGHT_PEDAL, NO_PEDAL };
    int pedalValue;
    PedalState currentPedalState, oldPedalState;


    /**
     * @brief getInputEvents gets the mouse relative deplacements.
     *
     * This method reads from file /dev/input/eventX and gets the mouse relative deplacements.
     */
    void getInputEvents();
    void manageEvent(const input_event &ev);
};

} // namespace misc

} // namespace component

} // namespace sofa

#endif
