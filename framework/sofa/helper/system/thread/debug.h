/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *
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
*                              SOFA :: Framework                              *
*                                                                             *
* Authors: The SOFA Team (see Authors.txt)                                    *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_HELPER_SYSTEM_THREAD_DEBUG_H
#define SOFA_HELPER_SYSTEM_THREAD_DEBUG_H

#include <vector>

#include <sofa/helper/system/thread/CTime.h>
#include <string.h>

#include <sofa/helper/helper.h>

namespace sofa
{

namespace helper
{

namespace system
{

namespace thread
{

enum SOFA_HELPER_API TraceLevel
{
    TRACE_DEBUG   = 0,
    TRACE_INFO    = 1,
    TRACE_ERROR   = 2,
    TRACE_WARNING = 3,
};

class SOFA_HELPER_API Trace
{
    static int mTraceLevel;
    static int mNbInstance;
public:
    Trace();

    static void setTraceLevel(int level);
    static void print(int level, const char *chaine);
};


class SOFA_HELPER_API TraceProfile
{
public:
    int index;
    char *name;
    int size;
    int *times;
    int sum;

    ctime_t beginTime;
    ctime_t endTime;

    TraceProfile(const char *name, int index, int size);
    ~TraceProfile();

    void addTime(int instant, int time);

    void begin();
    void end(int instant);
};



#ifdef TRACE_ENABLE

#define TRACE_LEVEl(level) { Trace::setTraceLevel(level); }
#define TRACE(level, chaine){ Trace::print((level), (char*)(chaine)); }

#else


#define TRACE_LEVEl(level) { }
#define TRACE(level, chaine){ }

#endif
} // namespace thread

} // namespace system

} // namespace helper

} // namespace sofa

#endif
