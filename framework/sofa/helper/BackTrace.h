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
#ifndef SOFA_HELPER_BACKTRACE_H
#define SOFA_HELPER_BACKTRACE_H

#include <sofa/helper/helper.h>

namespace sofa
{

namespace helper
{

class SOFA_HELPER_API BackTrace
{
public:
    /// Dump current backtrace to stderr.
    /// Currently only works on Linux. NOOP on other architectures.
    static void dump();

    /// Enable dump of backtrace when a signal is received.
    /// Useful to have information about crashes without starting a debugger (as it is not always easy to do, i.e. for parallel/distributed applications).
    /// Currently only works on Linux. NOOP on other architectures
    static void autodump();

protected:

    /// Callback for signals
    static void sig(int sig);
};

} // namespace helper

} // namespace sofa

#endif
