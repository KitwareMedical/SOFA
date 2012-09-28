/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 4      *
*                (c) 2006-2009 MGH, INRIA, USTL, UJF, CNRS                    *
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
*                               SOFA :: Tests                                 *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/

#include <sofa/helper/system/atomic.h>
#include <boost/test/auto_unit_test.hpp>

using sofa::helper::system::atomic;

BOOST_AUTO_TEST_CASE(dec_and_test_null)
{
    atomic<int> value(3);
    BOOST_CHECK_EQUAL(value.dec_and_test_null(), false);
    BOOST_CHECK_EQUAL(value, 2);
    BOOST_CHECK_EQUAL(value.dec_and_test_null(), false);
    BOOST_CHECK_EQUAL(value, 1);
    BOOST_CHECK_EQUAL(value.dec_and_test_null(), true);
    BOOST_CHECK_EQUAL(value, 0);
}

BOOST_AUTO_TEST_CASE(compare_and_swap)
{
    atomic<int> value(-1);
    BOOST_CHECK_EQUAL(value.compare_and_swap(-1, 10), -1);
    BOOST_CHECK_EQUAL(value, 10);

    BOOST_CHECK_EQUAL(value.compare_and_swap(5, 25), 10);
    BOOST_CHECK_EQUAL(value, 10);
}
