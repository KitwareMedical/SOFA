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
#ifndef SOFA_HELPER_KDTREE_H
#define SOFA_HELPER_KDTREE_H

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/helper/vector.h>
#include <set>


namespace sofa
{

namespace helper
{

/**
  This class implements classical kd tree for nearest neighbors search
  - the tree is rebuild from points by calling build(p)
  - N nearest points from point x (in terms of euclidean distance) are retrieved with getNClosest(distance/index_List , x , N)
  - Caching may be used to speed up retrieval: if dx< (d(n)-d(0))/2, then the closest point is in the n-1 cached points (updateCachedDistances is used to update the n-1 distances)
  see for instance: [zhang92] report and [simon96] thesis for more details
  **/


template<class Coord>
class kdTree
{
public:
    typedef typename Coord::value_type Real;
    enum { dim=Coord::total_size };
    typedef vector<Coord> VecCoord;

    typedef std::pair<Real,unsigned int> distanceToPoint;
    typedef std::set<distanceToPoint> distanceSet;
    typedef typename distanceSet::iterator distanceSetIt;
    typedef std::list<unsigned int> UIlist;


    typedef struct
    {
        unsigned char splitdir; // 0/1/2 -> x/y/z
        unsigned int left;	// index of the left node
        unsigned int right; // index of the right node
    } TREENODE;


    void build(const VecCoord& p);       // update tree (to be used whenever points p have changed)
    void getNClosest(distanceSet &cl, const Coord &x, const unsigned int n);  // get a order set of n distance/index pairs
    unsigned int getClosest(const Coord &x);
    void updateCachedDistances(distanceSet &cl, const Coord &x); // update n-1 distances and set last one to infinite (used in simon96 caching algorithm)

protected :

    const VecCoord* position;
    unsigned int N;
    vector< TREENODE > tree; unsigned int firstNode;

    unsigned int build(UIlist &list, unsigned char direction); // recursive function to build the kdtree
    void closest(distanceSet &cl, const Coord &x, const unsigned int &currentnode);     // recursive function to get closest points
    void closest(distanceToPoint &cl,const Coord &x, const unsigned int &currentnode);  // recursive function to get closest point
};



}
}

#endif
