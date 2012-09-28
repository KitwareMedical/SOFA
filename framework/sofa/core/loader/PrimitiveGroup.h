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

#ifndef SOFA_CORE_LOADER_PRIMITIVEGROUP_H_
#define SOFA_CORE_LOADER_PRIMITIVEGROUP_H_

#include <sofa/core/core.h>
#include <sofa/core/loader/Material.h>

namespace sofa
{

namespace core
{

namespace loader
{

class PrimitiveGroup
{
public:
    int p0, nbp;
    std::string materialName;
    std::string groupName;
    int materialId;
    inline friend std::ostream& operator << (std::ostream& out, const PrimitiveGroup &g)
    {
        out << g.groupName << " " << g.materialName << " " << g.materialId << " " << g.p0 << " " << g.nbp;
        return out;
    }
    inline friend std::istream& operator >> (std::istream& in, PrimitiveGroup &g)
    {
        in >> g.groupName >> g.materialName >> g.materialId >> g.p0 >> g.nbp;
        return in;
    }

    bool operator <(const PrimitiveGroup& p) const
    {
        return p0 < p.p0;
    }

    PrimitiveGroup() : p0(0), nbp(0), materialId(-1) {}
    PrimitiveGroup(int p0, int nbp, std::string materialName, std::string groupName, int materialId) : p0(p0), nbp(nbp), materialName(materialName), groupName(groupName), materialId(materialId) {}
};

} // namespace loader

} // namespace core

} // namespace sofa

#endif /* SOFA_CORE_LOADER_PRIMITIVEGROUP_H_ */
