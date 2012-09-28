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
#ifndef SOFA_COMPONENT_TOPOLOGY_MULTILEVELHEXAHEDRONSETTOPOLOGYCONTAINER_H
#define SOFA_COMPONENT_TOPOLOGY_MULTILEVELHEXAHEDRONSETTOPOLOGYCONTAINER_H

#include <sofa/component/topology/HexahedronSetTopologyContainer.h>
#include <sofa/component/topology/TopologyData.h>
#include <sofa/core/topology/Topology.h>
#include <sofa/defaulttype/Vec.h>
#include <set>

namespace sofa
{
namespace core
{
namespace topology
{
class TopologyChange;
}
}

namespace component
{

namespace topology
{
class MultilevelHexahedronSetTopologyModifier;

using sofa::defaulttype::Vec;
using sofa::core::topology::TopologyChange;

class SOFA_NON_UNIFORM_FEM_API MultilevelHexahedronSetTopologyContainer : public HexahedronSetTopologyContainer
{
    friend class MultilevelHexahedronSetTopologyModifier;

public:
    SOFA_CLASS(MultilevelHexahedronSetTopologyContainer,HexahedronSetTopologyContainer);

    typedef Vec<3,int>			Vec3i;
protected:
    MultilevelHexahedronSetTopologyContainer();

    MultilevelHexahedronSetTopologyContainer(const helper::vector< Hexahedron > &hexahedra);

    virtual ~MultilevelHexahedronSetTopologyContainer();
public:
    virtual void init();

    virtual void clear();

    void getHexaNeighbors(const unsigned int hexaId,
            helper::vector<unsigned int> &neighbors);

    void getHexaFaceNeighbors(const unsigned int hexaId,
            const unsigned int faceId,
            helper::vector<unsigned int> &neighbors);

    void getHexaVertexNeighbors(const unsigned int hexaId,
            const unsigned int vertexId,
            helper::vector<unsigned int> &neighbors);

    void addTopologyChangeFine(const TopologyChange *topologyChange)
    {
        m_changeListFine.push_back(topologyChange);
    }

    void resetTopologyChangeListFine()
    {
        for(std::list<const TopologyChange *>::iterator it = m_changeListFine.begin();
            it != m_changeListFine.end(); ++it)
        {
            delete (*it);
        }
        m_changeListFine.clear();
    }

    std::list<const TopologyChange *>::const_iterator beginChangeFine() const
    {
        return m_changeListFine.begin();
    }

    std::list<const TopologyChange *>::const_iterator endChangeFine() const
    {
        return m_changeListFine.end();
    }

    const std::list<const TopologyChange *>& getChangeListFine() const
    {
        return m_changeListFine;
    }

    int getLevel() const {return _level.getValue();}



    const Vec3i& getCoarseResolution() const { return _coarseResolution; }

    bool getHexaContainsPosition(const unsigned int hexaId, const defaulttype::Vector3& baryC) const;

    const Vec3i& getHexaIdxInCoarseRegularGrid(const unsigned int hexaId) const;
    int getHexaIdInCoarseRegularGrid(const unsigned int hexaId) const;

    const Vec3i& getHexaIdxInFineRegularGrid(const unsigned int hexaId) const;
    int getHexaIdInFineRegularGrid(const unsigned int hexaId) const;

    // gets a vector of fine hexahedra inside a specified coarse hexa
    int getHexaChildren(const unsigned int hexaId, helper::vector<unsigned int>& children) const;

    // gets a coarse hexa for a specified fine hexa
    int getHexaParent(const unsigned int hexaId) const;

    int getHexaInFineRegularGrid(const Vec3i& id) const;

    const std::set<Vec3i>& getHexaVoxels(const unsigned int hexaId) const;

    Data<int> _level;
    Data<Vec3i>	fineResolution;		// width, height, depth (number of hexa in each direction)
    Data<helper::vector<unsigned int> > hexaIndexInRegularGrid;

private:
    void setCoarseResolution(const Vec3i& res) { _coarseResolution = res; }

    void connectionToNodeAdjacency(const Vec3i& connection, std::map<unsigned int, unsigned int>& nodeMap) const;

    class SOFA_NON_UNIFORM_FEM_API Component
    {
    public:
        Component(const Vec3i& id, const std::set<Vec3i>& voxels);
        virtual ~Component();

        bool isEmpty() const;

        bool isConnected(const Component* other) const;
        bool getConnection(const Component* other, Vec3i& connection) const;
        bool merge(Component* other);

        void split(std::set<Component*>& newComponents);

        void clear();
        void removeVoxels(const std::set<Vec3i>& voxels);

        bool hasVoxel(const Vec3i& voxel) const;

        const Vec3i& getVoxelId() const {return _id;}

        bool isStronglyConnected() const;

        int getLevel() const;

        inline friend std::ostream& operator<< (std::ostream& out, const Component* /*t*/)
        {
            return out;
        }

        inline friend std::istream& operator>>(std::istream& in, Component* /*t*/)
        {
            return in;
        }

    private:
        Component(const Vec3i& id);
        bool isConnected(const std::set<Vec3i>&, const Vec3i&) const;

    public:
        Component*				_parent;
        std::set<Component*>	_children;
        std::set<Vec3i>			_voxels;

    private:
        Vec3i					_id;		// voxel id in the corresponding level
    };

private:


    std::list<const TopologyChange *>	m_changeListFine;

    HexahedronData<sofa::helper::vector<Component*> >		_coarseComponents;	// map between hexahedra and components - coarse
    HexahedronData<sofa::helper::vector<Component*> >		_fineComponents;	// map between hexahedra and components - fine

    // the fine mesh must be a regular grid - store its parameters here

    Vec3i	_coarseResolution;

    sofa::helper::vector<Component*>	_fineComponentInRegularGrid;
};

/** notifies change in the multilevel structure other than adding or removing coarse hexahedra */
class SOFA_NON_UNIFORM_FEM_API MultilevelModification : public core::topology::TopologyChange
{
public:
    static const int MULTILEVEL_MODIFICATION = core::topology::TOPOLOGYCHANGE_LASTID + 1;

    typedef Vec<3,int>	Vec3i;

    MultilevelModification(const sofa::helper::vector<unsigned int>& _tArray,
            const std::map<unsigned int, std::list<Vec3i> >& removedVoxels)
        : core::topology::TopologyChange((core::topology::TopologyChangeType) MULTILEVEL_MODIFICATION)
        , _modifiedHexahedraArray(_tArray)
        , _removedFineVoxels(removedVoxels)
    {}

    const sofa::helper::vector<unsigned int> &getArray() const
    {
        return _modifiedHexahedraArray;
    }

    const std::list<Vec3i> &getRemovedVoxels(const unsigned int hexaId) const
    {
        std::map<unsigned int, std::list<Vec3i> >::const_iterator it = _removedFineVoxels.find(hexaId);
        if(it != _removedFineVoxels.end())
            return it->second;
        else
            return __dummyList;
    }

    unsigned int getNbModifiedHexahedra() const
    {
        return _modifiedHexahedraArray.size();
    }

private:
    sofa::helper::vector<unsigned int>		_modifiedHexahedraArray;
    std::map<unsigned int, std::list<Vec3i> > _removedFineVoxels;

    const std::list<Vec3i>	__dummyList;
};

} // namespace topology

} // namespace component

} // namespace sofa

#endif
