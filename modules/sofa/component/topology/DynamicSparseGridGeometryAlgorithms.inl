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
#ifndef SOFA_COMPONENT_TOPOLOGY_DYNAMICSPARSEGRIDGEOMETRYALGORITHMS_INL
#define SOFA_COMPONENT_TOPOLOGY_DYNAMICSPARSEGRIDGEOMETRYALGORITHMS_INL

#include <sofa/component/topology/DynamicSparseGridGeometryAlgorithms.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/component/topology/CommonAlgorithms.h>
#include <sofa/component/container/MechanicalObject.inl>

namespace sofa
{

namespace component
{

namespace topology
{

template < class DataTypes >
void DynamicSparseGridGeometryAlgorithms<DataTypes>::init()
{
    HexahedronSetGeometryAlgorithms<DataTypes>::init();
    this->getContext()->get ( topoContainer );
    if ( !topoContainer )
    {
        serr << "buildTriangleMesh(). Error: can't find the mapping on the triangular topology." << sendl;
        exit(0);
    }
    this->getContext()->get( dof);
    if( !dof)
    {
        serr << "Can not find the dof" << sendl;
        return;
    }
}

template < class DataTypes >
HexaID DynamicSparseGridGeometryAlgorithms<DataTypes>::getTopoIndexFromRegularGridIndex ( unsigned int index, bool& existing )
{
    std::map< unsigned int, BaseMeshTopology::HexaID>::const_iterator it = topoContainer->idInRegularGrid2IndexInTopo.getValue().find( index);
    existing = !(it == topoContainer->idInRegularGrid2IndexInTopo.getValue().end());
    if( !existing)
    {
        //serr << "getTopoIndexFromRegularGridIndex(): Warning ! unexisting given index " << index << " !" << sendl;
        return 0;
    }
    return it->second;
}

template < class DataTypes >
unsigned int DynamicSparseGridGeometryAlgorithms<DataTypes>::getRegularGridIndexFromTopoIndex ( HexaID index )
{
    return topoContainer->idxInRegularGrid.getValue()[ index];
}

template < class DataTypes >
int DynamicSparseGridGeometryAlgorithms<DataTypes>::findNearestElementInRestPos(const Coord& pos, sofa::defaulttype::Vector3& baryC, Real& distance) const
{
    int index = -1;
    distance = 1e10;

    Vec3i resolution = topoContainer->resolution.getValue();
    const sofa::defaulttype::Vec3d& translation = dof->getTranslation();
    Vec3i currentIndex = Vec3i( (int)((pos[0] - translation[0]) / topoContainer->voxelSize.getValue()[0]), (int)((pos[1] - translation[1]) / topoContainer->voxelSize.getValue()[1]), (int)((pos[2] - translation[2]) / topoContainer->voxelSize.getValue()[2]));

//        std::cout << "Find Nearest : " << pos << " ; " << baryC << " distance " << distance << " : " << resolution << " : resolution " << currentIndex << " : currentIndex\n";
    // Projection sur la bbox si l'element est en dehors.
    if( currentIndex[0] < 0) currentIndex[0] = 0;
    if( currentIndex[1] < 0) currentIndex[1] = 0;
    if( currentIndex[2] < 0) currentIndex[2] = 0;
    if( currentIndex[0] > resolution[0]) currentIndex[0] = resolution[0];
    if( currentIndex[1] > resolution[1]) currentIndex[1] = resolution[1];
    if( currentIndex[2] > resolution[2]) currentIndex[2] = resolution[2];

    const std::map< unsigned int, BaseMeshTopology::HexaID>& regular2topo = topoContainer->idInRegularGrid2IndexInTopo.getValue();
    unsigned int regularGridIndex;
    std::map< unsigned int, BaseMeshTopology::HexaID>::const_iterator it;
    for( int k = 0; k < 3; k++)
    {
        if((((int)currentIndex[2])-1+k < 0) || (currentIndex[2]-1+k > resolution[2])) continue;
        for( int j = 0; j < 3; j++)
        {
            if((((int)currentIndex[1])-1+j < 0) || (currentIndex[1]-1+j > resolution[1])) continue;
            for( int i = 0; i < 3; i++)
            {
                if((((int)currentIndex[0])-1+i < 0) || (currentIndex[0]-1+i > resolution[0])) continue;
                regularGridIndex = (currentIndex[0]-1+i) + (currentIndex[1]-1+j)*resolution[0] + (currentIndex[2]-1+k)*resolution[0]*resolution[1];
//              std::cout << regularGridIndex << " Regular Grid Index\n";
                it = regular2topo.find( regularGridIndex);
                if( it != regular2topo.end())
                {
                    const Real d = this->computeElementRestDistanceMeasure(it->second, pos);
//                std::cout << "Distance : " << d << "\n";
                    if(d<distance)
                    {
                        distance = d;
                        index = it->second;
                    }
                }
            }
        }
    }
    if( index == -1)
    {
        // Dans le cas de projection ou autre.... il se peut que la zone ciblée ne contienne pas d'hexahedra, il faut alors tous les parcourrir.
        serr << "findNearestElementInRestPos(). Index not found" << sendl;
        //sout << "findNearestElementInRestPos(). Index not found => Search in all the hexahedra ! SLOW." << sendl;
        //sout << "pos: " << pos << ", currentIndex: " << currentIndex << ", et regular index: " << regularGridIndex << sendl;
        return HexahedronSetGeometryAlgorithms<DataTypes>::findNearestElementInRestPos( pos, baryC, distance);
    }

    distance = this->computeElementRestDistanceMeasure( index, pos);

    baryC = this->computeHexahedronRestBarycentricCoeficients(index, pos);

    return index;
}

} // namespace topology

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENTS_HexahedronSetTOPOLOGY_INL
