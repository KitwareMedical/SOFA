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
#ifndef SOFA_COMPONENT_ENGINE_DISTANCES_H
#define SOFA_COMPONENT_ENGINE_DISTANCES_H

#if !defined(__GNUC__) || (__GNUC__ > 3 || (_GNUC__ == 3 && __GNUC_MINOR__ > 3))
#pragma once
#endif

#include <sofa/core/DataEngine.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/component/topology/DynamicSparseGridTopologyContainer.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/helper/SVector.h>
#include <sofa/helper/set.h>
#include <sofa/helper/map.h>
#include <sofa/helper/OptionsGroup.h>

#define TYPE_GEODESIC 0
#define TYPE_HARMONIC 1
#define TYPE_STIFFNESS_DIFFUSION 2
#define TYPE_VORONOI 3
#define TYPE_HARMONIC_STIFFNESS 4

namespace sofa
{

namespace component
{

namespace topology
{

class HexahedronSetTopologyContainer;

class HexahedronSetTopologyModifier;

template < class DataTypes >
class DynamicSparseGridGeometryAlgorithms;

}

namespace engine
{

using helper::vector;
using helper::SVector;
using helper::set;
using std::map;
using std::string;
using core::behavior::MechanicalState;
using namespace sofa::component::topology;


/// This class can be overridden if needed for additionnal storage within template specializations.
template<class DataTypes>
class DistancesInternalData
{
public:
};

/**
 * This class computes distances between to set of mechanical objects.
 */
template <class DataTypes>
class Distances : public core::DataEngine
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(Distances,DataTypes),core::DataEngine);

    typedef std::pair< HexaID, double> Distance;
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename sofa::helper::vector< VecCoord > VecVecCoord;
    typedef SVector<SVector<double> > VVD;

protected:
    DistancesInternalData<DataTypes> data;
    friend class DistancesInternalData<DataTypes>;

    Distances ( DynamicSparseGridTopologyContainer* hexaTopoContainer, MechanicalState<DataTypes>* targetPointSet );

    virtual ~Distances() {}

public:
    Data<unsigned int> showMapIndex;
    Data<bool> showDistanceMap;
    Data<bool> showGoalDistanceMap;
    Data<double> showTextScaleFactor;
    Data<bool> showGradientMap;
    Data<double> showGradientsScaleFactor;
    Data<Coord> offset;
    Data<sofa::helper::OptionsGroup> distanceType;
    Data<bool> initTarget;
    Data<int> initTargetStep;
    Data<map<unsigned int, unsigned int> > zonesFramePair;
    Data<double> harmonicMaxValue;

    void init();

    void reinit();

    void update();

    /** \brief Compute the distance map depending ion the distance type.
    *
    * @param elt the point from which the distances are computed.
    * @param beginElts distance until we stop propagating.
    * @param distMax distance until we stop propagating.
    */
    void computeDistanceMap ( VecCoord beginElts = VecCoord(), const double& distMax = 0 );

    /** \brief Add a 'from' element and recompute the map of distances.
    *
    * @param elt the point from which the distances are computed.
    * @param beginElts distance until we stop propagating.
    * @param distMax distance until we stop propagating.
    */
    void addElt ( const Coord& elt, VecCoord beginElts = VecCoord(), const double& distMax = 0 );

    /** \brief Get the distance for a point set using the computed map.
    *
    * @param distances distance for each point of the topology.
    * @param gradients gradient of the distance for each point in the topology.
    * @param point the point from which the distances are computed.
    */
    void getDistances ( VVD& distances, VecVecCoord& gradients, const VecCoord& goals );

    void draw(const core::visual::VisualParams* vparams);

    /// Pre-construction check method called by ObjectFactory.
    ///
    /// This implementation read the object1 and object2 attributes and check
    /// if they are compatible with the input and output model types of this
    /// mapping.
    template<class T>
    static bool canCreate ( T*& obj, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg )
    {
        if ( arg->findObject ( arg->getAttribute ( "hexaContainerPath","../.." ) ) == NULL )
            context->serr << "Cannot create "<<className ( obj ) <<" as the hexas container is missing."<<context->sendl;
        if ( arg->findObject ( arg->getAttribute ( "targetPath",".." ) ) == NULL )
            context->serr << "Cannot create "<<className ( obj ) <<" as the target point set is missing."<<context->sendl;
        if ( dynamic_cast<DynamicSparseGridTopologyContainer*> ( arg->findObject ( arg->getAttribute ( "hexaContainerPath","../.." ) ) ) == NULL )
            return false;
        if ( dynamic_cast<MechanicalState<DataTypes>*> ( arg->findObject ( arg->getAttribute ( "targetPath",".." ) ) ) == NULL )
            return false;
        return true;
    }
    /// Construction method called by ObjectFactory.
    ///
    /// This implementation read the object1 and object2 attributes to
    /// find the input and output models of this mapping.
    template<class T>
    static typename T::SPtr create(T*, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg )
    {
        typename T::SPtr obj = sofa::core::objectmodel::New<T>(
                ( arg?dynamic_cast<DynamicSparseGridTopologyContainer*> ( arg->findObject ( arg->getAttribute ( "hexaContainerPath","../.." ) ) ) :NULL ),
                ( arg?dynamic_cast<MechanicalState<DataTypes>*> ( arg->findObject ( arg->getAttribute ( "targetPath",".." ) ) ) :NULL ) );

        if ( context ) context->addObject ( obj );

        if ( arg )
        {
            if ( arg->getAttribute ( "hexaContainerPath" ) )
            {
                obj->hexaContainerPath.setValue ( arg->getAttribute ( "hexaContainerPath" ) );
                arg->removeAttribute ( "hexaContainerPath" );
            }
            if ( arg->getAttribute ( "targetPath" ) )
            {
                obj->targetPath.setValue ( arg->getAttribute ( "targetPath" ) );
                arg->removeAttribute ( "targetPath" );
            }
            obj->parse ( arg );
        }

        return obj;
    }
    std::string getTemplateName() const
    {
        return templateName(this);
    }
    static std::string templateName(const Distances<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }


private:
    Data<string> fileDistance;
    Data<string> targetPath;
    MechanicalState<DataTypes>* target;

    Data<string> hexaContainerPath;
    DynamicSparseGridTopologyContainer* hexaContainer;
    DynamicSparseGridGeometryAlgorithms< DataTypes >* hexaGeoAlgo;
    const unsigned char * densityValues; // Density values
    const unsigned char * segmentIDData; // Density values

    VVD distanceMap; // distance for each hexa of the grid (topology order)

    /*************************/
    /*   Compute distances   */
    /*************************/
    inline void computeGeodesicalDistance ( const unsigned int& mapIndex, const VecCoord& beginElts, const bool& diffuseAccordingToStiffness, const double& distMax = 0 );
    // Store harmonic coords in the distanceMap structure of the class depending on the fixed values 'hfrom'
    inline void computeHarmonicCoords ( const unsigned int& mapIndex, const vector<HexaID>& hfrom, const bool& useStiffnessMap );
    inline void computeVoronoiDistances( const unsigned int& mapIndex, const VecCoord& beginElts, const double& distMax = 0 );


    /*************************/
    /*         Utils         */
    /*************************/
    inline void findCorrespondingHexas ( vector<HexaID>& hexas, const VecCoord& pointSet ); // Find indices from coord.
    inline void find1DCoord ( unsigned int& hexaID, const Coord& point );
    void getNeighbors ( const HexaID& hexaID, helper::set<HexaID>& neighbors ) const;
    void computeGradients ( const unsigned int mapIndex, vector<double>& distances, VecCoord& gradients, const vector<HexaID>& hexaGoal, const VecCoord& goals );
    inline void addContribution ( double& valueWrite, int& nbTest, const vector<double>& valueRead, const unsigned int& gridID, const int coeff );
    inline void addContribution ( double& valueWrite, int& nbTest, double*** valueRead, const int& x, const int& y, const int& z, const int coeff, const bool& useStiffnessMap );
};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_BUILD_MISC_ENGINE)
#ifndef SOFA_FLOAT
extern template class SOFA_MISC_ENGINE_API Distances<defaulttype::Vec3dTypes>;
//extern template class SOFA_MISC_ENGINE_API Distances<defaulttype::Rigid3dTypes>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
extern template class SOFA_MISC_ENGINE_API Distances<defaulttype::Vec3fTypes>;
//extern template class SOFA_MISC_ENGINE_API Distances<defaulttype::Rigid3fTypes>;
#endif //SOFA_DOUBLE
#endif

} // namespace engine

} // namespace component

} // namespace sofa

#endif
