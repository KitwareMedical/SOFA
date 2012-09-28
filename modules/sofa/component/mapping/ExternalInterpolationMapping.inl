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
#ifndef SOFA_COMPONENT_MAPPING_EXTERNALINTERPOLATIONMAPPING_INL
#define SOFA_COMPONENT_MAPPING_EXTERNALINTERPOLATIONMAPPING_INL

#include "ExternalInterpolationMapping.h"
#include <sofa/core/topology/BaseMeshTopology.h>

namespace sofa
{

namespace component
{

namespace mapping
{


template <class TIn, class TOut>
ExternalInterpolationMapping<TIn, TOut>::ExternalInterpolationMapping()
    : Inherit()
    , f_interpolationIndices( initData(&f_interpolationIndices, "InterpolationIndices", "Table that provides interpolation Indices"))
    , f_interpolationValues( initData(&f_interpolationValues, "InterpolationValues", "Table that provides interpolation Values"))
    , doNotMap(false)
{
}

template <class TIn, class TOut>
ExternalInterpolationMapping<TIn, TOut>::~ExternalInterpolationMapping()
{
}




// Handle topological changes
template <class TIn, class TOut>
void ExternalInterpolationMapping<TIn, TOut>::handleTopologyChange(core::topology::Topology* /*t*/)
{
    /*
     core::topology::BaseMeshTopology* topoFrom = this->fromModel->getContext()->getMeshTopology();
     if (t != topoFrom) return;
     std::list<const core::topology::TopologyChange *>::const_iterator itBegin=topoFrom->beginChange();
     std::list<const core::topology::TopologyChange *>::const_iterator itEnd=topoFrom->endChange();
     f_indices.beginEdit()->handleTopologyEvents(itBegin,itEnd,this->fromModel->getX()->size());
     f_indices.endEdit();
    */
}

template <class TIn, class TOut>
void ExternalInterpolationMapping<TIn, TOut>::init()
{
    // verification of the input table:

    const sofa::helper::vector<sofa::helper::vector< unsigned int > > table_indices = f_interpolationIndices.getValue();
    const sofa::helper::vector<sofa::helper::vector< Real > > table_values = f_interpolationValues.getValue();


    if(table_indices.size() != table_values.size())
    {
        serr<<"WARNING interpolationIndices and interpolationValues do not have the same size ! DoNotMap activated !!"<<sendl;
        doNotMap = true;

        for (unsigned int i=0; i<table_indices.size(); i++)
        {
            if (table_indices[i].size() != table_values[i].size() )
            {
                serr<<"WARNING interpolationIndices and interpolationValues do not have the same size for point "<< i <<" ! DoNotMap activated !!"<<sendl;
                doNotMap = true;
            }
            if (table_indices[i].size() == 0 )
            {
                serr<<"WARNING  no interpolation for point "<< i <<" ! DoNotMap activated !!"<<sendl;
                doNotMap = true;
            }

        }
    }

    this->Inherit::init();
}


template <class TIn, class TOut>
void ExternalInterpolationMapping<TIn, TOut>::apply( typename Out::VecCoord& out, const typename In::VecCoord& in )
{
    if(doNotMap)
        return;

    const sofa::helper::vector<sofa::helper::vector< unsigned int > > table_indices = f_interpolationIndices.getValue();
    const sofa::helper::vector<sofa::helper::vector< Real > > table_values = f_interpolationValues.getValue();


    out.resize(table_indices.size());
    for(unsigned int i = 0; i < out.size(); ++i)
    {
        out[i] = in[ table_indices[i][0] ]  * table_values[i][0];
        for(unsigned int j = 1; j < table_indices[i].size(); j++)
        {
            out[i] += in[ table_indices[i][j] ]  * table_values[i][j];
        }
    }
}

template <class TIn, class TOut>
void ExternalInterpolationMapping<TIn, TOut>::applyJ( typename Out::VecDeriv& out, const typename In::VecDeriv& in )
{
    if(doNotMap)
        return;

    const sofa::helper::vector<sofa::helper::vector< unsigned int > > table_indices = f_interpolationIndices.getValue();
    const sofa::helper::vector<sofa::helper::vector< Real > > table_values = f_interpolationValues.getValue();


    out.resize(table_indices.size());
    for(unsigned int i = 0; i < out.size(); ++i)
    {
        out[i] = in[ table_indices[i][0] ] * table_values[i][0];
        for(unsigned int j = 1; j < table_indices[i].size(); j++)
        {
            out[i] += in[ table_indices[i][j] ] * table_values[i][j];
        }

    }
}

template <class TIn, class TOut>
void ExternalInterpolationMapping<TIn, TOut>::applyJT( typename In::VecDeriv& out, const typename Out::VecDeriv& in )
{
    if(doNotMap)
        return;

    const sofa::helper::vector<sofa::helper::vector< unsigned int > > table_indices = f_interpolationIndices.getValue();
    const sofa::helper::vector<sofa::helper::vector< Real > > table_values = f_interpolationValues.getValue();


    for(unsigned int i = 0; i < in.size(); ++i)
    {
        for(unsigned int j = 0; j < table_indices[i].size(); j++)
        {
            out[table_indices[i][j]] += in[ i ] * table_values[i][j];
        }

    }

}

template <class TIn, class TOut>
void ExternalInterpolationMapping<TIn, TOut>::applyJT( typename In::MatrixDeriv& out, const typename Out::MatrixDeriv& in )
{
    using sofa::helper::vector;

    if(doNotMap)
        return;

    const vector< vector< unsigned int > > table_indices = f_interpolationIndices.getValue();
    const vector< vector< Real > > table_values = f_interpolationValues.getValue();

    typename Out::MatrixDeriv::RowConstIterator rowItEnd = in.end();

    unsigned int i = 0;

    for (typename Out::MatrixDeriv::RowConstIterator rowIt = in.begin(); rowIt != rowItEnd; ++rowIt)
    {
        typename Out::MatrixDeriv::ColConstIterator colIt = rowIt.begin();
        typename Out::MatrixDeriv::ColConstIterator colItEnd = rowIt.end();

        // Creates a constraints if the input constraint is not empty.
        if (colIt != colItEnd)
        {
            typename In::MatrixDeriv::RowIterator o = out.writeLine(rowIt.index());

            while (colIt != colItEnd)
            {
                const unsigned int indexIn = colIt.index();
                const OutDeriv data = colIt.val();

                const unsigned int tIndicesSize = table_indices[indexIn].size();

                for(unsigned int j = 0; j < tIndicesSize; j++)
                {
                    o.addCol(table_indices[indexIn][j] , data * table_values[indexIn][j] );
                }

                ++colIt;
            }
        }

        i++;
    }

    /*if(doNotMap)
    	return;

    const sofa::helper::vector<sofa::helper::vector< unsigned int > > table_indices = f_interpolationIndices.getValue();
    const sofa::helper::vector<sofa::helper::vector< Real > > table_values = f_interpolationValues.getValue();

    int offset = out.size();
    out.resize(offset+in.size());


    for(unsigned int i = 0; i < in.size(); ++i)
    {
        OutConstraintIterator itOut;
        std::pair< OutConstraintIterator, OutConstraintIterator > iter=in[i].data();

        for (itOut=iter.first;itOut!=iter.second;itOut++)
        {
            unsigned int indexIn = itOut->first;
            OutDeriv data = (OutDeriv) itOut->second;

            for(unsigned int j = 0; j < table_indices[i].size();j++){

                out[i+offset].add( table_indices[indexIn][j] , data*table_values[indexIn][j] );
            }
        }
    }*/
}

} // namespace mapping

} // namespace component

} // namespace sofa

#endif
