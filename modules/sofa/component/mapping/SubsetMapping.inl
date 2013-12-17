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
#ifndef SOFA_COMPONENT_MAPPING_SUBSETMAPPING_INL
#define SOFA_COMPONENT_MAPPING_SUBSETMAPPING_INL

#include "SubsetMapping.h"

#include <sofa/core/Mapping.inl>
#include <sofa/component/topology/TopologySubsetData.inl>

namespace sofa
{

namespace component
{

namespace mapping
{


template <class TIn, class TOut>
SubsetMapping<TIn, TOut>::SubsetMapping()
    : Inherit()
    , f_indices( initData(&f_indices, "indices", "list of input indices"))
	, f_definedSubsetSize( initData(&f_definedSubsetSize, -1, "definedSubsetSize", "Size of the wanted subset"))
    , f_first( initData(&f_first, -1, "first", "first index (use if indices are sequential)"))
    , f_last( initData(&f_last, -1, "last", "last index (use if indices are sequential)"))
    , f_radius( initData(&f_radius, (Real)1.0e-5, "radius", "search radius to find corresponding points in case no indices are given"))
    , f_handleTopologyChange( initData(&f_handleTopologyChange, true, "handleTopologyChange", "Enable support of topological changes for indices (disable if it is linked from SubsetTopologicalMapping::pointD2S)"))
    , matrixJ()
    , updateJ(false)
{
}

template <class TIn, class TOut>
SubsetMapping<TIn, TOut>::~SubsetMapping()
{
}


template <class TIn, class TOut>
void SubsetMapping<TIn, TOut>::clear(int reserve)
{
    IndexArray& indices = *f_indices.beginEdit();
    indices.clear();
    if (reserve) indices.reserve(reserve);
    f_indices.endEdit();
}

template <class TIn, class TOut>
int SubsetMapping<TIn, TOut>::addPoint(int index)
{
    IndexArray& indices = *f_indices.beginEdit();
    int i = indices.size();
    indices.push_back(index);
    f_indices.endEdit();
    return i;
}

// Handle topological changes
/*
template <class TIn, class TOut>
void SubsetMapping<TIn, TOut>::handleTopologyChange(core::topology::Topology* t)
{
    core::topology::BaseMeshTopology* topoFrom = this->fromModel->getContext()->getMeshTopology();
    if (t != topoFrom) return;

    std::list<const core::topology::TopologyChange *>::const_iterator itBegin=topoFrom->beginChange();
    std::list<const core::topology::TopologyChange *>::const_iterator itEnd=topoFrom->endChange();
    f_indices.beginEdit()->handleTopologyEvents(itBegin,itEnd,this->fromModel->getX()->size());
    f_indices.endEdit();
}
*/

template <class TIn, class TOut>
void SubsetMapping<TIn, TOut>::init()
{
    unsigned int inSize = this->fromModel->getX()->size();
    if (f_indices.getValue().empty() && f_first.getValue() != -1)
    {
        IndexArray& indices = *f_indices.beginEdit();
        unsigned int first = (unsigned int)f_first.getValue();
        unsigned int last = (unsigned int)f_last.getValue();
        if (first >= inSize)
            first = 0;
        if (last >= inSize)
            last = inSize-1;
        indices.resize(last-first+1);
        for (unsigned int i=0; i<indices.size(); ++i)
            indices[i] = first+i;
        f_indices.endEdit();
    }
    else if (f_indices.getValue().empty())
    {

        // We have to construct the correspondance index
        const InVecCoord& in   = *this->fromModel->getX();
        const OutVecCoord& out = *this->toModel->getX();
        IndexArray& indices = *f_indices.beginEdit();

		int defSize = f_definedSubsetSize.getValue();
		unsigned int unfound_card=0;
		if(defSize!=-1)
			indices.resize(defSize);
		else
			indices.resize(out.size());

        // searching for the first corresponding point in the 'from' model (there might be several ones).
        for (unsigned int i = 0; i < out.size(); ++i)
        {
            bool found = false;
            Real rmax = f_radius.getValue();
            for (unsigned int j = 0;  j < in.size(); ++j )
            {
                Real r = (Real)((out[i] - in[j]).norm());
                if ( r < rmax )
                {
                    indices[i-unfound_card] = j;
                    found = true;
                    rmax = r;
                }
            }
			
			if(!found)
			{
				if(defSize==-1)
				{
					sout<<"ERROR(SubsetMapping): point "<<i<<"="<<out[i]<<" not found in input model within a radius of "<<rmax<<"."<<sendl;
					indices[i] = -1;
				}
				else
					++unfound_card;
			}
        }
        f_indices.endEdit();
    }
    else
    {
        IndexArray& indices = *f_indices.beginEdit();
        for (unsigned int i=0; i<indices.size(); ++i)
        {
            if ((unsigned)indices[i] >= inSize)
            {
                serr << "ERROR(SubsetMapping): incorrect index "<<indices[i]<<" (input size "<<inSize<<")"<<sendl;
                indices.erase(indices.begin()+i);
                --i;
            }
        }
        f_indices.endEdit();
    }
    this->Inherit::init();

    topology = this->getContext()->getMeshTopology();

    // Initialize functions and parameters
    f_indices.createTopologicalEngine(topology);
    f_indices.registerTopologicalData();

    postInit();
}

template <class TIn, class TOut>
void SubsetMapping<TIn, TOut>::postInit()
{
}

template <class TIn, class TOut>
void SubsetMapping<TIn, TOut>::apply ( const core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, OutDataVecCoord& dOut, const InDataVecCoord& dIn )
{
    const IndexArray& indices = f_indices.getValue();

    const InVecCoord& in = dIn.getValue();
    OutVecCoord& out = *dOut.beginEdit();

    out.resize(indices.size());
    for(unsigned int i = 0; i < out.size(); ++i)
    {
		if(indices[i]>=0)
			out[i] = in[ indices[i] ];
    }

    dOut.endEdit();
}

template <class TIn, class TOut>
void SubsetMapping<TIn, TOut>::applyJ( const core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, OutDataVecDeriv& dOut, const InDataVecDeriv& dIn )
{
    const IndexArray& indices = f_indices.getValue();

    const InVecDeriv& in = dIn.getValue();
    OutVecDeriv& out = *dOut.beginEdit();

    out.resize(indices.size());
    for(unsigned int i = 0; i < out.size(); ++i)
    {
		if(indices[i]>=0)
			out[i] = in[ indices[i] ];
    }

    dOut.endEdit();
}

template <class TIn, class TOut>
void SubsetMapping<TIn, TOut>::applyJT ( const core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, InDataVecDeriv& dOut, const OutDataVecDeriv& dIn )
{
    const IndexArray& indices = f_indices.getValue();

    const OutVecDeriv& in = dIn.getValue();
    InVecDeriv& out = *dOut.beginEdit();

    if (indices.empty())
        return;

    if (!in.empty())
    {
        out.resize(this->fromModel->getSize());
    }

    for(unsigned int i = 0; i < in.size(); ++i)
    {
		if(indices[i]>=0)
			out[indices[i]] += in[ i ];
    }

    dOut.endEdit();
}

template <class TIn, class TOut>
void SubsetMapping<TIn, TOut>::applyJT ( const core::ConstraintParams * /*cparams*/ /* PARAMS FIRST */, InDataMatrixDeriv& dOut, const OutDataMatrixDeriv& dIn)
{
    const IndexArray& indices = f_indices.getValue();

    InMatrixDeriv& out = *dOut.beginEdit();
    const OutMatrixDeriv& in = dIn.getValue();

    typename Out::MatrixDeriv::RowConstIterator rowItEnd = in.end();

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
				if(indices[colIt.index()]>=0)
					o.addCol(indices[colIt.index()], colIt.val());
                ++colIt;
            }
        }
    }
    dOut.endEdit();
    //int offset = out.size();
    //out.resize(offset+in.size());

    //const IndexArray& indices = f_indices.getValue();
    //for(unsigned int i = 0; i < in.size(); ++i)
    //{
    //  OutConstraintIterator itOut;
    //  std::pair< OutConstraintIterator, OutConstraintIterator > iter=in[i].data();
    //
    //  for (itOut=iter.first;itOut!=iter.second;itOut++)
    //    {
    //      unsigned int indexIn = itOut->first;
    //      OutDeriv data = (OutDeriv) itOut->second;
    //      out[i+offset].add( indices[indexIn] , data );
    //    }
    //}
}

template<class TIn, class TOut>
const sofa::defaulttype::BaseMatrix* SubsetMapping<TIn, TOut>::getJ()
{
    if (matrixJ.get() == 0 || updateJ)
    {
        const OutVecCoord& out = *this->toModel->getX();
        const InVecCoord& in = *this->fromModel->getX();
        const IndexArray& indices = f_indices.getValue();
        assert(indices.size() == out.size());

        updateJ = false;
        if (matrixJ.get() == 0 ||
            matrixJ->rowBSize() != out.size() ||
            matrixJ->colBSize() != in.size())
        {
            matrixJ.reset(new MatrixType(out.size() * NOut, in.size() * NIn));
        }
        else
        {
            matrixJ->clear();
        }
        for (unsigned i = 0; i < indices.size(); ++i)
        {
			if(indices[i]>=0)
			{
				MBloc& block = *matrixJ->wbloc(i, indices[i], true);
				block.identity();
			}
        }
    }
    return matrixJ.get();
}

} // namespace mapping

} // namespace component

} // namespace sofa

#endif
