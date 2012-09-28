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
#include <sofa/component/linearsolver/EigenMatrixManipulator.h>
#include <sofa/core/visual/VisualParams.h>

namespace sofa
{

namespace component
{

namespace linearsolver
{

LLineManipulator& LLineManipulator::addCombination(unsigned int idxConstraint, SReal factor)
{
    _data.push_back(std::make_pair(idxConstraint, factor));
    return *this;
}

void LMatrixManipulator::init(const SparseMatrixEigen& L)
{
    const unsigned int numConstraint=L.rows();
    const unsigned int numDofs=L.cols();
    LMatrix.resize(numConstraint,SparseVectorEigen(numDofs));
    for (unsigned int i=0; i<LMatrix.size(); ++i) LMatrix[i].reserve(numDofs*3/10);
    for (int k=0; k<L.outerSize(); ++k)
    {
        for (SparseMatrixEigen::InnerIterator it(L,k); it; ++it)
        {
            const unsigned int row=it.row();
            const unsigned int col=it.col();
            const SReal value=it.value();
            LMatrix[row].insert(col)=value;
        }
    }
    for (unsigned int i=0; i<LMatrix.size(); ++i) LMatrix[i].finalize();
}



void LMatrixManipulator::buildLMatrix(const helper::vector<LLineManipulator> &lines, SparseMatrixEigen& matrix) const
{
    for (unsigned int l=0; l<lines.size(); ++l)
    {
        const LLineManipulator& lManip=lines[l];
        SparseVectorEigen vector;
        lManip.buildCombination(LMatrix,vector);
        matrix.startVec(l);
        for (SparseVectorEigen::InnerIterator it(vector); it; ++it)
        {
            matrix.insertBack(l,it.index())=it.value();
        }
    }
}


helper::vector< SparseVectorEigen > LMatrix;




}
}
}
