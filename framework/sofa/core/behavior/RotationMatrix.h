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
#ifndef SOFA_COMPONENT_LINEARSOLVER_ROTATIONMATRIX_H
#define SOFA_COMPONENT_LINEARSOLVER_ROTATIONMATRIX_H

#include <sofa/component/linearsolver/SparseMatrix.h>

namespace sofa
{

namespace component
{

namespace linearsolver
{

/// Direct linear solver based on Sparse LDL^T factorization, implemented with the CSPARSE library
template<class TReal>
class RotationMatrix : public defaulttype::BaseMatrix
{
public:
    typedef TReal Real;

    virtual unsigned int rowSize(void) const
    {
        return data.size()/3;
    }

    /// Number of columns
    virtual unsigned int colSize(void) const
    {
        return data.size()/3;
    }

    /// Read the value of the element at row i, column j (using 0-based indices)
    virtual SReal element(int i, int j) const
    {
        int bd = j-(i/3)*3;
        if ((bd<0) || (bd>2)) return 0.0 ;

        return data[i*3+bd];
    }

    /// Resize the matrix and reset all values to 0
    virtual void resize(int nbRow, int nbCol)
    {
        if (nbRow!=nbCol) return;
        data.resize(nbRow*3);
    }

    /// Reset all values to 0
    virtual void clear()
    {
        data.clear();
    }

    /// Write the value of the element at row i, column j (using 0-based indices)
    virtual void set(int i, int j, double v)
    {
        int bd = (i/3)*3;
        if ((j<bd) || (j>bd+2)) return;
        data[i*3+j-bd] = (Real)v;
    }

    /// Add v to the existing value of the element at row i, column j (using 0-based indices)
    virtual void add(int i, int j, double v)
    {
        int bd = (i/3)*3;
        if ((j<bd) || (j>bd+2)) return;

        data[i*3+j-bd] += (Real)v;
    }

    virtual helper::vector<Real> & getVector()
    {
        return data;
    }

    virtual void opMulV(defaulttype::BaseVector* result, const defaulttype::BaseVector* v) const
    {
        //Solve lv = R * lvR
        size_t k = 0,l = 0;
        while (k < data.size())
        {
            result->set(l+0,data[k + 0] * v->element(l+0) + data[k + 1] * v->element(l+1) + data[k + 2] * v->element(l+2));
            result->set(l+1,data[k + 3] * v->element(l+0) + data[k + 4] * v->element(l+1) + data[k + 5] * v->element(l+2));
            result->set(l+2,data[k + 6] * v->element(l+0) + data[k + 7] * v->element(l+1) + data[k + 8] * v->element(l+2));
            l+=3;
            k+=9;
        }
    }

    virtual void opMulTV(defaulttype::BaseVector* result, const defaulttype::BaseVector* v) const
    {
        size_t k = 0,l = 0;
        while (k < data.size())
        {
            result->set(l+0,data[k + 0] * v->element(l+0) + data[k + 3] * v->element(l+1) + data[k + 6] * v->element(l+2));
            result->set(l+1,data[k + 1] * v->element(l+0) + data[k + 4] * v->element(l+1) + data[k + 7] * v->element(l+2));
            result->set(l+2,data[k + 2] * v->element(l+0) + data[k + 5] * v->element(l+1) + data[k + 8] * v->element(l+2));
            l+=3;
            k+=9;
        }
    }

    /// multiply the transpose current matrix by m matrix and strore the result in m
    virtual void opMulTM(defaulttype::BaseMatrix * bresult,defaulttype::BaseMatrix * bm)
    {
        if (RotationMatrix<Real> * m = dynamic_cast<RotationMatrix<Real> * >(bm))
        {
            if (RotationMatrix<Real> * result = dynamic_cast<RotationMatrix<Real> * >(bresult))
            {
                Real tmp[9];
                unsigned datSz = data.size() < m->data.size() ? data.size() : m->data.size();
                unsigned minSz = datSz < result->data.size() ? datSz : result->data.size();

                for (size_t i=0; i<minSz; i+=9)
                {
                    tmp[0] = data[i+0] * m->data[i+0] + data[i+1] * m->data[i+1] + data[i+2] * m->data[i+2];
                    tmp[1] = data[i+0] * m->data[i+3] + data[i+1] * m->data[i+4] + data[i+2] * m->data[i+5];
                    tmp[2] = data[i+0] * m->data[i+6] + data[i+1] * m->data[i+7] + data[i+2] * m->data[i+8];

                    tmp[3] = data[i+3] * m->data[i+0] + data[i+4] * m->data[i+1] + data[i+5] * m->data[i+2];
                    tmp[4] = data[i+3] * m->data[i+3] + data[i+4] * m->data[i+4] + data[i+5] * m->data[i+5];
                    tmp[5] = data[i+3] * m->data[i+6] + data[i+4] * m->data[i+7] + data[i+5] * m->data[i+8];

                    tmp[6] = data[i+6] * m->data[i+0] + data[i+7] * m->data[i+1] + data[i+8] * m->data[i+2];
                    tmp[7] = data[i+6] * m->data[i+3] + data[i+7] * m->data[i+4] + data[i+8] * m->data[i+5];
                    tmp[8] = data[i+6] * m->data[i+6] + data[i+7] * m->data[i+7] + data[i+8] * m->data[i+8];

                    result->data[i+0] = tmp[0]; result->data[i+1] = tmp[1]; result->data[i+2] = tmp[2];
                    result->data[i+3] = tmp[3]; result->data[i+4] = tmp[4]; result->data[i+5] = tmp[5];
                    result->data[i+6] = tmp[6]; result->data[i+7] = tmp[7]; result->data[i+8] = tmp[8];
                }

                if (minSz < result->data.size())
                {
                    if (datSz<data.size())
                    {
                        for (size_t i=minSz; i<data.size(); i+=9)
                        {
                            result->data[i+0] = data[i+0]; result->data[i+1] = data[i+1]; result->data[i+2] = data[i+2];
                            result->data[i+3] = data[i+3]; result->data[i+4] = data[i+4]; result->data[i+5] = data[i+5];
                            result->data[i+6] = data[i+6]; result->data[i+7] = data[i+7]; result->data[i+8] = data[i+8];
                        }
                        minSz = data.size();
                    }
                    else if (datSz<m->data.size())
                    {
                        for (size_t i=datSz; i<m->data.size(); i+=9)
                        {
                            result->data[i+0] = m->data[i+0]; result->data[i+1] = m->data[i+1]; result->data[i+2] = m->data[i+2];
                            result->data[i+3] = m->data[i+3]; result->data[i+4] = m->data[i+4]; result->data[i+5] = m->data[i+5];
                            result->data[i+6] = m->data[i+6]; result->data[i+7] = m->data[i+7]; result->data[i+8] = m->data[i+8];
                        }
                        minSz = m->data.size();
                    }
                }

                if (minSz < result->data.size())
                {
                    for (size_t i=datSz; i<result->data.size(); i+=9)
                    {
                        result->data[i+0] = 1; result->data[i+1] = 0; result->data[i+2] = 0;
                        result->data[i+3] = 0; result->data[i+4] = 1; result->data[i+5] = 0;
                        result->data[i+6] = 0; result->data[i+7] = 0; result->data[i+8] = 1;
                    }
                }

                return;
            }
        }
        defaulttype::BaseMatrix::opMulTM(bresult,bm);
    }

    void rotateMatrix(defaulttype::BaseMatrix * mat,const defaulttype::BaseMatrix * Jmat)
    {
        if (mat!=Jmat) {
            mat->clear();
            mat->resize(Jmat->rowSize(),Jmat->colSize());
        }

        if (const SparseMatrix<float> * J = dynamic_cast<const SparseMatrix<float> * >(Jmat)) {
            for (typename sofa::component::linearsolver::SparseMatrix<float>::LineConstIterator jit1 = J->begin(), jend = J->end() ; jit1 != jend; jit1++)
            {
                int l = jit1->first;
                for (typename sofa::component::linearsolver::SparseMatrix<float>::LElementConstIterator i1 = jit1->second.begin(), jitend = jit1->second.end(); i1 != jitend;)
                {
                    int c = i1->first;
                    Real v0 = (Real)i1->second; i1++; if (i1==jitend) break;
                    Real v1 = (Real)i1->second; i1++; if (i1==jitend) break;
                    Real v2 = (Real)i1->second; i1++;
                    mat->set(l,c+0,v0 * data[(c+0)*3+0] + v1 * data[(c+1)*3+0] + v2 * data[(c+2)*3+0] );
                    mat->set(l,c+1,v0 * data[(c+0)*3+1] + v1 * data[(c+1)*3+1] + v2 * data[(c+2)*3+1] );
                    mat->set(l,c+2,v0 * data[(c+0)*3+2] + v1 * data[(c+1)*3+2] + v2 * data[(c+2)*3+2] );
                }
            }
        } else if (const SparseMatrix<double> * J = dynamic_cast<const SparseMatrix<double> * >(Jmat)) {
            for (typename sofa::component::linearsolver::SparseMatrix<double>::LineConstIterator jit1 = J->begin(), jend = J->end() ; jit1 != jend; jit1++)
            {
                int l = jit1->first;
                for (typename sofa::component::linearsolver::SparseMatrix<double>::LElementConstIterator i1 = jit1->second.begin(), jitend = jit1->second.end(); i1 != jitend;)
                {
                    int c = i1->first;
                    Real v0 = (Real)i1->second; i1++; if (i1==jitend) break;
                    Real v1 = (Real)i1->second; i1++; if (i1==jitend) break;
                    Real v2 = (Real)i1->second; i1++;
                    mat->set(l,c+0,v0 * data[(c+0)*3+0] + v1 * data[(c+1)*3+0] + v2 * data[(c+2)*3+0] );
                    mat->set(l,c+1,v0 * data[(c+0)*3+1] + v1 * data[(c+1)*3+1] + v2 * data[(c+2)*3+1] );
                    mat->set(l,c+2,v0 * data[(c+0)*3+2] + v1 * data[(c+1)*3+2] + v2 * data[(c+2)*3+2] );
                }
            }
        } else {
            std::cerr << "rotateMatrix for this kind of matrix is not implemented" << std::endl;
        }
    }

    friend std::ostream& operator << (std::ostream& out, const RotationMatrix<Real> & v )
    {
        std::cout.precision(4);
        out << "[";
        for (unsigned y=0; y<v.data.size(); y+=9)
        {
            for (int x=0; x<3; ++x)
            {
                out << "\n[" << std::fixed << v.data[y+x*3] << " " << std::fixed << v.data[y+x*3+1] << " " << std::fixed << v.data[y+x*3+2] << "]";
            }
        }
        out << "\n]";
        return out;
    }

    static const char* Name();

protected :
    helper::vector<Real> data;

};

template<>
inline const char* RotationMatrix<float>::Name() {
    return "RotationMatrixf";
}

template<>
inline const char* RotationMatrix<double>::Name() {
    return "RotationMatrixd";
}

} // namespace misc

} // namespace component

} // namespace sofa

#endif
