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
#ifndef SOFA_CORE_BEHAVIOR_MULTIVEC_H
#define SOFA_CORE_BEHAVIOR_MULTIVEC_H

#include <sofa/core/MultiVecId.h>
#include <sofa/core/behavior/BaseVectorOperations.h>
namespace sofa
{

namespace core
{

namespace behavior
{

/// Helper class providing a high-level view of underlying state vectors.
///
/// It is used to convert math-like operations to call to computation methods.
template<VecType vtype>
class TMultiVec
{
public:

    typedef TMultiVecId<vtype, V_WRITE> MyMultiVecId;
    typedef TMultiVecId<vtype, V_READ> ConstMyMultiVecId;
    typedef TMultiVecId<V_ALL, V_WRITE> AllMultiVecId;
    typedef TMultiVecId<V_ALL, V_READ> ConstAllMultiVecId;

protected:
    /// Solver who is using this vector
    BaseVectorOperations* vop;

    /// Identifier of this vector
    MyMultiVecId v;

    /// Flag indicating if this vector was dynamically allocated
    bool dynamic;

private:
    /// Copy-constructor is forbidden
    TMultiVec(const TMultiVec<vtype>& ) {}

public:
    /// Default constructor, to combine with the set method
    TMultiVec() : vop(NULL), v(MyMultiVecId::null()), dynamic(false)
    {}

    /// Refers to a state vector with the given ID (VecId::position(), VecId::velocity(), etc).
    TMultiVec( BaseVectorOperations* vop, MyMultiVecId v) : vop(vop), v(v), dynamic(false)
    {}

    /// Allocate a new temporary vector with the given type (sofa::core::V_COORD or sofa::core::V_DERIV).
    TMultiVec( BaseVectorOperations* vop) : vop(vop), v(MyMultiVecId::null()), dynamic(true)
    {
        BOOST_STATIC_ASSERT(vtype == V_COORD || vtype == V_DERIV);
        vop->v_alloc( v );
    }

    void set( BaseVectorOperations* vop )
    {
        // not yet allocated
        if( v.isNull() ) vop->v_alloc( v );
    }

    ~TMultiVec()
    {
        if (dynamic) vop->v_free(v);
    }

    /// Automatic conversion to the underlying VecId
    operator MyMultiVecId() {  return v;  }
    operator ConstMyMultiVecId() {  return v;  }
    operator AllMultiVecId() {  return v;  }
    operator ConstAllMultiVecId() {  return v;  }

    const MyMultiVecId& id() const { return v; }
    MyMultiVecId& id() { return v; }

    BaseVectorOperations* ops() { return vop; }
    void setOps(BaseVectorOperations* op) { vop = op; }

    /// v = 0
    void clear()
    {
        vop->v_clear(v);
    }

    /// v = a
    void eq(MyMultiVecId a)
    {
        vop->v_eq(v, a);
    }

    /// v = a*f
    void eq(MyMultiVecId a, double f)
    {
        vop->v_eq(v, a, f);
    }

    /// v += a*f
    void peq(AllMultiVecId a, double f=1.0)
    {
        vop->v_peq(v, a, f);
    }

    /// v *= f
    void teq(double f)
    {
        vop->v_teq(v, f);
    }

    /// v = a+b*f
    void eq(AllMultiVecId a, AllMultiVecId b, double f=1.0)
    {
        vop->v_op(v, a, b, f);
    }

    /// \return v.a
    double dot(MyMultiVecId a)
    {
        vop->v_dot(v, a);
        return vop->finish();
    }

    /// nullify values below given threshold
    void threshold( double threshold )
    {
        vop->v_threshold(v, threshold);
    }

    /// \return sqrt(v.v)
    double norm()
    {
        vop->v_dot(v, v);
        return sqrt( vop->finish() );
    }

    /// v = a
    void operator=(MyMultiVecId a)
    {
        eq(a);
    }

    /// v = a
    void operator=(const TMultiVec<vtype>& a)
    {
        eq(a.v);
    }

    /// v += a
    void operator+=(MyMultiVecId a)
    {
        peq(a);
    }

    /// v -= a
    void operator-=(MyMultiVecId a)
    {
        peq(a,-1);
    }

    /// v *= f
    void operator*=(double f)
    {
        teq(f);
    }

    /// v /= f
    void operator/=(double f)
    {
        teq(1.0/f);
    }

    /// return the scalar product dot(v,a)
    double operator*(MyMultiVecId a)
    {
        return dot(a);
    }

    friend std::ostream& operator << (std::ostream& out, const TMultiVec<vtype>& mv )
    {
        mv.vop->print(mv.v,out);
        return out;
    }
};

typedef TMultiVec<V_COORD> MultiVecCoord;
typedef TMultiVec<V_DERIV> MultiVecDeriv;
typedef TMultiVec<V_MATDERIV> MultiVecMatrixDeriv;

} // namespace behavior

} // namespace core

} // namespace sofa

#endif // SOFA_CORE_BEHAVIOR_MULTIVEC_H

