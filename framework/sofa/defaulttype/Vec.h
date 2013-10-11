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
#ifndef SOFA_DEFAULTTYPE_VEC_H
#define SOFA_DEFAULTTYPE_VEC_H

#include <sofa/helper/fixed_array.h>
#include <boost/static_assert.hpp>
#include <sofa/helper/rmath.h>
#include <sofa/defaulttype/DataTypeInfo.h>
#include <functional>
#include <limits>

#define EQUALITY_THRESHOLD 1e-6

namespace sofa
{

namespace defaulttype
{

enum NoInit { NOINIT }; ///< use when calling Vec or Mat constructor to skip initialization of values to 0

template <int N, typename real=float>
class Vec : public helper::fixed_array<real,N>
{
public:
    /// Compile-time constant specifying the number of scalars within this vector (equivalent to static_size and size() method)
    enum { total_size = N };
    /// Compile-time constant specifying the number of dimensions of space (equivalent to total_size here)
    enum { spatial_dimensions = N };

    /// Default constructor: sets all values to 0.
    Vec()
    {
        this->clear();
    }

    /// Fast constructor: no initialization
    explicit Vec(NoInit)
    {
    }

    /// Specific constructor for 1-element vectors.
    explicit Vec(real r1)
    {
        set( r1 );
    }

    /// Specific constructor for 1-element vectors.
    void operator=(real r1)
    {
        set( r1 );
    }

    /// Specific constructor for 2-elements vectors.
    Vec(real r1, real r2)
    {
        set( r1, r2 );
    }

    /// Specific constructor for 3-elements vectors.
    Vec(real r1, real r2, real r3)
    {
        set( r1, r2, r3 );
    }

    /// Specific constructor for 4-elements vectors.
    Vec(real r1, real r2, real r3, real r4)
    {
        set( r1, r2, r3, r4 );
    }

    /// Specific constructor for 5-elements vectors.
    Vec(real r1, real r2, real r3, real r4, real r5)
    {
        set( r1, r2, r3, r4, r5 );
    }

    /// Specific constructor for 6-elements vectors.
    Vec(real r1, real r2, real r3, real r4, real r5, real r6)
    {
        set( r1, r2, r3, r4, r5, r6 );
    }

    /// Specific constructor for 6-elements vectors.
    template<typename R, typename T>
    Vec( const Vec<3,R>& a , const Vec<3,T>& b )
    {
        set( a[0], a[1], a[2], b[0], b[1], b[2] );
    }

    /// Specific constructor for 8-elements vectors.
    Vec(real r1, real r2, real r3, real r4, real r5, real r6, real r7, real r8)
    {
        set( r1, r2, r3, r4, r5, r6, r7, r8 );
    }

    /// Specific constructor for 9-elements vectors.
    Vec(real r1, real r2, real r3, real r4, real r5, real r6, real r7, real r8, real r9)
    {
        set( r1, r2, r3, r4, r5, r6, r7, r8, r9 );
    }

    /// Specific constructor for 12-elements vectors.
    Vec(real r1, real r2, real r3, real r4, real r5, real r6, real r7, real r8, real r9, real r10, real r11, real r12)
    {
        set( r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12 );
    }

    /// Specific set for 1-element vectors.
    void set(real r1)
    {
        BOOST_STATIC_ASSERT(N==1);
        this->elems[0]=r1;
    }

    /// Specific set for 2-elements vectors.
    void set(real r1, real r2)
    {
        BOOST_STATIC_ASSERT(N == 2);
        this->elems[0]=r1;
        this->elems[1]=r2;
    }

    /// Specific set for 3-elements vectors.
    void set(real r1, real r2, real r3)
    {
        BOOST_STATIC_ASSERT(N == 3);
        this->elems[0]=r1;
        this->elems[1]=r2;
        this->elems[2]=r3;
    }

    /// Specific set for 4-elements vectors.
    void set(real r1, real r2, real r3, real r4)
    {
        BOOST_STATIC_ASSERT(N == 4);
        this->elems[0]=r1;
        this->elems[1]=r2;
        this->elems[2]=r3;
        this->elems[3]=r4;
    }

    /// Specific set for 5-elements vectors.
    void set(real r1, real r2, real r3, real r4, real r5)
    {
        BOOST_STATIC_ASSERT(N == 5);
        this->elems[0]=r1;
        this->elems[1]=r2;
        this->elems[2]=r3;
        this->elems[3]=r4;
        this->elems[4]=r5;
    }

    /// Specific set for 6-elements vectors.
    void set(real r1, real r2, real r3, real r4, real r5, real r6)
    {
        BOOST_STATIC_ASSERT(N == 6);
        this->elems[0]=r1;
        this->elems[1]=r2;
        this->elems[2]=r3;
        this->elems[3]=r4;
        this->elems[4]=r5;
        this->elems[5]=r6;
    }

    /// Specific set for 8-elements vectors.
    void set(real r1, real r2, real r3, real r4, real r5, real r6, real r7, real r8)
    {
        BOOST_STATIC_ASSERT(N == 8);
        this->elems[0]=r1;
        this->elems[1]=r2;
        this->elems[2]=r3;
        this->elems[3]=r4;
        this->elems[4]=r5;
        this->elems[5]=r6;
        this->elems[6]=r7;
        this->elems[7]=r8;
    }

    /// Specific set for 9-elements vectors.
    void set(real r1, real r2, real r3, real r4, real r5, real r6, real r7, real r8, real r9)
    {
        BOOST_STATIC_ASSERT(N == 9);
        this->elems[0]=r1;
        this->elems[1]=r2;
        this->elems[2]=r3;
        this->elems[3]=r4;
        this->elems[4]=r5;
        this->elems[5]=r6;
        this->elems[6]=r7;
        this->elems[7]=r8;
        this->elems[8]=r9;
    }

    /// Specific set for 12-elements vectors.
    void set(real r1, real r2, real r3, real r4, real r5, real r6, real r7, real r8, real r9, real r10, real r11, real r12)
    {
        BOOST_STATIC_ASSERT(N == 12);
        this->elems[0]=r1;
        this->elems[1]=r2;
        this->elems[2]=r3;
        this->elems[3]=r4;
        this->elems[4]=r5;
        this->elems[5]=r6;
        this->elems[6]=r7;
        this->elems[7]=r8;
        this->elems[8]=r9;
        this->elems[9]=r10;
        this->elems[10]=r11;
        this->elems[11]=r12;
    }

    /// Specific set from a different size vector (given default value and ignored outside entries)
    template<int N2, class real2>
    void set(const Vec<N2,real2>& v, real defaultvalue=0)
    {
        int maxN = std::min( N, N2 );
        for(int i=0; i<maxN; i++)
            this->elems[i] = v[i];
        for(int i=maxN; i<N ; i++)
            this->elems[i] = defaultvalue;
    }


    /// Constructor from an N-1 elements vector and an additional value (added at the end).
    Vec(const Vec<N-1,real>& v, real r1)
    {
        BOOST_STATIC_ASSERT(N > 1);
        set( v, r1 );
    }

    template<int N2, typename real2>
    Vec(const Vec<N2,real2>& v)
    {
        BOOST_STATIC_ASSERT(N > 1);
        set( v, 0 );
    }

    template<typename real2>
    Vec(const Vec<N, real2>& p)
    {
        //std::copy(p.begin(), p.end(), this->begin());
        for(int i=0; i<N; i++)
            this->elems[i] = (real)p(i);
    }

    /// Constructor from an array of values.
    template<typename real2>
    explicit Vec(const real2* p)
    {
        //std::copy(p, p+N, this->begin());
        for(int i=0; i<N; i++)
            this->elems[i] = (real)p[i];
    }

    /// Special access to first element.
    real& x()
    {
        BOOST_STATIC_ASSERT(N >= 1);
        return this->elems[0];
    }
    /// Special access to second element.
    real& y()
    {
        BOOST_STATIC_ASSERT(N >= 2);
        return this->elems[1];
    }
    /// Special access to third element.
    real& z()
    {
        BOOST_STATIC_ASSERT(N >= 3);
        return this->elems[2];
    }
    /// Special access to fourth element.
    real& w()
    {
        BOOST_STATIC_ASSERT(N >= 4);
        return this->elems[3];
    }

    /// Special const access to first element.
    const real& x() const
    {
        BOOST_STATIC_ASSERT(N >= 1);
        return this->elems[0];
    }
    /// Special const access to second element.
    const real& y() const
    {
        BOOST_STATIC_ASSERT(N >= 2);
        return this->elems[1];
    }
    /// Special const access to third element.
    const real& z() const
    {
        BOOST_STATIC_ASSERT(N >= 3);
        return this->elems[2];
    }
    /// Special const access to fourth element.
    const real& w() const
    {
        BOOST_STATIC_ASSERT(N >= 4);
        return this->elems[3];
    }

    /// Assignment operator from an array of values.
    template<typename real2>
    void operator=(const real2* p)
    {
        //std::copy(p, p+N, this->begin());
        for(int i=0; i<N; i++)
            this->elems[i] = (real)p[i];
    }

    /// Assignment from a vector with different dimensions.
    template<int M, typename real2>
    void operator=(const Vec<M,real2>& v)
    {
        //std::copy(v.begin(), v.begin()+(N>M?M:N), this->begin());
        for(int i=0; i<(N>M?M:N); i++)
            this->elems[i] = (real)v(i);
    }

    /// Sets every element to 0.
    void clear()
    {
        this->assign(real());
    }

    /// Sets every element to r.
    void fill(real r)
    {
        this->assign(r);
    }

    // Access to i-th element.
    // Already in fixed_array
    //real& operator[](int i)
    //{
    //    return this->elems[i];
    //}

    // Access to i-th element.
    // Already in fixed_array
    /// Const access to i-th element.
    //const real& operator[](int i) const
    //{
    //    return this->elems[i];
    //}

    /// Access to i-th element.
    real& operator()(int i)
    {
        return this->elems[i];
    }

    /// Const access to i-th element.
    const real& operator()(int i) const
    {
        return this->elems[i];
    }

    /// Cast into a const array of values.
    /// CHANGE(Jeremie A.): removed it as it confuses some compilers. Use ptr() or data() instead
    //operator const real*() const
    //{
    //    return this->elems;
    //}

    /// Cast into an array of values.
    /// CHANGE(Jeremie A.): removed it as it confuses some compilers. Use ptr() or data() instead
    //operator real*()
    //{
    //    return this->elems;
    //}

    /// Cast into a const array of values.
    const real* ptr() const
    {
        return this->elems;
    }

    /// Cast into an array of values.
    real* ptr()
    {
        return this->elems;
    }

    // LINEAR ALGEBRA

    // BUG (J.A. 12/31/2010): gcc 4.0 does not support templated
    // operators that are restricted to scalar types using static_assert.
    // So for now we are defining them as templated method, and the
    // operators then simply call them with the right types.

    /// Multiplication by a scalar f.
    template<class real2>
    Vec<N,real> mulscalar(real2 f) const
    {
        BOOST_STATIC_ASSERT(DataTypeInfo<real2>::ValidInfo && DataTypeInfo<real2>::Size==1);
        Vec<N,real> r(NOINIT);
        for (int i=0; i<N; i++)
            r[i] = this->elems[i]*(real)f;
        return r;
    }

    Vec<N,real> operator*(         float     f) const {  return mulscalar(f);  }
    Vec<N,real> operator*(         double    f) const {  return mulscalar(f);  }
    Vec<N,real> operator*(         int       f) const {  return mulscalar(f);  }
    Vec<N,real> operator*(unsigned int       f) const {  return mulscalar(f);  }
    Vec<N,real> operator*(         long      f) const {  return mulscalar(f);  }
    Vec<N,real> operator*(unsigned long      f) const {  return mulscalar(f);  }
    Vec<N,real> operator*(         long long f) const {  return mulscalar(f);  }
    Vec<N,real> operator*(unsigned long long f) const {  return mulscalar(f);  }

    /// In-place multiplication by a scalar f.
    template<class real2>
    void eqmulscalar(real2 f)
    {
        BOOST_STATIC_ASSERT(DataTypeInfo<real2>::ValidInfo && DataTypeInfo<real2>::Size==1);
        for (int i=0; i<N; i++)
            this->elems[i]*=(real)f;
    }

    void operator*=(         float     f) {  eqmulscalar(f);  }
    void operator*=(         double    f) {  eqmulscalar(f);  }
    void operator*=(         int       f) {  eqmulscalar(f);  }
    void operator*=(unsigned int       f) {  eqmulscalar(f);  }
    void operator*=(         long      f) {  eqmulscalar(f);  }
    void operator*=(unsigned long      f) {  eqmulscalar(f);  }
    void operator*=(         long long f) {  eqmulscalar(f);  }
    void operator*=(unsigned long long f) {  eqmulscalar(f);  }

    /// Division by a scalar f.
    template<class real2>
    Vec<N,real> divscalar(real2 f) const
    {
        BOOST_STATIC_ASSERT(DataTypeInfo<real2>::ValidInfo && DataTypeInfo<real2>::Size==1);
        Vec<N,real> r(NOINIT);
        for (int i=0; i<N; i++)
            r[i] = this->elems[i]/(real)f;
        return r;
    }

    Vec<N,real> operator/(         float     f) const {  return divscalar(f);  }
    Vec<N,real> operator/(         double    f) const {  return divscalar(f);  }
    Vec<N,real> operator/(         int       f) const {  return divscalar(f);  }
    Vec<N,real> operator/(unsigned int       f) const {  return divscalar(f);  }
    Vec<N,real> operator/(         long      f) const {  return divscalar(f);  }
    Vec<N,real> operator/(unsigned long      f) const {  return divscalar(f);  }
    Vec<N,real> operator/(         long long f) const {  return divscalar(f);  }
    Vec<N,real> operator/(unsigned long long f) const {  return divscalar(f);  }

    /// In-place division by a scalar f.
    template<class real2>
    void eqdivscalar(real2 f)
    {
        BOOST_STATIC_ASSERT(DataTypeInfo<real2>::ValidInfo && DataTypeInfo<real2>::Size==1);
        for (int i=0; i<N; i++)
            this->elems[i]/=(real)f;
    }

    void operator/=(         float     f) {  eqdivscalar(f);  }
    void operator/=(         double    f) {  eqdivscalar(f);  }
    void operator/=(         int       f) {  eqdivscalar(f);  }
    void operator/=(unsigned int       f) {  eqdivscalar(f);  }
    void operator/=(         long      f) {  eqdivscalar(f);  }
    void operator/=(unsigned long      f) {  eqdivscalar(f);  }
    void operator/=(         long long f) {  eqdivscalar(f);  }
    void operator/=(unsigned long long f) {  eqdivscalar(f);  }

    /// Dot product.
    template<class real2>
    real operator*(const Vec<N,real2>& v) const
    {
        real r = (real)(this->elems[0]*v[0]);
        for (int i=1; i<N; i++)
            r += (real)(this->elems[i]*v[i]);
        return r;
    }

    /// linear product.
    template<class real2>
    Vec<N,real> linearProduct(const Vec<N,real2>& v) const
    {
        Vec<N,real> r(NOINIT);
        for (int i=0; i<N; i++)
            r[i]=this->elems[i]*(real)v[i];
        return r;
    }


    /// linear division.
    template<class real2>
    Vec<N,real> linearDivision(const Vec<N,real2>& v) const
    {
        Vec<N,real> r(NOINIT);
        for (int i=0; i<N; i++)
            r[i]=this->elems[i]/(real)v[i];
        return r;
    }

    /// Vector addition.
    template<class real2>
    Vec<N,real> operator+(const Vec<N,real2>& v) const
    {
        Vec<N,real> r(NOINIT);
        for (int i=0; i<N; i++)
            r[i]=this->elems[i]+(real)v[i];
        return r;
    }

    /// In-place vector addition.
    template<class real2>
    void operator+=(const Vec<N,real2>& v)
    {
        for (int i=0; i<N; i++)
            this->elems[i]+=(real)v[i];
    }

    /// Vector subtraction.
    template<class real2>
    Vec<N,real> operator-(const Vec<N,real2>& v) const
    {
        Vec<N,real> r(NOINIT);
        for (int i=0; i<N; i++)
            r[i]=this->elems[i]-(real)v[i];
        return r;
    }

    /// In-place vector subtraction.
    template<class real2>
    void operator-=(const Vec<N,real2>& v)
    {
        for (int i=0; i<N; i++)
            this->elems[i]-=(real)v[i];
    }

    /// Vector negation.
    Vec<N,real> operator-() const
    {
        Vec<N,real> r(NOINIT);
        for (int i=0; i<N; i++)
            r[i]=-this->elems[i];
        return r;
    }

    /// Squared norm.
    real norm2() const
    {
        real r = this->elems[0]*this->elems[0];
        for (int i=1; i<N; i++)
            r += this->elems[i]*this->elems[i];
        return r;
    }

    /// Euclidean norm.
    real norm() const
    {
        return helper::rsqrt(norm2());
    }

    /// Normalize the vector with the given norm
    void normalizeWithNorm(real norm, real threshold=std::numeric_limits<real>::epsilon())
    {
        if (norm>threshold)
            for (int i=0; i<N; i++)
                this->elems[i]/=norm;
    }

    /// Normalize the vector.
    void normalize(real threshold=std::numeric_limits<real>::epsilon())
    {
        normalizeWithNorm(norm(),threshold);
    }

    /// return true iff norm()==1
    bool isNormalized( real threshold=std::numeric_limits<real>::epsilon()*(real)10 ) const { return helper::rabs<real>( norm2()-(real)1 ) <= threshold; }

    template<typename R>
    Vec cross( const Vec<3,R>& b ) const
    {
        BOOST_STATIC_ASSERT(N == 3);
        return Vec(
                (real)((*this)[1]*b[2] - (*this)[2]*b[1]),
                (real)((*this)[2]*b[0] - (*this)[0]*b[2]),
                (real)((*this)[0]*b[1] - (*this)[1]*b[0])
                );
    }


    /// sum of all elements of the vector
    real sum() const
    {
        real sum = 0.0;
        for (int i=0; i<N; i++)
            sum += this->elems[i];
        return sum;
    }


    /// @name Tests operators
    /// @{

    bool operator==(const Vec& b) const
    {
        for (int i=0; i<N; i++)
            if ( fabs( (float)(this->elems[i] - b[i]) ) > EQUALITY_THRESHOLD ) return false;
        return true;
    }

    bool operator!=(const Vec& b) const
    {
        for (int i=0; i<N; i++)
            if ( fabs( (float)(this->elems[i] - b[i]) ) > EQUALITY_THRESHOLD ) return true;
        return false;
    }

    /// @}
};

/// Same as Vec except the values are not initialized by default
template <int N, typename real=float>
class VecNoInit : public Vec<N,real>
{
public:
    VecNoInit()
        : Vec<N,real>(NOINIT)
    {
    }

    /// Assignment operator from an array of values.
    template<typename real2>
    void operator=(const real2* p)
    {
        this->Vec<N,real>::operator=(p);
    }

    /// Assignment from a vector with different dimensions.
    template<int M, typename real2>
    void operator=(const Vec<M,real2>& v)
    {
        this->Vec<N,real>::operator=(v);
    }

    /// Scalar vector multiplication operator.
    friend Vec<N,real> operator*(real r, const Vec<N,real>& v)
    {
        return v*r;
    }
};

/// Read from an input stream
template<int N,typename Real>
std::istream& operator >> ( std::istream& in, Vec<N,Real>& v )
{
    for( int i=0; i<N; ++i )
        in>>v[i];
    return in;
}

/// Write to an output stream
template<int N,typename Real>
std::ostream& operator << ( std::ostream& out, const Vec<N,Real>& v )
{
    for( int i=0; i<N-1; ++i )
        out<<v[i]<<" ";
    out<<v[N-1];
    return out;
}

/// Cross product for 3-elements vectors.
template<typename real1, typename real2 >
inline Vec<3,real1> cross(const Vec<3,real1>& a, const Vec<3,real2>& b)
{
    return Vec<3,real1>((real1)(a.y()*b.z() - a.z()*b.y()),
            (real1)(a.z()*b.x() - a.x()*b.z()),
            (real1)(a.x()*b.y() - a.y()*b.x()));
}

/// Cross product for 2-elements vectors.
template <typename real1, typename real2>
real1 cross(const defaulttype::Vec<2,real1>& a, const defaulttype::Vec<2,real2>& b )
{
    return (real1)(a[0]*b[1] - a[1]*b[0]);
}

/// Dot product (alias for operator*)
template<int N,typename real>
inline real dot(const Vec<N,real>& a, const Vec<N,real>& b)
{
    return a*b;
}

/// multiplication with a scalar \returns a*V
template <int N, typename real>
Vec<N,real> operator*(const real& a, const Vec<N,real>& V)
{
    return V * a;
}

typedef Vec<1,float> Vec1f;
typedef Vec<1,double> Vec1d;
typedef Vec<1,int> Vec1i;
typedef Vec<1,unsigned> Vec1u;


typedef Vec<2,float> Vec2f;
typedef Vec<2,double> Vec2d;
typedef Vec<2,int> Vec2i;
typedef Vec<2,unsigned> Vec2u;


typedef Vec<3,float> Vec3f;
typedef Vec<3,double> Vec3d;
typedef Vec<3,int> Vec3i;
typedef Vec<3,unsigned> Vec3u;


typedef Vec<4,float> Vec4f;
typedef Vec<4,double> Vec4d;
typedef Vec<4,int> Vec4i;
typedef Vec<4,unsigned> Vec4u;


typedef Vec<6,float> Vec6f;
typedef Vec<6,double> Vec6d;
typedef Vec<6,int> Vec6i;
typedef Vec<6,unsigned> Vec6u;


#ifdef SOFA_FLOAT
typedef Vec1f Vector1; ///< alias
typedef Vec2f Vector2; ///< alias
typedef Vec3f Vector3; ///< alias
typedef Vec4f Vector4; ///< alias
typedef Vec6f Vector6; ///< alias
#else
typedef Vec1d Vector1; ///< alias
typedef Vec2d Vector2; ///< alias
typedef Vec3d Vector3; ///< alias
typedef Vec4d Vector4; ///< alias
typedef Vec6d Vector6; ///< alias
#endif

} // namespace defaulttype

} // namespace sofa

// Specialization of the defaulttype::DataTypeInfo type traits template

namespace sofa
{

namespace defaulttype
{

template<int N, typename real>
struct DataTypeInfo< sofa::defaulttype::Vec<N,real> > : public FixedArrayTypeInfo<sofa::defaulttype::Vec<N,real> >
{
    static std::string name() { std::ostringstream o; o << "Vec<" << N << "," << DataTypeName<real>::name() << ">"; return o.str(); }
};

template<int N, typename real>
struct DataTypeInfo< sofa::defaulttype::VecNoInit<N,real> > : public FixedArrayTypeInfo<sofa::defaulttype::VecNoInit<N,real> >
{
    static std::string name() { std::ostringstream o; o << "VecNoInit<" << N << "," << DataTypeName<real>::name() << ">"; return o.str(); }
};

// The next line hides all those methods from the doxygen documentation
/// \cond TEMPLATE_OVERRIDES

template<> struct DataTypeName<defaulttype::Vec1f> { static const char* name() { return "Vec1f"; } };
template<> struct DataTypeName<defaulttype::Vec1d> { static const char* name() { return "Vec1d"; } };
template<> struct DataTypeName<defaulttype::Vec2f> { static const char* name() { return "Vec2f"; } };
template<> struct DataTypeName<defaulttype::Vec2d> { static const char* name() { return "Vec2d"; } };
template<> struct DataTypeName<defaulttype::Vec3f> { static const char* name() { return "Vec3f"; } };
template<> struct DataTypeName<defaulttype::Vec3d> { static const char* name() { return "Vec3d"; } };
template<> struct DataTypeName<defaulttype::Vec4f> { static const char* name() { return "Vec4f"; } };
template<> struct DataTypeName<defaulttype::Vec4d> { static const char* name() { return "Vec4d"; } };
template<> struct DataTypeName<defaulttype::Vec6f> { static const char* name() { return "Vec6f"; } };
template<> struct DataTypeName<defaulttype::Vec6d> { static const char* name() { return "Vec6d"; } };

/// \endcond

} // namespace defaulttype

} // namespace sofa

// Specialization of the std comparison function, to use Vec as std::map key
namespace std
{

// template <>
template<int N, class T>
struct less< sofa::defaulttype::Vec<N,T> > : public binary_function< sofa::defaulttype::Vec<N,T>,  sofa::defaulttype::Vec<N,T>, bool>
{
    bool operator()(const  sofa::defaulttype::Vec<N,T>& x, const  sofa::defaulttype::Vec<N,T>& y) const
    {
        //std::cerr<<"specialized std::less, x = "<<x<<", y = "<<y<<std::endl;
        for( unsigned i=0; i<N; ++i )
        {
            if( x[i]<y[i] )
                return true;
            else if( y[i]<x[i] )
                return false;
        }
        return false;
    }
};

} // namespace std

#endif

