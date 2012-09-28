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
#ifndef SOFA_HELPER_SVECTOR_H
#define SOFA_HELPER_SVECTOR_H

#include <sofa/helper/vector.h>

namespace sofa
{

namespace helper
{

//======================================================================
/// Same as helper::vector, + delimitors on serialization
//======================================================================
template<
class T
>
class SVector: public helper::vector<T, helper::CPUMemoryManager<T> >
{
public:
    typedef helper::CPUMemoryManager<T>  Alloc;
    /// size_type
    typedef typename helper::vector<T,Alloc>::size_type size_type;
    /// reference to a value (read-write)
    //typedef typename helper::vector<T,Alloc>::reference reference;
    /// const reference to a value (read only)
    //typedef typename helper::vector<T,Alloc>::const_reference const_reference;

    /// Basic onstructor
    SVector() : helper::vector<T,Alloc>() {}
    /// Constructor
    SVector(size_type n, const T& value): helper::vector<T,Alloc>(n,value) {}
    /// Constructor
    SVector(int n, const T& value): helper::vector<T,Alloc>(n,value) {}
    /// Constructor
    SVector(long n, const T& value): helper::vector<T,Alloc>(n,value) {}
    /// Constructor
    explicit SVector(size_type n): helper::vector<T,Alloc>(n) {}
    /// Constructor
    SVector(const helper::vector<T, Alloc>& x): helper::vector<T,Alloc>(x) {}
    /// Constructor
    vector<T, Alloc>& operator=(const helper::vector<T, Alloc>& x)
    {
        helper::vector<T,Alloc>::operator = (x);
        return (*this);
    }

#ifdef __STL_MEMBER_TEMPLATES
    /// Constructor
    template <class InputIterator>
    SVector(InputIterator first, InputIterator last): helper::vector<T,Alloc>(first,last) {}
#else /* __STL_MEMBER_TEMPLATES */
    /// Constructor
    SVector(typename SVector<T>::const_iterator first, typename SVector<T>::const_iterator last): helper::vector<T,Alloc>(first,last) {}
#endif /* __STL_MEMBER_TEMPLATES */


    std::ostream& write ( std::ostream& os ) const
    {
        if ( !this->empty() )
        {
            typename SVector<T>::const_iterator i = this->begin();
            os << "[ " << *i;
            ++i;
            for ( ; i!=this->end(); ++i )
                os << ", " << *i;
            os << "]";

        }
        return os;
    }

    std::istream& read ( std::istream& in )
    {
        T t;
        this->clear();
        char c;
        in >> c;
        if ( c != '[' )
        {
            std::cerr << "SVector::read : Bad begin character : " << c << ", expected  [" << std::endl;
            return in;
        }
        c = ',';
        while( !in.eof() && c == ',')
        {
            in >> t;
            this->push_back ( t );
            in >> c;
        }
        if ( c != ']' )
        {
            std::cerr << "SVector::read : Bad end character : " << c << ", expected  ]" << std::endl;
            return in;
        }
        return in;
    }

/// Output stream
    inline friend std::ostream& operator<< ( std::ostream& os, const SVector<T>& vec )
    {
        return vec.write(os);
    }

/// Input stream
    inline friend std::istream& operator>> ( std::istream& in, SVector<T>& vec )
    {
        return vec.read(in);
    }

};

} // namespace helper

} // namespace sofa

#endif


