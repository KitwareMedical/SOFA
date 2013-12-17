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
#ifndef SOFA_CORE_OBJECTMODEL_DATAFIELD_H
#define SOFA_CORE_OBJECTMODEL_DATAFIELD_H

#if !defined(__GNUC__) || (__GNUC__ > 3 || (_GNUC__ == 3 && __GNUC_MINOR__ > 3))
#pragma once
#endif

#include <sofa/core/core.h>
#include <sofa/core/objectmodel/BaseData.h>
#include <sofa/helper/accessor.h>
#include <boost/shared_ptr.hpp>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>

namespace sofa
{

namespace core
{

namespace objectmodel
{

/**
 *  \brief Abstract templated data, readable and writable from/to a string.
 *
 */
template < class T >
class TData : public BaseData
{
public:
    typedef T value_type;

    /// @name Class reflection system
    /// @{
    typedef TClass<TData<T>,BaseData> MyClass;
    static const MyClass* GetClass() { return MyClass::get(); }
    virtual const BaseClass* getClass() const
    { return GetClass(); }

    static std::string templateName(const TData<T>* = NULL)
    {
        T* ptr = NULL;
        return BaseData::typeName(ptr);
    }
    /// @}

    explicit TData(const BaseInitData& init)
        : BaseData(init), parentData(initLink("parentSameType", "Linked Data in case it stores exactly the same type of Data, and efficient copies can be made (by value or by sharing pointers with Copy-on-Write)"))
    {
    }

    TData( const char* helpMsg=0, bool isDisplayed=true, bool isReadOnly=false)
        : BaseData(helpMsg, isDisplayed, isReadOnly), parentData(initLink("parentSameType", "Linked Data in case it stores exactly the same type of Data, and efficient copies can be made (by value or by sharing pointers with Copy-on-Write)"))
    {
    }

    virtual ~TData()
    {}

    inline void printValue(std::ostream& out) const;
    inline std::string getValueString() const;
    inline std::string getValueTypeString() const; // { return std::string(typeid(m_value).name()); }

    /// Get info about the value type of the associated variable
    virtual const sofa::defaulttype::AbstractTypeInfo* getValueTypeInfo() const
    {
        return sofa::defaulttype::VirtualTypeInfo<T>::get();
    }

    virtual const T& virtualGetValue() const = 0;
    virtual void virtualSetValue(const T& v) = 0;
    virtual void virtualSetLink(const BaseData& bd) = 0;
    virtual T* virtualBeginEdit() = 0;
    virtual void virtualEndEdit() = 0;

    /// Get current value as a void pointer (use getValueTypeInfo to find how to access it)
    virtual const void* getValueVoidPtr() const
    {
        return &(virtualGetValue());
    }

    /// Begin edit current value as a void pointer (use getValueTypeInfo to find how to access it)
    virtual void* beginEditVoidPtr()
    {
        return virtualBeginEdit();
    }

    /// End edit current value as a void pointer (use getValueTypeInfo to find how to access it)
    virtual void endEditVoidPtr()
    {
        virtualEndEdit();
    }

    /** Try to read argument value from an input stream.
    Return false if failed
     */
    virtual bool read( const std::string& s )
    {
        if (s.empty())
            return false;
        //serr<<"Field::read "<<s.c_str()<<sendl;
        std::istringstream istr( s.c_str() );
        istr >> *virtualBeginEdit();
        virtualEndEdit();
        if( istr.fail() )
        {
            return false;
        }
        else
        {
            return true;
        }
    }

    virtual bool isCounterValid() const {return true;}

    bool copyValue(const TData<T>* parent)
    {
        virtualSetValue(parent->virtualGetValue());
        return true;
    }

    virtual bool copyValue(const BaseData* parent)
    {
        const TData<T>* p = dynamic_cast<const TData<T>*>(parent);
        if (p)
        {
            virtualSetValue(p->virtualGetValue());
            return true;
        }
        return BaseData::copyValue(parent);
    }


    virtual bool validParent(BaseData* parent)
    {
        if (dynamic_cast<TData<T>*>(parent))
            return true;
        return BaseData::validParent(parent);
    }

protected:

    BaseLink::InitLink<TData<T> >
    initLink(const char* name, const char* help)
    {
        return BaseLink::InitLink<TData<T> >(this, name, help);
    }

    void doSetParent(BaseData* parent)
    {
        parentData.set(dynamic_cast<TData<T>*>(parent));
        BaseData::doSetParent(parent);
    }

    bool updateFromParentValue(const BaseData* parent)
    {
        if (parent == parentData.get())
        {
            //virtualSetValue(parentData->virtualGetValue());
            virtualSetLink(*parentData.get());
            return true;
        }
        else
            return BaseData::updateFromParentValue(parent);
    }

    SingleLink<TData<T>,TData<T>, BaseLink::FLAG_DATALINK|BaseLink::FLAG_DUPLICATE> parentData;
};

template <class T, bool COW>
class DataValue;

template <class T>
class DataValue<T, false>
{
    T data;
public:

    DataValue()
        : data(T())// BUGFIX (Jeremie A.): Force initialization of basic types to 0 (bool, int, float, etc).
    {
    }

    explicit DataValue(const T &value)
        : data(value)
    {
    }

    DataValue(const DataValue& dc)
        : data(dc.getValue())
    {
    }

    DataValue& operator=(const DataValue& dc )
    {
        data = dc.getValue();
        return *this;
    }

    T* beginEdit() { return &data; }
    void endEdit() {}
    const T& getValue() const { return data; }
    void setValue(const T& value)
    {
        data = value;
    }
    void release()
    {
    }
//    T& value() { return data; }
};


template <class T>
class DataValue<T, true>
{
    boost::shared_ptr<T> ptr;
public:

    DataValue()
        : ptr(new T(T())) // BUGFIX (Jeremie A.): Force initialization of basic types to 0 (bool, int, float, etc).
    {
    }

    explicit DataValue(const T& value)
        : ptr(new T(value))
    {
    }

    DataValue(const DataValue& dc)
        : ptr(dc.ptr)
    {
    }

    ~DataValue()
    {
    }

    DataValue& operator=(const DataValue& dc )
    {
        //avoid self reference
        if(&dc != this)
        {
            ptr = dc.ptr;
        }

        return *this;
    }

    T* beginEdit()
    {
        if(!ptr.unique())
        {
            ptr.reset(new T(*ptr));
        }
        return ptr.get();
    }

    void endEdit()
    {
    }

    const T& getValue() const
    {
        return *ptr;
    }

    void setValue(const T& value)
    {
        if(!ptr.unique())
        {
            ptr.reset(new T(value));
        }
        else
        {
            *ptr = value;
        }
    }

    void release()
    {
        ptr.reset();
    }

//    T& value()
//    {
//        return *beginEdit();
//    }
};



/**
 *  \brief Container of data, readable and writable from/to a string.
 *
 */
template < class T = void* >
class Data : public TData<T>
{
public:

    /// @name Class reflection system
    /// @{
    typedef TClass<Data<T>, TData<T> > MyClass;
    static const MyClass* GetClass() { return MyClass::get(); }
    virtual const BaseClass* getClass() const
    { return GetClass(); }

    static std::string templateName(const Data<T>* = NULL)
    {
        T* ptr = NULL;
        return BaseData::typeName(ptr);
    }
    /// @}

    /// @name Construction / destruction
    /// @{

    /// This internal class is used by the initData() methods to store initialization parameters of a Data
    class InitData : public BaseData::BaseInitData
    {
    public:
        InitData() : value(T()) {}
        InitData(const T& v) : value(v) {}
        InitData(const BaseData::BaseInitData& i) : BaseData::BaseInitData(i), value(T()) {}

        T value;
    };

    /** Constructor
        this constructor should be used through the initData() methods
     */
    explicit Data(const BaseData::BaseInitData& init)
        : TData<T>(init)
        , shared(NULL)
    {
    }

    /** Constructor
        this constructor should be used through the initData() methods
     */
    explicit Data(const InitData& init)
        : TData<T>(init)
        , m_values()
        , shared(NULL)
    {
        m_values[DDGNode::currentAspect()] = ValueType(init.value);
    }

    /** Constructor
    \param helpMsg help on the field
     */
    Data( const char* helpMsg=0, bool isDisplayed=true, bool isReadOnly=false)
        : TData<T>(helpMsg, isDisplayed, isReadOnly)
        , m_values()
        , shared(NULL)
    {
        ValueType val;
        m_values.assign(val);
    }

    /** Constructor
    \param value default value
    \param helpMsg help on the field
     */
    Data( const T& value, const char* helpMsg=0, bool isDisplayed=true, bool isReadOnly=false)
        : TData<T>(helpMsg, isDisplayed, isReadOnly)
        , m_values()
        , shared(NULL)
    {
        m_values[DDGNode::currentAspect()] = ValueType(value);
    }

    virtual ~Data()
    {}

    /// @}

    /// @name Simple edition and retrieval API
    /// @{

    inline T* beginEdit(const core::ExecParams* params = 0)
    {
        size_t aspect = DDGNode::currentAspect(params);
        this->updateIfDirty(params);
        ++this->m_counters[aspect];
        this->m_isSets[aspect] = true;
        BaseData::setDirtyOutputs(params);
        return m_values[aspect].beginEdit();
    }

    inline void endEdit(const core::ExecParams* params = 0)
    {
        m_values[DDGNode::currentAspect(params)].endEdit();
    }

    inline void setValue(const T& value)
    {
        *beginEdit() = value;
        endEdit();
    }

    inline void setValue(const core::ExecParams* params, const T& value)
    {
        *beginEdit(params) = value;
        endEdit(params);
    }

    inline const T& getValue(const core::ExecParams* params = 0) const
    {
        this->updateIfDirty(params);
        return m_values[DDGNode::currentAspect(params)].getValue();
    }

    void copyAspect(int destAspect, int srcAspect)
    {
        m_values[destAspect] = m_values[srcAspect];
        BaseData::copyAspect(destAspect, srcAspect);
    }

    void releaseAspect(int aspect)
    {
        m_values[aspect].release();
    }
    /// @}

    /// @name Virtual edition and retrieval API (for generic TData parent API, deprecated)
    /// @{

    virtual const T& virtualGetValue() const { return getValue(); }
    virtual void virtualSetValue(const T& v) { setValue(v); }

    virtual void virtualSetLink(const BaseData& bd)
    {
        const Data<T>* d = dynamic_cast< const Data<T>* >(&bd);
        if (d)
        {
            size_t aspect = DDGNode::currentAspect();
            this->m_values[aspect] = d->m_values[aspect];
            //FIX: update counter
            ++this->m_counters[aspect];
            this->m_isSets[aspect] = true;
            BaseData::setDirtyOutputs();
        }
    }

    virtual T* virtualBeginEdit() { return beginEdit(); }
    virtual void virtualEndEdit() { endEdit(); }


    /// @}

    inline friend std::ostream & operator << (std::ostream &out, const Data& df)
    {
        out<<df.getValue();
        return out;
    }

    inline bool operator ==( const T& value ) const
    {
        return getValue()==value;
    }

    inline bool operator !=( const T& value ) const
    {
        return getValue()!=value;
    }

    inline void operator =( const T& value )
    {
        this->setValue(value);
    }

protected:

    typedef DataValue<T, sofa::defaulttype::DataTypeInfo<T>::CopyOnWrite> ValueType;

    /// Value
    helper::fixed_array<ValueType, SOFA_DATA_MAX_ASPECTS> m_values;

public:
    mutable void* shared;

private:
    Data(const Data& );
    Data& operator=(const Data& );
};

/// Specialization for reading strings
template<>
bool TData<std::string>::read( const std::string& str );


/// Specialization for reading booleans
template<>
bool TData<bool>::read( const std::string& str );


/// General case for printing default value
template<class T>
inline
void TData<T>::printValue( std::ostream& out) const
{
    out << virtualGetValue() << " ";
}

/// General case for printing default value
template<class T>
inline
std::string TData<T>::getValueString() const
{
    std::ostringstream out;
    out << virtualGetValue();
    return out.str();
}

template<class T>
inline
std::string TData<T>::getValueTypeString() const
{
    return BaseData::typeName(&virtualGetValue());
}


#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_BUILD_CORE)

extern template class SOFA_CORE_API TData< std::string >;
extern template class SOFA_CORE_API Data< std::string >;
extern template class SOFA_CORE_API TData< bool >;
extern template class SOFA_CORE_API Data< bool >;

#endif

} // namespace objectmodel

} // namespace core

// Overload helper::ReadAccessor and helper::WriteAccessor

namespace helper
{

template<class T>
class ReadAccessor< core::objectmodel::Data<T> > : public ReadAccessor<T>
{
public:
    typedef ReadAccessor<T> Inherit;
    typedef core::objectmodel::Data<T> data_container_type;
    typedef T container_type;

protected:
    const data_container_type* data;
public:
    ReadAccessor(const data_container_type& d) : Inherit(d.getValue()), data(&d) {}
    ReadAccessor(const data_container_type* d) : Inherit(d->getValue()), data(d) {}
    ReadAccessor(const core::ExecParams* params, const data_container_type& d) : Inherit(d.getValue(params)), data(&d) {}
    ReadAccessor(const core::ExecParams* params, const data_container_type* d) : Inherit(d->getValue(params)), data(d) {}
};

template<class T>
class WriteAccessor< core::objectmodel::Data<T> > : public WriteAccessor<T>
{
public:
    typedef WriteAccessor<T> Inherit;
    typedef core::objectmodel::Data<T> data_container_type;
    typedef T container_type;

	// these are forbidden (until c++11 move semantics) as they break
	// RAII encapsulation. the reference member 'data' prevents them
	// anyways, but the intent is more obvious like this.
	WriteAccessor(const WriteAccessor& );
	WriteAccessor& operator=(const WriteAccessor& );

protected:
    data_container_type& data;
    const core::ExecParams* dparams;

public:
    WriteAccessor(data_container_type& d) : Inherit(*d.beginEdit()), data(d), dparams(NULL) {}
    WriteAccessor(data_container_type* d) : Inherit(*d->beginEdit()), data(*d), dparams(NULL) {}
    WriteAccessor(const core::ExecParams* params, data_container_type& d) : Inherit(*d.beginEdit(params)), data(d), dparams(params) {}
    WriteAccessor(const core::ExecParams* params, data_container_type* d) : Inherit(*d->beginEdit(params)), data(*d), dparams(params) {}
    ~WriteAccessor() { if (dparams) data.endEdit(dparams); else data.endEdit(); }
};

/// Easy syntax for getting write access to a Data using operator ->. Example: write(someFlagData)->setFlagValue(true);
template<class T>
inline WriteAccessor<core::objectmodel::Data<T> > write(core::objectmodel::Data<T>& data) { return WriteAccessor<core::objectmodel::Data<T> >(data); }


} // namespace helper

// the Data class is used everywhere
using core::objectmodel::Data;

} // namespace sofa

#endif
