/***************************************************************************
 *  include/stxxl/bits/common/simple_vector.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2008, 2011 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_COMMON_SIMPLE_VECTOR_HEADER
#define STXXL_COMMON_SIMPLE_VECTOR_HEADER

#include <algorithm>
#include <cstring>
#include <cassert>
#include <stxxl/bits/noncopyable.h>
#include <stxxl/bits/verbose.h>
#include <stxxl/bits/common/utils.h>

STXXL_BEGIN_NAMESPACE

//! \addtogroup support
//! \{

/*!
 * Simpler non-growing vector without initialization.
 *
 * simple_vector can be used a replacement for std::vector when only a
 * non-growing array is needed. The advantages of simple_vector are that it
 * does not initilize memory for POD types (faster), allows simpler inlines and
 * is less error prone to copying and other problems..
 */
template <typename ValueType>
class simple_vector : private noncopyable
{
public:
    typedef ValueType value_type;
    typedef size_t size_type;

protected:
    //! size of allocated memory
    size_type m_size;

    //! pointer to allocated memory area
    value_type* m_array;

public:
    // *** simple pointer iterators

    typedef value_type* iterator;
    typedef const value_type* const_iterator;
    typedef value_type& reference;
    typedef const value_type& const_reference;

public:
    //! allocate empty simple vector
    simple_vector()
        : m_size(0), m_array(NULL)
    { }
    //! allocate vector's memory
    simple_vector(size_type sz)
        : m_size(sz), m_array(NULL)
    {
        if (m_size > 0)
            m_array = new value_type[m_size];
    }
    //! swap vector with another one
    void swap(simple_vector& obj)
    {
        std::swap(m_size, obj.m_size);
        std::swap(m_array, obj.m_array);
    }
    //! delete vector
    ~simple_vector()
    {
        delete[] m_array;
    }
    //! return iterator to beginning of vector
    iterator data()
    {
        return m_array;
    }
    //! return iterator to beginning of vector
    const_iterator data() const
    {
        return m_array;
    }
    //! return mutable iterator to first element
    iterator begin()
    {
        return m_array;
    }
    //! return constant iterator to first element
    const_iterator begin() const
    {
        return m_array;
    }
    //! return constant iterator to first element
    const_iterator cbegin() const
    {
        return begin();
    }
    //! return mutable iterator beyond last element
    iterator end()
    {
        return m_array + m_size;
    }
    //! return constant iterator beyond last element
    const_iterator end() const
    {
        return m_array + m_size;
    }
    //! return constant iterator beyond last element
    const_iterator cend() const
    {
        return end();
    }
    //! return number of items in vector
    size_type size() const
    {
        return m_size;
    }
    //! return the i-th position of the vector
    reference operator [] (size_type i)
    {
        assert(i < m_size);
        return *(begin() + i);
    }
    //! return constant reference to the i-th position of the vector
    const_reference operator [] (size_type i) const
    {
        assert(i < m_size);
        return *(begin() + i);
    }
    //! resize the array to contain exactly newsize items
    void resize(size_type newsize)
    {
        if (m_array)
        {
            STXXL_MSG("Warning: resizing non-empty simple_vector");
            value_type* tmp = m_array;
            m_array = new value_type[newsize];
            memcpy((void*)m_array, (void*)tmp,
                   sizeof(value_type) * STXXL_MIN(m_size, newsize));
            delete[] tmp;
            m_size = newsize;
        }
        else
        {
            m_array = new value_type[newsize];
            m_size = newsize;
        }
    }
    //! Zero the whole array content.
    void memzero()
    {
        memset(m_array, 0, m_size * sizeof(value_type));
    }
};

// \}

STXXL_END_NAMESPACE

namespace std {

template <class ValueType>
void swap(stxxl::simple_vector<ValueType>& a,
          stxxl::simple_vector<ValueType>& b)
{
    a.swap(b);
}

} // namespace std

#endif // !STXXL_COMMON_SIMPLE_VECTOR_HEADER
