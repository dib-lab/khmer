/***************************************************************************
 *  include/stxxl/bits/common/swap_vector.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2014 Thomas Keh <thomas.keh@student.kit.edu>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_COMMON_SWAP_VECTOR_HEADER
#define STXXL_COMMON_SWAP_VECTOR_HEADER

#include <algorithm>
#include <cassert>
#include <stxxl/bits/verbose.h>
#include <stxxl/bits/noncopyable.h>

STXXL_BEGIN_NAMESPACE

//! \addtogroup support
//! \{

/*!
 * Vector that avoids copying of ValueType objects in push_back() (here:
 * swap_back()) and resize() operations.  Values are swapped with
 * default-constructed instances instead.  Important: A template spezialization
 * for std::swap(ValueType&,ValueType&) must be provided. Make shure the swap
 * implementation is located above these lines.
 */
template <typename ValueType>
class swap_vector : private noncopyable
{
public:
    typedef ValueType value_type;
    typedef size_t size_type;

protected:
    //! size of vector
    size_type m_size;

    //! size of allocated memory
    size_type m_capacity;

    //! pointer to allocated memory area
    value_type* m_array;

public:
    // *** simple pointer iterators

    typedef value_type* iterator;
    typedef const value_type* const_iterator;
    typedef value_type& reference;
    typedef const value_type& const_reference;

public:
    //! Create an empty vector.
    swap_vector()
        : m_size(0), m_capacity(0), m_array(NULL)
    { }
    //! Create a vector with the spezified size.
    swap_vector(size_type size)
        : m_size(size), m_capacity(size), m_array(NULL)
    {
        if (m_size > 0)
            m_array = new value_type[m_size];
    }
    //! Create a vector with the spezified size and reserve (possibly more)
    //! space.
    swap_vector(size_type size, size_type capacity)
        : m_size(size), m_capacity(std::max(size, capacity)), m_array(NULL)
    {
        if (m_capacity > 0)
            m_array = new value_type[m_capacity];
    }
    //! Swap the vector with another one.
    void swap(swap_vector& obj)
    {
        using std::swap;
        swap(m_size, obj.m_size);
        swap(m_capacity, obj.m_capacity);
        swap(m_array, obj.m_array);
    }
    //! Delete the vector.
    ~swap_vector()
    {
        delete[] m_array;
    }
    //! Return the vector size.
    size_type size() const
    {
        return m_size;
    }
    //! Return the vector size.
    bool empty() const
    {
        return (m_size == 0);
    }
    //! Return the size of the underlaying array.
    size_type capacity() const
    {
        return m_capacity;
    }
    //! Return iterator to the beginning of vector.
    iterator data()
    {
        return m_array;
    }
    //! Return iterator to the beginning of vector.
    const_iterator data() const
    {
        return m_array;
    }
    //! Return mutable iterator to the first element.
    iterator begin()
    {
        return m_array;
    }
    //! Return constant iterator to the first element.
    const_iterator begin() const
    {
        return m_array;
    }
    //! Return constant iterator to the first element.
    const_iterator cbegin() const
    {
        return begin();
    }
    //! Return mutable iterator beyond the last element.
    iterator end()
    {
        return m_array + m_size;
    }
    //! Return constant iterator beyond the last element.
    const_iterator end() const
    {
        return m_array + m_size;
    }
    //! Return constant iterator beyond the last element.
    const_iterator cend() const
    {
        return end();
    }
    //! Return the i-th position of the vector.
    reference operator [] (size_type i)
    {
        assert(i < m_size);
        return *(begin() + i);
    }
    //! Return constant reference to the i-th position of the vector.
    const_reference operator [] (size_type i) const
    {
        assert(i < m_size);
        return *(begin() + i);
    }
    //! Return reference to first element.
    reference front()
    {
        assert(m_size > 0);
        return *(m_array);
    }
    //! Return constant reference to first element.
    const_reference front() const
    {
        assert(m_size > 0);
        return *(m_array);
    }
    //! Return reference to last element.
    reference back()
    {
        assert(m_size > 0);
        return *(m_array + m_size - 1);
    }
    //! Return reference to last element.
    const_reference back() const
    {
        assert(m_size > 0);
        return *(m_array + m_size - 1);
    }
    //! Resize the underlaying array to contain at least newsize items.
    void reserve(size_type newsize)
    {
        if (newsize > m_capacity)
        {
            if (m_array)
            {
                value_type* tmp = m_array;
                m_array = new value_type[newsize];
                for (unsigned i = 0; i < m_size; ++i)
                {
                    using std::swap;
                    swap(tmp[i], m_array[i]);
                }
                delete[] tmp;
            }
            else
            {
                m_array = new value_type[newsize];
            }
            m_capacity = newsize;
        }
    }
    //! Resize the vector to contain at least newsize items.
    void resize(size_type newsize)
    {
        reserve(newsize);
        m_size = newsize;
    }
    //! Create a new value_type object at the end of the vector and then swap
    //! it with val.
    void swap_back(reference val)
    {
        if (m_size + 1 > m_capacity)
        {
            reserve(std::max((size_type)2, 2 * m_capacity));
        }
        using std::swap;
        swap(m_array[m_size], val);
        ++m_size;
    }
    //! Clear the vector. The capacity doesn't change.
    void clear()
    {
        m_size = 0;
    }
    //! Erase the element at the given position by swapping it to the and and
    //! then reducing the vector size.
    iterator erase(iterator position)
    {
        return erase(position, position + 1);
    }
    //! Erase the elements at in the range [begin,last) by swapping them to the
    //! and and then reducing the vector size.
    iterator erase(iterator first, iterator last)
    {
        assert(first >= begin());
        assert(last <= end());
        if (last < m_array + m_size - 1)
        {
            // iteratively swap forward the elements behind last
            iterator f = first;
            iterator l = last;
            while (l < end())
            {
                using std::swap;
                swap(*f, *l);
                ++f;
                ++l;
            }
        }
        m_size -= (last - first);
        return first;
    }
};

/*!
 * Transforms the range [first,last) into a range with all the elements for
 * which pred returns true removed, and returns an iterator to the new end of
 * that range.
 *
 * The function is compatible to std::remove_if, but uses std::swap instead of
 * copy-assignment (resp. move-assignment in C++11).
 */
template <class ForwardIterator, class UnaryPredicate>
ForwardIterator swap_remove_if(ForwardIterator first, ForwardIterator last, UnaryPredicate pred)
{
    ForwardIterator result = first;
    while (first != last)
    {
        if (!pred(*first))
        {
            using std::swap;
            swap(*first, *result);
            ++result;
        }
        ++first;
    }
    return result;
}

// \}

STXXL_END_NAMESPACE

namespace std {

template <class ValueType>
void swap(stxxl::swap_vector<ValueType>& a,
          stxxl::swap_vector<ValueType>& b)
{
    a.swap(b);
}

} // namespace std

#endif // !STXXL_COMMON_SWAP_VECTOR_HEADER
