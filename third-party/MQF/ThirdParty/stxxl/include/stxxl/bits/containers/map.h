/***************************************************************************
 *  include/stxxl/bits/containers/map.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2006 Roman Dementiev <dementiev@ira.uka.de>
 *  Copyright (C) 2008, 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_CONTAINERS_MAP_HEADER
#define STXXL_CONTAINERS_MAP_HEADER

#include <stxxl/bits/noncopyable.h>
#include <stxxl/bits/containers/btree/btree.h>

STXXL_BEGIN_NAMESPACE

namespace btree {

template <class KeyType,
          class DataType,
          class CompareType,
          unsigned LogNodeSize,
          unsigned LogLeafSize,
          class PDAllocStrategy
          >
class btree;

} // namespace btree

//! \addtogroup stlcont
//! \{

//! External associative container (map). \n
//! <b> Introduction </b> to map container: see \ref tutorial_map tutorial. \n
//! <b> Design and Internals </b> of map container: see \ref design_map
//!
//! \tparam KeyType key type (POD with no references to internal memory)
//! \tparam DataType data type (POD with no references to internal memory)
//! \tparam CompareType comparison type used to determine
//! whether a key is smaller than another one.
//! If CompareType()(x,y) is true, then x is smaller than y.
//! CompareType must also provide a static \c max_value method, that returns
//! a value of type KeyType that is
//! larger than any key stored in map : i.e. for all \b x in map holds
//! CompareType()(x,CompareType::max_value())
//!
//! <BR>
//! Example: :
//! \verbatim
//! struct CmpIntGreater
//! {
//!   bool operator () (const int & a, const int & b) const { return a>b; }
//!   static int max_value() { return std::numeric_limits<int>::min(); }
//! };
//! \endverbatim
//! Another example:
//! \verbatim
//! struct CmpIntLess
//! {
//!   bool operator () (const int & a, const int & b) const { return a<b; }
//!   static int max_value() const  { return std::numeric_limits<int>::max(); }
//! };
//! \endverbatim
//! Note that CompareType must define a strict weak ordering.
//! (<A HREF="http://www.sgi.com/tech/stl/StrictWeakOrdering.html">see what it is</A>)
//! \tparam RawNodeSize size of internal nodes of map in bytes (btree implementation).
//! \tparam RawLeafSize size of leaves of map in bytes (btree implementation).
//! \tparam PDAllocStrategy parallel disk allocation strategy (\c stxxl::SR is recommended and default)
//!
template <class KeyType,
          class DataType,
          class CompareType,
          unsigned RawNodeSize = 16* 1024,      // 16 KBytes default
          unsigned RawLeafSize = 128* 1024,     // 128 KBytes default
          class PDAllocStrategy = stxxl::SR
          >
class map : private noncopyable
{
    typedef btree::btree<KeyType, DataType, CompareType, RawNodeSize, RawLeafSize, PDAllocStrategy> impl_type;

    impl_type impl;

public:
    typedef typename impl_type::node_block_type node_block_type;
    typedef typename impl_type::leaf_block_type leaf_block_type;

    typedef typename impl_type::key_type key_type;
    typedef typename impl_type::data_type data_type;
    typedef typename impl_type::data_type mapped_type;
    typedef typename impl_type::value_type value_type;
    typedef typename impl_type::key_compare key_compare;
    typedef typename impl_type::value_compare value_compare;
    typedef typename impl_type::pointer pointer;
    typedef typename impl_type::const_pointer const_pointer;
    typedef typename impl_type::reference reference;
    typedef typename impl_type::const_reference const_reference;
    typedef typename impl_type::size_type size_type;
    typedef typename impl_type::difference_type difference_type;
    typedef typename impl_type::iterator iterator;
    typedef typename impl_type::const_iterator const_iterator;
    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    //! \name Iterators
    //! \{

    iterator begin() { return impl.begin(); }
    iterator end() { return impl.end(); }
    const_iterator begin() const { return impl.begin(); }
    const_iterator end() const { return impl.end(); }
    const_iterator cbegin() const { return begin(); }
    const_iterator cend() const { return end(); }

    reverse_iterator rbegin()
    {
        return reverse_iterator(end());
    }
    const_reverse_iterator rbegin() const
    {
        return const_reverse_iterator(end());
    }
    const_reverse_iterator crbegin() const
    {
        return const_reverse_iterator(end());
    }
    reverse_iterator rend()
    {
        return reverse_iterator(begin());
    }
    const_reverse_iterator rend() const
    {
        return const_reverse_iterator(begin());
    }
    const_reverse_iterator crend() const
    {
        return const_reverse_iterator(begin());
    }

    //! \}

    //! \name Capacity
    //! \{

    size_type size() const { return impl.size(); }
    size_type max_size() const { return impl.max_size(); }
    bool empty() const { return impl.empty(); }

    //! \}

    //! \name Observers
    //! \{

    key_compare key_comp() const { return impl.key_comp(); }
    value_compare value_comp() const { return impl.value_comp(); }

    //! \}

    //! \name Constructors/Destructors
    //! \{

    //! A constructor
    //! \param node_cache_size_in_bytes size of node cache in bytes (btree implementation)
    //! \param leaf_cache_size_in_bytes size of leaf cache in bytes (btree implementation)
    map(unsigned_type node_cache_size_in_bytes,
        unsigned_type leaf_cache_size_in_bytes
        ) : impl(node_cache_size_in_bytes, leaf_cache_size_in_bytes)
    { }

    //! A constructor
    //! \param c_ comparator object
    //! \param node_cache_size_in_bytes size of node cache in bytes (btree implementation)
    //! \param leaf_cache_size_in_bytes size of leaf cache in bytes (btree implementation)
    map(const key_compare& c_,
        unsigned_type node_cache_size_in_bytes,
        unsigned_type leaf_cache_size_in_bytes
        ) : impl(c_, node_cache_size_in_bytes, leaf_cache_size_in_bytes)
    { }

    //! Constructs a map from a given input range
    //! \param b beginning of the range
    //! \param e end of the range
    //! \param node_cache_size_in_bytes size of node cache in bytes (btree implementation)
    //! \param leaf_cache_size_in_bytes size of leaf cache in bytes (btree implementation)
    //! \param range_sorted if \c true than the constructor assumes that the range is sorted
    //! and performs a fast bottom-up bulk construction of the map (btree implementation)
    //! \param node_fill_factor node fill factor in [0,1] for bulk construction
    //! \param leaf_fill_factor leaf fill factor in [0,1] for bulk construction
    template <class InputIterator>
    map(InputIterator b,
        InputIterator e,
        unsigned_type node_cache_size_in_bytes,
        unsigned_type leaf_cache_size_in_bytes,
        bool range_sorted = false,
        double node_fill_factor = 0.75,
        double leaf_fill_factor = 0.6
        ) : impl(b, e, node_cache_size_in_bytes, leaf_cache_size_in_bytes,
                 range_sorted, node_fill_factor, leaf_fill_factor)
    { }

    //! Constructs a map from a given input range
    //! \param b beginning of the range
    //! \param e end of the range
    //! \param c_ comparator object
    //! \param node_cache_size_in_bytes size of node cache in bytes (btree implementation)
    //! \param leaf_cache_size_in_bytes size of leaf cache in bytes (btree implementation)
    //! \param range_sorted if \c true than the constructor assumes that the range is sorted
    //! and performs a fast bottom-up bulk construction of the map (btree implementation)
    //! \param node_fill_factor node fill factor in [0,1] for bulk construction
    //! \param leaf_fill_factor leaf fill factor in [0,1] for bulk construction
    template <class InputIterator>
    map(InputIterator b,
        InputIterator e,
        const key_compare& c_,
        unsigned_type node_cache_size_in_bytes,
        unsigned_type leaf_cache_size_in_bytes,
        bool range_sorted = false,
        double node_fill_factor = 0.75,
        double leaf_fill_factor = 0.6
        ) : impl(b, e, c_, node_cache_size_in_bytes, leaf_cache_size_in_bytes,
                 range_sorted, node_fill_factor, leaf_fill_factor)
    { }

    //! \}

    //! \name Modifiers
    //! \{

    void swap(map& obj) { std::swap(impl, obj.impl); }
    std::pair<iterator, bool> insert(const value_type& x)
    {
        return impl.insert(x);
    }
    iterator insert(iterator pos, const value_type& x)
    {
        return impl.insert(pos, x);
    }
    template <class InputIterator>
    void insert(InputIterator b, InputIterator e)
    {
        impl.insert(b, e);
    }
    void erase(iterator pos)
    {
        impl.erase(pos);
    }
    size_type erase(const key_type& k)
    {
        return impl.erase(k);
    }
    void erase(iterator first, iterator last)
    {
        impl.erase(first, last);
    }
    void clear()
    {
        impl.clear();
    }

    //! \}

    //! \name Operations
    //! \{

    iterator find(const key_type& k)
    {
        return impl.find(k);
    }
    const_iterator find(const key_type& k) const
    {
        return impl.find(k);
    }
    size_type count(const key_type& k)
    {
        return impl.count(k);
    }
    iterator lower_bound(const key_type& k)
    {
        return impl.lower_bound(k);
    }
    const_iterator lower_bound(const key_type& k) const
    {
        return impl.lower_bound(k);
    }
    iterator upper_bound(const key_type& k)
    {
        return impl.upper_bound(k);
    }
    const_iterator upper_bound(const key_type& k) const
    {
        return impl.upper_bound(k);
    }
    std::pair<iterator, iterator> equal_range(const key_type& k)
    {
        return impl.equal_range(k);
    }
    std::pair<const_iterator, const_iterator> equal_range(const key_type& k) const
    {
        return impl.equal_range(k);
    }

    //! \}

    //! \name Operators
    //! \{

    data_type& operator [] (const key_type& k)
    {
        return impl[k];
    }

    //! \}

    //! \name Miscellaneous
    //! \{

    //! Enables leaf prefetching during scanning
    void enable_prefetching()
    {
        impl.enable_prefetching();
    }

    //! Disables leaf prefetching during scanning
    void disable_prefetching()
    {
        impl.disable_prefetching();
    }

    //! Returns the status of leaf prefetching during scanning
    bool prefetching_enabled()
    {
        return impl.prefetching_enabled();
    }

    //! Prints cache statistics
    void print_statistics(std::ostream& o) const
    {
        impl.print_statistics(o);
    }

    //! Resets cache statistics
    void reset_statistics()
    {
        impl.reset_statistics();
    }

    //! \}

    //////////////////////////////////////////////////
    template <class KeyType_,
              class DataType_,
              class CompareType_,
              unsigned RawNodeSize_,
              unsigned RawLeafSize_,
              class PDAllocStrategy_>
    friend bool operator == (const map<KeyType_, DataType_, CompareType_, RawNodeSize_, RawLeafSize_, PDAllocStrategy_>& a,
                             const map<KeyType_, DataType_, CompareType_, RawNodeSize_, RawLeafSize_, PDAllocStrategy_>& b);
    //////////////////////////////////////////////////
    template <class KeyType_,
              class DataType_,
              class CompareType_,
              unsigned RawNodeSize_,
              unsigned RawLeafSize_,
              class PDAllocStrategy_>
    friend bool operator < (const map<KeyType_, DataType_, CompareType_, RawNodeSize_, RawLeafSize_, PDAllocStrategy_>& a,
                            const map<KeyType_, DataType_, CompareType_, RawNodeSize_, RawLeafSize_, PDAllocStrategy_>& b);
    //////////////////////////////////////////////////
    template <class KeyType_,
              class DataType_,
              class CompareType_,
              unsigned RawNodeSize_,
              unsigned RawLeafSize_,
              class PDAllocStrategy_>
    friend bool operator > (const map<KeyType_, DataType_, CompareType_, RawNodeSize_, RawLeafSize_, PDAllocStrategy_>& a,
                            const map<KeyType_, DataType_, CompareType_, RawNodeSize_, RawLeafSize_, PDAllocStrategy_>& b);
    //////////////////////////////////////////////////
    template <class KeyType_,
              class DataType_,
              class CompareType_,
              unsigned RawNodeSize_,
              unsigned RawLeafSize_,
              class PDAllocStrategy_>
    friend bool operator != (const map<KeyType_, DataType_, CompareType_, RawNodeSize_, RawLeafSize_, PDAllocStrategy_>& a,
                             const map<KeyType_, DataType_, CompareType_, RawNodeSize_, RawLeafSize_, PDAllocStrategy_>& b);
    //////////////////////////////////////////////////
    template <class KeyType_,
              class DataType_,
              class CompareType_,
              unsigned RawNodeSize_,
              unsigned RawLeafSize_,
              class PDAllocStrategy_>
    friend bool operator <= (const map<KeyType_, DataType_, CompareType_, RawNodeSize_, RawLeafSize_, PDAllocStrategy_>& a,
                             const map<KeyType_, DataType_, CompareType_, RawNodeSize_, RawLeafSize_, PDAllocStrategy_>& b);
    //////////////////////////////////////////////////
    template <class KeyType_,
              class DataType_,
              class CompareType_,
              unsigned RawNodeSize_,
              unsigned RawLeafSize_,
              class PDAllocStrategy_>
    friend bool operator >= (const map<KeyType_, DataType_, CompareType_, RawNodeSize_, RawLeafSize_, PDAllocStrategy_>& a,
                             const map<KeyType_, DataType_, CompareType_, RawNodeSize_, RawLeafSize_, PDAllocStrategy_>& b);
    //////////////////////////////////////////////////
};

template <class KeyType,
          class DataType,
          class CompareType,
          unsigned RawNodeSize,
          unsigned RawLeafSize,
          class PDAllocStrategy
          >
inline bool operator == (const map<KeyType, DataType, CompareType, RawNodeSize, RawLeafSize, PDAllocStrategy>& a,
                         const map<KeyType, DataType, CompareType, RawNodeSize, RawLeafSize, PDAllocStrategy>& b)
{
    return a.impl == b.impl;
}

template <class KeyType,
          class DataType,
          class CompareType,
          unsigned RawNodeSize,
          unsigned RawLeafSize,
          class PDAllocStrategy
          >
inline bool operator < (const map<KeyType, DataType, CompareType, RawNodeSize, RawLeafSize, PDAllocStrategy>& a,
                        const map<KeyType, DataType, CompareType, RawNodeSize, RawLeafSize, PDAllocStrategy>& b)
{
    return a.impl < b.impl;
}

template <class KeyType,
          class DataType,
          class CompareType,
          unsigned RawNodeSize,
          unsigned RawLeafSize,
          class PDAllocStrategy
          >
inline bool operator > (const map<KeyType, DataType, CompareType, RawNodeSize, RawLeafSize, PDAllocStrategy>& a,
                        const map<KeyType, DataType, CompareType, RawNodeSize, RawLeafSize, PDAllocStrategy>& b)
{
    return a.impl > b.impl;
}

template <class KeyType,
          class DataType,
          class CompareType,
          unsigned RawNodeSize,
          unsigned RawLeafSize,
          class PDAllocStrategy
          >
inline bool operator != (const map<KeyType, DataType, CompareType, RawNodeSize, RawLeafSize, PDAllocStrategy>& a,
                         const map<KeyType, DataType, CompareType, RawNodeSize, RawLeafSize, PDAllocStrategy>& b)
{
    return a.impl != b.impl;
}

template <class KeyType,
          class DataType,
          class CompareType,
          unsigned RawNodeSize,
          unsigned RawLeafSize,
          class PDAllocStrategy
          >
inline bool operator <= (const map<KeyType, DataType, CompareType, RawNodeSize, RawLeafSize, PDAllocStrategy>& a,
                         const map<KeyType, DataType, CompareType, RawNodeSize, RawLeafSize, PDAllocStrategy>& b)
{
    return a.impl <= b.impl;
}

template <class KeyType,
          class DataType,
          class CompareType,
          unsigned RawNodeSize,
          unsigned RawLeafSize,
          class PDAllocStrategy
          >
inline bool operator >= (const map<KeyType, DataType, CompareType, RawNodeSize, RawLeafSize, PDAllocStrategy>& a,
                         const map<KeyType, DataType, CompareType, RawNodeSize, RawLeafSize, PDAllocStrategy>& b)
{
    return a.impl >= b.impl;
}

//! \}

STXXL_END_NAMESPACE

namespace std {

template <class KeyType,
          class DataType,
          class CompareType,
          unsigned RawNodeSize,
          unsigned RawLeafSize,
          class PDAllocStrategy
          >
void swap(stxxl::map<KeyType, DataType, CompareType, RawNodeSize, RawLeafSize, PDAllocStrategy>& a,
          stxxl::map<KeyType, DataType, CompareType, RawNodeSize, RawLeafSize, PDAllocStrategy>& b
          )
{
    a.swap(b);
}

} // namespace std

#endif // !STXXL_CONTAINERS_MAP_HEADER
