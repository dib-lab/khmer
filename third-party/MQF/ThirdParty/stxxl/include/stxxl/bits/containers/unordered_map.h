/***************************************************************************
 *  include/stxxl/bits/containers/unordered_map.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2008 Markus Westphal <marwes@users.sourceforge.net>
 *  Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_CONTAINERS_UNORDERED_MAP_HEADER
#define STXXL_CONTAINERS_UNORDERED_MAP_HEADER

#include <stxxl/bits/noncopyable.h>
#include <stxxl/bits/containers/hash_map/hash_map.h>

STXXL_BEGIN_NAMESPACE

namespace hash_map {

template <
    class KeyType,
    class DataType,
    class HashType,
    class CompareType,
    unsigned SubBlockSize,
    unsigned SubBlocksPerBlock,
    class Alloc
    >
class hash_map;

} // namespace hash_map

//! \addtogroup stlcont
//! \{

/*!
 * An external memory implementation of the STL unordered_map container, which
 * is based on an external memory hash map. For more information see \ref
 * tutorial_unordered_map.
 *
 * \tparam KeyType the key type
 * \tparam MappedType the mapped type associated with a key
 * \tparam HashType a hash functional
 * \tparam CompareType a less comparison relation for KeyType
 * \tparam SubBlockSize the raw size of a subblock (caching granularity)
 * (default: 8192)
 * \tparam SubBlocksPerBlock the number of subblocks per external block
 * (default: 256 -> 2MB blocks)
 * \tparam AllocType allocator for internal-memory buffer
 */
template <
    class KeyType,
    class MappedType,
    class HashType,
    class CompareType,
    unsigned SubBlockSize = 8* 1024,
    unsigned SubBlocksPerBlock = 256,
    class AllocType = std::allocator<std::pair<const KeyType, MappedType> >
    >
class unordered_map : private noncopyable
{
    typedef hash_map::hash_map<KeyType, MappedType, HashType, CompareType,
                               SubBlockSize, SubBlocksPerBlock, AllocType> impl_type;

    impl_type impl;

public:
    //! \name Types
    //! \{

    //! the first template parameter (Key)
    typedef typename impl_type::key_type key_type;
    //! the second template parameter (T)
    typedef typename impl_type::mapped_type mapped_type;
    //! pair<const key_type,mapped_type>
    typedef typename impl_type::value_type value_type;
    //! the third template parameter (HashType)
    typedef typename impl_type::hasher hasher;
    //! the fourth template parameter (CompareType) (!!! not: equality compare)
    typedef typename impl_type::key_compare key_compare;
    //! the fifth template parameter (AllocType)
    typedef AllocType allocator_type;

    typedef typename impl_type::reference reference;
    typedef typename impl_type::const_reference const_reference;
    typedef typename impl_type::pointer pointer;
    typedef typename impl_type::const_pointer const_pointer;

    typedef typename impl_type::external_size_type size_type;
    typedef typename impl_type::difference_type difference_type;

    typedef typename impl_type::external_size_type external_size_type;
    typedef typename impl_type::internal_size_type internal_size_type;

    typedef typename impl_type::iterator iterator;
    typedef typename impl_type::const_iterator const_iterator;

    //! constructed equality predicate for key
    typedef typename impl_type::key_equal key_equal;

    //! \}

    //! \name Constructors
    //! \{

    /*!
     * Construct a new hash-map
     *
     * \param n initial number of buckets
     * \param hf hash-function
     * \param cmp comparator-object
     * \param buffer_size size of internal-memory buffer in bytes
     * \param a allocation-strategory for internal-memory buffer
     */
    unordered_map(internal_size_type n = 0,
                  const hasher& hf = hasher(),
                  const key_compare& cmp = key_compare(),
                  internal_size_type buffer_size = 100*1024*1024,
                  const allocator_type& a = allocator_type())
        : impl(n, hf, cmp, buffer_size, a)
    { }

    /*!
     * Construct a new hash-map and insert all values in the range [begin,end)
     *
     * \param begin beginning of the range
     * \param end end of the range
     * \param mem_to_sort internal memory that may be used for
     * bulk-construction (not to be confused with the buffer-memory)
     * \param n initial number of buckets
     * \param hf hash-function
     * \param cmp comparator-object
     * \param buffer_size size of internal-memory buffer in bytes
     * \param a allocation-strategory for internal-memory buffer
     */
    template <class InputIterator>
    unordered_map(InputIterator begin, InputIterator end,
                  internal_size_type mem_to_sort = 256*1024*1024,
                  internal_size_type n = 0,
                  const hasher& hf = hasher(),
                  const key_compare& cmp = key_compare(),
                  internal_size_type buffer_size = 100*1024*1024,
                  const allocator_type& a = allocator_type())
        : impl(begin, end, mem_to_sort, n, hf, cmp, buffer_size, a)
    { }

    //! \}

    //! \name Size and Capacity
    //! \{

    //! Number of values currently stored. Note: If the correct number is
    //! currently unknown (because **oblivous-methods** were used), external
    //! memory will be scanned.
    external_size_type size() const
    {
        return impl.size();
    }

    //! The hash-map may store up to this number of values
    external_size_type max_size() const
    {
        return impl.max_size();
    }

    //! Check if container is empty, see size() about oblivious-methods.
    bool empty() const
    {
        return impl.empty();
    }

    //! \}

    //! \name Iterators
    //! \{

    //! iterator pointing to the beginnning of the hash-map
    iterator begin()
    {
        return impl.begin();
    }

    //! iterator pointing to the end of the hash-map (iterator-type as
    //! template-parameter)
    iterator end()
    {
        return impl.end();
    }

    //! iterator pointing to the beginnning of the hash-map
    const_iterator begin() const
    {
        return impl.begin();
    }

    //! iterator pointing to the end of the hash-map (iterator-type as
    //! template-parameter)
    const_iterator end() const
    {
        return impl.end();
    }

    //! \}

    //! \name Lookup and Element Access
    //! \{

    //! Convenience operator to quickly insert or find values. Use with caution
    //! since using this operator will check external-memory.
    mapped_type& operator [] (const key_type& key)
    {
        return impl[key];
    }

    //! Look up value by key. Non-const access.
    //! \param key key for value to look up
    iterator find(const key_type& key)
    {
        return impl.find(key);
    }

    //! Look up value by key. Const access.
    //! \param key key for value to look up
    const_iterator find(const key_type& key) const
    {
        return impl.find(key);
    }

    //! Number of values with given key
    //! \param key key for value to look up
    //! \return 0 or 1 depending on the presence of a value with the given key
    external_size_type count(const key_type& key) const
    {
        return impl.count(key);
    }

    //! Finds a range containing all values with given key. Non-const access
    //! \param key key to look for#
    //! \return range may be empty or contains exactly one value
    std::pair<iterator, iterator>
    equal_range(const key_type& key)
    {
        return impl.equal_range(key);
    }

    //! Finds a range containing all values with given key. Const access
    //! \param key key to look for#
    //! \return range may be empty or contains exactly one value
    std::pair<const_iterator, const_iterator>
    equal_range(const key_type& key) const
    {
        return impl.equal_range(key);
    }

    //! \}

    //! \name Modifiers: Insert
    //! \{

    /*!
     * Insert a new value if no value with the same key is already present;
     * external memory must therefore be accessed
     *
     * \param value what to insert
     * \return a tuple whose second part is true iff the value was actually
     * added (no value with the same key present); the first part is an
     * iterator pointing to the newly inserted or already stored value
     */
    std::pair<iterator, bool> insert(const value_type& value)
    {
        return impl.insert(value);
    }

    //! Insert a value; external memory is not accessed so that another value
    //! with the same key may be overwritten
    //! \param value what to insert
    //! \return iterator pointing to the inserted value
    iterator insert_oblivious(const value_type& value)
    {
        return impl.insert_oblivious(value);
    }

    //! Bulk-insert of values in the range [f, l)
    //! \param first beginning of the range
    //! \param last end of the range
    //! \param mem internal memory that may be used (note: this memory will be
    //! used additionally to the buffer). The more the better
    template <class InputIterator>
    void insert(InputIterator first, InputIterator last, internal_size_type mem)
    {
        impl.insert(first, last, mem);
    }

    //! \}

    //! \name Modifiers: Erase
    //! \{

    //! Erase value by iterator
    //! \param it iterator pointing to the value to erase
    void erase(const_iterator it)
    {
        impl.erase(it);
    }

    //! Erase value by key; check external memory
    //! \param key key of value to erase
    //! \return number of values actually erased (0 or 1)
    external_size_type erase(const key_type& key)
    {
        return impl.erase(key);
    }

    //! Erase value by key but without looking at external memory
    //! \param key key for value to release
    void erase_oblivious(const key_type& key)
    {
        impl.erase_oblivious(key);
    }

    //! Reset hash-map: erase all values, invalidate all iterators
    void clear()
    {
        impl.clear();
    }

    //! Exchange stored values with another hash-map
    //! \param obj hash-map to swap values with
    void swap(unordered_map& obj)
    {
        std::swap(impl, obj.impl);
    }

    //! \}

    //! \name Bucket Interface
    //! \{

    //! Number of buckets
    internal_size_type bucket_count() const
    {
        return impl.bucket_count();
    }

    //! Maximum number of buckets
    internal_size_type max_bucket_count() const
    {
        return impl.max_bucket_count();
    }

    // bucket_size()?

    //! Bucket-index for values with given key.
    internal_size_type bucket(const key_type& k) const
    {
        return impl.bucket_index(k);
    }

    //! \}

    //! \name Hash Policy
    //! \{

    //! Average number of (sub)blocks occupied by a bucket.
    float load_factor() const
    {
        return impl.load_factor();
    }

    //! Set desired load-factor
    float opt_load_factor() const
    {
        return impl.opt_load_factor();
    }

    //! Set desired load-factor
    void opt_load_factor(float z)
    {
        impl.opt_load_factor(z);
    }

    //! Rehash with (at least) n buckets
    void rehash(internal_size_type n)
    {
        impl.rehash(n);
    }

    //! \}

    //! \name Observers
    //! \{

    //! Hash-function used by this hash-map
    hasher hash_function() const
    {
        return impl.hash_function();
    }

    //! Strict-weak-ordering used by this hash-map
    key_compare key_comp() const
    {
        return impl.key_cmp();
    }

    //! Constructed equality predicate used by this hash-map
    key_equal key_eq() const
    {
        return impl.key_eq();
    }

    //! Get node memory allocator
    allocator_type get_allocator() const
    {
        return impl.get_allocator();
    }

    //! \}

    //! \name Internal Memory Buffer Policy
    //! \{

    //! Number of bytes occupied by buffer
    internal_size_type buffer_size() const
    {
        return impl.buffer_size();
    }

    //! Maximum buffer size in byte
    internal_size_type max_buffer_size() const
    {
        return impl.max_buffer_size();
    }

    //! Set maximum buffer size
    //! \param buffer_size new size in byte
    void max_buffer_size(internal_size_type buffer_size)
    {
        impl.max_buffer_size(buffer_size);
    }

    //! \}

    //! \name Statistics
    //! \{

    //! Reset hash-map statistics
    void reset_statistics()
    {
        impl.reset_statistics();
    }

    //! Print short general statistics to output stream
    void print_statistics(std::ostream& o = std::cout) const
    {
        impl.print_statistics(o);
    }

    //! Even more statistics: Number of buckets, number of values, buffer-size,
    //! values per bucket
    void print_load_statistics(std::ostream& o = std::cout) const
    {
        impl.print_load_statistics(o);
    }

    // \}
};

//! \}

STXXL_END_NAMESPACE

namespace std {

template <
    class KeyType,
    class MappedType,
    class HashType,
    class CompareType,
    unsigned SubBlockSize,
    unsigned SubBlocksPerBlock,
    class AllocType
    >
void swap(stxxl::unordered_map<KeyType, MappedType, HashType, CompareType,
                               SubBlockSize, SubBlocksPerBlock, AllocType>& a,
          stxxl::unordered_map<KeyType, MappedType, HashType, CompareType,
                               SubBlockSize, SubBlocksPerBlock, AllocType>& b
          )
{
    a.swap(b);
}

} // namespace std

#endif // !STXXL_CONTAINERS_UNORDERED_MAP_HEADER
