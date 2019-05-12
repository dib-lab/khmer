/***************************************************************************
 *  include/stxxl/bits/containers/hash_map/hash_map.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2007 Markus Westphal <marwes@users.sourceforge.net>
 *  Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_CONTAINERS_HASH_MAP_HASH_MAP_HEADER
#define STXXL_CONTAINERS_HASH_MAP_HASH_MAP_HEADER

#include <functional>

#include <stxxl/bits/noncopyable.h>
#include <stxxl/bits/namespace.h>
#include <stxxl/bits/mng/block_manager.h>
#include <stxxl/bits/common/tuple.h>
#include <stxxl/bits/stream/stream.h>
#include <stxxl/bits/stream/sort_stream.h>

#include <stxxl/bits/containers/hash_map/iterator.h>
#include <stxxl/bits/containers/hash_map/iterator_map.h>
#include <stxxl/bits/containers/hash_map/block_cache.h>
#include <stxxl/bits/containers/hash_map/util.h>

STXXL_BEGIN_NAMESPACE

#define STXXL_VERBOSE_HASH_MAP(m) \
    STXXL_VERBOSE1("hash_map[" << static_cast<const void*>(this) << "]::" << m)

//! External memory hash-map
namespace hash_map {

/*!
 * Main implementation of external memory hash map.
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
template <class KeyType,
          class MappedType,
          class HashType,
          class KeyCompareType,
          unsigned SubBlockSize = 4*1024,
          unsigned SubBlocksPerBlock = 256,
          class AllocatorType = std::allocator<std::pair<const KeyType, MappedType> >
          >
class hash_map : private noncopyable
{
protected:
    typedef hash_map<KeyType, MappedType, HashType, KeyCompareType,
                     SubBlockSize, SubBlocksPerBlock> self_type;

public:
    //! type of the keys being used
    typedef KeyType key_type;
    //! type of the data to be stored
    typedef MappedType mapped_type;
    //! actually store (key-data)-pairs
    typedef std::pair<KeyType, MappedType> value_type;
    //! type for value-references
    typedef value_type& reference;
    //! type for constant value-references
    typedef const value_type& const_reference;
    //! pointer to type of keys
    typedef value_type* pointer;
    //! const pointer to type of keys
    typedef value_type const* const_pointer;

    typedef stxxl::external_size_type external_size_type;
    typedef stxxl::internal_size_type internal_size_type;
    typedef stxxl::int64 difference_type;

    //! type of (mother) hash-function
    typedef HashType hasher;
    //! functor that imposes a ordering on keys (but see _lt())
    typedef KeyCompareType key_compare;
    //! allocator template type
    typedef AllocatorType allocator_type;

    typedef hash_map_iterator<self_type> iterator;
    typedef hash_map_const_iterator<self_type> const_iterator;

    //! subblock- and block-size in bytes
    enum {
        block_raw_size = SubBlocksPerBlock * SubBlockSize,
        subblock_raw_size = SubBlockSize
    };

    //! Subblock-size as number of elements, block-size as number of subblocks
    enum {
        subblocks_per_block = SubBlocksPerBlock,
        subblock_size = SubBlockSize / sizeof(value_type)
    };

    //! a subblock consists of subblock_size values
    typedef typed_block<subblock_raw_size, value_type> subblock_type;
    //! a block consists of block_size subblocks
    typedef typed_block<block_raw_size, subblock_type> block_type;

    //! block-identifier for subblocks
    typedef typename subblock_type::bid_type subblock_bid_type;
    //! block-identifier for blocks
    typedef typename block_type::bid_type bid_type;
    //! container for block-bids
    typedef std::vector<bid_type> bid_container_type;
    //! iterator for block-bids
    typedef typename bid_container_type::iterator bid_iterator_type;

    enum source_type { src_internal, src_external, src_unknown };

    //! nodes for internal-memory buffer
    typedef node<value_type> node_type;
    //! buckets
    typedef bucket<node_type> bucket_type;

    typedef std::vector<bucket_type> buckets_container_type;

    //! for tracking active iterators
    typedef iterator_map<self_type> iterator_map_type;

    typedef block_cache<block_type> block_cache_type;

    typedef buffered_reader<block_cache_type, bid_iterator_type> reader_type;

    typedef typename allocator_type::template rebind<node_type>::other node_allocator_type;

protected:
    //! user supplied mother hash-function
    hasher hash_;
    //! user supplied strict-weak-ordering for keys
    key_compare cmp_;
    //! array of bucket
    buckets_container_type buckets_;
    //! blocks-ids of allocated blocks
    bid_container_type bids_;
    //! size of internal-memory buffer in number of entries
    internal_size_type buffer_size_;
    //! maximum size for internal-memory buffer
    internal_size_type max_buffer_size_;
    //! keeps track of all active iterators
    iterator_map_type iterator_map_;

    mutable block_cache_type block_cache_;
    //! used to allocate new nodes for internal buffer
    node_allocator_type node_allocator_;
    //! false if the total-number of values is correct (false) or true if
    //! estimated (true); see *oblivious_-methods
    mutable bool oblivious_;
    //! (estimated) number of values
    mutable external_size_type num_total_;
    //! desired load factor after rehashing
    float opt_load_factor_;

public:
    /*!
     * Construct a new hash-map
     * \param n initial number of buckets
     * \param hf hash-function
     * \param cmp comparator-object
     * \param buffer_size size of internal-memory buffer in bytes
     * \param a allocation-strategory for internal-memory buffer
     */
    hash_map(internal_size_type n = 0,
             const hasher& hf = hasher(),
             const key_compare& cmp = key_compare(),
             internal_size_type buffer_size = 128*1024*1024,
             const allocator_type& a = allocator_type())
        : hash_(hf),
          cmp_(cmp),
          buckets_(n),
          bids_(0),
          buffer_size_(0),
          iterator_map_(this),
          block_cache_(tuning::get_instance()->blockcache_size),
          node_allocator_(a),
          oblivious_(false),
          num_total_(0),
          opt_load_factor_(0.875)
    {
        max_buffer_size_ = buffer_size / sizeof(node_type);
    }

    /*!
     * Construct a new hash-map and insert all values in the range [f,l)
     *
     * \param begin beginning of the range
     * \param end end of the range
     * \param mem_to_sort internal memory that may be used for bulk-construction (not
     * to be confused with the buffer-memory)
     * \param n initial number of buckets
     * \param hf hash-function
     * \param cmp comparator-object
     * \param buffer_size size of internal-memory buffer in bytes
     * \param a allocation-strategory for internal-memory buffer
     */
    template <class InputIterator>
    hash_map(InputIterator begin, InputIterator end,
             internal_size_type mem_to_sort = 256*1024*1024,
             internal_size_type n = 0,
             const hasher& hf = hasher(),
             const key_compare& cmp = key_compare(),
             internal_size_type buffer_size = 128*1024*1024,
             const allocator_type& a = allocator_type())
        : hash_(hf),
          cmp_(cmp),
          buckets_(n),                 // insert will determine a good size
          bids_(0),
          buffer_size_(0),
          iterator_map_(this),
          block_cache_(tuning::get_instance()->blockcache_size),
          node_allocator_(a),
          oblivious_(false),
          num_total_(0),
          opt_load_factor_(0.875)
    {
        max_buffer_size_ = buffer_size / sizeof(node_type);
        insert(begin, end, mem_to_sort);
    }

    ~hash_map()
    {
        clear();
    }

public:
    //! Hash-function used by this hash-map
    hasher hash_function() const
    { return hash_; }

    //! Strict-weak-ordering used by this hash-map
    key_compare key_cmp() const
    { return cmp_; }

    //! Get node memory allocator
    allocator_type get_allocator() const
    { return node_allocator_; }

protected:
    /*!
     * After using *oblivious_-methods only an estimate for the total number of
     * elements can be given.  This method accesses external memory to
     * calculate the exact number.
     */
    void _make_conscious()
    { /* const */                           //! TODO: make const again
        if (!oblivious_)
            return;

        typedef HashedValuesStream<self_type, reader_type> values_stream_type;

        // this will start prefetching automatically
        reader_type reader(bids_.begin(), bids_.end(), block_cache_);
        values_stream_type values(buckets_.begin(), buckets_.end(),
                                  reader, bids_.begin(), *this);

        num_total_ = 0;
        while (!values.empty())
        {
            ++num_total_;
            ++values;
        }
        oblivious_ = false;
    }

public:
    //! Number of values currently stored. Note: If the correct number is
    //! currently unknown (because *_oblivous-methods were used), external
    //! memory will be scanned.
    external_size_type size() const
    {
        if (oblivious_)
            ((self_type*)this)->_make_conscious();
        return num_total_;
    }

    //! The hash-map may store up to this number of values
    external_size_type max_size() const
    {
        return std::numeric_limits<external_size_type>::max();
    }

    //! Check if container is empty.
    bool empty() const
    {
        return size() == 0;
    }

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
        if (buckets_.size() == 0)
            _rebuild_buckets(128);

        internal_size_type i_bucket = _bkt_num(value.first);
        bucket_type& bucket = buckets_[i_bucket];
        node_type* node = _find_key_internal(bucket, value.first);

        // found value in internal memory
        if (node && _eq(node->value_.first, value.first))
        {
            bool old_deleted = node->deleted();
            if (old_deleted)
            {
                node->set_deleted(false);
                node->value_ = value;
                ++num_total_;
            }
            return std::pair<iterator, bool>(
                iterator(this, i_bucket, node,
                         0, src_internal, false, value.first), old_deleted);
        }
        // search external memory ...
        else
        {
            tuple<external_size_type, value_type> result
                = _find_key_external(bucket, value.first);

            external_size_type i_external = result.first;
            value_type ext_value = result.second;

            // ... if found, return iterator pointing to external position ...
            if (i_external < bucket.n_external_ && _eq(ext_value.first, value.first))
            {
                return std::pair<iterator, bool>(
                    iterator(this, i_bucket, node,
                             i_external, src_external, true, value.first), false);
            }
            // ... otherwise create a new buffer-node to add the value
            else
            {
                ++num_total_;
                node_type* new_node =
                    node
                    ? node->set_next(_new_node(value, node->next(), false))
                    : (bucket.list_ = _new_node(value, bucket.list_, false));

                iterator it(this, i_bucket, new_node,
                            0, src_internal, false, value.first);

                ++buffer_size_;
                if (buffer_size_ >= max_buffer_size_)
                    _rebuild_buckets();                 // will fix it as well

                return std::pair<iterator, bool>(it, true);
            }
        }
    }

    //! Insert a value; external memory is not accessed so that another value
    //! with the same key may be overwritten
    //! \param value what to insert
    //! \return iterator pointing to the inserted value
    iterator insert_oblivious(const value_type& value)
    {
        internal_size_type i_bucket = _bkt_num(value.first);
        bucket_type& bucket = buckets_[i_bucket];
        node_type* node = _find_key_internal(bucket, value.first);

        // found value in internal memory
        if (node && _eq(node->value_.first, value.first))
        {
            if (node->deleted())
                ++num_total_;

            node->set_deleted(false);
            node->value_ = value;
            return iterator(this, i_bucket, node,
                            0, src_internal, false, value.first);
        }
        // not found; ignore external memory and add a new node to the
        // internal-memory buffer
        else
        {
            oblivious_ = true;
            ++num_total_;
            node_type* new_node =
                node
                ? node->set_next(_new_node(value, node->next(), false))
                : (bucket.list_ = _new_node(value, bucket.list_, false));

            // there may be some iterators that reference the newly inserted
            // value in external memory these need to be fixed (make them point
            // to new_node)
            iterator_map_.fix_iterators_2int(i_bucket, value.first, new_node);

            iterator it(this, i_bucket, new_node,
                        0, src_internal, false, value.first);

            ++buffer_size_;
            if (buffer_size_ >= max_buffer_size_)
                _rebuild_buckets();

            return it;
        }
    }

    //! Erase value by iterator
    //! \param it iterator pointing to the value to erase
    void erase(const_iterator it)
    {
        --num_total_;
        bucket_type& bucket = buckets_[it.i_bucket_];

        if (it.source_ == src_internal)
        {
            it.node_->set_deleted(true);
            iterator_map_.fix_iterators_2end(it.i_bucket_, it.key_);
        }
        else {
            // find biggest value < iterator's value
            node_type* node = _find_key_internal(bucket, it.key_);
            assert(!node || !_eq(node->value_.first, it.key_));

            // add delete-node to buffer
            if (node)
                node->set_next(_new_node(value_type(it.key_, mapped_type()), node->next(), true));
            else
                bucket.list_ = _new_node(value_type(it.key_, mapped_type()), bucket.list_, true);

            iterator_map_.fix_iterators_2end(it.i_bucket_, it.key_);

            ++buffer_size_;
            if (buffer_size_ >= max_buffer_size_)
                _rebuild_buckets();
        }
    }

    //! Erase value by key; check external memory
    //! \param key key of value to erase
    //! \return number of values actually erased (0 or 1)
    external_size_type erase(const key_type& key)
    {
        internal_size_type i_bucket = _bkt_num(key);
        bucket_type& bucket = buckets_[i_bucket];
        node_type* node = _find_key_internal(bucket, key);

        // found in internal memory
        if (node && _eq(node->value_.first, key))
        {
            if (!node->deleted())
            {
                node->set_deleted(true);
                --num_total_;
                iterator_map_.fix_iterators_2end(i_bucket, key);
                return 1;
            }
            else
                return 0;               // already deleted
        }
        // check external memory
        else
        {
            tuple<external_size_type, value_type> result
                = _find_key_external(bucket, key);

            external_size_type i_external = result.first;
            value_type ext_value = result.second;

            // found in external memory; add delete-node
            if (i_external < bucket.n_external_ && _eq(ext_value.first, key))
            {
                --num_total_;

                if (node)
                    node->set_next(_new_node(value_type(key, mapped_type()), node->next(), true));
                else
                    bucket.list_ = _new_node(value_type(key, mapped_type()), bucket.list_, true);

                iterator_map_.fix_iterators_2end(i_bucket, key);

                ++buffer_size_;
                if (buffer_size_ >= max_buffer_size_)
                    _rebuild_buckets();

                return 1;
            }
            // no value with given key
            else
                return 0;
        }
    }

    //! Erase value by key but without looking at external memory
    //! \param key key for value to release
    void erase_oblivious(const key_type& key)
    {
        internal_size_type i_bucket = _bkt_num(key);
        bucket_type& bucket = buckets_[i_bucket];
        node_type* node = _find_key_internal(bucket, key);

        // found value in internal-memory
        if (node && _eq(node->value_.first, key))
        {
            if (!node->deleted())
            {
                --num_total_;
                node->set_deleted(true);
                iterator_map_.fix_iterators_2end(i_bucket, key);
            }
        }
        // not found; ignore external memory and add delete-node
        else
        {
            oblivious_ = true;
            --num_total_;

            if (node)
                node->set_next(_new_node(value_type(key, mapped_type()), node->next(), true));
            else
                bucket.list_ = _new_node(value_type(key, mapped_type()), bucket.list_, true);

            iterator_map_.fix_iterators_2end(i_bucket, key);

            ++buffer_size_;
            if (buffer_size_ >= max_buffer_size_)
                _rebuild_buckets();
        }
    }

    //! Reset hash-map: erase all values, invalidate all iterators
    void clear()
    {
        STXXL_VERBOSE_HASH_MAP("clear()");

        iterator_map_.fix_iterators_all2end();
        block_cache_.flush();
        block_cache_.clear();

        // reset buckets and release buffer-memory
        for (internal_size_type i_bucket = 0;
             i_bucket < buckets_.size(); i_bucket++)
        {
            _erase_nodes(buckets_[i_bucket].list_, NULL);
            buckets_[i_bucket] = bucket_type();
        }
        oblivious_ = false;
        num_total_ = 0;
        buffer_size_ = 0;

        // free external memory
        block_manager* bm = block_manager::get_instance();
        bm->delete_blocks(bids_.begin(), bids_.end());
        bids_.clear();
    }

    //! Exchange stored values with another hash-map
    //! \param obj hash-map to swap values with
    void swap(self_type& obj)
    {
        std::swap(buckets_, obj.buckets_);
        std::swap(bids_, obj.bids_);

        std::swap(oblivious_, obj.oblivious_);
        std::swap(num_total_, obj.num_total_);

        std::swap(node_allocator_, obj.node_allocator_);

        std::swap(hash_, obj.hash_);
        std::swap(cmp_, obj.cmp_);

        std::swap(buffer_size_, obj.buffer_size_);
        std::swap(max_buffer_size_, obj.max_buffer_size_);

        std::swap(opt_load_factor_, obj.opt_load_factor_);

        std::swap(iterator_map_, obj.iterator_map_);

        std::swap(block_cache_, obj.block_cache_);
    }

protected:
    // find statistics
    mutable external_size_type n_subblocks_loaded;
    mutable external_size_type n_found_internal;
    mutable external_size_type n_found_external;
    mutable external_size_type n_not_found;

public:
    //! Reset hash-map statistics
    void reset_statistics()
    {
        block_cache_.reset_statistics();
        n_subblocks_loaded = n_found_external = n_found_internal = n_not_found = 0;
    }

    //! Print short general statistics to output stream
    void print_statistics(std::ostream& o = std::cout) const
    {
        o << "Find-statistics:" << std::endl;
        o << "  Found internal     : " << n_found_internal << std::endl;
        o << "  Found external     : " << n_found_external << std::endl;
        o << "  Not found          : " << n_not_found << std::endl;
        o << "  Subblocks searched : " << n_subblocks_loaded << std::endl;

        iterator_map_.print_statistics(o);
        block_cache_.print_statistics(o);
    }

    //! Look up value by key. Non-const access.
    //! \param key key for value to look up
    iterator find(const key_type& key)
    {
        if (buffer_size_ + 1 >= max_buffer_size_)       // (*)
            _rebuild_buckets();

        internal_size_type i_bucket = _bkt_num(key);
        bucket_type& bucket = buckets_[i_bucket];
        node_type* node = _find_key_internal(bucket, key);

        // found in internal-memory buffer
        if (node && _eq(node->value_.first, key)) {
            n_found_internal++;
            if (node->deleted())
                return this->_end<iterator>();
            else
                return iterator(this, i_bucket, node, 0, src_internal, false, key);
        }
        // search external elements
        else {
            tuple<external_size_type, value_type> result
                = _find_key_external(bucket, key);

            external_size_type i_external = result.first;
            value_type value = result.second;

            // found in external memory
            if (i_external < bucket.n_external_ && _eq(value.first, key)) {
                n_found_external++;

                // we ultimately expect the user to de-reference the returned
                // iterator to change its value (non-const!).  to prevent an
                // additional disk-access, we create a new node in the
                // internal-memory buffer overwriting the external value.
                // note: by checking and rebuilding (if neccessary) in (*) we
                // made sure that the new node will fit into the buffer and no
                // rebuild is neccessary here.
                node_type* new_node =
                    node
                    ? node->set_next(_new_node(value, node->next(), false))
                    : (bucket.list_ = _new_node(value, bucket.list_, false));

                ++buffer_size_;

                iterator_map_.fix_iterators_2int(i_bucket, value.first, new_node);

                return iterator(this, i_bucket, new_node, i_external + 1, src_internal, true, key);
            }
            // not found in external memory
            else {
                n_not_found++;
                return this->_end<iterator>();
            }
        }
    }

    //! Look up value by key. Const access.
    //! \param key key for value to look up
    const_iterator find(const key_type& key) const
    {
        internal_size_type i_bucket = _bkt_num(key);
        const bucket_type& bucket = buckets_[i_bucket];
        node_type* node = _find_key_internal(bucket, key);

        // found in internal-memory buffer
        if (node && _eq(node->value_.first, key)) {
            n_found_internal++;
            if (node->deleted())
                return this->_end<const_iterator>();
            else
                return const_iterator((self_type*)this, i_bucket, node, 0, src_internal, false, key);
        }
        // search external elements
        else {
            tuple<external_size_type, value_type> result
                = _find_key_external(bucket, key);

            external_size_type i_external = result.first;
            value_type value = result.second;

            // found in external memory
            if (i_external < bucket.n_external_ && _eq(value.first, key)) {
                n_found_external++;
                return const_iterator((self_type*)this, i_bucket, node, i_external, src_external, true, key);
            }
            // not found in external memory
            else {
                n_not_found++;
                return this->_end<const_iterator>();
            }
        }
    }

    //! Number of values with given key
    //! \param k key for value to look up
    //! \return 0 or 1 depending on the presence of a value with the given key
    external_size_type count(const key_type& k) const
    {
        const_iterator cit = find(k);
        return (cit == end()) ? 0 : 1;
    }

    //! Finds a range containing all values with given key. Non-const access
    //! \param key key to look for#
    //! \return range may be empty or contains exactly one value
    std::pair<iterator, iterator> equal_range(const key_type& key)
    {
        iterator it = find(key);
        return std::pair<iterator, iterator>(it, it);
    }

    //! Finds a range containing all values with given key. Const access
    //! \param key key to look for#
    //! \return range may be empty or contains exactly one value
    std::pair<const_iterator, const_iterator> equal_range(const key_type& key) const
    {
        const_iterator cit = find(key);
        return std::pair<const_iterator, const_iterator>(cit, cit);
    }

    //! Convenience operator to quickly insert or find values. Use with caution
    //! since using this operator will check external-memory.
    mapped_type& operator [] (const key_type& key)
    {
        if (buffer_size_ + 1 >= max_buffer_size_)       // (*)
            _rebuild_buckets();

        internal_size_type i_bucket = _bkt_num(key);
        bucket_type& bucket = buckets_[i_bucket];
        node_type* node = _find_key_internal(bucket, key);

        // found in internal-memory buffer
        if (node && _eq(node->value_.first, key)) {
            if (node->deleted()) {
                node->set_deleted(false);
                node->value_.second = mapped_type();
                ++num_total_;
            }
            return node->value_.second;
        }
        // search external elements
        else {
            tuple<external_size_type, value_type> result
                = _find_key_external(bucket, key);

            external_size_type i_external = result.first;
            value_type found_value = result.second;

            value_type buffer_value =
                (i_external < bucket.n_external_ && _eq(found_value.first, key))
                ? found_value
                : value_type(key, mapped_type());

            // add a new node to the buffer. this new node's value overwrites
            // the external value if it was found and otherwise is set to (key,
            // mapped_type())
            node_type* new_node =
                node
                ? node->set_next(_new_node(buffer_value, node->next(), false))
                : (bucket.list_ = _new_node(buffer_value, bucket.list_, false));

            ++buffer_size_;
            // note that we already checked the buffer-size in (*)

            return new_node->value_.second;
        }
    }

    //! Number of buckets
    internal_size_type bucket_count() const
    { return buckets_.size(); }

    //! Maximum number of buckets
    internal_size_type max_bucket_count() const
    { return (internal_size_type)(max_size() / subblock_size); }

    //! Bucket-index for values with given key.
    internal_size_type bucket_index(const key_type& k) const
    { return _bkt_num(k); }

public:
    //! Average number of (sub)blocks occupied by a bucket.
    float load_factor() const
    { return (float)num_total_ / ((float)subblock_size * (float)buckets_.size()); }

    //! Get desired load-factor
    float opt_load_factor() const { return opt_load_factor_; }

    //! Set desired load-factor
    void opt_load_factor(float z)
    {
        opt_load_factor_ = z;
        if (load_factor() > opt_load_factor_)
            _rebuild_buckets();
    }

    //! Rehash with (at least) n buckets
    void rehash(internal_size_type n = 0)
    {
        _rebuild_buckets(n);
    }

    //! Number of bytes occupied by buffer
    internal_size_type buffer_size() const
    {
        // buffer-size internally stored as number of nodes
        return buffer_size_ * sizeof(node_type);
    }

    //! Maximum buffer size in byte
    internal_size_type max_buffer_size() const
    {
        return max_buffer_size_ * sizeof(node_type);
    }

    //! Set maximum buffer size
    //! \param buffer_size new size in byte
    void max_buffer_size(internal_size_type buffer_size)
    {
        max_buffer_size_ = buffer_size / sizeof(node_type);
        if (buffer_size_ >= max_buffer_size_)
            _rebuild_buckets();
    }

protected:
    //! iterator pointing to the beginnning of the hash-map
    template <class Iterator>
    Iterator _begin() const
    {
        self_type* non_const_this = (self_type*)this;

        if (buckets_.size() == 0)
            return _end<Iterator>();

        // correct key will be set by find_next()
        Iterator it(non_const_this, 0, buckets_[0].list_,
                    0, src_unknown, true, key_type());
        it.find_next();

        return it;
    }

    //! iterator pointing to the end of the hash-map (iterator-type as
    //! template-parameter)
    template <class Iterator>
    Iterator _end() const
    {
        self_type* non_const_this = (self_type*)this;
        return Iterator(non_const_this);
    }

public:
    //! Returns an iterator pointing to the beginning of the hash-map
    iterator begin() { return _begin<iterator>(); }

    //! Returns a const_interator pointing to the beginning of the hash-map
    const_iterator begin() const { return _begin<const_iterator>(); }

    //! Returns an iterator pointing to the end of the hash-map
    iterator end() { return _end<iterator>(); }

    //! Returns a const_iterator pointing to the end of the hash-map
    const_iterator end() const { return _end<const_iterator>(); }

protected:
    //! Allocate a new buffer-node
    node_type * _get_node()
    {
        return node_allocator_.allocate(1);
    }

    //! Free given node
    void _put_node(node_type* node)
    {
        node_allocator_.deallocate(node, 1);
    }

    //! Allocate a new buffer-node and initialize with given value, node and
    //! deleted-flag
    node_type * _new_node(const value_type& value, node_type* nxt, bool del)
    {
        node_type* node = _get_node();
        node->value_ = value;
        node->set_next(nxt);
        node->set_deleted(del);
        return node;
    }

    //! Free nodes in range [first, last). If last is NULL all nodes will be
    //! freed.
    void _erase_nodes(node_type* first, node_type* last)
    {
        node_type* curr = first;
        while (curr != last)
        {
            node_type* next = curr->next();
            _put_node(curr);
            curr = next;
        }
    }

    //! Bucket-index for values with given key
    internal_size_type _bkt_num(const key_type& key) const
    {
        return _bkt_num(key, buckets_.size());
    }

    /*!
     * Bucket-index for values with given key. The total number of buckets has
     * to be specified as well.  The bucket number is determined by \f$
     * bucket_num = (hash/max_hash)*n_buckets \f$ max_hash is in fact 2^63-1
     * (internal_size_type=uint64 (or uint32)) but we rather divide by 2^64, so
     * we can use plain integer arithmetic easily (there should be only a small
     * difference): this way we must only calculate the upper 64 bits of the
     * product hash*n_buckets and we're done. See
     * http://www.cs.uaf.edu/~cs301/notes/Chapter5/node5.html
    */
    internal_size_type _bkt_num(const key_type& key, internal_size_type n) const
    {
        //! TODO maybe specialize double arithmetic to integer. the old code
        //! was faulty -tb.
        return (internal_size_type)(
            (double)n * ((double)hash_(key) / (double)std::numeric_limits<internal_size_type>::max())
            );
    }

    /*!
     * Locate the given key in the internal-memory chained list.  If the key is
     * not present, the node with the biggest key smaller than the given key is
     * returned.  Note that the returned value may be zero: either because the
     * chained list is empty or because the given key is smaller than all other
     * keys in the chained list.
     */
    node_type*
    _find_key_internal(const bucket_type& bucket, const key_type& key) const
    {
        node_type* old = NULL;
        for (node_type* curr = bucket.list_;
             curr && _leq(curr->value_.first, key);
             curr = curr->next())
        {
            old = curr;
        }
        return old;
    }

    /*!
     * Search for key in external part of bucket. Return value is (i_external,
     * value), where i_ext = bucket._num_external if key could not be found.
     */
    tuple<external_size_type, value_type>
    _find_key_external(const bucket_type& bucket, const key_type& key) const
    {
        subblock_type* subblock;

        // number of subblocks occupied by bucket
        internal_size_type n_subblocks = (internal_size_type)(
            bucket.n_external_ / subblock_size
            );
        if (bucket.n_external_ % subblock_size != 0)
            n_subblocks++;

        for (internal_size_type i_subblock = 0;
             i_subblock < n_subblocks; i_subblock++)
        {
            subblock = _load_subblock(bucket, i_subblock);
            // number of values in i-th subblock
            internal_size_type n_values =
                (i_subblock + 1 < n_subblocks)
                ? (internal_size_type)subblock_size
                : (internal_size_type)(
                    bucket.n_external_ - i_subblock * subblock_size
                    );

            //! TODO: replace with bucket.n_external_ % subblock_size

            // biggest key in current subblock still too small => next subblock
            if (_lt((*subblock)[n_values - 1].first, key))
                continue;

            // binary search in current subblock
            internal_size_type i_lower = 0, i_upper = n_values;
            while (i_lower + 1 != i_upper)
            {
                internal_size_type i_middle = (i_lower + i_upper) / 2;
                if (_leq((*subblock)[i_middle].first, key))
                    i_lower = i_middle;
                else
                    i_upper = i_middle;
            }

            value_type value = (*subblock)[i_lower];

            if (_eq(value.first, key))
                return tuple<external_size_type, value_type>
                           (i_subblock * subblock_size + i_lower, value);
            else
                return tuple<external_size_type, value_type>
                           (bucket.n_external_, value_type());
        }

        return tuple<external_size_type, value_type>
                   (bucket.n_external_, value_type());
    }

    /*!
     * Load the given bucket's i-th subblock.
     * Since a bucket may be spread over several blocks, we must
     * 1. determine in which block the requested subblock is located
     * 2. at which position within the obove-mentioned block the questioned subblock is located
     */
    subblock_type*
    _load_subblock(const bucket_type& bucket, internal_size_type which_subblock) const
    {
        n_subblocks_loaded++;

        // index of the requested subblock counted from the very beginning of
        // the bucket's first block
        external_size_type i_abs_subblock = bucket.i_subblock_ + which_subblock;

        /* 1. */
        bid_type bid = bids_[bucket.i_block_ + (internal_size_type)(i_abs_subblock / subblocks_per_block)];
        /* 2. */
        internal_size_type i_subblock_within = (internal_size_type)(i_abs_subblock % subblocks_per_block);

        return block_cache_.get_subblock(bid, i_subblock_within);
    }

    typedef HashedValue<self_type> hashed_value_type;

    //! Functor to extracts the actual value from a HashedValue-struct
    struct HashedValueExtractor
    {
        value_type& operator () (hashed_value_type& hvalue)
        { return hvalue.value_; }
    };

    /*!
     * Will return from its input-stream all values that are to be stored in
     * the given bucket.  Those values must appear in consecutive order
     * beginning with the input-stream's current value.
     */
    template <class InputStream, class ValueExtractor>
    struct HashingStream
    {
        typedef typename InputStream::value_type value_type;

        self_type* map_;
        InputStream& input_;
        internal_size_type i_bucket_;
        external_size_type bucket_size_;
        value_type value_;
        bool empty_;
        ValueExtractor vextract_;

        HashingStream(InputStream& input, internal_size_type i_bucket,
                      ValueExtractor vextract, self_type* map)
            : map_(map),
              input_(input),
              i_bucket_(i_bucket),
              bucket_size_(0),
              vextract_(vextract)
        {
            empty_ = find_next();
        }

        const value_type& operator * () { return value_; }

        bool empty() const { return empty_; }

        void operator ++ ()
        {
            ++input_;
            empty_ = find_next();
        }

        bool find_next()
        {
            if (input_.empty())
                return true;
            value_ = *input_;
            if (map_->_bkt_num(vextract_(value_).first) != i_bucket_)
                return true;

            ++bucket_size_;
            return false;
        }
    };

    /*  Rebuild hash-map. The desired number of buckets may be supplied. */
    void _rebuild_buckets(internal_size_type n_desired = 0)
    {
        STXXL_VERBOSE_HASH_MAP("_rebuild_buckets()");

        typedef buffered_writer<block_type, bid_container_type> writer_type;
        typedef HashedValuesStream<self_type, reader_type> values_stream_type;
        typedef HashingStream<values_stream_type, HashedValueExtractor> hashing_stream_type;

        const int_type write_buffer_size = config::get_instance()->disks_number() * 4;

        // determine new number of buckets from desired load_factor ...
        internal_size_type n_new;
        n_new = (internal_size_type)ceil((double)num_total_ / ((double)subblock_size * (double)opt_load_factor()));

        // ... but give the user the chance to request even more buckets
        if (n_desired > n_new)
            n_new = std::min<internal_size_type>(n_desired, max_bucket_count());

        // allocate new buckets and bids
        buckets_container_type old_buckets(n_new);
        std::swap(buckets_, old_buckets);

        bid_container_type old_bids;
        std::swap(bids_, old_bids);

        // read stored values in consecutive order

        // use new to control point of destruction (see below)
        reader_type* reader
            = new reader_type(old_bids.begin(), old_bids.end(), block_cache_);

        values_stream_type values_stream(old_buckets.begin(), old_buckets.end(),
                                         *reader, old_bids.begin(), *this);

        writer_type writer(&bids_, write_buffer_size, write_buffer_size / 2);

        // re-distribute values among new buckets.

        // this makes use of the fact that if value1 preceeds value2 before
        // resizing, value1 will preceed value2 after resizing as well (uniform
        // rehashing)
        num_total_ = 0;
        for (internal_size_type i_bucket = 0;
             i_bucket < buckets_.size(); i_bucket++)
        {
            buckets_[i_bucket] = bucket_type();
            buckets_[i_bucket].i_block_ = writer.i_block();
            buckets_[i_bucket].i_subblock_ = writer.i_subblock();

            hashing_stream_type hasher(values_stream, i_bucket, HashedValueExtractor(), this);
            external_size_type i_ext = 0;
            while (!hasher.empty())
            {
                const hashed_value_type& hvalue = *hasher;
                iterator_map_.fix_iterators_2ext(hvalue.i_bucket_, hvalue.value_.first, i_bucket, i_ext);

                writer.append(hvalue.value_);
                ++hasher;
                ++i_ext;
            }

            writer.finish_subblock();
            buckets_[i_bucket].n_external_ = hasher.bucket_size_;
            num_total_ += hasher.bucket_size_;
        }
        writer.flush();
        // reader must be deleted before deleting old_bids because its
        // destructor will dereference the bid-iterator
        delete reader;
        block_cache_.clear();

        // get rid of old blocks and buckets
        block_manager* bm = stxxl::block_manager::get_instance();
        bm->delete_blocks(old_bids.begin(), old_bids.end());

        for (internal_size_type i_bucket = 0;
             i_bucket < old_buckets.size(); i_bucket++)
        {
            _erase_nodes(old_buckets[i_bucket].list_, NULL);
            old_buckets[i_bucket] = bucket_type();
        }

        buffer_size_ = 0;
        oblivious_ = false;
    }

    /*!
     * Stream for filtering duplicates. Used to eliminate duplicated values
     * when bulk-inserting Note: input has to be sorted, so that duplicates
     * will occure in row
     */
    template <class InputStream>
    struct UniqueValueStream
    {
        typedef typename InputStream::value_type value_type;
        self_type& map_;
        InputStream& in_;

        UniqueValueStream(InputStream& input, self_type& map)
            : map_(map), in_(input)
        { }

        bool empty() const { return in_.empty(); }

        const value_type& operator * () { return *in_; }

        void operator ++ ()
        {
            value_type v_old = *in_;
            ++in_;
            while (!in_.empty() && v_old.first == (*in_).first)
                ++in_;
        }
    };

    template <class InputStream>
    struct AddHashStream
    {
        //! (hash,value)
        typedef std::pair<internal_size_type, typename InputStream::value_type> value_type;
        self_type& map_;
        InputStream& in_;

        AddHashStream(InputStream& input, self_type& map)
            : map_(map), in_(input)
        { }

        bool empty() const { return in_.empty(); }

        value_type operator * ()
        { return value_type(map_.hash_((*in_).first), *in_); }

        void operator ++ () { ++in_; }
    };

    /*!
     * Extracts the value-part (ignoring the hashvalue); required by
     * HashingStream (see above)
     */
    struct StripHashFunctor
    {
        const value_type& operator () (std::pair<internal_size_type, value_type>& v)
        { return v.second; }
    };

    /*!
     * Comparator object for values as required by stxxl::sort. Sorting is done
     * lexicographically by <hash-value, key> Note: the hash-value has already
     * been computed.
     */
    struct Cmp : public std::binary_function<
                     std::pair<internal_size_type, value_type>,
                     std::pair<internal_size_type, value_type>, bool
                     >
    {
        self_type& map_;
        Cmp(self_type& map) : map_(map) { }

        bool operator () (const std::pair<internal_size_type, value_type>& a,
                          const std::pair<internal_size_type, value_type>& b) const
        {
            return (a.first < b.first) ||
                   ((a.first == b.first) && map_.cmp_(a.second.first, b.second.first));
        }
        std::pair<internal_size_type, value_type> min_value() const
        {
            return std::pair<internal_size_type, value_type>(
                std::numeric_limits<internal_size_type>::min(),
                value_type(map_.cmp_.min_value(), mapped_type())
                );
        }
        std::pair<internal_size_type, value_type> max_value() const
        {
            return std::pair<internal_size_type, value_type>(
                std::numeric_limits<internal_size_type>::max(),
                value_type(map_.cmp_.max_value(), mapped_type())
                );
        }
    };

public:
    //! Bulk-insert of values in the range [f, l)
    //! \param f beginning of the range
    //! \param l end of the range
    //! \param mem internal memory that may be used (note: this memory will be used additionally to the buffer). The more the better
    template <class InputIterator>
    void insert(InputIterator f, InputIterator l, internal_size_type mem)
    {
        //! values already stored in the hashtable ("old values")
        typedef HashedValuesStream<self_type, reader_type> old_values_stream;
        //! old values, that are to be stored in a certain (new) bucket
        typedef HashingStream<old_values_stream, HashedValueExtractor> old_hashing_stream;

        //! values to insert ("new values")
        typedef typename stxxl::stream::streamify_traits<InputIterator>::stream_type input_stream;

        //! new values with added hash: (hash, (key, mapped))
        typedef AddHashStream<input_stream> new_values_stream;
        //! new values sorted by <hash-value, key>
        typedef stxxl::stream::sort<new_values_stream, Cmp> new_sorted_values_stream;
        //! new values sorted by <hash-value, key> with duplicates eliminated
        typedef UniqueValueStream<new_sorted_values_stream> new_unique_values_stream;
        //! new values, that are to be stored in a certain bucket
        typedef HashingStream<new_unique_values_stream, StripHashFunctor> new_hashing_stream;

        typedef buffered_writer<block_type, bid_container_type> writer_type;

        int_type write_buffer_size = config::get_instance()->disks_number() * 2;

        // calculate new number of buckets
        external_size_type num_total_new = num_total_ + (l - f);         // estimated number of elements
        external_size_type n_buckets_new = (external_size_type)ceil((double)num_total_new / ((double)subblock_size * (double)opt_load_factor()));
        if (n_buckets_new > max_bucket_count())
            n_buckets_new = max_bucket_count();

        STXXL_VERBOSE_HASH_MAP("insert() items=" << (l - f) << " buckets_new=" << n_buckets_new);

        // prepare new buckets and bids
        buckets_container_type old_buckets((internal_size_type)n_buckets_new);
        std::swap(buckets_, old_buckets);
        // writer will allocate new blocks as necessary
        bid_container_type old_bids;
        std::swap(bids_, old_bids);

        // already stored values ("old values")
        reader_type* reader = new reader_type(old_bids.begin(), old_bids.end(),
                                              block_cache_);
        old_values_stream old_values(old_buckets.begin(), old_buckets.end(),
                                     *reader, old_bids.begin(), *this);

        // values to insert ("new values")
        input_stream input = stxxl::stream::streamify(f, l);
        new_values_stream new_values(input, *this);
        new_sorted_values_stream new_sorted_values(new_values, Cmp(*this), mem);
        new_unique_values_stream new_unique_values(new_sorted_values, *this);

        writer_type writer(&bids_, write_buffer_size, write_buffer_size / 2);

        num_total_ = 0;
        for (internal_size_type i_bucket = 0; i_bucket < buckets_.size(); i_bucket++)
        {
            buckets_[i_bucket] = bucket_type();
            buckets_[i_bucket].i_block_ = writer.i_block();
            buckets_[i_bucket].i_subblock_ = writer.i_subblock();

            old_hashing_stream old_hasher(old_values, i_bucket, HashedValueExtractor(), this);
            new_hashing_stream new_hasher(new_unique_values, i_bucket, StripHashFunctor(), this);
            internal_size_type bucket_size = 0;

            // more old and new values for the current bucket => choose smallest
            while (!old_hasher.empty() && !new_hasher.empty())
            {
                internal_size_type old_hash = hash_((*old_hasher).value_.first);
                internal_size_type new_hash = (*new_hasher).first;
                key_type old_key = (*old_hasher).value_.first;
                key_type new_key = (*new_hasher).second.first;

                // old value wins
                if ((old_hash < new_hash) || (old_hash == new_hash && cmp_(old_key, new_key)))                // (_lt((*old_hasher)._value.first, (*new_hasher).second.first))
                {
                    const hashed_value_type& hvalue = *old_hasher;
                    iterator_map_.fix_iterators_2ext(hvalue.i_bucket_, hvalue.value_.first, i_bucket, bucket_size);
                    writer.append(hvalue.value_);
                    ++old_hasher;
                }
                // new value smaller or equal => new value wins
                else
                {
                    if (_eq(old_key, new_key))
                    {
                        const hashed_value_type& hvalue = *old_hasher;
                        iterator_map_.fix_iterators_2ext(hvalue.i_bucket_, hvalue.value_.first, i_bucket, bucket_size);
                        ++old_hasher;
                    }
                    writer.append((*new_hasher).second);
                    ++new_hasher;
                }
                ++bucket_size;
            }
            // no more new values for the current bucket
            while (!old_hasher.empty())
            {
                const hashed_value_type& hvalue = *old_hasher;
                iterator_map_.fix_iterators_2ext(hvalue.i_bucket_, hvalue.value_.first, i_bucket, bucket_size);
                writer.append(hvalue.value_);
                ++old_hasher;
                ++bucket_size;
            }
            // no more old values for the current bucket
            while (!new_hasher.empty())
            {
                writer.append((*new_hasher).second);
                ++new_hasher;
                ++bucket_size;
            }

            writer.finish_subblock();
            buckets_[i_bucket].n_external_ = bucket_size;
            num_total_ += bucket_size;
        }
        writer.flush();
        delete reader;
        block_cache_.clear();

        // release old blocks
        block_manager* bm = stxxl::block_manager::get_instance();
        bm->delete_blocks(old_bids.begin(), old_bids.end());

        // free nodes in old bucket lists
        for (internal_size_type i_bucket = 0;
             i_bucket < old_buckets.size(); i_bucket++)
        {
            _erase_nodes(old_buckets[i_bucket].list_, NULL);
            old_buckets[i_bucket] = bucket_type();
        }

        buffer_size_ = 0;
        oblivious_ = false;
    }

protected:
    /* 1 iff a <  b
       The comparison is done lexicographically by (hash-value, key)
    */
    bool _lt(const key_type& a, const key_type& b) const
    {
        internal_size_type hash_a = hash_(a);
        internal_size_type hash_b = hash_(b);

        return (hash_a < hash_b) ||
               ((hash_a == hash_b) && cmp_(a, b));
    }

    //! true iff a >  b
    bool _gt(const key_type& a, const key_type& b) const { return _lt(b, a); }
    //! true iff a <= b
    bool _leq(const key_type& a, const key_type& b) const { return !_gt(a, b); }
    //! true iff a >= b
    bool _geq(const key_type& a, const key_type& b) const { return !_lt(a, b); }

    //! true iff a == b. note: it is mandatory that equal keys yield equal
    //! hash-values => hashing not neccessary for equality-testing.
    bool _eq(const key_type& a, const key_type& b) const
    { return !cmp_(a, b) && !cmp_(b, a); }

    friend class hash_map_iterator_base<self_type>;
    friend class hash_map_iterator<self_type>;
    friend class hash_map_const_iterator<self_type>;
    friend class iterator_map<self_type>;
    friend class block_cache<block_type>;
    friend struct HashedValuesStream<self_type, reader_type>;

#if 1
    void _dump_external()
    {
        reader_type reader(bids_.begin(), bids_.end(), &block_cache_);

        for (internal_size_type i_block = 0; i_block < bids_.size(); i_block++) {
            std::cout << "block " << i_block << ":\n";

            for (internal_size_type i_subblock = 0; i_subblock < subblocks_per_block; i_subblock++) {
                std::cout << "  subblock " << i_subblock << ":\n    ";

                for (external_size_type i_element = 0; i_element < subblocks_per_block; i_element++) {
                    std::cout << reader.const_value().first << ", ";
                    ++reader;
                }
                std::cout << std::endl;
            }
        }
    }

    void _dump_buckets()
    {
        reader_type reader(bids_.begin(), bids_.end(), &block_cache_);

        std::cout << "number of buckets: " << buckets_.size() << std::endl;
        for (internal_size_type i_bucket = 0; i_bucket < buckets_.size(); i_bucket++) {
            const bucket_type& bucket = buckets_[i_bucket];
            reader.skip_to(bids_.begin() + bucket.i_block_, bucket.i_subblock_);

            std::cout << "  bucket " << i_bucket << ": block=" << bucket.i_block_ << ", subblock=" << bucket.i_subblock_ << ", external=" << bucket.n_external_ << std::endl;

            node_type* node = bucket.list_;
            std::cout << "     internal_list=";
            while (node) {
                std::cout << node->value_.first << " (del=" << node->deleted() << "), ";
                node = node->next();
            }
            std::cout << std::endl;

            std::cout << "     external=";
            for (external_size_type i_element = 0; i_element < bucket.n_external_; i_element++) {
                std::cout << reader.const_value().first << ", ";
                ++reader;
            }
            std::cout << std::endl;
        }
    }

    void _dump_bucket_statistics()
    {
        std::cout << "number of buckets: " << buckets_.size() << std::endl;
        for (internal_size_type i_bucket = 0; i_bucket < buckets_.size(); i_bucket++) {
            const bucket_type& bucket = buckets_[i_bucket];
            std::cout << "  bucket " << i_bucket << ": block=" << bucket.i_block_ << ", subblock=" << bucket.i_subblock_ << ", external=" << bucket.n_external_ << ", list=" << bucket.list_ << std::endl;
        }
    }
#endif

public:
    //! Construct an equality predicate from the comparison operator
    struct equal_to : public std::binary_function<key_type, key_type, bool>
    {
        //! reference to hash_map
        const self_type& m_map;

        //! constructor requires reference to hash_map
        equal_to(const self_type& map) : m_map(map) { }

        //! return whether the arguments compare equal (x==y).
        bool operator () (const key_type& x, const key_type& y) const
        {
            return m_map._eq(x, y);
        }

        //! C++11 required type
        typedef key_type first_argument_type;
        //! C++11 required type
        typedef key_type second_argument_type;
        //! C++11 required type
        typedef bool result_type;
    };

    //! Type of constructed equality predicate
    typedef equal_to key_equal;

    //! Constructed equality predicate used by this hash-map
    key_equal key_eq() const
    {
        return equal_to(*this);
    }

public:
    //! Even more statistics: Number of buckets, number of values, buffer-size,
    //! values per bucket
    void print_load_statistics(std::ostream& o = std::cout) const
    {
        external_size_type sum_external = 0;
        external_size_type square_sum_external = 0;
        external_size_type max_external = 0;

        for (internal_size_type i_bucket = 0; i_bucket < buckets_.size(); i_bucket++)
        {
            const bucket_type& b = buckets_[i_bucket];

            sum_external += b.n_external_;
            square_sum_external += b.n_external_ * b.n_external_;
            if (b.n_external_ > max_external)
                max_external = b.n_external_;
        }

        double avg_external = (double)sum_external / (double)buckets_.size();
        double std_external = sqrt(((double)square_sum_external / (double)buckets_.size()) - (avg_external * avg_external));

        o << "Bucket count         : " << buckets_.size() << std::endl;
        o << "Values total         : " << num_total_ << std::endl;
        o << "Values buffered      : " << buffer_size_ << std::endl;
        o << "Max Buffer-Size      : " << max_buffer_size_ << std::endl;
        o << "Max external/bucket  : " << max_external << std::endl;
        o << "Avg external/bucket  : " << avg_external << std::endl;
        o << "Std external/bucket  : " << std_external << std::endl;
        o << "Load-factor          : " << load_factor() << std::endl;
        o << "Blocks allocated     : " << bids_.size() << " => " << (bids_.size() * block_type::raw_size) << " bytes" << std::endl;
        o << "Bytes per value      : " << ((double)(bids_.size() * block_type::raw_size) / (double)num_total_) << std::endl;
    }
};     /* end of class hash_map */

} // namespace hash_map

STXXL_END_NAMESPACE

namespace std {

template <class KeyType, class MappedType, class HashType, class KeyCompareType,
          unsigned SubBlockSize, unsigned SubBlocksPerBlock, class AllocType>
void swap(stxxl::hash_map::hash_map<KeyType, MappedType, HashType, KeyCompareType,
                                    SubBlockSize, SubBlocksPerBlock, AllocType>& a,
          stxxl::hash_map::hash_map<KeyType, MappedType, HashType, KeyCompareType,
                                    SubBlockSize, SubBlocksPerBlock, AllocType>& b)
{
    if (&a != &b)
        a.swap(b);
}

} // namespace std

#endif // !STXXL_CONTAINERS_HASH_MAP_HASH_MAP_HEADER
