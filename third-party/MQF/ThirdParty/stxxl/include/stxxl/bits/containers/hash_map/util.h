/***************************************************************************
 *  include/stxxl/bits/containers/hash_map/util.h
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

#ifndef STXXL_CONTAINERS_HASH_MAP_UTIL_HEADER
#define STXXL_CONTAINERS_HASH_MAP_UTIL_HEADER
#define STXXL_CONTAINERS_HASHMAP__UTIL_H

#include <stxxl/bits/mng/block_manager.h>
#include <stxxl/bits/mng/buf_writer.h>

#include <stxxl/bits/containers/hash_map/tuning.h>
#include <stxxl/bits/containers/hash_map/block_cache.h>

STXXL_BEGIN_NAMESPACE

namespace hash_map {

// For internal memory chaining: struct to compose next-pointer and delete-flag
// share the same memory: the lowest bit is occupied by the del-flag.
template <class ValueType>
struct node
{
    node<ValueType>* next_and_del_;
    ValueType value_;

    //! check if the next node is deleted.
    bool deleted()
    {
        return ((int_type)next_and_del_ & 0x01) == 1;
    }
    //! change deleted flag on the next node
    bool set_deleted(bool d)
    {
        next_and_del_ = (node<ValueType>*)(((int_type)next_and_del_ & ~0x01) | (int_type)d);
        return d;
    }

    //! return the next node, without the "next" flag.
    node<ValueType> * next()
    {
        return (node<ValueType>*)((int_type)next_and_del_ & ~0x01);
    }
    //! change the "next" value of next node pointer
    node<ValueType> * set_next(node<ValueType>* n)
    {
        next_and_del_ = (node<ValueType>*)(((int_type)next_and_del_ & 0x01) | (int_type)n);
        return n;
    }
};

template <class NodeType>
struct bucket
{
    //! entry point to the chain in internal memory
    NodeType* list_;

    //! number of elements in external memory
    external_size_type n_external_;

    //! index of first block's bid (to be used as index for hash_map's
    //! bids_-array
    internal_size_type i_block_;

    //! index of first subblock
    internal_size_type i_subblock_;

    bucket()
        : list_(NULL),
          n_external_(0),
          i_block_(0),
          i_subblock_(0)
    { }

    bucket(NodeType* list, external_size_type n_external,
           internal_size_type i_block, internal_size_type i_subblock)
        : list_(list),
          n_external_(n_external),
          i_block_(i_block),
          i_subblock_(i_subblock)
    { }
};

//! Used to scan external memory with prefetching.
template <class CacheType, class BidIterator>
class buffered_reader : private noncopyable
{
public:
    typedef CacheType cache_type;
    typedef BidIterator bid_iterator;

    typedef typename cache_type::block_type block_type;
    typedef typename block_type::value_type subblock_type;
    typedef typename subblock_type::value_type value_type;

    typedef typename bid_iterator::value_type bid_type;

    enum { block_size = block_type::size, subblock_size = subblock_type::size };

private:
    //! index within current block
    unsigned_type i_value_;
    //! points to the beginning of the block-sequence
    bid_iterator begin_bid_;
    //! points to the current block
    bid_iterator curr_bid_;
    //! points to the end of the block-sequence
    bid_iterator end_bid_;
    //! points to the next block to prefetch
    bid_iterator pref_bid_;

    //! shared block-cache
    cache_type& cache_;

    //! true if prefetching enabled
    bool prefetch_;
    //! pages, which are read at once from disk, consist of this many blocks
    unsigned_type page_size_;
    //! number of pages to prefetch
    unsigned_type prefetch_pages_;

    //! current block dirty ?
    bool dirty_;
    //! current subblock
    subblock_type* subblock_;

public:
    //! Create a new buffered reader to read the blocks in [seq_begin, seq_end)
    //! \param seq_begin First block's bid
    //! \param seq_end Last block's bid
    //! \param cache Block-cache used for prefetching
    //! \param i_subblock Start reading from this subblock
    //! \param prefetch Enable/Disable prefetching
    buffered_reader(bid_iterator seq_begin, bid_iterator seq_end,
                    cache_type& cache,
                    internal_size_type i_subblock = 0, bool prefetch = true)
        : i_value_(0),
          begin_bid_(seq_begin),
          curr_bid_(seq_begin),
          end_bid_(seq_end),
          cache_(cache),
          prefetch_(false),
          page_size_(tuning::get_instance()->prefetch_page_size),
          prefetch_pages_(tuning::get_instance()->prefetch_pages),
          dirty_(false),
          subblock_(NULL)
    {
        if (seq_begin == seq_end)
            return;

        if (prefetch)
            enable_prefetching();

        // will (amongst other things) set subblock_ and retain current block
        skip_to(seq_begin, i_subblock);
    }

    ~buffered_reader()
    {
        if (curr_bid_ != end_bid_)
            cache_.release_block(*curr_bid_);
    }

    void enable_prefetching()
    {
        if (prefetch_)
            return;

        prefetch_ = true;
        pref_bid_ = curr_bid_;
        // start prefetching page_size*prefetch_pages blocks beginning with current one
        for (unsigned_type i = 0; i < page_size_ * prefetch_pages_; i++)
        {
            if (pref_bid_ == end_bid_)
                break;

            cache_.prefetch_block(*pref_bid_);
            ++pref_bid_;
        }
    }

    //! Get const-reference to current value.
    const value_type & const_value()
    {
        return (*subblock_)[i_value_ % subblock_size];
    }

    //! Get reference to current value. The current value's block's dirty flag
    //! will be set.
    value_type & value()
    {
        if (!dirty_) {
            cache_.make_dirty(*curr_bid_);
            dirty_ = true;
        }

        return (*subblock_)[i_value_ % subblock_size];
    }

    //! Advance to the next value
    //! \return false if last value has been reached, otherwise true.
    bool operator ++ ()
    {
        if (curr_bid_ == end_bid_)
            return false;

        // same block
        if (i_value_ + 1 < block_size * subblock_size)
        {
            i_value_++;
        }
        // entered new block
        else
        {
            cache_.release_block(*curr_bid_);

            i_value_ = 0;
            dirty_ = false;
            ++curr_bid_;

            if (curr_bid_ == end_bid_)
                return false;

            cache_.retain_block(*curr_bid_);

            // if a complete page has been consumed, prefetch the next one
            if (prefetch_ && (curr_bid_ - begin_bid_) % page_size_ == 0)
            {
                for (unsigned i = 0; i < page_size_; i++)
                {
                    if (pref_bid_ == end_bid_)
                        break;
                    cache_.prefetch_block(*pref_bid_);
                    ++pref_bid_;
                }
            }
        }

        // entered new subblock
        if (i_value_ % subblock_size == 0)
        {
            subblock_ = cache_.get_subblock(*curr_bid_, i_value_ / subblock_size);
        }

        return true;
    }

    //! Skip remaining values of the current subblock.
    void next_subblock()
    {
        i_value_ = (i_value_ / subblock_size + 1) * subblock_size - 1;
        operator ++ ();         // takes care of prefetching etc
    }

    //! Continue reading at given block and subblock.
    void skip_to(bid_iterator bid, internal_size_type i_subblock)
    {
        if (curr_bid_ == end_bid_)
            return;

        if (bid != curr_bid_)
            dirty_ = false;

        cache_.release_block(*curr_bid_);

        if (bid == end_bid_)
            return;

        // skip to block
        while (curr_bid_ != bid) {
            ++curr_bid_;

            if (prefetch_ && (curr_bid_ - begin_bid_) % page_size_ == 0)
            {
                for (unsigned i = 0; i < page_size_; i++)
                {
                    if (pref_bid_ == end_bid_)
                        break;
                    cache_.prefetch_block(*pref_bid_);
                    ++pref_bid_;
                }
            }
        }
        // skip to subblock
        i_value_ = i_subblock * subblock_size;
        subblock_ = cache_.get_subblock(*curr_bid_, i_subblock);
        cache_.retain_block(*curr_bid_);
    }
};

//! Buffered writing of values. New Blocks are allocated as needed.
template <class BlockType, class BidContainer>
class buffered_writer : private noncopyable
{
public:
    typedef BlockType block_type;
    typedef BidContainer bid_container_type;

    typedef typename block_type::value_type subblock_type;
    typedef typename subblock_type::value_type value_type;

    typedef stxxl::buffered_writer<block_type> writer_type;

    enum {
        block_size = block_type::size,
        subblock_size = subblock_type::size
    };

private:
    //! buffered writer
    writer_type writer_;
    //! current buffer-block
    block_type* block_;

    //! sequence of allocated blocks (to be expanded as needed)
    bid_container_type* bids_;

    //! current block's index
    unsigned_type i_block_;
    //! current value's index in the range of [0..\#values per block[
    unsigned_type i_value_;
    //! number of blocks to allocate in a row
    unsigned_type increase_;

public:
    //! Create a new buffered writer.
    //! \param c write values to these blocks (c holds the bids)
    //! \param buffer_size Number of write-buffers to use
    //! \param batch_size bulk buffered writing
    buffered_writer(bid_container_type* c,
                    int_type buffer_size, int_type batch_size)
        : writer_(buffer_size, batch_size),
          bids_(c),
          i_block_(0),
          i_value_(0),
          increase_(config::get_instance()->disks_number() * 3)
    {
        block_ = writer_.get_free_block();
    }

    ~buffered_writer()
    {
        flush();
    }

    //! Write all values from given stream.
    template <class StreamType>
    void append_from_stream(StreamType& stream)
    {
        while (!stream.empty())
        {
            append(*stream);
            ++stream;
        }
    }

    //! Write given value.
    void append(const value_type& value)
    {
        internal_size_type i_subblock = (i_value_ / subblock_size);
        (*block_)[i_subblock][i_value_ % subblock_size] = value;

        if (i_value_ + 1 < block_size * subblock_size)
            i_value_++;
        // reached end of a block
        else
        {
            i_value_ = 0;

            // allocate new blocks if neccessary ...
            if (i_block_ == bids_->size())
            {
                bids_->resize(bids_->size() + increase_);
                block_manager* bm = stxxl::block_manager::get_instance();
                bm->new_blocks(striping(), bids_->end() - increase_, bids_->end());
            }
            // ... and write current block
            block_ = writer_.write(block_, (*bids_)[i_block_]);

            i_block_++;
        }
    }

    //! Continue writing at the beginning of the next subblock. TODO more
    //! efficient
    void finish_subblock()
    {
        i_value_ = (i_value_ / subblock_size + 1) * subblock_size - 1;
        append(value_type());           // writing and allocating blocks etc
    }

    //! Flushes not yet written blocks.
    void flush()
    {
        i_value_ = 0;
        if (i_block_ == bids_->size())
        {
            bids_->resize(bids_->size() + increase_);
            block_manager* bm = stxxl::block_manager::get_instance();
            bm->new_blocks(striping(), bids_->end() - increase_, bids_->end());
        }
        block_ = writer_.write(block_, (*bids_)[i_block_]);
        i_block_++;

        writer_.flush();
    }

    //! Index of current block.
    internal_size_type i_block() { return i_block_; }

    //! Index of current subblock.
    internal_size_type i_subblock() { return i_value_ / subblock_size; }
};

/*!
 * Additional information about a stored value:
 * - the bucket in which it can be found
 * - where it is currently stored (intern or extern)
 * - the buffer-node
 * - the position in external memory
 */
template <class HashMap>
struct HashedValue
{
    typedef HashMap hash_map_type;
    typedef typename hash_map_type::value_type value_type;
    typedef typename hash_map_type::source_type source_type;
    typedef typename hash_map_type::node_type node_type;

    typedef typename hash_map_type::internal_size_type internal_size_type;
    typedef typename hash_map_type::external_size_type external_size_type;

    value_type value_;
    internal_size_type i_bucket_;
    source_type source_;
    node_type* node_;
    external_size_type i_external_;

    HashedValue()
        : i_bucket_(internal_size_type(-1))
    { }

    HashedValue(const value_type& value, internal_size_type i_bucket,
                source_type src, node_type* node, external_size_type i_external)
        : value_(value),
          i_bucket_(i_bucket),
          source_(src),
          node_(node),
          i_external_(i_external)
    { }
};

/*!
 * Stream interface for all value-pairs currently stored in the map. Returned
 * values are HashedValue-objects (actual value enriched with information on
 * where the value can be found (bucket-number, internal, external)).  Values,
 * marked as deleted in internal-memory, are not returned; for modified values
 * only the one in internal memory is returned.
*/
template <class HashMap, class Reader>
struct HashedValuesStream
{
    typedef HashMap hash_map_type;
    typedef HashedValue<HashMap> value_type;

    typedef typename hash_map_type::node_type node_type;
    typedef typename hash_map_type::bid_container_type::iterator bid_iterator;
    typedef typename hash_map_type::buckets_container_type::iterator bucket_iterator;

    typedef typename hash_map_type::internal_size_type internal_size_type;
    typedef typename hash_map_type::external_size_type external_size_type;

    hash_map_type& map_;
    Reader& reader_;
    bucket_iterator curr_bucket_;
    bucket_iterator end_bucket_;
    bid_iterator begin_bid_;
    internal_size_type i_bucket_;
    node_type* node_;
    external_size_type i_external_;
    value_type value_;

    HashedValuesStream(bucket_iterator begin_bucket, bucket_iterator end_bucket,
                       Reader& reader, bid_iterator begin_bid,
                       hash_map_type& map)
        : map_(map),
          reader_(reader),
          curr_bucket_(begin_bucket),
          end_bucket_(end_bucket),
          begin_bid_(begin_bid),
          i_bucket_(0),
          node_(curr_bucket_ != end_bucket_ ? curr_bucket_->list_ : NULL),
          i_external_(0)
    {
        if (!empty())
            value_ = find_next();
    }

    const value_type& operator * () { return value_; }

    bool empty() const { return curr_bucket_ == end_bucket_; }

    void operator ++ ()
    {
        if (value_.source_ == hash_map_type::src_internal)
            node_ = node_->next();
        else
        {
            ++reader_;
            ++i_external_;
        }
        value_ = find_next();
    }

    value_type find_next()
    {
        while (true)
        {
            // internal and external elements available
            while (node_ && i_external_ < curr_bucket_->n_external_)
            {
                if (map_._leq(node_->value_.first, reader_.const_value().first))
                {
                    if (map_._eq(node_->value_.first, reader_.const_value().first))
                    {
                        ++reader_;
                        ++i_external_;
                    }

                    if (!node_->deleted())
                        return value_type(node_->value_, i_bucket_, hash_map_type::src_internal, node_, i_external_);
                    else
                        node_ = node_->next();
                }
                else
                    return value_type(reader_.const_value(), i_bucket_, hash_map_type::src_external, node_, i_external_);
            }
            // only internal elements left
            while (node_)
            {
                if (!node_->deleted())
                    return value_type(node_->value_, i_bucket_, hash_map_type::src_internal, node_, i_external_);
                else
                    node_ = node_->next();
            }
            // only external elements left
            while (i_external_ < curr_bucket_->n_external_)
                return value_type(reader_.const_value(), i_bucket_, hash_map_type::src_external, node_, i_external_);

            // if we made it to this point there are obviously no more values in the current bucket
            // let's try the next one (outer while-loop!)
            ++curr_bucket_;
            ++i_bucket_;
            if (curr_bucket_ == end_bucket_)
                return value_type();

            node_ = curr_bucket_->list_;
            i_external_ = 0;
            reader_.skip_to(begin_bid_ + curr_bucket_->i_block_, curr_bucket_->i_subblock_);
        }
    }
};

} // namespace hash_map

STXXL_END_NAMESPACE

#endif // !STXXL_CONTAINERS_HASH_MAP_UTIL_HEADER
