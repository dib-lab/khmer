/***************************************************************************
 *  include/stxxl/bits/containers/hash_map/iterator.h
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

#ifndef STXXL_CONTAINERS_HASH_MAP_ITERATOR_HEADER
#define STXXL_CONTAINERS_HASH_MAP_ITERATOR_HEADER

#include <stxxl/bits/namespace.h>
#include <stxxl/bits/mng/block_manager.h>

#include <stxxl/bits/containers/hash_map/util.h>

STXXL_BEGIN_NAMESPACE

namespace hash_map {

template <class HashMap>
class iterator_map;
template <class HashMap>
class hash_map_iterator;
template <class HashMap>
class hash_map_const_iterator;
template <class HashMap>
class block_cache;

template <class HashMap>
class hash_map_iterator_base
{
public:
    friend class iterator_map<HashMap>;
    friend void HashMap::erase(hash_map_const_iterator<HashMap> it);

    typedef HashMap hash_map_type;
    typedef typename hash_map_type::internal_size_type internal_size_type;
    typedef typename hash_map_type::external_size_type external_size_type;
    typedef typename hash_map_type::value_type value_type;
    typedef typename hash_map_type::key_type key_type;
    typedef typename hash_map_type::reference reference;
    typedef typename hash_map_type::const_reference const_reference;
    typedef typename hash_map_type::node_type node_type;
    typedef typename hash_map_type::bucket_type bucket_type;
    typedef typename hash_map_type::bid_iterator_type bid_iterator_type;
    typedef typename hash_map_type::source_type source_type;

    typedef buffered_reader<typename hash_map_type::block_cache_type, bid_iterator_type> reader_type;

    typedef std::forward_iterator_tag iterator_category;

protected:
    HashMap* map_;
    reader_type* reader_;
    //! true if prefetching enabled; false by default, will be set to true when
    //! incrementing (see find_next())
    bool prefetch_;
    //! index of current bucket
    internal_size_type i_bucket_;
    //! source of current value: external or internal
    source_type source_;
    //! current (source=internal) or old (src=external) internal node
    node_type* node_;
    //! position of current (source=external) or next (source=internal)
    //! external value (see _ext_valid)
    external_size_type i_external_;
    //! key of current value
    key_type key_;
    /*! true if i_external points to the current or next external value

      example: iterator was created by hash_map::find() and the value was found
      in internal memory

      => iterator pointing to internal node is created and location of next
      external value is unknown (_ext_valid == false)

      => when incrementing the iterator the external values will be scanned
      from the beginning of the bucket to find the valid external index
    */
    bool ext_valid_;
    //! true if iterator equals end()
    bool end_;

public:
    //! Construct a new iterator
    hash_map_iterator_base(HashMap* map, internal_size_type i_bucket, node_type* node,
                           external_size_type i_external, source_type source,
                           bool ext_valid, key_type key)
        : map_(map),
          reader_(NULL),
          prefetch_(false),
          i_bucket_(i_bucket),
          source_(source),
          node_(node),
          i_external_(i_external),
          key_(key),
          ext_valid_(ext_valid),
          end_(false)
    {
        STXXL_VERBOSE3("hash_map_iterator_base parameter construct addr=" << this);
        map_->iterator_map_.register_iterator(*this);
    }

    //! Construct a new iterator pointing to the end of the given hash-map.
    hash_map_iterator_base(hash_map_type* map)
        : map_(map),
          reader_(NULL),
          prefetch_(false),
          i_bucket_(0),
          source_(hash_map_type::src_unknown),
          node_(NULL),
          i_external_(0),
          ext_valid_(false),
          end_(true)
    { }

    //! Construct a new iterator from an existing one
    hash_map_iterator_base(const hash_map_iterator_base& obj)
        : map_(obj.map_),
          reader_(NULL),
          prefetch_(obj.prefetch_),
          i_bucket_(obj.i_bucket_),
          source_(obj.source_),
          node_(obj.node_),
          i_external_(obj.i_external_),
          key_(obj.key_),
          ext_valid_(obj.ext_valid_),
          end_(obj.end_)
    {
        STXXL_VERBOSE3("hash_map_iterator_base constr from" << (&obj) << " to " << this);

        if (!end_ && map_)
            map_->iterator_map_.register_iterator(*this);
    }

    //! Assignment operator
    hash_map_iterator_base& operator = (const hash_map_iterator_base& obj)
    {
        STXXL_VERBOSE3("hash_map_iterator_base copy from" << (&obj) << " to " << this);

        if (&obj != this)
        {
            if (map_ && !end_)
                map_->iterator_map_.unregister_iterator(*this);

            reset_reader();

            map_ = obj.map_;
            i_bucket_ = obj.i_bucket_;
            node_ = obj.node_;
            source_ = obj.source_;
            i_external_ = obj.i_external_;
            ext_valid_ = obj.ext_valid_;
            prefetch_ = obj.prefetch_;
            end_ = obj.end_;
            key_ = obj.key_;

            if (map_ && !end_)
                map_->iterator_map_.register_iterator(*this);
        }
        return *this;
    }

    //! Two iterators are equal if the point to the same value in the same map
    bool operator == (const hash_map_iterator_base& obj) const
    {
        if (end_ && obj.end_)
            return true;

        if (map_ != obj.map_ ||
            i_bucket_ != obj.i_bucket_ ||
            source_ != obj.source_)
            return false;

        if (source_ == hash_map_type::src_internal)
            return node_ == obj.node_;
        else
            return i_external_ == obj.i_external_;
    }

    bool operator != (const hash_map_iterator_base& obj) const
    {
        return ! operator == (obj);
    }

protected:
    //! Initialize reader object to scan external values
    void init_reader()
    {
        const bucket_type& bucket = map_->buckets_[i_bucket_];

        bid_iterator_type begin = map_->bids_.begin() + bucket.i_block_;
        bid_iterator_type end = map_->bids_.end();

        reader_ = new reader_type(begin, end, map_->block_cache_,
                                  bucket.i_subblock_, prefetch_);

        // external value's index already known
        if (ext_valid_)
        {
            // TODO: speed this up (go directly to i_external_
            for (external_size_type i = 0; i < i_external_; i++)
                ++(*reader_);
        }
        // otw lookup external value.
        // case I: no internal value => first external value is the desired one
        else if (node_ == NULL)
        {
            i_external_ = 0;
            ext_valid_ = true;
        }
        // case II: search for smallest external value > internal value
        else
        {
            i_external_ = 0;
            while (i_external_ < bucket.n_external_)
            {
                if (map_->_gt(reader_->const_value().first, node_->value_.first))
                    break;

                ++(*reader_);
                ++i_external_;
            }
            // note: i_external==num_external just means that there was no
            // external value > internal value (which is perfectly OK)
            ext_valid_ = true;
        }
    }

    //! Reset reader-object
    void reset_reader()
    {
        if (reader_) {
            delete reader_;
            reader_ = NULL;
        }
    }

public:
    //! Advance iterator to the next value
    //! The next value is determined in the following way
    //! - if there are remaining internal or external values in the current
    //!   bucket, choose the smallest among them, that is not marked as deleted
    //! - otherwise continue with the next bucket
    void find_next(bool start_prefetching = false)
    {
        // invariant: current external value is always > current internal value
        assert(!end_);

        internal_size_type i_bucket_old = i_bucket_;
        bucket_type bucket = map_->buckets_[i_bucket_];

        if (reader_ == NULL)
            init_reader();

        // when incremented once, more increments are likely to follow;
        // therefore start prefetching
        if (start_prefetching && !prefetch_)
        {
            reader_->enable_prefetching();
            prefetch_ = true;
        }

        // determine starting-points for comparision, which are given by:
        // - tmp_node: smallest internal value > old value (tmp_node may be NULL)
        // - reader_: smallest external value > old value (external value may not exists)
        node_type* tmp_node = (node_) ? node_ : bucket.list_;
        if (source_ == hash_map_type::src_external)
        {
            while (tmp_node && map_->_leq(tmp_node->value_.first, key_))
                tmp_node = tmp_node->next();

            ++i_external_;
            ++(*reader_);
        }
        else if (source_ == hash_map_type::src_internal)
            tmp_node = node_->next();
        // else (source unknown): tmp_node and reader_ already point to the
        // correct values

        while (true) {
            // internal and external values available
            while (tmp_node && i_external_ < bucket.n_external_)
            {
                // internal value less or equal external value => internal wins
                if (map_->_leq(tmp_node->value_.first, reader_->const_value().first))
                {
                    node_ = tmp_node;
                    if (map_->_eq(node_->value_.first, reader_->const_value().first))
                    {
                        ++i_external_;
                        ++(*reader_);
                    }

                    if (!node_->deleted())
                    {
                        key_ = node_->value_.first;
                        source_ = hash_map_type::src_internal;
                        goto end_search;       // just this once - I promise...
                    }
                    else
                        // continue search if internal value flaged as deleted
                        tmp_node = tmp_node->next();
                }
                // otherwise external wins
                else
                {
                    key_ = reader_->const_value().first;
                    source_ = hash_map_type::src_external;
                    goto end_search;
                }
            }
            // only external values left
            if (i_external_ < bucket.n_external_)
            {
                key_ = reader_->const_value().first;
                source_ = hash_map_type::src_external;
                goto end_search;
            }
            // only internal values left
            while (tmp_node)
            {
                node_ = tmp_node;
                if (!node_->deleted())
                {
                    key_ = node_->value_.first;
                    source_ = hash_map_type::src_internal;
                    goto end_search;
                }
                else
                    tmp_node = tmp_node->next();        // continue search
            }

            // at this point there are obviously no more values in the current
            // bucket let's try the next one (outer while-loop!)
            i_bucket_++;
            if (i_bucket_ == map_->buckets_.size())
            {
                end_ = true;
                reset_reader();
                goto end_search;
            }
            else
            {
                bucket = map_->buckets_[i_bucket_];
                i_external_ = 0;
                tmp_node = bucket.list_;
                node_ = NULL;
                reader_->skip_to(map_->bids_.begin() + bucket.i_block_, bucket.i_subblock_);
            }
        }

end_search:
        if (end_)
        {
            this->map_->iterator_map_.unregister_iterator(*this, i_bucket_old);
        }
        else if (i_bucket_old != i_bucket_)
        {
            this->map_->iterator_map_.unregister_iterator(*this, i_bucket_old);
            this->map_->iterator_map_.register_iterator(*this, i_bucket_);
        }
    }

    virtual ~hash_map_iterator_base()
    {
        STXXL_VERBOSE3("hash_map_iterator_base deconst " << this);

        if (map_ && !end_)
            map_->iterator_map_.unregister_iterator(*this);
        reset_reader();
    }
};

template <class HashMap>
class hash_map_iterator : public hash_map_iterator_base<HashMap>
{
public:
    typedef HashMap hash_map_type;
    typedef typename hash_map_type::internal_size_type internal_size_type;
    typedef typename hash_map_type::external_size_type external_size_type;
    typedef typename hash_map_type::value_type value_type;
    typedef typename hash_map_type::key_type key_type;
    typedef typename hash_map_type::reference reference;
    typedef typename hash_map_type::const_reference const_reference;
    typedef typename hash_map_type::pointer pointer;
    typedef typename hash_map_type::const_pointer const_pointer;
    typedef typename hash_map_type::node_type node_type;
    typedef typename hash_map_type::bucket_type bucket_type;
    typedef typename hash_map_type::bid_iterator_type bid_iterator_type;
    typedef typename hash_map_type::source_type source_type;

    typedef buffered_reader<typename hash_map_type::block_cache_type,
                            bid_iterator_type> reader_type;

    typedef std::forward_iterator_tag iterator_category;

    typedef stxxl::hash_map::hash_map_iterator_base<hash_map_type> base_type;
    typedef stxxl::hash_map::hash_map_const_iterator<hash_map_type> hash_map_const_iterator;

public:
    hash_map_iterator(hash_map_type* map, internal_size_type i_bucket,
                      node_type* node, external_size_type i_external,
                      source_type source, bool ext_valid, key_type key)
        : base_type(map, i_bucket, node, i_external, source, ext_valid, key)
    { }

    hash_map_iterator()
        : base_type(NULL)
    { }

    hash_map_iterator(hash_map_type* map)
        : base_type(map)
    { }

    hash_map_iterator(const hash_map_iterator& obj)
        : base_type(obj)
    { }

    hash_map_iterator& operator = (const hash_map_iterator& obj)
    {
        base_type::operator = (obj);
        return *this;
    }

    bool operator == (const hash_map_iterator& obj) const
    {
        return base_type::operator == (obj);
    }

    bool operator == (const hash_map_const_iterator& obj) const
    {
        return base_type::operator == (obj);
    }

    bool operator != (const hash_map_iterator& obj) const
    {
        return base_type::operator != (obj);
    }

    bool operator != (const hash_map_const_iterator& obj) const
    {
        return base_type::operator != (obj);
    }

    //! Return reference to current value. If source is external, mark the
    //! value's block as dirty
    reference operator * ()
    {
        if (this->source_ == hash_map_type::src_internal)
        {
            return this->node_->value_;
        }
        else
        {
            if (this->reader_ == NULL)
                base_type::init_reader();

            return this->reader_->value();
        }
    }

    //! Return reference to current value. If source is external, mark the
    //! value's block as dirty
    pointer operator -> ()
    {
        return &operator * ();
    }

    //! Increment iterator
    hash_map_iterator<hash_map_type>& operator ++ ()
    {
        base_type::find_next(true);
        return *this;
    }
};

template <class HashMap>
class hash_map_const_iterator : public hash_map_iterator_base<HashMap>
{
public:
    typedef HashMap hash_map_type;
    typedef typename hash_map_type::internal_size_type internal_size_type;
    typedef typename hash_map_type::external_size_type external_size_type;
    typedef typename hash_map_type::value_type value_type;
    typedef typename hash_map_type::key_type key_type;
    typedef typename hash_map_type::reference reference;
    typedef typename hash_map_type::const_reference const_reference;
    typedef typename hash_map_type::pointer pointer;
    typedef typename hash_map_type::const_pointer const_pointer;
    typedef typename hash_map_type::node_type node_type;
    typedef typename hash_map_type::bucket_type bucket_type;
    typedef typename hash_map_type::bid_iterator_type bid_iterator_type;
    typedef typename hash_map_type::source_type source_type;

    typedef buffered_reader<typename hash_map_type::block_cache_type,
                            bid_iterator_type> reader_type;

    typedef std::forward_iterator_tag iterator_category;

    typedef stxxl::hash_map::hash_map_iterator_base<hash_map_type> base_type;
    typedef stxxl::hash_map::hash_map_iterator<hash_map_type> hash_map_iterator;

public:
    hash_map_const_iterator(hash_map_type* map, internal_size_type i_bucket,
                            node_type* node, external_size_type i_external,
                            source_type source, bool ext_valid, key_type key)
        : base_type(map, i_bucket, node, i_external, source, ext_valid, key)
    { }

    hash_map_const_iterator()
        : base_type(NULL)
    { }

    hash_map_const_iterator(hash_map_type* map)
        : base_type(map)
    { }

    hash_map_const_iterator(const hash_map_iterator& obj)
        : base_type(obj)
    { }

    hash_map_const_iterator(const hash_map_const_iterator& obj)
        : base_type(obj)
    { }

    hash_map_const_iterator& operator = (const hash_map_const_iterator& obj)
    {
        base_type::operator = (obj);
        return *this;
    }

    bool operator == (const hash_map_const_iterator& obj) const
    {
        return base_type::operator == (obj);
    }

    bool operator == (const hash_map_iterator& obj) const
    {
        return base_type::operator == (obj);
    }

    bool operator != (const hash_map_const_iterator& obj) const
    {
        return base_type::operator != (obj);
    }

    bool operator != (const hash_map_iterator& obj) const
    {
        return base_type::operator != (obj);
    }

    //! Return const-reference to current value
    const_reference operator * ()
    {
        if (this->source_ == hash_map_type::src_internal)
        {
            return this->node_->value_;
        }
        else
        {
            if (this->reader_ == NULL)
                base_type::init_reader();

            return this->reader_->const_value();
        }
    }

    //! Increment iterator
    hash_map_const_iterator<hash_map_type>& operator ++ ()
    {
        base_type::find_next(true);
        return *this;
    }
};

} // namespace hash_map

STXXL_END_NAMESPACE

#endif // !STXXL_CONTAINERS_HASH_MAP_ITERATOR_HEADER
