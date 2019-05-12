/***************************************************************************
 *  include/stxxl/bits/containers/hash_map/iterator_map.h
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

#ifndef STXXL_CONTAINERS_HASH_MAP_ITERATOR_MAP_HEADER
#define STXXL_CONTAINERS_HASH_MAP_ITERATOR_MAP_HEADER

#include <map>

#include <stxxl/bits/noncopyable.h>
#include <stxxl/bits/compat/hash_map.h>
#include <stxxl/bits/containers/hash_map/iterator.h>

STXXL_BEGIN_NAMESPACE

namespace hash_map {

template <class HashMap>
class iterator_map : private noncopyable
{
public:
    typedef HashMap hash_map_type;
    typedef typename hash_map_type::node_type node_type;
    typedef typename hash_map_type::source_type source_type;
    typedef typename hash_map_type::key_type key_type;

    typedef typename hash_map_type::internal_size_type internal_size_type;
    typedef typename hash_map_type::external_size_type external_size_type;

    typedef hash_map_iterator_base<hash_map_type> iterator_base;

private:
#if 0
    typedef std::multimap<internal_size_type, iterator_base*> multimap_type;
#else
    struct hasher
    {
        size_t operator () (const internal_size_type& key) const
        {
            return longhash1(key);
        }
#if STXXL_MSVC
        bool operator () (const internal_size_type& a, const internal_size_type& b) const
        {
            return (a < b);
        }
        enum
        {                                       // parameters for hash table
            bucket_size = 4,                    // 0 < bucket_size
            min_buckets = 8                     // min_buckets = 2 ^^ N, 0 < N
        };
#endif
    };
    // store iterators by bucket-index
    typedef typename compat_hash_multimap<
            internal_size_type, iterator_base*, hasher
            >::result multimap_type;
#endif

    //! bucket-index and pointer to iterator_base
    typedef typename multimap_type::value_type pair_type;
    typedef typename multimap_type::iterator mmiterator_type;
    typedef typename multimap_type::const_iterator const_mmiterator_type;

    hash_map_type* map_;
    multimap_type it_map_;

public:
    iterator_map(hash_map_type* map)
        : map_(map)
    { }

    ~iterator_map()
    {
        it_map_.clear();
    }

    void register_iterator(iterator_base& it)
    {
        register_iterator(it, it.i_bucket_);
    }

    void register_iterator(iterator_base& it, internal_size_type i_bucket)
    {
        STXXL_VERBOSE2("hash_map::iterator_map register_iterator addr=" << &it << " bucket=" << i_bucket);
        it_map_.insert(pair_type(i_bucket, &it));
    }

    void unregister_iterator(iterator_base& it)
    {
        unregister_iterator(it, it.i_bucket_);
    }

    void unregister_iterator(iterator_base& it, internal_size_type i_bucket)
    {
        STXXL_VERBOSE2("hash_map::iterator_map unregister_iterator addr=" << &it << " bucket=" << i_bucket);

        std::pair<mmiterator_type, mmiterator_type> range
            = it_map_.equal_range(i_bucket);

        assert(range.first != range.second);

        for (mmiterator_type i = range.first; i != range.second; ++i)
        {
            if (i->second == &it)
            {
                it_map_.erase(i);
                return;
            }
        }

        throw std::runtime_error("unregister_iterator Panic in hash_map::iterator_map, can not find and unregister iterator");
    }

    //! Update iterators with given key and bucket and make them point to the
    //! specified location in external memory (will be called during
    //! re-hashing)
    void fix_iterators_2ext(internal_size_type i_bucket_old, const key_type& key,
                            internal_size_type i_bucket_new, external_size_type i_ext)
    {
        STXXL_VERBOSE2("hash_map::iterator_map fix_iterators_2ext i_bucket=" << i_bucket_old << " new_i_ext=" << i_ext);

        std::vector<iterator_base*> its2fix;
        _find(i_bucket_old, its2fix);

        for (typename std::vector<iterator_base*>::iterator
             it2fix = its2fix.begin(); it2fix != its2fix.end(); ++it2fix)
        {
            if (!map_->_eq(key, (**it2fix).key_))
                continue;

            if (i_bucket_old != i_bucket_new)
            {
                unregister_iterator(**it2fix);
                register_iterator(**it2fix, i_bucket_new);
            }

            (**it2fix).i_bucket_ = i_bucket_new;
            (**it2fix).node_ = NULL;
            (**it2fix).i_external_ = i_ext;
            (**it2fix).source_ = hash_map_type::src_external;
            // external position is now known (i_ext) and therefore valid
            (**it2fix).ext_valid_ = true;
            (**it2fix).reset_reader();
            (**it2fix).reader_ = NULL;
        }
    }

    //! Update iterators with given key and bucket and make them point to the
    //! specified node in internal memory (will be called by insert_oblivious)
    void fix_iterators_2int(internal_size_type i_bucket, const key_type& key, node_type* node)
    {
        STXXL_VERBOSE2("hash_map::iterator_map fix_iterators_2int i_bucket=" << i_bucket << " node=" << node);

        std::vector<iterator_base*> its2fix;
        _find(i_bucket, its2fix);

        for (typename std::vector<iterator_base*>::iterator
             it2fix = its2fix.begin(); it2fix != its2fix.end(); ++it2fix)
        {
            if (!map_->_eq((**it2fix).key_, key))
                continue;

            assert((**it2fix).source_ == hash_map_type::src_external);

            (**it2fix).source_ = hash_map_type::src_internal;
            (**it2fix).node_ = node;
            (**it2fix).i_external_++;
            if ((**it2fix).reader_)
                (**it2fix).reader_->operator ++ ();
        }
    }

    //! Update iterators with given key and bucket and make them point to the
    //! end of the hash-map (called by erase and erase_oblivious)
    void fix_iterators_2end(internal_size_type i_bucket, const key_type& key)
    {
        STXXL_VERBOSE2("hash_map::iterator_map fix_iterators_2end i_bucket=" << i_bucket);

        std::vector<iterator_base*> its2fix;
        _find(i_bucket, its2fix);

        for (typename std::vector<iterator_base*>::iterator
             it2fix = its2fix.begin(); it2fix != its2fix.end(); ++it2fix)
        {
            if (!map_->_eq(key, (**it2fix).key_))
                continue;

            (**it2fix).end_ = true;
            (**it2fix).reset_reader();
            unregister_iterator(**it2fix);
        }
    }

    //! Update all iterators and make them point to the end of the hash-map
    //! (used by clear())
    void fix_iterators_all2end()
    {
        for (mmiterator_type it2fix = it_map_.begin();
             it2fix != it_map_.end(); ++it2fix)
        {
            (*it2fix).second->end_ = true;
            (*it2fix).second->reset_reader();
        }
        it_map_.clear();
    }

private:
    //! Find all iterators registered with given bucket and add them to outc
    template <class OutputContainer>
    void _find(internal_size_type i_bucket, OutputContainer& outc)
    {
        std::pair<mmiterator_type, mmiterator_type> range
            = it_map_.equal_range(i_bucket);

        for (mmiterator_type i = range.first; i != range.second; ++i)
            outc.push_back((*i).second);
    }

    // changes hash_map pointer in all contained iterators
    void change_hash_map_pointers(hash_map_type* map)
    {
        for (mmiterator_type it = it_map_.begin(); it != it_map_.end(); ++it)
            ((*it).second)->map_ = map;
    }

public:
    void swap(iterator_map<HashMap>& obj)
    {
        std::swap(it_map_, obj.it_map_);
        std::swap(map_, obj.map_);

        change_hash_map_pointers(map_);
        obj.change_hash_map_pointers(obj.map_);
    }

    void print_statistics(std::ostream& o = std::cout) const
    {
        o << "Registered iterators: " << it_map_.size() << "\n";

        for (const_mmiterator_type i = it_map_.begin(); i != it_map_.end(); ++i)
        {
            o << "  Address=" << i->second
              << ", Bucket=" << i->second->i_bucket_
              << ", Node=" << i->second->node_
              << ", i_ext=" << i->second->i_external_
              << ", "
              << ((i->second->source_ == hash_map_type::src_external)
                ? "external" : "internal") << std::endl;
        }
    }
};

} // namespace hash_map

STXXL_END_NAMESPACE

namespace std {

template <class HashMapType>
void swap(stxxl::hash_map::iterator_map<HashMapType>& a,
          stxxl::hash_map::iterator_map<HashMapType>& b)
{
    if (&a != &b)
        a.swap(b);
}

} // namespace std

#endif // !STXXL_CONTAINERS_HASH_MAP_ITERATOR_MAP_HEADER
