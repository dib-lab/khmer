/***************************************************************************
 *  include/stxxl/bits/containers/btree/iterator_map.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2006 Roman Dementiev <dementiev@ira.uka.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_CONTAINERS_BTREE_ITERATOR_MAP_HEADER
#define STXXL_CONTAINERS_BTREE_ITERATOR_MAP_HEADER

#include <map>

#include <stxxl/bits/noncopyable.h>
#include <stxxl/bits/containers/btree/iterator.h>
#include <stxxl/bits/common/error_handling.h>

STXXL_BEGIN_NAMESPACE

namespace btree {

template <class BTreeType>
class iterator_map : private noncopyable
{
public:
    typedef BTreeType btree_type;
    typedef typename btree_type::leaf_bid_type bid_type;
    typedef btree_iterator_base<btree_type> iterator_base;

private:
    struct Key
    {
        bid_type bid;
        unsigned_type pos;
        Key() { }
        Key(const bid_type& b, unsigned_type p)
            : bid(b), pos(p) { }
    };

    struct bid_comp
    {
        bool operator () (const bid_type& a, const bid_type& b) const
        {
            return (a.storage < b.storage) || (a.storage == b.storage && a.offset < b.offset);
        }
    };
    struct KeyCmp
    {
        bid_comp BIDComp;
        bool operator () (const Key& a, const Key& b) const
        {
            return BIDComp(a.bid, b.bid) || (a.bid == b.bid && a.pos < b.pos);
        }
    };

    typedef std::multimap<Key, iterator_base*, KeyCmp> multimap_type;

    multimap_type m_it2addr;
    btree_type* m_btree;

    typedef typename multimap_type::value_type pair_type;
    typedef typename multimap_type::iterator mmiterator_type;
    typedef typename multimap_type::const_iterator mmconst_iterator_type;

    //! changes btree pointer in all contained iterators
    void change_btree_pointers(btree_type* b)
    {
        for (mmconst_iterator_type it = m_it2addr.begin();
             it != m_it2addr.end(); ++it)
        {
            (it->second)->btree = b;
        }
    }

public:
    iterator_map(btree_type* b)
        : m_btree(b)
    { }

    void register_iterator(iterator_base& it)
    {
        STXXL_VERBOSE2("btree::iterator_map register_iterator addr=" << &it <<
                       " BID: " << it.bid << " POS: " << it.pos);
        m_it2addr.insert(pair_type(Key(it.bid, it.pos), &it));
    }
    void unregister_iterator(iterator_base& it)
    {
        STXXL_VERBOSE2("btree::iterator_map unregister_iterator addr=" << &it <<
                       " BID: " << it.bid << " POS: " << it.pos);
        assert(!m_it2addr.empty());
        Key key(it.bid, it.pos);
        std::pair<mmiterator_type, mmiterator_type> range =
            m_it2addr.equal_range(key);

        assert(range.first != range.second);

        mmiterator_type i = range.first;
        for ( ; i != range.second; ++i)
        {
            assert(it.bid == i->first.bid);
            assert(it.pos == i->first.pos);

            if (i->second == &it)
            {
                // found it
                m_it2addr.erase(i);
                return;
            }
        }

        STXXL_THROW2(std::runtime_error, "btree::iterator_map::unregister_iterator", "Panic in btree::iterator_map, can not find and unregister iterator");
    }
    template <class OutputContainer>
    void find(const bid_type& bid,
              unsigned_type first_pos,
              unsigned_type last_pos,
              OutputContainer& out)
    {
        Key firstkey(bid, first_pos);
        Key lastkey(bid, last_pos);
        mmconst_iterator_type begin = m_it2addr.lower_bound(firstkey);
        mmconst_iterator_type end = m_it2addr.upper_bound(lastkey);

        for (mmconst_iterator_type i = begin;
             i != end; ++i)
        {
            assert(bid == i->first.bid);
            out.push_back(i->second);
        }
    }

    virtual ~iterator_map()
    {
        for (mmconst_iterator_type it = m_it2addr.begin();
             it != m_it2addr.end(); ++it)
        {
            it->second->make_invalid();
        }
    }

    void swap(iterator_map& obj)
    {
        std::swap(m_it2addr, obj.m_it2addr);
        change_btree_pointers(m_btree);
        obj.change_btree_pointers(obj.m_btree);
    }
};

} // namespace btree

STXXL_END_NAMESPACE

namespace std {

template <class BTreeType>
void swap(stxxl::btree::iterator_map<BTreeType>& a,
          stxxl::btree::iterator_map<BTreeType>& b)
{
    a.swap(b);
}

} // namespace std

#endif // !STXXL_CONTAINERS_BTREE_ITERATOR_MAP_HEADER
