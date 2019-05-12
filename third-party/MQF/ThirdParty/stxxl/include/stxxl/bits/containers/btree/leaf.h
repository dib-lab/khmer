/***************************************************************************
 *  include/stxxl/bits/containers/btree/leaf.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2006 Roman Dementiev <dementiev@ira.uka.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_CONTAINERS_BTREE_LEAF_HEADER
#define STXXL_CONTAINERS_BTREE_LEAF_HEADER

#include <stxxl/bits/containers/btree/iterator.h>
#include <stxxl/bits/containers/btree/node_cache.h>

STXXL_BEGIN_NAMESPACE

namespace btree {

template <class NodeType, class BTreeType>
class node_cache;

template <class KeyType, class DataType, class KeyCmp, unsigned RawSize, class BTreeType>
class normal_leaf : private noncopyable
{
public:
    typedef normal_leaf<KeyType, DataType, KeyCmp, RawSize, BTreeType> self_type;

    friend class node_cache<self_type, BTreeType>;

    typedef KeyType key_type;
    typedef DataType data_type;
    typedef KeyCmp key_compare;
    typedef std::pair<key_type, data_type> value_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;

    enum {
        raw_size = RawSize
    };

    typedef BID<raw_size> bid_type;

    struct metainfo_type
    {
        bid_type me, pred, succ;
        unsigned cur_size;
    };

    typedef typed_block<raw_size, value_type, 0, metainfo_type> block_type;
    enum {
        nelements = block_type::size - 1,
        max_size = nelements,
        min_size = nelements / 2
    };

    typedef BTreeType btree_type;
    typedef typename btree_type::size_type size_type;
    typedef btree_iterator_base<btree_type> iterator_base;
    typedef btree_iterator<btree_type> iterator;
    typedef btree_const_iterator<btree_type> const_iterator;

    typedef node_cache<normal_leaf, btree_type> leaf_cache_type;

public:
    struct value_compare : public std::binary_function<value_type, value_type, bool>
    {
        key_compare comp;

        value_compare(key_compare c) : comp(c) { }

        bool operator () (const value_type& x, const value_type& y) const
        {
            return comp(x.first, y.first);
        }
    };

private:
    block_type* m_block;
    btree_type* m_btree;

    key_compare m_cmp;
    value_compare m_vcmp;

    void split(std::pair<key_type, bid_type>& splitter)
    {
        bid_type new_bid;
        m_btree->m_leaf_cache.get_new_node(new_bid);                         // new (left) leaf
        normal_leaf* new_leaf = m_btree->m_leaf_cache.get_node(new_bid, true);
        assert(new_leaf);

        // update links
        new_leaf->succ() = my_bid();
        normal_leaf* pred_leaf = NULL;
        if (pred().valid())
        {
            new_leaf->pred() = pred();
            pred_leaf = m_btree->m_leaf_cache.get_node(pred());
            assert(pred_leaf);
            assert(m_vcmp(pred_leaf->back(), front()));
            pred_leaf->succ() = new_bid;
        }
        pred() = new_bid;

        typedef std::vector<iterator_base*> iterators2fix_type;
        iterators2fix_type iterators2fix;
        m_btree->m_iterator_map.find(my_bid(), 0, size(), iterators2fix);

        const unsigned end_of_smaller_part = size() / 2;

        splitter.first = ((*m_block)[end_of_smaller_part - 1]).first;
        splitter.second = new_bid;

        const unsigned old_size = size();
        // copy the smaller part
        std::copy(m_block->begin(), m_block->begin() + end_of_smaller_part,
                  new_leaf->m_block->begin());
        new_leaf->m_block->info.cur_size = end_of_smaller_part;
        // copy the larger part
        std::copy(m_block->begin() + end_of_smaller_part,
                  m_block->begin() + old_size, m_block->begin());
        m_block->info.cur_size = old_size - end_of_smaller_part;
        assert(size() + new_leaf->size() == old_size);

        // fix iterators
        for (typename iterators2fix_type::iterator it2fix = iterators2fix.begin();
             it2fix != iterators2fix.end(); ++it2fix)
        {
            m_btree->m_iterator_map.unregister_iterator(**it2fix);

            if ((*it2fix)->pos < end_of_smaller_part)     // belongs to the smaller part
                (*it2fix)->bid = new_bid;

            else
                (*it2fix)->pos -= end_of_smaller_part;

            m_btree->m_iterator_map.register_iterator(**it2fix);
        }

        STXXL_VERBOSE1("btree::normal_leaf split leaf " << this
                                                        << " splitter: " << splitter.first);

#if STXXL_VERBOSE_LEVEL >= 1
        if (pred_leaf)
        {
            STXXL_VERBOSE1("btree::normal_leaf pred_part.smallest    = " << pred_leaf->front().first);
            STXXL_VERBOSE1("btree::normal_leaf pred_part.largest     = " << pred_leaf->back().first);
        }
#endif
        STXXL_VERBOSE1("btree::normal_leaf smaller_part.smallest = " << new_leaf->front().first);
        STXXL_VERBOSE1("btree::normal_leaf smaller_part.largest  = " << new_leaf->back().first);
        STXXL_VERBOSE1("btree::normal_leaf larger_part.smallest  = " << front().first);
        STXXL_VERBOSE1("btree::normal_leaf larger_part.largest   = " << back().first);

        m_btree->m_leaf_cache.unfix_node(new_bid);
    }

public:
    virtual ~normal_leaf()
    {
        delete m_block;
    }

    normal_leaf(btree_type* btree,
                key_compare cmp)
        : m_block(new block_type),
          m_btree(btree),
          m_cmp(cmp),
          m_vcmp(cmp)
    {
        assert(min_nelements() >= 2);
        assert(2 * min_nelements() - 1 <= max_nelements());
        assert(max_nelements() <= nelements);
        assert(unsigned(block_type::size) >= nelements + 1);                       // extra space for an overflow
    }

    bool overflows() const { return m_block->info.cur_size > max_nelements(); }
    bool underflows() const { return m_block->info.cur_size < min_nelements(); }

    static unsigned max_nelements() { return max_size; }
    static unsigned min_nelements() { return min_size; }

    bid_type & succ()
    {
        return m_block->info.succ;
    }
    bid_type & pred()
    {
        return m_block->info.pred;
    }

    const bid_type & succ() const
    {
        return m_block->info.succ;
    }
    const bid_type & pred() const
    {
        return m_block->info.pred;
    }

    /*
       template <class InputIterator>
       normal_leaf(InputIterator begin_, InputIterator end_,
            btree_type * btree,
            key_compare cmp):
            m_block(new block_type),
            m_btree(btree),
            m_cmp(cmp),
            m_vcmp(cmp)
       {
            assert(min_nelements() >=2);
            assert(2*min_nelements() - 1 <= max_nelements());
            assert(max_nelements() <= nelements);
            assert(unsigned(block_type::size) >= nelements +1); // extra space for an overflow element

            unsigned new_size = end_ - begin_;
            assert(new_size <= max_nelements());
            assert(new_size >= min_nelements());

            std::copy(begin_,end_,m_block->begin());
            assert(stxxl::is_sorted(m_block->begin(),m_block->begin() + new_size, m_vcmp));
            m_block->info.cur_size = new_size;
       }*/

    unsigned size() const
    {
        return m_block->info.cur_size;
    }

    const bid_type & my_bid() const
    {
        return m_block->info.me;
    }

    void save()
    {
        request_ptr req = m_block->write(my_bid());
        req->wait();
    }

    request_ptr load(const bid_type& bid)
    {
        request_ptr req = m_block->read(bid);
        req->wait();
        assert(bid == my_bid());
        return req;
    }

    request_ptr prefetch(const bid_type& bid)
    {
        return m_block->read(bid);
    }

    void init(const bid_type& my_bid_)
    {
        m_block->info.me = my_bid_;
        m_block->info.succ = bid_type();
        m_block->info.pred = bid_type();
        m_block->info.cur_size = 0;
    }

    reference operator [] (unsigned_type i)
    {
        return (*m_block)[i];
    }

    const_reference operator [] (unsigned_type i) const
    {
        return (*m_block)[i];
    }

    reference back()
    {
        return (*m_block)[size() - 1];
    }

    reference front()
    {
        return *(m_block->begin());
    }

    const_reference back() const
    {
        return (*m_block)[size() - 1];
    }

    const_reference front() const
    {
        return *(m_block->begin());
    }

    void dump()
    {
        STXXL_VERBOSE2("Dump of leaf " << this);
        for (unsigned i = 0; i < size(); ++i)
            STXXL_VERBOSE2((*this)[i].first << " " << (*this)[i].second);
    }

    std::pair<iterator, bool> insert(
        const value_type& x,
        std::pair<key_type, bid_type>& splitter)
    {
        assert(size() <= max_nelements());
        splitter.first = key_compare::max_value();

        typename block_type::iterator it =
            std::lower_bound(m_block->begin(), m_block->begin() + size(), x, m_vcmp);

        if (!(m_vcmp(*it, x) || m_vcmp(x, *it)) && it != (m_block->begin() + size()))                    // *it == x
        {
            // already exists
            return std::pair<iterator, bool>(
                iterator(m_btree, my_bid(), unsigned(it - m_block->begin())),
                false);
        }

        typename block_type::iterator cur = m_block->begin() + size() - 1;

        for ( ; cur >= it; --cur)
            *(cur + 1) = *cur;
        // move elements to make space for the new element

        *it = x;

        std::vector<iterator_base*> iterators2fix;
        m_btree->m_iterator_map.find(my_bid(), unsigned(it - m_block->begin()), size(), iterators2fix);
        typename std::vector<iterator_base*>::iterator it2fix = iterators2fix.begin();
        for ( ; it2fix != iterators2fix.end(); ++it2fix)
        {
            m_btree->m_iterator_map.unregister_iterator(**it2fix);
            ++((*it2fix)->pos);                          // fixing iterators
            m_btree->m_iterator_map.register_iterator(**it2fix);
        }

        ++(m_block->info.cur_size);

        std::pair<iterator, bool> result(iterator(m_btree, my_bid(), unsigned(it - m_block->begin())), true);

        if (size() <= max_nelements())
        {
            // no overflow
            dump();
            return result;
        }

        // overflow

        split(splitter);

        return result;
    }

    iterator begin()
    {
        return iterator(m_btree, my_bid(), 0);
    }

    const_iterator begin() const
    {
        return const_iterator(m_btree, my_bid(), 0);
    }

    iterator end()
    {
        return iterator(m_btree, my_bid(), size());
    }

    void increment_iterator(iterator_base& it) const
    {
        assert(it.bid == my_bid());
        assert(it.pos != size());

        m_btree->m_iterator_map.unregister_iterator(it);

        ++(it.pos);
        if (it.pos == size() && succ().valid())
        {
            // run to the end of the leaf
            STXXL_VERBOSE1("btree::normal_leaf jumping to the next block");
            it.pos = 0;
            it.bid = succ();
        }
        // increment of pos from 0 to 1
        else if (it.pos == 1 && m_btree->m_prefetching_enabled)
        {
            // prefetch the succ leaf
            if (succ().valid())
                m_btree->m_leaf_cache.prefetch_node(succ());
        }
        m_btree->m_iterator_map.register_iterator(it);
    }

    void decrement_iterator(iterator_base& it) const
    {
        assert(it.bid == my_bid());

        m_btree->m_iterator_map.unregister_iterator(it);

        if (it.pos == 0)
        {
            assert(pred().valid());

            it.bid = pred();
            normal_leaf const* pred_leaf = m_btree->m_leaf_cache.get_const_node(pred(), true);
            assert(pred_leaf);
            it.pos = pred_leaf->size() - 1;

            // prefetch the pred leaf of pred_leaf
            if (m_btree->m_prefetching_enabled && pred_leaf->pred().valid())
                m_btree->m_leaf_cache.prefetch_node(pred_leaf->pred());

            m_btree->m_leaf_cache.unfix_node(pred());
        }
        else
            --it.pos;

        m_btree->m_iterator_map.register_iterator(it);
    }

    iterator find(const key_type& k)
    {
        value_type search_val(k, data_type());
        typename block_type::iterator lb =
            std::lower_bound(m_block->begin(), m_block->begin() + size(), search_val, m_vcmp);
        if (lb == m_block->begin() + size() || lb->first != k)
            return m_btree->end();

        return iterator(m_btree, my_bid(), unsigned(lb - m_block->begin()));
    }

    const_iterator find(const key_type& k) const
    {
        value_type search_val(k, data_type());
        typename block_type::iterator lb =
            std::lower_bound(m_block->begin(), m_block->begin() + size(), search_val, m_vcmp);
        if (lb == m_block->begin() + size() || lb->first != k)
            return m_btree->end();

        return const_iterator(m_btree, my_bid(), unsigned(lb - m_block->begin()));
    }

    iterator lower_bound(const key_type& k)
    {
        value_type search_val(k, data_type());

        typename block_type::iterator lb =
            std::lower_bound(m_block->begin(), m_block->begin() + size(), search_val, m_vcmp);

        // lower_bound is in the succ block
        if (lb == m_block->begin() + size() && succ().valid())
        {
            return iterator(m_btree, succ(), 0);
        }

        return iterator(m_btree, my_bid(), unsigned(lb - m_block->begin()));
    }

    const_iterator lower_bound(const key_type& k) const
    {
        value_type search_val(k, data_type());
        typename block_type::iterator lb =
            std::lower_bound(m_block->begin(), m_block->begin() + size(), search_val, m_vcmp);

        // lower_bound is in the succ block
        if (lb == m_block->begin() + size() && succ().valid())
        {
            return iterator(m_btree, succ(), 0);
        }

        return const_iterator(m_btree, my_bid(), unsigned(lb - m_block->begin()));
    }

    iterator upper_bound(const key_type& k)
    {
        value_type search_val(k, data_type());
        typename block_type::iterator lb =
            std::upper_bound(m_block->begin(), m_block->begin() + size(), search_val, m_vcmp);

        // upper_bound is in the succ block
        if (lb == m_block->begin() + size() && succ().valid())
        {
            return iterator(m_btree, succ(), 0);
        }

        return iterator(m_btree, my_bid(), unsigned(lb - m_block->begin()));
    }

    const_iterator upper_bound(const key_type& k) const
    {
        value_type search_val(k, data_type());
        typename block_type::iterator lb =
            std::upper_bound(m_block->begin(), m_block->begin() + size(), search_val, m_vcmp);

        // upper_bound is in the succ block
        if (lb == m_block->begin() + size() && succ().valid())
        {
            return const_iterator(m_btree, succ(), 0);
        }

        return const_iterator(m_btree, my_bid(), unsigned(lb - m_block->begin()));
    }

    size_type erase(const key_type& k)
    {
        value_type search_val(k, data_type());
        typename block_type::iterator it =
            std::lower_bound(m_block->begin(), m_block->begin() + size(), search_val, m_vcmp);

        if (it == m_block->begin() + size() || it->first != k)
            return 0;
        // no such element

        // move elements one position left
        std::copy(it + 1, m_block->begin() + size(), it);

        std::vector<iterator_base*> iterators2fix;
        m_btree->m_iterator_map.find(my_bid(), unsigned(it + 1 - m_block->begin()), size(), iterators2fix);
        typename std::vector<iterator_base*>::iterator it2fix = iterators2fix.begin();
        for ( ; it2fix != iterators2fix.end(); ++it2fix)
        {
            STXXL_VERBOSE2("btree::normal_leaf updating iterator " << (*it2fix) << " (pos--)");
            m_btree->m_iterator_map.unregister_iterator(**it2fix);
            --((*it2fix)->pos);                          // fixing iterators
            m_btree->m_iterator_map.register_iterator(**it2fix);
        }

        --(m_block->info.cur_size);

        return 1;
    }

    void fuse(const normal_leaf& src)
    {
        STXXL_VERBOSE1("btree::normal_leaf Fusing");
        assert(m_vcmp(src.back(), front()));
        const unsigned src_size = src.size();

        typename block_type::iterator cur = m_block->begin() + size() - 1;
        typename block_type::const_iterator begin = m_block->begin();

        for ( ; cur >= begin; --cur)
            *(cur + src_size) = *cur;
        // move elements to make space for Src elements

        // copy Src to *this leaf
        std::copy(src.m_block->begin(), src.m_block->begin() + src_size, m_block->begin());

        std::vector<iterator_base*> iterators2fix;
        m_btree->m_iterator_map.find(my_bid(), 0, size(), iterators2fix);
        typename std::vector<iterator_base*>::iterator it2fix = iterators2fix.begin();
        for ( ; it2fix != iterators2fix.end(); ++it2fix)
        {
            STXXL_VERBOSE2("btree::normal_leaf updating iterator " << (*it2fix) <<
                           " (pos+" << src_size << ")");
            m_btree->m_iterator_map.unregister_iterator(**it2fix);
            ((*it2fix)->pos) += src_size;                           // fixing iterators
            m_btree->m_iterator_map.register_iterator(**it2fix);
        }

        iterators2fix.clear();
        m_btree->m_iterator_map.find(src.my_bid(), 0, src_size, iterators2fix);
        for (it2fix = iterators2fix.begin(); it2fix != iterators2fix.end(); ++it2fix)
        {
            STXXL_VERBOSE2("btree::normal_leaf updating iterator " << (*it2fix) <<
                           " (bid=" << my_bid() << ")");
            m_btree->m_iterator_map.unregister_iterator(**it2fix);
            ((*it2fix)->bid) = my_bid();                             // fixing iterators
            m_btree->m_iterator_map.register_iterator(**it2fix);
        }

        m_block->info.cur_size += src_size;

        // update links
        pred() = src.pred();
        if (pred().valid())
        {                         // update successor link
            normal_leaf* new_pred = m_btree->m_leaf_cache.get_node(pred());
            assert(new_pred);
            new_pred->succ() = my_bid();
        }
    }

    key_type balance(normal_leaf& left)
    {
        STXXL_VERBOSE1("btree::normal_leaf Balancing leaves with bids " <<
                       left.my_bid() << " and " << my_bid());
        const unsigned total_size = left.size() + size();
        unsigned new_left_size = total_size / 2;
        assert(new_left_size <= left.max_nelements());
        assert(new_left_size >= left.min_nelements());
        unsigned new_right_size = total_size - new_left_size;
        assert(new_right_size <= max_nelements());
        assert(new_right_size >= min_nelements());

        assert(m_vcmp(left.back(), front()) || size() == 0);

        if (new_left_size < left.size())
        {
            // #elements to move from left to *this
            const unsigned nEl2Move = left.size() - new_left_size;

            typename block_type::iterator cur = m_block->begin() + size() - 1;
            typename block_type::const_iterator begin = m_block->begin();

            for ( ; cur >= begin; --cur)
                *(cur + nEl2Move) = *cur;
            // move elements to make space for Src elements

            // copy Left to *this leaf
            std::copy(left.m_block->begin() + new_left_size,
                      left.m_block->begin() + left.size(), m_block->begin());

            std::vector<iterator_base*> iterators2fix1;
            std::vector<iterator_base*> iterators2fix2;
            m_btree->m_iterator_map.find(my_bid(), 0, size(), iterators2fix1);
            m_btree->m_iterator_map.find(left.my_bid(), new_left_size, left.size(), iterators2fix2);

            typename std::vector<iterator_base*>::iterator it2fix = iterators2fix1.begin();
            for ( ; it2fix != iterators2fix1.end(); ++it2fix)
            {
                STXXL_VERBOSE2("btree::normal_leaf updating iterator " << (*it2fix) <<
                               " (pos+" << nEl2Move << ")");
                m_btree->m_iterator_map.unregister_iterator(**it2fix);
                ((*it2fix)->pos) += nEl2Move;                               // fixing iterators
                m_btree->m_iterator_map.register_iterator(**it2fix);
            }

            it2fix = iterators2fix2.begin();
            for ( ; it2fix != iterators2fix2.end(); ++it2fix)
            {
                STXXL_VERBOSE2("btree::normal_leaf updating iterator " << (*it2fix) <<
                               " (pos-" << new_left_size << " bid=" << my_bid() << ")");
                m_btree->m_iterator_map.unregister_iterator(**it2fix);
                ((*it2fix)->bid) = my_bid();                                     // fixing iterators
                ((*it2fix)->pos) -= new_left_size;                               // fixing iterators
                m_btree->m_iterator_map.register_iterator(**it2fix);
            }
        }
        else
        {
            assert(new_right_size < size());

            const unsigned nEl2Move = size() - new_right_size;                            // #elements to move from *this to Left

            // copy *this to Left
            std::copy(m_block->begin(),
                      m_block->begin() + nEl2Move, left.m_block->begin() + left.size());
            // move elements in *this
            std::copy(m_block->begin() + nEl2Move,
                      m_block->begin() + size(), m_block->begin());

            std::vector<iterator_base*> iterators2fix1;
            std::vector<iterator_base*> iterators2fix2;
            m_btree->m_iterator_map.find(my_bid(), nEl2Move, size(), iterators2fix1);
            m_btree->m_iterator_map.find(my_bid(), 0, nEl2Move - 1, iterators2fix2);

            typename std::vector<iterator_base*>::iterator it2fix = iterators2fix1.begin();
            for ( ; it2fix != iterators2fix1.end(); ++it2fix)
            {
                STXXL_VERBOSE2("btree::normal_leaf updating iterator " << (*it2fix) <<
                               " (pos-" << nEl2Move << ")");
                m_btree->m_iterator_map.unregister_iterator(**it2fix);
                ((*it2fix)->pos) -= nEl2Move;                                 // fixing iterators
                m_btree->m_iterator_map.register_iterator(**it2fix);
            }

            it2fix = iterators2fix2.begin();
            for ( ; it2fix != iterators2fix2.end(); ++it2fix)
            {
                STXXL_VERBOSE2("btree::normal_leaf updating iterator " << (*it2fix) <<
                               " (pos+" << left.size() << " bid=" << left.my_bid() << ")");
                m_btree->m_iterator_map.unregister_iterator(**it2fix);
                ((*it2fix)->bid) = left.my_bid();                                 // fixing iterators
                ((*it2fix)->pos) += left.size();                                  // fixing iterators
                m_btree->m_iterator_map.register_iterator(**it2fix);
            }
        }

        m_block->info.cur_size = new_right_size;                             // update size
        left.m_block->info.cur_size = new_left_size;                         // update size

        return left.back().first;
    }

    void push_back(const value_type& x)
    {
        (*this)[size()] = x;
        ++(m_block->info.cur_size);
    }
};

} // namespace btree

STXXL_END_NAMESPACE

#endif // !STXXL_CONTAINERS_BTREE_LEAF_HEADER
