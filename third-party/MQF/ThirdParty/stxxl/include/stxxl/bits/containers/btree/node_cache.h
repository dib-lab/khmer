/***************************************************************************
 *  include/stxxl/bits/containers/btree/node_cache.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2006 Roman Dementiev <dementiev@ira.uka.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_CONTAINERS_BTREE_NODE_CACHE_HEADER
#define STXXL_CONTAINERS_BTREE_NODE_CACHE_HEADER

#include <stxxl/bits/config.h>
#include <stxxl/bits/compat/hash_map.h>
#include <stxxl/bits/io/request.h>
#include <stxxl/bits/mng/block_manager.h>
#include <stxxl/bits/mng/typed_block.h>
#include <stxxl/bits/containers/pager.h>
#include <stxxl/bits/common/error_handling.h>

STXXL_BEGIN_NAMESPACE

#define STXXL_BTREE_CACHE_VERBOSE STXXL_VERBOSE2

// TODO:  speedup BID2node_ access using search result iterator in the methods

namespace btree {

template <class NodeType, class BTreeType>
class node_cache : private noncopyable
{
public:
    typedef BTreeType btree_type;
    typedef NodeType node_type;
    typedef typename node_type::block_type block_type;
    typedef typename node_type::bid_type bid_type;
    typedef typename btree_type::key_compare key_compare;

    typedef typename btree_type::alloc_strategy_type alloc_strategy_type;
    typedef stxxl::lru_pager<> pager_type;

private:
    btree_type* m_btree;
    key_compare m_cmp;

/*
        struct bid_comp
        {
            bool operator ()  (const bid_type & a, const bid_type & b) const
            {
                return (a.storage < b.storage) || ( a.storage == b.storage && a.offset < b.offset);
            }
        };
*/

    struct bid_hash
    {
        size_t operator () (const bid_type& bid) const
        {
            size_t result =
                longhash1(bid.offset + reinterpret_cast<uint64>(bid.storage));
            return result;
        }
#if STXXL_MSVC
        bool operator () (const bid_type& a, const bid_type& b) const
        {
            return (a.storage < b.storage) || (a.storage == b.storage && a.offset < b.offset);
        }
        enum
        {                                       // parameters for hash table
            bucket_size = 4,                    // 0 < bucket_size
            min_buckets = 8                     // min_buckets = 2 ^^ N, 0 < N
        };
#endif
    };

    std::vector<node_type*> m_nodes;
    std::vector<request_ptr> m_reqs;
    std::vector<bool> m_fixed;
    std::vector<bool> m_dirty;
    std::vector<int_type> m_free_nodes;
    typedef typename compat_hash_map<bid_type, int_type, bid_hash>::result hash_map_type;

    //typedef std::map<bid_type,int_type,bid_comp> BID2node_type;
    typedef hash_map_type bid2node_type;

    bid2node_type m_bid2node;
    pager_type m_pager;
    block_manager* m_bm;
    alloc_strategy_type m_alloc_strategy;

    int64 n_found;
    int64 n_not_found;
    int64 n_created;
    int64 n_deleted;
    int64 n_read;
    int64 n_written;
    int64 n_clean_forced;

    // changes btree pointer in all contained iterators
    void change_btree_pointers(btree_type* b)
    {
        for (typename std::vector<node_type*>::const_iterator it = m_nodes.begin();
             it != m_nodes.end(); ++it)
        {
            (*it)->m_btree = b;
        }
    }

public:
    node_cache(unsigned_type cache_size_in_bytes,
               btree_type* btree,
               key_compare cmp)
        : m_btree(btree),
          m_cmp(cmp),
          m_bm(block_manager::get_instance()),
          n_found(0),
          n_not_found(0),
          n_created(0),
          n_deleted(0),
          n_read(0),
          n_written(0),
          n_clean_forced(0)
    {
        const unsigned_type nnodes = cache_size_in_bytes / block_type::raw_size;
        STXXL_BTREE_CACHE_VERBOSE("btree::node_cache constructor nodes=" << nnodes);
        if (nnodes < 3)
        {
            STXXL_THROW2(std::runtime_error, "btree::node_cache::node_cache", "Too few memory for a node cache (<3)");
        }
        m_nodes.reserve(nnodes);
        m_reqs.resize(nnodes);
        m_free_nodes.reserve(nnodes);
        m_fixed.resize(nnodes, false);
        m_dirty.resize(nnodes, true);
        for (unsigned_type i = 0; i < nnodes; ++i)
        {
            m_nodes.push_back(new node_type(m_btree, m_cmp));
            m_free_nodes.push_back(i);
        }

        pager_type tmp_pager(nnodes);
        std::swap(m_pager, tmp_pager);
    }

    unsigned_type size() const
    {
        return m_nodes.size();
    }

    // returns the number of fixed pages
    unsigned_type nfixed() const
    {
        typename bid2node_type::const_iterator i = m_bid2node.begin();
        typename bid2node_type::const_iterator end = m_bid2node.end();
        unsigned_type cnt = 0;
        for ( ; i != end; ++i)
        {
            if (m_fixed[(*i).second])
                ++cnt;
        }
        return cnt;
    }

    ~node_cache()
    {
        STXXL_BTREE_CACHE_VERBOSE("btree::node_cache destructor addr=" << this);
        typename bid2node_type::const_iterator i = m_bid2node.begin();
        typename bid2node_type::const_iterator end = m_bid2node.end();
        for ( ; i != end; ++i)
        {
            const unsigned_type p = (*i).second;
            if (m_reqs[p].valid())
                m_reqs[p]->wait();

            if (m_dirty[p])
                m_nodes[p]->save();
        }
        for (unsigned_type i = 0; i < size(); ++i)
            delete m_nodes[i];
    }

    node_type * get_new_node(bid_type& new_bid)
    {
        ++n_created;

        if (m_free_nodes.empty())
        {
            // need to kick a node
            int_type node2kick;
            unsigned_type i = 0;
            const unsigned_type max_tries = size() + 1;
            do
            {
                ++i;
                node2kick = m_pager.kick();
                if (i == max_tries)
                {
                    STXXL_ERRMSG(
                        "The node cache is too small, no node can be kicked out (all nodes are fixed) !");
                    STXXL_ERRMSG("Returning NULL node.");
                    return NULL;
                }
                m_pager.hit(node2kick);
            } while (m_fixed[node2kick]);

            if (m_reqs[node2kick].valid())
                m_reqs[node2kick]->wait();

            node_type& node = *(m_nodes[node2kick]);

            if (m_dirty[node2kick])
            {
                node.save();
                ++n_written;
            }
            else
                ++n_clean_forced;

            //reqs_[node2kick] = request_ptr(); // reset request

            assert(m_bid2node.find(node.my_bid()) != m_bid2node.end());
            m_bid2node.erase(node.my_bid());
            m_bm->new_block(m_alloc_strategy, new_bid);

            m_bid2node[new_bid] = node2kick;

            node.init(new_bid);

            m_dirty[node2kick] = true;

            assert(size() == m_bid2node.size() + m_free_nodes.size());

            STXXL_BTREE_CACHE_VERBOSE("btree::node_cache get_new_node, need to kick node " << node2kick);

            return &node;
        }

        int_type free_node = m_free_nodes.back();
        m_free_nodes.pop_back();
        assert(m_fixed[free_node] == false);

        m_bm->new_block(m_alloc_strategy, new_bid);
        m_bid2node[new_bid] = free_node;
        node_type& node = *(m_nodes[free_node]);
        node.init(new_bid);

        // assert(!(reqs_[free_node].valid()));

        m_pager.hit(free_node);

        m_dirty[free_node] = true;

        assert(size() == m_bid2node.size() + m_free_nodes.size());

        STXXL_BTREE_CACHE_VERBOSE("btree::node_cache get_new_node, free node " << free_node << "available");

        return &node;
    }

    node_type * get_node(const bid_type& bid, bool fix = false)
    {
        typename bid2node_type::const_iterator it = m_bid2node.find(bid);
        ++n_read;

        if (it != m_bid2node.end())
        {
            // the node is in cache
            const int_type nodeindex = it->second;
            STXXL_BTREE_CACHE_VERBOSE("btree::node_cache get_node, the node " << nodeindex << "is in cache , fix=" << fix);
            m_fixed[nodeindex] = fix;
            m_pager.hit(nodeindex);
            m_dirty[nodeindex] = true;

            if (m_reqs[nodeindex].valid() && !m_reqs[nodeindex]->poll())
                m_reqs[nodeindex]->wait();

            ++n_found;
            return m_nodes[nodeindex];
        }

        ++n_not_found;

        // the node is not in cache
        if (m_free_nodes.empty())
        {
            // need to kick a node
            int_type node2kick;
            unsigned_type i = 0;
            const unsigned_type max_tries = size() + 1;
            do
            {
                ++i;
                node2kick = m_pager.kick();
                if (i == max_tries)
                {
                    STXXL_ERRMSG(
                        "The node cache is too small, no node can be kicked out (all nodes are fixed) !");
                    STXXL_ERRMSG("Returning NULL node.");
                    return NULL;
                }
                m_pager.hit(node2kick);
            } while (m_fixed[node2kick]);

            if (m_reqs[node2kick].valid())
                m_reqs[node2kick]->wait();

            node_type& node = *(m_nodes[node2kick]);

            if (m_dirty[node2kick])
            {
                node.save();
                ++n_written;
            }
            else
                ++n_clean_forced;

            m_bid2node.erase(node.my_bid());

            m_reqs[node2kick] = node.load(bid);
            m_bid2node[bid] = node2kick;

            m_fixed[node2kick] = fix;

            m_dirty[node2kick] = true;

            assert(size() == m_bid2node.size() + m_free_nodes.size());

            STXXL_BTREE_CACHE_VERBOSE("btree::node_cache get_node, need to kick node" << node2kick << " fix=" << fix);

            return &node;
        }

        int_type free_node = m_free_nodes.back();
        m_free_nodes.pop_back();
        assert(m_fixed[free_node] == false);

        node_type& node = *(m_nodes[free_node]);
        m_reqs[free_node] = node.load(bid);
        m_bid2node[bid] = free_node;

        m_pager.hit(free_node);

        m_fixed[free_node] = fix;

        m_dirty[free_node] = true;

        assert(size() == m_bid2node.size() + m_free_nodes.size());

        STXXL_BTREE_CACHE_VERBOSE("btree::node_cache get_node, free node " << free_node << "available, fix=" << fix);

        return &node;
    }

    node_type const * get_const_node(const bid_type& bid, bool fix = false)
    {
        typename bid2node_type::const_iterator it = m_bid2node.find(bid);
        ++n_read;

        if (it != m_bid2node.end())
        {
            // the node is in cache
            const int_type nodeindex = it->second;
            STXXL_BTREE_CACHE_VERBOSE("btree::node_cache get_node, the node " << nodeindex << "is in cache , fix=" << fix);
            m_fixed[nodeindex] = fix;
            m_pager.hit(nodeindex);

            if (m_reqs[nodeindex].valid() && !m_reqs[nodeindex]->poll())
                m_reqs[nodeindex]->wait();

            ++n_found;
            return m_nodes[nodeindex];
        }

        ++n_not_found;

        // the node is not in cache
        if (m_free_nodes.empty())
        {
            // need to kick a node
            int_type node2kick;
            unsigned_type i = 0;
            const unsigned_type max_tries = size() + 1;
            do
            {
                ++i;
                node2kick = m_pager.kick();
                if (i == max_tries)
                {
                    STXXL_ERRMSG(
                        "The node cache is too small, no node can be kicked out (all nodes are fixed) !");
                    STXXL_ERRMSG("Returning NULL node.");
                    return NULL;
                }
                m_pager.hit(node2kick);
            } while (m_fixed[node2kick]);

            if (m_reqs[node2kick].valid())
                m_reqs[node2kick]->wait();

            node_type& node = *(m_nodes[node2kick]);
            if (m_dirty[node2kick])
            {
                node.save();
                ++n_written;
            }
            else
                ++n_clean_forced;

            m_bid2node.erase(node.my_bid());

            m_reqs[node2kick] = node.load(bid);
            m_bid2node[bid] = node2kick;

            m_fixed[node2kick] = fix;

            m_dirty[node2kick] = false;

            assert(size() == m_bid2node.size() + m_free_nodes.size());

            STXXL_BTREE_CACHE_VERBOSE("btree::node_cache get_node, need to kick node" << node2kick << " fix=" << fix);

            return &node;
        }

        int_type free_node = m_free_nodes.back();
        m_free_nodes.pop_back();
        assert(m_fixed[free_node] == false);

        node_type& node = *(m_nodes[free_node]);
        m_reqs[free_node] = node.load(bid);
        m_bid2node[bid] = free_node;

        m_pager.hit(free_node);

        m_fixed[free_node] = fix;

        m_dirty[free_node] = false;

        assert(size() == m_bid2node.size() + m_free_nodes.size());

        STXXL_BTREE_CACHE_VERBOSE("btree::node_cache get_node, free node " << free_node << "available, fix=" << fix);

        return &node;
    }

    void delete_node(const bid_type& bid)
    {
        typename bid2node_type::const_iterator it = m_bid2node.find(bid);
        try
        {
            if (it != m_bid2node.end())
            {
                // the node is in the cache
                const int_type nodeindex = it->second;
                STXXL_BTREE_CACHE_VERBOSE("btree::node_cache delete_node " << nodeindex << " from cache.");
                if (m_reqs[nodeindex].valid())
                    m_reqs[nodeindex]->wait();

                //reqs_[nodeindex] = request_ptr(); // reset request
                m_free_nodes.push_back(nodeindex);
                m_bid2node.erase(bid);
                m_fixed[nodeindex] = false;
            }
            ++n_deleted;
        } catch (const io_error& ex)
        {
            m_bm->delete_block(bid);
            throw io_error(ex.what());
        }
        m_bm->delete_block(bid);
    }

    void prefetch_node(const bid_type& bid)
    {
        if (m_bid2node.find(bid) != m_bid2node.end())
            return;

        // the node is not in cache
        if (m_free_nodes.empty())
        {
            // need to kick a node
            int_type node2kick;
            unsigned_type i = 0;
            const unsigned_type max_tries = size() + 1;
            do
            {
                ++i;
                node2kick = m_pager.kick();
                if (i == max_tries)
                {
                    STXXL_ERRMSG(
                        "The node cache is too small, no node can be kicked out (all nodes are fixed) !");
                    STXXL_ERRMSG("Returning NULL node.");
                    return;
                }
                m_pager.hit(node2kick);
            } while (m_fixed[node2kick]);

            if (m_reqs[node2kick].valid())
                m_reqs[node2kick]->wait();

            node_type& node = *(m_nodes[node2kick]);

            if (m_dirty[node2kick])
            {
                node.save();
                ++n_written;
            }
            else
                ++n_clean_forced;

            m_bid2node.erase(node.my_bid());

            m_reqs[node2kick] = node.prefetch(bid);
            m_bid2node[bid] = node2kick;

            m_fixed[node2kick] = false;

            m_dirty[node2kick] = false;

            assert(size() == m_bid2node.size() + m_free_nodes.size());

            STXXL_BTREE_CACHE_VERBOSE("btree::node_cache prefetch_node, need to kick node" << node2kick << " ");

            return;
        }

        int_type free_node = m_free_nodes.back();
        m_free_nodes.pop_back();
        assert(m_fixed[free_node] == false);

        node_type& node = *(m_nodes[free_node]);
        m_reqs[free_node] = node.prefetch(bid);
        m_bid2node[bid] = free_node;

        m_pager.hit(free_node);

        m_fixed[free_node] = false;

        m_dirty[free_node] = false;

        assert(size() == m_bid2node.size() + m_free_nodes.size());

        STXXL_BTREE_CACHE_VERBOSE("btree::node_cache prefetch_node, free node " << free_node << "available");

        return;
    }

    void unfix_node(const bid_type& bid)
    {
        assert(m_bid2node.find(bid) != m_bid2node.end());
        m_fixed[m_bid2node[bid]] = false;
        STXXL_BTREE_CACHE_VERBOSE("btree::node_cache unfix_node,  node " << m_bid2node[bid]);
    }

    void swap(node_cache& obj)
    {
        std::swap(m_cmp, obj.m_cmp);
        std::swap(m_nodes, obj.m_nodes);
        std::swap(m_reqs, obj.m_reqs);
        change_btree_pointers(m_btree);
        obj.change_btree_pointers(obj.m_btree);
        std::swap(m_fixed, obj.m_fixed);
        std::swap(m_free_nodes, obj.m_free_nodes);
        std::swap(m_bid2node, obj.m_bid2node);
        std::swap(m_pager, obj.m_pager);
        std::swap(m_alloc_strategy, obj.m_alloc_strategy);
        std::swap(n_found, obj.n_found);
        std::swap(n_not_found, obj.n_found);
        std::swap(n_created, obj.n_created);
        std::swap(n_deleted, obj.n_deleted);
        std::swap(n_read, obj.n_read);
        std::swap(n_written, obj.n_written);
        std::swap(n_clean_forced, obj.n_clean_forced);
    }

    void print_statistics(std::ostream& o) const
    {
        if (n_read)
            o << "Found blocks                      : " << n_found << " (" <<
                100. * double(n_found) / double(n_read) << "%)" << std::endl;
        else
            o << "Found blocks                      : " << n_found << " (" <<
                100 << "%)" << std::endl;

        o << "Not found blocks                  : " << n_not_found << std::endl;
        o << "Created in the cache blocks       : " << n_created << std::endl;
        o << "Deleted blocks                    : " << n_deleted << std::endl;
        o << "Read blocks                       : " << n_read << std::endl;
        o << "Written blocks                    : " << n_written << std::endl;
        o << "Clean blocks forced from the cache: " << n_clean_forced << std::endl;
    }
    void reset_statistics()
    {
        n_found = 0;
        n_not_found = 0;
        n_created = 0;
        n_deleted = 0;
        n_read = 0;
        n_written = 0;
        n_clean_forced = 0;
    }
};

} // namespace btree

STXXL_END_NAMESPACE

namespace std {

template <class NodeType, class BTreeType>
void swap(stxxl::btree::node_cache<NodeType, BTreeType>& a,
          stxxl::btree::node_cache<NodeType, BTreeType>& b)
{
    a.swap(b);
}

} // namespace std

#endif // !STXXL_CONTAINERS_BTREE_NODE_CACHE_HEADER
