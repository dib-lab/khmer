/***************************************************************************
 *  include/stxxl/bits/algo/run_cursor.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2003 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_ALGO_RUN_CURSOR_HEADER
#define STXXL_ALGO_RUN_CURSOR_HEADER

#include <cstdlib>
#include <stxxl/bits/common/types.h>

STXXL_BEGIN_NAMESPACE

template <typename BlockType>
struct run_cursor
{
    unsigned_type pos;
    BlockType* buffer;

    run_cursor() : pos(0), buffer(NULL) { }

    inline typename BlockType::const_reference current() const
    {
        return (*buffer)[pos];
    }
    inline void operator ++ ()
    {
        ++pos;
    }
};

#ifdef STXXL_SORT_SINGLE_PREFETCHER

template <typename MustBeVoid = void>
struct have_prefetcher
{
    static void* untyped_prefetcher;
};

#endif

template <typename BlockType,
          typename PrefetcherType>
struct run_cursor2 : public run_cursor<BlockType>
#ifdef STXXL_SORT_SINGLE_PREFETCHER
                     , public have_prefetcher<>
#endif
{
    typedef BlockType block_type;
    typedef PrefetcherType prefetcher_type;
    typedef run_cursor2<block_type, prefetcher_type> _Self;
    typedef typename block_type::value_type value_type;

    using run_cursor<block_type>::pos;
    using run_cursor<block_type>::buffer;

#ifdef STXXL_SORT_SINGLE_PREFETCHER
    static prefetcher_type* const prefetcher()     // sorry, a hack
    {
        return reinterpret_cast<prefetcher_type*>(untyped_prefetcher);
    }
    static void set_prefetcher(prefetcher_type* pfptr)
    {
        untyped_prefetcher = pfptr;
    }
    run_cursor2() { }
#else
    prefetcher_type* prefetcher_;
    prefetcher_type* & prefetcher()  // sorry, a hack
    {
        return prefetcher_;
    }

    run_cursor2(prefetcher_type* p = NULL) : prefetcher_(p) { }
#endif

    inline bool empty() const
    {
        return (pos >= block_type::size);
    }
    inline void operator ++ ()
    {
        assert(!empty());
        ++pos;
        if (UNLIKELY(pos >= block_type::size))
        {
            if (prefetcher()->block_consumed(buffer))
                pos = 0;
        }
    }
    inline void make_inf()
    {
        pos = block_type::size;
    }
};

#ifdef STXXL_SORT_SINGLE_PREFETCHER
template <typename MustBeVoid>
void* have_prefetcher<MustBeVoid>::untyped_prefetcher = NULL;
#endif

#if 0
template <typename block_type>
struct run_cursor_cmp
{
    typedef run_cursor<block_type> cursor_type;

    inline bool operator () (const cursor_type& a, const cursor_type& b)        // greater or equal
    {
        return !((*a.buffer)[a.pos] < (*b.buffer)[b.pos]);
    }
};
#endif

STXXL_END_NAMESPACE

#endif // !STXXL_ALGO_RUN_CURSOR_HEADER
// vim: et:ts=4:sw=4
