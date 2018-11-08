/***************************************************************************
 *  include/stxxl/bits/parallel/equally_split.h
 *
 *  Function to split a sequence into parts of almost equal size.
 *  Extracted from MCSTL - http://algo2.iti.uni-karlsruhe.de/singler/mcstl/
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2007 Johannes Singler <singler@ira.uka.de>
 *  Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_PARALLEL_EQUALLY_SPLIT_HEADER
#define STXXL_PARALLEL_EQUALLY_SPLIT_HEADER

#include <stxxl/bits/namespace.h>
#include <stxxl/bits/parallel/base.h>
#include <stxxl/bits/parallel/types.h>

STXXL_BEGIN_NAMESPACE

namespace parallel {

/*!
 * Split a sequence into parts of almost equal size.
 *
 * The resulting sequence s of length p+1 contains the splitting positions when
 * splitting the range [0,n) into parts of almost equal size (plus minus 1).
 * The first entry is 0, the last one n. There may result empty parts.
 *
 * \param n Number of elements
 * \param p Number of parts
 * \param s Splitters
 * \returns End of splitter sequence, i. e. \c s+p+1
 */
template <typename DiffType, typename DiffTypeOutputIterator>
DiffTypeOutputIterator equally_split(DiffType n, thread_index_t p,
                                     DiffTypeOutputIterator s)
{
    DiffType chunk_length = n / p, split = n % p, start = 0;
    for (thread_index_t i = 0; i < p; i++)
    {
        *s++ = start;
        start += ((DiffType)i < split) ? (chunk_length + 1) : chunk_length;
    }
    *s++ = n;

    return s;
}

} // namespace parallel

STXXL_END_NAMESPACE

#endif // !STXXL_PARALLEL_EQUALLY_SPLIT_HEADER
