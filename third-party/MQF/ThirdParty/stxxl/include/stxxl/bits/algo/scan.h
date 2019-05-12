/***************************************************************************
 *  include/stxxl/bits/algo/scan.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002-2004 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2008, 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_ALGO_SCAN_HEADER
#define STXXL_ALGO_SCAN_HEADER

#include <stxxl/bits/namespace.h>
#include <stxxl/bits/mng/config.h>
#include <stxxl/bits/mng/buf_istream.h>
#include <stxxl/bits/mng/buf_ostream.h>

STXXL_BEGIN_NAMESPACE

//! \addtogroup stlalgo
//! \{

/*!
 * External equivalent of std::for_each, see \ref design_algo_foreach.
 *
 * stxxl::for_each applies the function object \c functor to each element in
 * the range [first, last); \c functor's return value, if any, is
 * ignored. Applications are performed in forward order, i.e. from first to
 * last. stxxl::for_each returns the function object after it has been applied
 * to each element.  To overlap I/O and computation \c nbuffers used (a value
 * at least \a D is recommended). The size of the buffers is derived from the
 * container that is pointed by the iterators.
 *
 * \remark The implementation exploits STXXL buffered streams (computation and I/O overlapped).
 *
 * \param begin object of model of \c ext_random_access_iterator concept
 * \param end object of model of \c ext_random_access_iterator concept
 * \param functor function object of model of \c std::UnaryFunction concept
 * \param nbuffers number of buffers (blocks) for internal use (should be at least 2*D )
 * \return function object \c functor after it has been applied to the each element of the given range
 *
 * \warning nested stxxl::for_each are not supported
 */
template <typename ExtIterator, typename UnaryFunction>
UnaryFunction for_each(ExtIterator begin, ExtIterator end,
                       UnaryFunction functor, int_type nbuffers = 0)
{
    if (begin == end)
        return functor;

    typedef typename ExtIterator::value_type value_type;

    typedef buf_istream<
            typename ExtIterator::block_type,
            typename ExtIterator::bids_container_iterator
            > buf_istream_type;

    begin.flush();     // flush container

    if (nbuffers == 0)
        nbuffers = 2 * config::get_instance()->disks_number();

    // create prefetching stream,
    buf_istream_type in(begin.bid(), end.bid() + ((end.block_offset()) ? 1 : 0), nbuffers);

    ExtIterator cur = begin - begin.block_offset();

    // leave part of the block before begin untouched (e.g. copy)
    for ( ; cur != begin; ++cur)
    {
        value_type tmp;
        in >> tmp;
    }

    // apply functor to the range [begin,end)
    for ( ; cur != end; ++cur)
    {
        value_type tmp;
        in >> tmp;
        functor(tmp);
    }

    // leave part of the block after end untouched
    if (end.block_offset())
    {
        ExtIterator last_block_end = end - end.block_offset() + ExtIterator::block_type::size;
        for ( ; cur != last_block_end; ++cur)
        {
            value_type tmp;
            in >> tmp;
        }
    }

    return functor;
}

/*!
 * External equivalent of std::for_each (mutating), see \ref design_algo_foreachm
 *
 * stxxl::for_each_m applies the function object \c functor to each element in
 * the range [first, last); \c functor's return value, if any, is
 * ignored. Applications are performed in forward order, i.e. from first to
 * last. stxxl::for_each_m returns the function object after it has been
 * applied to each element. To overlap I/O and computation \c nbuffers are used
 * (a value at least \a 2D is recommended). The size of the buffers is derived
 * from the container that is pointed by the iterators.
 *
 * \remark The implementation exploits STXXL buffered streams (computation and
 * I/O overlapped)
 *
 * \param begin object of model of \c ext_random_access_iterator concept
 * \param end object of model of \c ext_random_access_iterator concept
 * \param functor object of model of \c std::UnaryFunction concept
 * \param nbuffers number of buffers (blocks) for internal use (should be at least 2*D )
 * \return function object \c functor after it has been applied to the each element of the given range
 *
 * \warning nested stxxl::for_each_m are not supported
 */
template <typename ExtIterator, typename UnaryFunction>
UnaryFunction for_each_m(ExtIterator begin, ExtIterator end,
                         UnaryFunction functor, int_type nbuffers = 0)
{
    if (begin == end)
        return functor;

    typedef typename ExtIterator::value_type value_type;

    typedef buf_istream<
            typename ExtIterator::block_type,
            typename ExtIterator::bids_container_iterator
            > buf_istream_type;

    typedef buf_ostream<
            typename ExtIterator::block_type,
            typename ExtIterator::bids_container_iterator
            > buf_ostream_type;

    begin.flush();     // flush container

    if (nbuffers == 0)
        nbuffers = 2 * config::get_instance()->disks_number();

    // create prefetching stream,
    buf_istream_type in(begin.bid(), end.bid() + ((end.block_offset()) ? 1 : 0), nbuffers / 2);
    // create buffered write stream for blocks
    buf_ostream_type out(begin.bid(), nbuffers / 2);
    // REMARK: these two streams do I/O while
    //         functor is being computed (overlapping for free)

    ExtIterator cur = begin - begin.block_offset();

    // leave part of the block before begin untouched (e.g. copy)
    for ( ; cur != begin; ++cur)
    {
        value_type tmp;
        in >> tmp;
        out << tmp;
    }

    // apply functor to the range [begin,end)
    for ( ; cur != end; ++cur)
    {
        value_type tmp;
        in >> tmp;
        functor(tmp);
        out << tmp;
    }

    // leave part of the block after end untouched
    if (end.block_offset())
    {
        ExtIterator _last_block_end = end - end.block_offset() + ExtIterator::block_type::size;
        for ( ; cur != _last_block_end; ++cur)
        {
            value_type tmp;
            in >> tmp;
            out << tmp;
        }
    }

    return functor;
}

/*!
 * External equivalent of std::generate, see \ref design_algo_generate.
 *
 * Generate assigns the result of invoking \c generator, a function object that
 * takes no arguments, to each element in the range [first, last). To overlap
 * I/O and computation \c nbuffers are used (a value at least \a D is
 * recommended). The size of the buffers is derived from the container that is
 * pointed by the iterators.
 *
 * \remark The implementation exploits STXXL buffered streams (computation and
 * I/O overlapped).
 *
 * \param begin object of model of \c ext_random_access_iterator concept
 * \param end object of model of \c ext_random_access_iterator concept
 * \param generator function object of model of \c std::generator concept
 * \param nbuffers number of buffers (blocks) for internal use (should be at least 2*D, or zero for automaticl 2*D)
 */
template <typename ExtIterator, typename Generator>
void generate(ExtIterator begin, ExtIterator end,
              Generator generator, int_type nbuffers = 0)
{
    typedef typename ExtIterator::block_type block_type;

    typedef buf_ostream<
            block_type, typename ExtIterator::bids_container_iterator
            > buf_ostream_type;

    while (begin.block_offset())    //  go to the beginning of the block
    //  of the external vector
    {
        if (begin == end)
            return;

        *begin = generator();
        ++begin;
    }

    begin.flush();     // flush container

    if (nbuffers == 0)
        nbuffers = 2 * config::get_instance()->disks_number();

    // create buffered write stream for blocks
    buf_ostream_type outstream(begin.bid(), nbuffers);

    assert(begin.block_offset() == 0);

    // delay calling block_externally_updated() until the block is
    // completely filled (and written out) in outstream
    typename ExtIterator::const_iterator prev_block = begin;

    while (end != begin)
    {
        if (begin.block_offset() == 0) {
            if (prev_block != begin) {
                prev_block.block_externally_updated();
                prev_block = begin;
            }
        }

        *outstream = generator();
        ++begin;
        ++outstream;
    }

    typename ExtIterator::const_iterator out = begin;

    while (out.block_offset())    // filling the rest of the block
    {
        *outstream = *out;
        ++out;
        ++outstream;
    }

    if (prev_block != out)
        prev_block.block_externally_updated();

    begin.flush();
}

/*!
 * External equivalent of std::find, see \ref design_algo_find.
 *
 * Returns the first iterator \a i in the range [first, last) such that <tt>*i
 * == value</tt>. Returns last if no such iterator exists.  To overlap I/O and
 * computation \c nbuffers are used (a value at least \a D is recommended). The
 * size of the buffers is derived from the container that is pointed by the
 * iterators.
 *
 * \remark The implementation exploits STXXL buffered streams (computation and
 * I/O overlapped).
 *
 * \param begin object of model of \c ext_random_access_iterator concept
 * \param end object of model of \c ext_random_access_iterator concept
 * \param value value that is equality comparable to the ExtIterator's value type
 * \param nbuffers number of buffers (blocks) for internal use (should be at least 2*D)
 * \return first iterator \c i in the range [begin,end) such that *( \c i ) == \c value, if no
 *         such exists then \c end
 */
template <typename ExtIterator, typename EqualityComparable>
ExtIterator find(ExtIterator begin, ExtIterator end,
                 const EqualityComparable& value, int_type nbuffers = 0)
{
    if (begin == end)
        return end;

    typedef buf_istream<
            typename ExtIterator::block_type,
            typename ExtIterator::bids_container_iterator
            > buf_istream_type;

    begin.flush();     // flush container

    if (nbuffers == 0)
        nbuffers = 2 * config::get_instance()->disks_number();

    // create prefetching stream,
    buf_istream_type in(begin.bid(), end.bid() + ((end.block_offset()) ? 1 : 0), nbuffers);

    ExtIterator cur = begin - begin.block_offset();

    // skip part of the block before begin untouched
    for ( ; cur != begin; ++cur)
        ++in;

    // search in the the range [begin,end)
    for ( ; cur != end; ++cur)
    {
        typename ExtIterator::value_type tmp;
        in >> tmp;
        if (tmp == value)
            return cur;
    }

    return cur;
}

//! \}

STXXL_END_NAMESPACE

#endif // !STXXL_ALGO_SCAN_HEADER
// vim: et:ts=4:sw=4
