/***************************************************************************
 *  include/stxxl/bits/algo/random_shuffle.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2007 Manuel Krings
 *  Copyright (C) 2007 Markus Westphal <mail@markuswestphal.de>
 *  Copyright (C) 2009, 2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_ALGO_RANDOM_SHUFFLE_HEADER
#define STXXL_ALGO_RANDOM_SHUFFLE_HEADER

// TODO: improve main memory consumption in recursion
//        (free stacks buffers)
// TODO: shuffle small input in internal memory

#include <stxxl/bits/stream/stream.h>
#include <stxxl/scan>
#include <stxxl/stack>

STXXL_BEGIN_NAMESPACE

//! \addtogroup stlalgo
//! \{

//! External equivalent of std::random_shuffle
//! \param first begin of the range to shuffle
//! \param last end of the range to shuffle
//! \param rand random number generator object (functor)
//! \param M number of bytes for internal use
//! \param AS parallel disk allocation strategy
//!
//! - BlockSize size of the block to use for external memory data structures
//! - PageSize page size in blocks to use for external memory data structures
template <typename ExtIterator,
          typename RandomNumberGenerator,
          unsigned BlockSize,
          unsigned PageSize,
          typename AllocStrategy>
void random_shuffle(ExtIterator first,
                    ExtIterator last,
                    RandomNumberGenerator& rand,
                    unsigned_type M,
                    AllocStrategy AS = STXXL_DEFAULT_ALLOC_STRATEGY())
{
    STXXL_UNUSED(AS);  // FIXME: Why is this not being used?
    typedef typename ExtIterator::value_type value_type;
    typedef typename STACK_GENERATOR<
            value_type, external, grow_shrink2, PageSize,
            BlockSize, void, 0, AllocStrategy
            >::result stack_type;
    typedef typename stack_type::block_type block_type;

    STXXL_VERBOSE1("random_shuffle: Plain Version");
    STXXL_STATIC_ASSERT(int(BlockSize) < 0 && "This implementation was never tested. Please report to the stxxl developers if you have an ExtIterator that works with this implementation.");

    int64 n = last - first; // the number of input elements

    // make sure we have at least 6 blocks + 1 page
    if (M < 6 * BlockSize + PageSize * BlockSize) {
        STXXL_ERRMSG("random_shuffle: insufficient memory, " << M << " bytes supplied,");
        M = 6 * BlockSize + PageSize * BlockSize;
        STXXL_ERRMSG("random_shuffle: increasing to " << M << " bytes (6 blocks + 1 page)");
    }

    int_type k = M / (3 * BlockSize); // number of buckets

    int64 i, j, size = 0;

    value_type* temp_array;
    typedef typename VECTOR_GENERATOR<
            value_type, PageSize, 4, BlockSize, AllocStrategy
            >::result temp_vector_type;
    temp_vector_type* temp_vector;

    STXXL_VERBOSE1("random_shuffle: " << M / BlockSize - k << " write buffers for " << k << " buckets");
    read_write_pool<block_type> pool(0, M / BlockSize - k);  // no read buffers and M/B-k write buffers

    stack_type** buckets;

    // create and put buckets into container
    buckets = new stack_type*[k];
    for (j = 0; j < k; j++)
        buckets[j] = new stack_type(pool, 0);

    ///// Reading input /////////////////////
    typedef typename stream::streamify_traits<ExtIterator>::stream_type input_stream;
    input_stream in = stream::streamify(first, last);

    // distribute input into random buckets
    int_type random_bucket = 0;
    for (i = 0; i < n; ++i) {
        random_bucket = rand(k);
        buckets[random_bucket]->push(*in); // reading the current input element
        ++in;                              // go to the next input element
    }

    ///// Processing //////////////////////
    // resize buffers
    pool.resize_write(0);
    pool.resize_prefetch(PageSize);

    unsigned_type space_left = M - k * BlockSize -
                               PageSize * BlockSize; // remaining int space
    ExtIterator Writer = first;
    ExtIterator it = first;

    for (i = 0; i < k; i++) {
        STXXL_VERBOSE1("random_shuffle: bucket no " << i << " contains " << buckets[i]->size() << " elements");
    }

    // shuffle each bucket
    for (i = 0; i < k; i++) {
        buckets[i]->set_prefetch_aggr(PageSize);
        size = buckets[i]->size();

        // does the bucket fit into memory?
        if (size * sizeof(value_type) < space_left) {
            STXXL_VERBOSE1("random_shuffle: no recursion");

            // copy bucket into temp. array
            temp_array = new value_type[size];
            for (j = 0; j < size; j++) {
                temp_array[j] = buckets[i]->top();
                buckets[i]->pop();
            }

            // shuffle
            potentially_parallel::
            random_shuffle(temp_array, temp_array + size, rand);

            // write back
            for (j = 0; j < size; j++) {
                *Writer = temp_array[j];
                ++Writer;
            }

            // free memory
            delete[] temp_array;
        }
        else {
            STXXL_VERBOSE1("random_shuffle: recursion");

            // copy bucket into temp. stxxl::vector
            temp_vector = new temp_vector_type(size);

            for (j = 0; j < size; j++) {
                (*temp_vector)[j] = buckets[i]->top();
                buckets[i]->pop();
            }

            pool.resize_prefetch(0);
            space_left += PageSize * BlockSize;
            STXXL_VERBOSE1("random_shuffle: Space left: " << space_left);

            // recursive shuffle
            stxxl::random_shuffle(temp_vector->begin(),
                                  temp_vector->end(), rand, space_left);

            pool.resize_prefetch(PageSize);

            // write back
            for (j = 0; j < size; j++) {
                *Writer = (*temp_vector)[j];
                ++Writer;
            }

            // free memory
            delete temp_vector;
        }

        // free bucket
        delete buckets[i];
        space_left += BlockSize;
    }

    delete[] buckets;
}

//! External equivalent of std::random_shuffle (specialization for stxxl::vector)
//! \param first begin of the range to shuffle
//! \param last end of the range to shuffle
//! \param rand random number generator object (functor)
//! \param M number of bytes for internal use
template <typename Type, typename AllocStrategy, typename SizeType, typename DiffType,
          unsigned BlockSize, typename PageType, unsigned PageSize, typename RandomNumberGenerator>
void random_shuffle(
    stxxl::vector_iterator<Type, AllocStrategy, SizeType, DiffType,
                           BlockSize, PageType, PageSize> first,
    stxxl::vector_iterator<Type, AllocStrategy, SizeType, DiffType,
                           BlockSize, PageType, PageSize> last,
    RandomNumberGenerator& rand,
    unsigned_type M)
{
    typedef stxxl::vector_iterator<Type, AllocStrategy, SizeType, DiffType, BlockSize, PageType, PageSize> ExtIterator;
    typedef typename ExtIterator::value_type value_type;
    typedef typename ExtIterator::bids_container_iterator bids_container_iterator;
    typedef typename stxxl::STACK_GENERATOR<value_type, stxxl::external,
                                            stxxl::grow_shrink2, PageSize, BlockSize>::result stack_type;
    typedef typename stack_type::block_type block_type;

    STXXL_VERBOSE1("random_shuffle: Vector Version");

    // make sure we have at least 6 blocks + 1 page
    if (M < 6 * BlockSize + PageSize * BlockSize) {
        STXXL_ERRMSG("random_shuffle: insufficient memory, " << M << " bytes supplied,");
        M = 6 * BlockSize + PageSize * BlockSize;
        STXXL_ERRMSG("random_shuffle: increasing to " << M << " bytes (6 blocks + 1 page)");
    }

    stxxl::int64 n = last - first;    // the number of input elements
    int_type k = M / (3 * BlockSize); // number of buckets

    stxxl::int64 i, j, size = 0;

    value_type* temp_array;
    typedef typename stxxl::VECTOR_GENERATOR<
            value_type, PageSize, 4, BlockSize, AllocStrategy
            >::result temp_vector_type;
    temp_vector_type* temp_vector;

    // no read buffers and M/B-k write buffers
    stxxl::read_write_pool<block_type> pool(0, M / BlockSize - k);

    stack_type** buckets;

    // create and put buckets into container
    buckets = new stack_type*[k];
    for (j = 0; j < k; j++)
        buckets[j] = new stack_type(pool, 0);

    typedef buf_istream<block_type, bids_container_iterator> buf_istream_type;
    typedef buf_ostream<block_type, bids_container_iterator> buf_ostream_type;

    first.flush();     // flush container

    // create prefetching stream,
    buf_istream_type in(first.bid(), last.bid() + ((last.block_offset()) ? 1 : 0), 2);
    // create buffered write stream for blocks
    buf_ostream_type out(first.bid(), 2);

    ExtIterator _cur = first - first.block_offset();

    // leave part of the block before _begin untouched (e.g. copy)
    for ( ; _cur != first; ++_cur)
    {
        typename ExtIterator::value_type tmp;
        in >> tmp;
        out << tmp;
    }

    ///// Reading input /////////////////////

    // distribute input into random buckets
    int_type random_bucket = 0;
    for (i = 0; i < n; ++i, ++_cur) {
        random_bucket = rand((unsigned)k);
        typename ExtIterator::value_type tmp;
        in >> tmp;
        buckets[random_bucket]->push(tmp); // reading the current input element
    }

    ///// Processing //////////////////////
    // resize buffers
    pool.resize_write(0);
    pool.resize_prefetch(PageSize);

    // remaining int space
    unsigned_type space_left = M - k * BlockSize - PageSize * BlockSize;

    for (i = 0; i < k; i++) {
        STXXL_VERBOSE1("random_shuffle: bucket no " << i << " contains " << buckets[i]->size() << " elements");
    }

    // shuffle each bucket
    for (i = 0; i < k; i++) {
        buckets[i]->set_prefetch_aggr(PageSize);
        size = buckets[i]->size();

        // does the bucket fit into memory?
        if (size * sizeof(value_type) < space_left) {
            STXXL_VERBOSE1("random_shuffle: no recursion");

            // copy bucket into temp. array
            temp_array = new value_type[(size_t)size];
            for (j = 0; j < size; j++) {
                temp_array[j] = buckets[i]->top();
                buckets[i]->pop();
            }

            // shuffle
            potentially_parallel::
            random_shuffle(temp_array, temp_array + size, rand);

            // write back
            for (j = 0; j < size; j++) {
                typename ExtIterator::value_type tmp;
                tmp = temp_array[j];
                out << tmp;
            }

            // free memory
            delete[] temp_array;
        }
        else {
            STXXL_VERBOSE1("random_shuffle: recursion");
            // copy bucket into temp. stxxl::vector
            temp_vector = new temp_vector_type(size);

            for (j = 0; j < size; j++) {
                (*temp_vector)[j] = buckets[i]->top();
                buckets[i]->pop();
            }

            pool.resize_prefetch(0);
            space_left += PageSize * BlockSize;

            STXXL_VERBOSE1("random_shuffle: Space left: " << space_left);

            // recursive shuffle
            stxxl::random_shuffle(temp_vector->begin(),
                                  temp_vector->end(), rand, space_left);

            pool.resize_prefetch(PageSize);

            // write back
            for (j = 0; j < size; j++) {
                typename ExtIterator::value_type tmp;
                tmp = (*temp_vector)[j];
                out << tmp;
            }

            // free memory
            delete temp_vector;
        }

        // free bucket
        delete buckets[i];
        space_left += BlockSize;
    }

    delete[] buckets;

    // leave part of the block after _end untouched
    if (last.block_offset())
    {
        ExtIterator last_block_end = last + (block_type::size - last.block_offset());
        for ( ; _cur != last_block_end; ++_cur)
        {
            typename ExtIterator::value_type tmp;
            in >> tmp;
            out << tmp;
        }
    }
}

//! External equivalent of std::random_shuffle (specialization for stxxl::vector)
//! \param first begin of the range to shuffle
//! \param last end of the range to shuffle
//! \param M number of bytes for internal use
template <typename Type, typename AllocStrategy, typename SizeType, typename DiffType,
          unsigned BlockSize, typename PageType, unsigned PageSize>
inline
void random_shuffle(
    stxxl::vector_iterator<Type, AllocStrategy, SizeType, DiffType,
                           BlockSize, PageType, PageSize> first,
    stxxl::vector_iterator<Type, AllocStrategy, SizeType, DiffType,
                           BlockSize, PageType, PageSize> last,
    unsigned_type M)
{
    stxxl::random_number<> rand;
    stxxl::random_shuffle(first, last, rand, M);
}

//! \}

STXXL_END_NAMESPACE

#endif // !STXXL_ALGO_RANDOM_SHUFFLE_HEADER
// vim: et:ts=4:sw=4
