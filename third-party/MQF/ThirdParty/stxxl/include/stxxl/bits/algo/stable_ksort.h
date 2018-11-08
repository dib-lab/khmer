/***************************************************************************
 *  include/stxxl/bits/algo/stable_ksort.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2003 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2008 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_ALGO_STABLE_KSORT_HEADER
#define STXXL_ALGO_STABLE_KSORT_HEADER

// it is a first try: distribution sort without sampling
// I rework the stable_ksort when I would have a time

#include <stxxl/bits/mng/block_manager.h>
#include <stxxl/bits/mng/buf_istream.h>
#include <stxxl/bits/mng/buf_ostream.h>
#include <stxxl/bits/common/simple_vector.h>
#include <stxxl/bits/algo/intksort.h>
#include <stxxl/bits/algo/sort_base.h>
#include <stxxl/bits/common/utils.h>

#ifndef STXXL_VERBOSE_STABLE_KSORT
#define STXXL_VERBOSE_STABLE_KSORT STXXL_VERBOSE1
#endif

STXXL_BEGIN_NAMESPACE

//! \addtogroup stlalgo
//! \{

/*! \internal
 */
namespace stable_ksort_local {

template <class Type, class TypeKey>
void classify_block(Type* begin, Type* end, TypeKey*& out,
                    int_type* bucket, typename Type::key_type offset,
                    unsigned shift)
{
    for (Type* p = begin; p < end; p++, out++)      // count & create references
    {
        out->ptr = p;
        typename Type::key_type key = p->key();
        int_type ibucket = (int_type)((key - offset) >> shift);
        out->key = key;
        bucket[ibucket]++;
    }
}

template <typename Type>
struct type_key
{
    typedef typename Type::key_type key_type;
    key_type key;
    Type* ptr;

    type_key() { }
    type_key(key_type k, Type* p) : key(k), ptr(p)
    { }
};

template <typename Type>
bool operator < (const type_key<Type>& a, const type_key<Type>& b)
{
    return a.key < b.key;
}

template <typename Type>
bool operator > (const type_key<Type>& a, const type_key<Type>& b)
{
    return a.key > b.key;
}

template <typename BIDType, typename AllocStrategy>
class bid_sequence
{
public:
    typedef BIDType bid_type;
    typedef bid_type& reference;
    typedef AllocStrategy alloc_strategy;
    typedef typename simple_vector<bid_type>::size_type size_type;
    typedef typename simple_vector<bid_type>::iterator iterator;

protected:
    simple_vector<bid_type>* bids;
    alloc_strategy alloc_strategy_;

public:
    bid_sequence() : bids(NULL) { }
    bid_sequence(size_type size_)
    {
        bids = new simple_vector<bid_type>(size_);
        block_manager* mng = block_manager::get_instance();
        mng->new_blocks(alloc_strategy_, bids->begin(), bids->end());
    }
    void init(size_type size_)
    {
        bids = new simple_vector<bid_type>(size_);
        block_manager* mng = block_manager::get_instance();
        mng->new_blocks(alloc_strategy_, bids->begin(), bids->end());
    }
    reference operator [] (size_type i)
    {
        size_type size_ = size();                 // cache size in a register
        if (i < size_)
            return *(bids->begin() + i);

        block_manager* mng = block_manager::get_instance();
        simple_vector<bid_type>* larger_bids = new simple_vector<bid_type>((i + 1) * 2);
        std::copy(bids->begin(), bids->end(), larger_bids->begin());
        mng->new_blocks(alloc_strategy_, larger_bids->begin() + size_, larger_bids->end());
        delete bids;
        bids = larger_bids;
        return *(larger_bids->begin() + i);
    }
    size_type size() { return bids->size(); }
    iterator begin() { return bids->begin(); }
    ~bid_sequence()
    {
        block_manager::get_instance()->delete_blocks(bids->begin(), bids->end());
        delete bids;
    }
};

template <typename ExtIterator>
void distribute(
    bid_sequence<typename ExtIterator::vector_type::block_type::bid_type,
                 typename ExtIterator::vector_type::alloc_strategy_type>* bucket_bids,
    int64* bucket_sizes,
    const int_type nbuckets,
    const int_type lognbuckets,
    ExtIterator first,
    ExtIterator last,
    const int_type nread_buffers,
    const int_type nwrite_buffers)
{
    typedef typename ExtIterator::vector_type::value_type value_type;
    typedef typename value_type::key_type key_type;
    typedef typename ExtIterator::block_type block_type;
    typedef typename ExtIterator::bids_container_iterator bids_container_iterator;

    typedef buf_istream<block_type, bids_container_iterator> buf_istream_type;

    int_type i = 0;

    buf_istream_type in(first.bid(), last.bid() + ((first.block_offset()) ? 1 : 0),
                        nread_buffers);

    buffered_writer<block_type> out(
        nbuckets + nwrite_buffers,
        nwrite_buffers);

    unsigned_type* bucket_block_offsets = new unsigned_type[nbuckets];
    unsigned_type* bucket_iblock = new unsigned_type[nbuckets];
    block_type** bucket_blocks = new block_type*[nbuckets];

    std::fill(bucket_sizes, bucket_sizes + nbuckets, 0);
    std::fill(bucket_iblock, bucket_iblock + nbuckets, 0);
    std::fill(bucket_block_offsets, bucket_block_offsets + nbuckets, 0);

    for (i = 0; i < nbuckets; i++)
        bucket_blocks[i] = out.get_free_block();

    ExtIterator cur = first - first.block_offset();

    // skip part of the block before first untouched
    for ( ; cur != first; cur++)
        ++in;

    const int_type shift = sizeof(key_type) * 8 - lognbuckets;
    // search in the the range [_begin,_end)
    STXXL_VERBOSE_STABLE_KSORT("Shift by: " << shift << " bits, lognbuckets: " << lognbuckets);
    for ( ; cur != last; cur++)
    {
        key_type cur_key = in.current().key();
        int_type ibucket = (int_type)(cur_key >> shift);

        int_type block_offset = bucket_block_offsets[ibucket];
        in >> (bucket_blocks[ibucket]->elem[block_offset++]);
        if (block_offset == block_type::size)
        {
            block_offset = 0;
            int_type iblock = bucket_iblock[ibucket]++;
            bucket_blocks[ibucket] = out.write(bucket_blocks[ibucket], bucket_bids[ibucket][iblock]);
        }
        bucket_block_offsets[ibucket] = block_offset;
    }
    for (i = 0; i < nbuckets; i++)
    {
        if (bucket_block_offsets[i])
        {
            out.write(bucket_blocks[i], bucket_bids[i][bucket_iblock[i]]);
        }
        bucket_sizes[i] = int64(block_type::size) * bucket_iblock[i] +
                          bucket_block_offsets[i];
        STXXL_VERBOSE_STABLE_KSORT("Bucket " << i << " has size " << bucket_sizes[i] <<
                                   ", estimated size: " << ((last - first) / int64(nbuckets)));
    }

    delete[] bucket_blocks;
    delete[] bucket_block_offsets;
    delete[] bucket_iblock;
}

} // namespace stable_ksort_local

//! Sort records with integer keys
//! \param first object of model of \c ext_random_access_iterator concept
//! \param last object of model of \c ext_random_access_iterator concept
//! \param M amount of memory for internal use (in bytes)
//! \remark Elements must provide a method key() which returns the integer key.
//! \remark Not yet fully implemented, it assumes that the keys are uniformly
//! distributed between [0,std::numeric_limits<key_type>::max().
template <typename ExtIterator>
void stable_ksort(ExtIterator first, ExtIterator last, unsigned_type M)
{
    STXXL_MSG("Warning: stable_ksort is not yet fully implemented, it assumes that the keys are uniformly distributed between [0,std::numeric_limits<key_type>::max()]");
    typedef typename ExtIterator::vector_type::value_type value_type;
    typedef typename value_type::key_type key_type;
    typedef typename ExtIterator::block_type block_type;
    typedef typename ExtIterator::bids_container_iterator bids_container_iterator;
    typedef typename block_type::bid_type bid_type;
    typedef typename ExtIterator::vector_type::alloc_strategy_type alloc_strategy;
    typedef stable_ksort_local::bid_sequence<bid_type, alloc_strategy> bucket_bids_type;
    typedef stable_ksort_local::type_key<value_type> type_key_;

    first.flush();     // flush container

    double begin = timestamp();

    unsigned_type i = 0;
    config* cfg = config::get_instance();
    const unsigned_type m = M / block_type::raw_size;
    assert(2 * block_type::raw_size <= M);
    const unsigned_type write_buffers_multiple = 2;
    const unsigned_type read_buffers_multiple = 2;
    const unsigned_type ndisks = cfg->disks_number();
    const unsigned_type min_num_read_write_buffers = (write_buffers_multiple + read_buffers_multiple) * ndisks;
    const unsigned_type nmaxbuckets = m - min_num_read_write_buffers;
    const unsigned int lognbuckets = ilog2_floor(nmaxbuckets);
    const unsigned_type nbuckets = unsigned_type(1) << lognbuckets;
    const unsigned_type est_bucket_size = (unsigned_type)div_ceil((last - first) / nbuckets, block_type::size);      //in blocks

    if (m < min_num_read_write_buffers + 2 || nbuckets < 2) {
        STXXL_ERRMSG("stxxl::stable_ksort: Not enough memory. Blocks available: " << m <<
                     ", required for r/w buffers: " << min_num_read_write_buffers <<
                     ", required for buckets: 2, nbuckets: " << nbuckets);
        throw bad_parameter("stxxl::stable_ksort(): INSUFFICIENT MEMORY provided, please increase parameter 'M'");
    }
    STXXL_VERBOSE_STABLE_KSORT("Elements to sort: " << (last - first));
    STXXL_VERBOSE_STABLE_KSORT("Number of buckets has to be reduced from " << nmaxbuckets << " to " << nbuckets);
    const unsigned_type nread_buffers = (m - nbuckets) * read_buffers_multiple / (read_buffers_multiple + write_buffers_multiple);
    const unsigned_type nwrite_buffers = (m - nbuckets) * write_buffers_multiple / (read_buffers_multiple + write_buffers_multiple);

    STXXL_VERBOSE_STABLE_KSORT("Read buffers in distribution phase: " << nread_buffers);
    STXXL_VERBOSE_STABLE_KSORT("Write buffers in distribution phase: " << nwrite_buffers);

    bucket_bids_type* bucket_bids = new bucket_bids_type[nbuckets];
    for (i = 0; i < nbuckets; ++i)
        bucket_bids[i].init(est_bucket_size);

    int64* bucket_sizes = new int64[nbuckets];

    disk_queues::get_instance()->set_priority_op(request_queue::WRITE);

    stable_ksort_local::distribute(
        bucket_bids,
        bucket_sizes,
        nbuckets,
        lognbuckets,
        first,
        last,
        nread_buffers,
        nwrite_buffers);

    double dist_end = timestamp(), end;
    double io_wait_after_d = stats::get_instance()->get_io_wait_time();

    {
        // sort buckets
        unsigned_type write_buffers_multiple_bs = 2;
        unsigned_type max_bucket_size_bl = (m - write_buffers_multiple_bs * ndisks) / 2; // in number of blocks
        int64 max_bucket_size_rec = int64(max_bucket_size_bl) * block_type::size;        // in number of records
        int64 max_bucket_size_act = 0;                                                   // actual max bucket size
        // establish output stream

        for (i = 0; i < nbuckets; i++)
        {
            max_bucket_size_act = STXXL_MAX(bucket_sizes[i], max_bucket_size_act);
            if (bucket_sizes[i] > max_bucket_size_rec)
            {
                STXXL_ERRMSG("Bucket " << i << " is too large: " << bucket_sizes[i] <<
                             " records, maximum: " << max_bucket_size_rec);
                STXXL_ERRMSG("Recursion on buckets is not yet implemented, aborting.");
                abort();
            }
        }
        // here we can increase write_buffers_multiple_b knowing max(bucket_sizes[i])
        // ... and decrease max_bucket_size_bl
        const int_type max_bucket_size_act_bl = (int_type)div_ceil(max_bucket_size_act, block_type::size);
        STXXL_VERBOSE_STABLE_KSORT("Reducing required number of required blocks per bucket from " <<
                                   max_bucket_size_bl << " to " << max_bucket_size_act_bl);
        max_bucket_size_rec = max_bucket_size_act;
        max_bucket_size_bl = max_bucket_size_act_bl;
        const unsigned_type nwrite_buffers_bs = m - 2 * max_bucket_size_bl;
        STXXL_VERBOSE_STABLE_KSORT("Write buffers in bucket sorting phase: " << nwrite_buffers_bs);

        typedef buf_ostream<block_type, bids_container_iterator> buf_ostream_type;
        buf_ostream_type out(first.bid(), nwrite_buffers_bs);

        disk_queues::get_instance()->set_priority_op(request_queue::READ);

        if (first.block_offset())
        {
            // has to skip part of the first block
            block_type* block = new block_type;
            request_ptr req;
            req = block->read(*first.bid());
            req->wait();

            for (i = 0; i < first.block_offset(); i++)
            {
                out << block->elem[i];
            }
            delete block;
        }
        block_type* blocks1 = new block_type[max_bucket_size_bl];
        block_type* blocks2 = new block_type[max_bucket_size_bl];
        request_ptr* reqs1 = new request_ptr[max_bucket_size_bl];
        request_ptr* reqs2 = new request_ptr[max_bucket_size_bl];
        type_key_* refs1 = new type_key_[(size_t)max_bucket_size_rec];
        type_key_* refs2 = new type_key_[(size_t)max_bucket_size_rec];

        // submit reading first 2 buckets (Peter's scheme)
        unsigned_type nbucket_blocks = (unsigned_type)div_ceil(bucket_sizes[0], block_type::size);
        for (i = 0; i < nbucket_blocks; i++)
            reqs1[i] = blocks1[i].read(bucket_bids[0][i]);

        nbucket_blocks = (unsigned_type)div_ceil(bucket_sizes[1], block_type::size);
        for (i = 0; i < nbucket_blocks; i++)
            reqs2[i] = blocks2[i].read(bucket_bids[1][i]);

        key_type offset = 0;
        const unsigned log_k1 = STXXL_MAX<unsigned>(ilog2_ceil(max_bucket_size_rec * sizeof(type_key_) / STXXL_L2_SIZE), 1);
        unsigned_type k1 = unsigned_type(1) << log_k1;
        int_type* bucket1 = new int_type[k1];

        const unsigned int shift = (unsigned int)(sizeof(key_type) * 8 - lognbuckets);
        const unsigned int shift1 = shift - log_k1;

        STXXL_VERBOSE_STABLE_KSORT("Classifying " << nbuckets << " buckets, max size:" << max_bucket_size_rec <<
                                   " block size:" << block_type::size << " log_k1:" << log_k1);

        for (unsigned_type k = 0; k < nbuckets; k++)
        {
            nbucket_blocks = (unsigned_type)div_ceil(bucket_sizes[k], block_type::size);
            const unsigned log_k1_k = STXXL_MAX<unsigned>(ilog2_ceil(bucket_sizes[k] * sizeof(type_key_) / STXXL_L2_SIZE), 1);
            assert(log_k1_k <= log_k1);
            k1 = (unsigned_type)(1) << log_k1_k;
            std::fill(bucket1, bucket1 + k1, 0);

            STXXL_VERBOSE_STABLE_KSORT("Classifying bucket " << k << " size:" << bucket_sizes[k] <<
                                       " blocks:" << nbucket_blocks << " log_k1:" << log_k1_k);
            // classify first nbucket_blocks-1 blocks, they are full
            type_key_* ref_ptr = refs1;
            key_type offset1 = offset + (key_type(1) << key_type(shift)) * key_type(k);
            for (i = 0; i < nbucket_blocks - 1; i++)
            {
                reqs1[i]->wait();
                stable_ksort_local::classify_block(blocks1[i].begin(), blocks1[i].end(),
                                                   ref_ptr, bucket1, offset1, shift1 /*,k1*/);
            }
            // last block might be non-full
            const unsigned_type last_block_size =
                (unsigned_type)(bucket_sizes[k] - (nbucket_blocks - 1) * block_type::size);
            reqs1[i]->wait();

            //STXXL_MSG("block_type::size: "<<block_type::size<<" last_block_size:"<<last_block_size);

            classify_block(blocks1[i].begin(), blocks1[i].begin() + last_block_size,
                           ref_ptr, bucket1, offset1, shift1);

            exclusive_prefix_sum(bucket1, k1);
            classify(refs1, refs1 + bucket_sizes[k], refs2, bucket1, offset1, shift1);

            type_key_* c = refs2;
            type_key_* d = refs1;
            for (i = 0; i < k1; i++)
            {
                type_key_* cEnd = refs2 + bucket1[i];
                type_key_* dEnd = refs1 + bucket1[i];

                const unsigned log_k2 = ilog2_floor(bucket1[i]) - 1;        // adaptive bucket size
                const unsigned_type k2 = unsigned_type(1) << log_k2;
                int_type* bucket2 = new int_type[k2];
                const unsigned shift2 = shift1 - log_k2;

                // STXXL_MSG("Sorting bucket "<<k<<":"<<i);
                l1sort(c, cEnd, d, bucket2, k2,
                       offset1 + (key_type(1) << key_type(shift1)) * key_type(i),
                       shift2);

                // write out all
                for (type_key_* p = d; p < dEnd; p++)
                    out << (*(p->ptr));

                delete[] bucket2;
                c = cEnd;
                d = dEnd;
            }
            // submit next read
            const unsigned_type bucket2submit = k + 2;
            if (bucket2submit < nbuckets)
            {
                nbucket_blocks = (unsigned_type)div_ceil(bucket_sizes[bucket2submit], block_type::size);
                for (i = 0; i < nbucket_blocks; i++)
                    reqs1[i] = blocks1[i].read(bucket_bids[bucket2submit][i]);
            }

            std::swap(blocks1, blocks2);
            std::swap(reqs1, reqs2);
        }

        delete[] bucket1;
        delete[] refs1;
        delete[] refs2;
        delete[] blocks1;
        delete[] blocks2;
        delete[] reqs1;
        delete[] reqs2;
        delete[] bucket_bids;
        delete[] bucket_sizes;

        if (last.block_offset())
        {
            // has to skip part of the first block
            block_type* block = new block_type;
            request_ptr req = block->read(*last.bid());
            req->wait();

            for (i = last.block_offset(); i < block_type::size; i++)
            {
                out << block->elem[i];
            }
            delete block;
        }

        end = timestamp();
    }

    STXXL_VERBOSE("Elapsed time        : " << end - begin << " s. Distribution time: " <<
                  dist_end - begin << " s");
    STXXL_VERBOSE("Time in I/O wait(ds): " << io_wait_after_d << " s");
    STXXL_VERBOSE(*stats::get_instance());
}

//! \}

STXXL_END_NAMESPACE

#endif // !STXXL_ALGO_STABLE_KSORT_HEADER
