// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_POOL_SORTER_H
#define SEQAN_HEADER_POOL_SORTER_H

namespace SEQAN_NAMESPACE_MAIN
{

    template < typename TValue, typename Compare >
    struct MergeStreamComparer : public std::binary_function < PageBucket<TValue>,
                                                               PageBucket<TValue>,
                                                               bool>
    {
        Compare C;
        MergeStreamComparer(Compare &tmpC): C(tmpC) { }
        inline bool operator()(const PageBucket<TValue> &a,
                               const PageBucket<TValue> &b) const {
            return C(static_cast<const TValue&>(*(a.cur)),
                     static_cast<const TValue&>(*(b.cur))) > 0;
        }
    };

    template < typename TCompare >
    struct AdaptorCompare2Less :
        public std::binary_function <
            typename TCompare::first_argument_type,
            typename TCompare::second_argument_type,
            bool >
    {
        TCompare const & C;
        AdaptorCompare2Less(TCompare const & tmpC): C(tmpC) { }
        inline bool operator() (
            typename TCompare::first_argument_type const &a,
            typename TCompare::second_argument_type const &b) const
        {
            return C(a, b) < 0;
        }
    };

/*!
 * @class SorterConfigSize
 * @extends SorterSpec
 * @headerfile <seqan/pipe.h>
 *
 * @brief Configuration of Sorter.
 *
 * @signature template <typename TCompare, typename TSize[, typename TFile>
 *            struct SorterConfigSize;
 *
 * @tparam TCompare The compare function (see STL's <tt>binary_function</tt>).
 * @tparam TSize    The Sorter's size type.
 * @tparam TFile    The underlying File type.  <tt>File&lt;&gt;</tt>
 *
 * The requirement on TCompare are as follows: let <tt>comp</tt> be an object of type <tt>TCompare</tt>.  <tt>comp(a,
 * b)</tt> should return a value less, equal, or greater than 0 if <tt>a &lt; b</tt>, <tt>a == b</tt>, or <tt>a &gt;
 * b</tt>.
 *
 * @see SorterConfig
 */

    template < typename TCompare,
               typename TSize,
               typename TFile = File<> >
    struct SorterConfigSize {
        typedef TCompare    Compare;
        typedef TSize        SizeType;
        typedef TFile        File;
    };

/*!
 * @class SorterConfig
 * @extends SorterSpec
 * @headerfile <seqan/pipe.h>
 * @brief Configuration of Sorter.
 *
 * @signature template <typename TCompare[, typename TFile]>
 *            struct SorterConfig;
 *
 * @tparam TCompare The compare function (see STL's <tt>binary_function</tt>).
 * @tparam TFile The underlying File type, defaults to <tt>File&lt;&gt;</tt>.
 *
 * The requirement on TCompare are as follows: let <tt>comp</tt> be an object of type <tt>TCompare</tt>.  <tt>comp(a,
 * b)</tt> should return a value less, equal, or greater than 0 if <tt>a &lt; b</tt>, <tt>a == b</tt>, or <tt>a &gt;
 * b</tt>.
 *
 * Using this configuration spec., the Sorter's size type is <tt>Size&lt;TFile&gt;::Type</tt>.  To use a custom size
 * type @link SorterConfigSize @endlink should be used.
 *
 * @see SorterConfigSize
 */

    template < typename TCompare, typename TFile = File<> >
    struct SorterConfig {
        typedef TCompare                    Compare;
        typedef typename Size<TFile>::Type    SizeType;
        typedef TFile                        File;
    };

/*!
 * @class SorterSpec
 * @extends Pool
 * @headerfile <seqan/pipe.h>
 * @brief Sorts all elements using a custom compare function.
 *
 * @signature template <typename TValue, typename TConfig>
 *            struct Pool<TValue, SorterSpec<TConfig > >;
 *
 * @tparam TConfig Configuration Spec.  Defines compare function, size type, and file type.
 * @tparam TValue  The value type, that is the type of the stream elements.
 *
 * The Pool's input/output type is <tt>TValue</tt> and the size type is determined by the <tt>TConfig</tt>.
 *
 * If there exists an explicit function mapping input elements to their destined positions in the output
 * stream, @link MapperSpec @endlink should be preferred.
 *
 * @section Example
 *
 * @include demos/pipe/pool_sorter.cpp
 *
 * @code{.console}
 * Sorted elements:
 * -1
 * 3
 * 10
 * @endcode
 */

    template < typename TConfig >
    struct SorterSpec {
        typedef TConfig Config;
    };

    template < typename TValue,
               typename TConfig >
    struct HandlerArgs< Pool< TValue, SorterSpec<TConfig> > >
    {
        typedef typename TConfig::Compare Type;
    };


    //////////////////////////////////////////////////////////////////////////////
    // cache bucket based synchronous multiway merge
    struct ReadSorterSpec_;
    typedef Tag<ReadSorterSpec_> ReadSorterSpec;

    template <typename TValue, typename TPoolSpec>
    struct Handler<Pool<TValue, TPoolSpec>, ReadSorterSpec>
    {
        typedef Pool<TValue, TPoolSpec>                 TPool;
        typedef typename TPool::TBuffer                    TBuffer;
        typedef typename TPoolSpec::Config::Compare     TCompare;
        typedef PageBucketExtended<TValue>              TPageBucket;

        typedef MergeStreamComparer<TValue, TCompare>   TStreamComparer;
/*        typedef std::priority_queue <
            TPageBucket,
            std::vector<TPageBucket>,
            MergeStreamComparer<Type, Compare> >        PQueue;
*/
        typedef PriorityType<TPageBucket, TStreamComparer> TPrioQueue;

        TPool       &pool;
        TBuffer     bucketBuffer;
        TPrioQueue  pqueue;

        Handler(TPool &_pool):
            pool(_pool),
            pqueue(TStreamComparer(_pool.handlerArgs)) { }

        ~Handler() {
            cancel();
        }

        struct insertBucket : public std::unary_function<TPageBucket,void>
        {
            Handler &me;
            insertBucket(Handler &_me): me(_me) {}

            inline void operator() (TPageBucket &pb) const
            {
                pb.pageNo = length(me.pqueue);
                readBucket(pb, pb.pageNo, me.pool.pageSize, me.pool.dataSize(pb.pageNo), me.pool.file);
                push(me.pqueue, pb);
            }
        };

        bool begin()
        {
            // 1. initially fill priority queue
//            pqueue.reserve(pool.pages);
            equiDistantDistribution(
                bucketBuffer, pool.bucketBufferSize, *this,
                pool._size, pool.pageSize,
                insertBucket(*this));
            return true;
        }

        inline TValue const & front() const
        {
            return *(top(pqueue).cur);
        }

        inline void pop(TValue &Ref_)
        {
            TPageBucket &pb = top(pqueue);
            SEQAN_ASSERT_LEQ(pb.cur, pb.end);

            Ref_ = *pb.cur;
            if (++pb.cur == pb.end)
            {
                // bucket is empty, we have to fetch the next bucket
                if (!readBucket(pb, pb.pageNo, pool.pageSize, pool.dataSize(pb.pageNo), pool.file)) {
                    seqan::pop(pqueue);
                    return;
                }
            }
            adjustTop(pqueue);
        }

        inline void pop()
        {
            TPageBucket &pb = top(pqueue);
            SEQAN_ASSERT_LEQ(pb.cur, pb.end);

            if (++pb.cur == pb.end)
                // bucket is empty, we have to fetch the next bucket
                if (!readBucket(pb, pb.pageNo, pool.pageSize, pool.dataSize(pb.pageNo), pool.file)) {
                    seqan::pop(pqueue);
                    return;
                }
            adjustTop(pqueue);
        }

        inline bool eof() const
        {
            return empty(pqueue);
        }

        inline void end()
        {
            cancel();
        }

        void cancel()
        {
            clear(pqueue);
            freePage(bucketBuffer, *this);
        }

        inline void process() {}

    };


    template <typename TValue, typename TPoolSpec>
    struct BufferHandler< Pool<TValue, TPoolSpec>, ReadSorterSpec >
    {
        typedef Pool<TValue, TPoolSpec>                 TPool;
        typedef typename TPool::Buffer                  TBuffer;
        typedef typename TPool::Spec::Compare           TCompare;

        typedef PageBucketExtended<TValue>              TPageBucket;

        typedef MergeStreamComparer<TValue, TCompare>   StreamComparer;
        typedef std::priority_queue <
            TPageBucket,
            std::vector<TPageBucket>,
            MergeStreamComparer<TValue, TCompare> >     TPrioQueue;

        TPool       &pool;
        unsigned    mergeBufferSize;
        TBuffer     bucketBuffer, mergeBuffer, tmpBuffer;
        TPrioQueue  pqueue;

        BufferHandler(TPool &_pool):
            pool(_pool),
            mergeBufferSize(_pool.pageSize),
            pqueue(StreamComparer(_pool.handlerData.compare)) { }

        BufferHandler(TPool &_pool, unsigned _requestedBufferSize):
            pool(_pool),
            mergeBufferSize(_min(_pool.size(), _requestedBufferSize)),
            pqueue(StreamComparer(_pool.handlerData.compare)) { }

        ~BufferHandler() {
            cancel();
        }

        struct insertBucket : public std::unary_function<TPageBucket, void> {
            BufferHandler &me;
            insertBucket(BufferHandler &_me): me(_me) {}

            inline void operator() (TPageBucket &pb) const
            {
                pb.pageNo = length(pqueue);
                readBucket(pb, pb.pageNo, me.pool.pageSize, me.pool.dataSize(pb.pageNo), pool.file);
                push(pqueue, pb);
            }
        };

        inline TBuffer & first() {
            // 1. initially fill priority queue
//            pqueue.reserve(pool.pages);
            equiDistantDistribution(
                bucketBuffer, pool.bucketBufferSize, *this,
                pool._size, pool.pageSize,
                insertBucket(*this));
            allocPage(mergeBuffer, mergeBufferSize, *this);
            return merge();
        }

        inline TBuffer & next() {
            return merge();
        }

        inline void end() {
            cancel();
        }

        void cancel()
        {
            clear(pqueue);
            freePage(mergeBuffer, *this);
            freePage(bucketBuffer, *this);
        }

        inline void process() {}

    private:

        TBuffer & merge()
        {
            // 2. merge streams into mergeBuffer

            typename TPrioQueue::size_type pqsize = length(pqueue);
            if (!pqsize) {
                resize(tmpBuffer, 0);
                return tmpBuffer;
            }

            if (pqsize == 1)
            {
                // only one stream left => retrieve what's left in stream es
                TPageBucket &pb = pqueue.top();

                // flush stream
                if (pb.cur != pb.end) {
                    tmpBuffer.begin = pb.cur;
                    tmpBuffer.end = pb.end;
                    pb.cur = pb.end;
                    if (pb.pageOfs == pool.dataSize(pb.pageNo))
                        pop(pqueue);
                    return tmpBuffer;
                }

                // read directly from disk
                pb.begin = mergeBuffer.begin;
                pb.end = pb.cur = pb.begin + pageSize(mergeBuffer);

                resize(mergeBuffer, readBucket(pb, pb.pageNo, pool.pageSize, pool.dataSize(pb.pageNo)));
                if (pb.pageOfs == pool.dataSize(pb.pageNo))
                    pop(pqueue);
            }
            else
            {
                for(TValue *cur = mergeBuffer.begin; cur != mergeBuffer.end; ++cur) {
                    TPageBucket &pb = pqueue.top();
                    *cur = *pb.cur;
                    if (++pb.cur == pb.end) {
                        // bucket is empty, we have to fetch the next bucket
                        if (!readBucket(pb, pb.pageNo, pool.pageSize, pool.dataSize(pb.pageNo), pool.file)) {
                            pop(pqueue);
                            // queue contains only one stream
                            // => we return what we have merged
                            if (--pqsize == 1) {
                                resize(cur - mergeBuffer.begin + 1);
                                return mergeBuffer;
                            }
                        }
                    }
                    adjustTop(pqueue);
                }
                resize(mergeBuffer, pageSize(mergeBuffer));
            }

            return mergeBuffer;
        }
    };

    template < typename TValue,
               typename TConfig >
    inline Buffer<TValue, PageFrame<typename TConfig::File, Dynamic> > & processBuffer(
        Buffer<TValue, PageFrame<typename TConfig::File, Dynamic> > &buf,
        BufferHandler< Pool< TValue, SorterSpec<TConfig> >, WriteFileSpec > &me)
    {
        AdaptorCompare2Less<typename TConfig::Compare> cmp(me.pool.handlerArgs);
        std::sort(buf.begin, buf.end, cmp);
        return buf;
    }


    //////////////////////////////////////////////////////////////////////////////
    // character and buffer based handler definitions
    template < typename TValue,
               typename TConfig >
    inline Buffer<TValue> & processBuffer(
        Buffer<TValue> &buf,
        BufferHandler< Pool< TValue, SorterSpec<TConfig> >, MemorySpec > &me)
    {
        AdaptorCompare2Less<typename TConfig::Compare> cmp(me.pool.handlerArgs);
        std::sort(buf.begin, buf.end, cmp);
        return buf;
    }

    template < typename TValue,
               typename TConfig >
    struct BufReadHandler< Pool< TValue, SorterSpec<TConfig> > >
    {
        typedef BufferHandler< Bundle2<
            BufferHandler< Pool< TValue, SorterSpec<TConfig> >, MemorySpec >,
            BufferHandler< Pool< TValue, SorterSpec<TConfig> >, ReadSorterSpec >
        >, MultiplexSpec > Type;
    };

    template < typename TValue,
               typename TConfig >
    struct ReadHandler< Pool< TValue, SorterSpec<TConfig> > >
    {
        typedef Handler< Bundle2<
            Handler< BufferHandler < Pool< TValue, SorterSpec<TConfig> >, MemorySpec >, AdapterSpec >,
            Handler<                 Pool< TValue, SorterSpec<TConfig> >, ReadSorterSpec >
        >, MultiplexSpec > Type;
    };

}

#endif
