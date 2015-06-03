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
// Thread-safe queue
// ==========================================================================
// This queue (almost) lock-free. It only uses spin-locks when the queue must
// be resized or for blocking if the queue is full or empty.
//
// If the queue is often empty or full, the suspendable queue should be used.

#ifndef SEQAN_PARALLEL_PARALLEL_QUEUE_H_
#define SEQAN_PARALLEL_PARALLEL_QUEUE_H_

namespace seqan {

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class ConcurrentQueue
// ----------------------------------------------------------------------------
/*!
 * @class ConcurrentQueue Concurrent Queue
 * @headerfile <seqan/parallel.h>
 * @brief Thread-safe queue for multiple producers and multiple consumers.
 *
 * @signature template <typename TValue, typename TSpec>
 *            class ConcurrentQueue;
 *
 * @tparam TValue Element type of the queue.
 * @tparam TSpec  Tag for further specializing the Concurrent Queue. Default is <tt>void</tt>.
 *
 * The Concurrent Queue is a thread-safe FIFO queue that supports multiple producers and multiple consumers (MPMC).
 * Elements are enqueued via @link ConcurrentQueue#appendValue @endlink and dequeued with @link
 * ConcurrentQueue#tryPopFront @endlink or @link ConcurrentQueue#popFront @endlink.
 * Depending on the expansion tag of appendValue it can grow dynamically or have a fixed size.
 *
 * The implementation is lock-free and uses a @Class.AllocString@ as ring buffer.
 *
 * @section Examples
 *
 * Simple example for a single producer single consumer (SPSC) dynamic queue.
 *
 * @include demos/parallel/queue_example.cpp
 *
 * The output is:
 *
 * @include demos/parallel/queue_example.cpp.stdout
 */

template <typename TValue, typename TSpec = void>
class ConcurrentQueue
{
public:
    typedef typename Host<ConcurrentQueue>::Type                TString;
    typedef typename Size<TString>::Type                        TSize;
    typedef typename Atomic<TSize>::Type                        TAtomicSize;
    typedef typename If<typename IsSameType<TSpec, Limit>::Type,
                        Serial,
                        typename DefaultParallelSpec<ConcurrentQueue>::Type
                     >::Type                                    TLockTag;

    TString                 data;
    mutable ReadWriteLock   lock;           char pad1[SEQAN_CACHE_LINE_SIZE - sizeof(ReadWriteLock)];
    size_t                  readerCount;
    size_t                  writerCount;

    TAtomicSize headPos;                    char pad4[SEQAN_CACHE_LINE_SIZE - sizeof(TAtomicSize)];
    TAtomicSize headReadPos;                char pad5[SEQAN_CACHE_LINE_SIZE - sizeof(TAtomicSize)];
    TAtomicSize tailPos;                    char pad6[SEQAN_CACHE_LINE_SIZE - sizeof(TAtomicSize)];
    TAtomicSize tailWritePos;               char pad7[SEQAN_CACHE_LINE_SIZE - sizeof(TAtomicSize)];
    TAtomicSize roundSize;                  char pad8[SEQAN_CACHE_LINE_SIZE - sizeof(TAtomicSize)];
    Atomic<bool>::Type virgin;              char pad9[SEQAN_CACHE_LINE_SIZE - sizeof(Atomic < bool > ::Type)];

    ConcurrentQueue() :
        readerCount(0),
        writerCount(0),
        headPos(0),
        headReadPos(0),
        tailPos(0),
        tailWritePos(0),
        roundSize(0),
        virgin(true)
    {}

    // you can set the initial capacity here
    explicit
    ConcurrentQueue(TSize initCapacity) :
        readerCount(0),
        writerCount(0),
        headPos(0),
        headReadPos(0),
        tailPos(0),
        tailWritePos(0),
        virgin(true)
    {
        reserve(data, initCapacity + 1, Exact());
        roundSize = (TSize)1 << (log2(capacity(data) - 1) + 1);
    }

    explicit
    ConcurrentQueue(TString & data) :
        data(data),
        readerCount(0),
        writerCount(0),
        headPos(0),
        headReadPos(0),
        tailPos(length(data)),
        tailWritePos(length(data)),
        virgin(true)
    {
        roundSize = (TSize)1 << (log2(capacity(data) - 1) + 1);
    }

    ~ConcurrentQueue()
    {
        SEQAN_ASSERT_EQ(tailPos, tailWritePos);
        SEQAN_ASSERT_EQ(headPos, headReadPos);
        SEQAN_ASSERT(empty(lock));
        SEQAN_ASSERT_EQ(writerCount, 0u);

        TSize mask = roundSize - 1;
        headPos &= mask;
        tailPos &= mask;

        // wait for all pending readers to finish
        while (readerCount != 0)
        {}

        typename Iterator<TString, Standard>::Type arrayBegin = begin(data, Standard());

        if (headPos <= tailPos)
        {
            arrayDestruct(arrayBegin + headPos, arrayBegin + tailPos);
        }
        else
        {
            arrayDestruct(arrayBegin, arrayBegin + tailPos);
            arrayDestruct(arrayBegin + headPos, arrayBegin + capacity(data));
        }
        _setLength(data, 0);
    }

private:
    ConcurrentQueue(ConcurrentQueue const &);
    void operator=(ConcurrentQueue const &);
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction DefaultParallelSpec
// ----------------------------------------------------------------------------

template <typename TValue>
struct DefaultParallelSpec<ConcurrentQueue<TValue, Serial> >
{
    typedef Serial Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct Value<ConcurrentQueue<TValue, TSpec> >
{
    typedef TValue Type;
};

template <typename TValue, typename TSpec>
struct Host<ConcurrentQueue<TValue, TSpec> >
{
    typedef String<TValue> Type;
};

template <typename TValue, typename TSpec>
struct Size<ConcurrentQueue<TValue, TSpec> >:
    Size<Host<ConcurrentQueue<TValue, TSpec> > >
{};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function lockReading() / unlockReading()
// ----------------------------------------------------------------------------

/*!
 * @fn ConcurrentQueue#lockReading
 * @brief Register a reader.
 *
 * @signature void lockReading(queue);
 *
 * @param[in] queue The queue to register a reader at.
 *
 * The destructor of the queue will spinlock until all readers are deregistered.
 */

template <typename TValue, typename TSpec>
inline void
lockReading(ConcurrentQueue<TValue, TSpec> & me)
{
    atomicInc(me.readerCount);
}

/*!
 * @fn ConcurrentQueue#unlockReading
 * @brief Deregister a reader.
 *
 * @signature void unlockReading(queue);
 *
 * @param[in] queue The queue to deregister a reader from.
 *
 * The destructor of the queue will spinlock until all readers are deregistered.
 */

template <typename TValue, typename TSpec>
inline void
unlockReading(ConcurrentQueue<TValue, TSpec> & me)
{
    atomicDec(me.readerCount);
}

// ----------------------------------------------------------------------------
// Function lockWriting() / unlockWriting()
// ----------------------------------------------------------------------------

/*!
 * @fn ConcurrentQueue#lockWriting
 * @brief Register a writer.
 *
 * @signature void lockWriting(queue);
 *
 * @param[in] queue The queue to register a writer at.
 *
 */

template <typename TValue, typename TSpec>
inline void
lockWriting(ConcurrentQueue<TValue, TSpec> & me)
{
    atomicInc(me.writerCount);
}

/*!
 * @fn ConcurrentQueue#unlockWriting
 * @brief Deregister a writer.
 *
 * @signature void unlockWriting(queue);
 *
 * @param[in] queue The queue to deregister a writer from.
 *
 */

template <typename TValue, typename TSpec>
inline void
unlockWriting(ConcurrentQueue<TValue, TSpec> & me)
{
    atomicDec(me.writerCount);
}

// ----------------------------------------------------------------------------
// Function _cyclicInc()
// ----------------------------------------------------------------------------

template <typename TValue>
inline TValue
_cyclicInc(TValue value, TValue modulo, TValue roundSize)
{
    // invariants:
    //   - roundSize is a power of 2
    //   - (value % roundSize) is in [0, modulo)
    //
    // return the next greater value that fulfills the invariants
    if ((++value & (roundSize - 1)) >= modulo)
        value += roundSize - modulo;
    return value;
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

/*!
 * @fn ConcurrentQueue#empty
 * @brief Returns whether a queue is empty.
 *
 * @signature bool empty(queue);
 *
 * @param[in] queue The queue to query.
 *
 * @return bool Whether or not the queue is empty.
 */

template <typename TValue, typename TSpec>
inline bool
empty(ConcurrentQueue<TValue, TSpec> const & me)
{
    typedef ConcurrentQueue<TValue, TSpec>  TQueue;
    typedef typename TQueue::TLockTag       TLockTag;

    ScopedWriteLock<ReadWriteLock, TLockTag> writeLock(me.lock);
    return me.headPos == me.tailPos;
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

/*!
 * @fn ConcurrentQueue#length
 * @brief Returns the size of a queue.
 *
 * @signature TSize length(queue);
 *
 * @param[in] queue The queue to query for its size.
 *
 * @return TSize The number of elements in the queue.
 */

template <typename TValue, typename TSpec>
inline typename Size<ConcurrentQueue<TValue, TSpec> >::Type
length(ConcurrentQueue<TValue, TSpec> const & me)
{
    typedef ConcurrentQueue<TValue, TSpec>  TQueue;
    typedef typename TQueue::TLockTag       TLockTag;
    typedef typename Size<TQueue>::Type     TSize;

    ScopedWriteLock<ReadWriteLock, TLockTag> writeLock(me.lock);
    TSize mask = me.roundSize - 1;
    if ((me.headPos & mask) <= (me.tailPos & mask))
        return me.tailPos - me.headPos;
    else
        return me.tailPos - me.headPos - (me.roundSize - capacity(me.data));
}

// ----------------------------------------------------------------------------
// Function capacity()
// ----------------------------------------------------------------------------

/*!
 * @fn ConcurrentQueue#capacity
 * @brief Returns the capacity of a queue.
 *
 * @signature TSize capacity(queue);
 *
 * @param[in] queue The queue to query for its capacity.
 *
 * @return TSize Returns the capacity of the queue.
 *
 * The capacity is the number of elements that can be enqueued at the same time without reallocating memory.
 */

template <typename TValue, typename TSpec>
inline typename Size<ConcurrentQueue<TValue, TSpec> >::Type
capacity(ConcurrentQueue<TValue, TSpec> const & me)
{
    typedef ConcurrentQueue<TValue, TSpec>  TQueue;
    typedef typename TQueue::TLockTag       TLockTag;

    ScopedReadLock<ReadWriteLock, TLockTag> writeLock(me.lock);
    return capacity(me.data) - 1;
}

// ----------------------------------------------------------------------------
// Function tryPopFront()
// ----------------------------------------------------------------------------

/*!
 * @fn ConcurrentQueue#tryPopFront
 * @headerfile <seqan/parallel.h>
 * @brief Try to dequeue a value from a queue.
 *
 * @signature bool tryPopFront(result, queue[, parallelTag]);
 *
 *
 * @param[in,out] queue       A queue.
 * @param[out]    result      The dequeued value (if available).
 * @param[in]     parallelTag The concurrency scheme. If multiple threads dequeue values concurrently this tag must be
 *                            @link ParallelismTags#Parallel @endlink. The more efficient @link ParallelismTags#Serial
 *                            @endlink tag can only be used if one thread calls <tt>popFront</tt> at a time.
 *                            Default is @link ParallelismTags#Parallel @endlink.
 * @return        bool        Returns <tt>true</tt> if a value could be dequeued and <tt>false</tt> otherwise.
 */

//
//  [  ?  ]  [  4  ]  [  3  ]  [  8  ]  [  0  ]  [  x  ]  [  ?  ]
//                       |                          ^
//                       v                          |
//             head            headRead   tail  tailWrite
//
// empty = (head == tail)
// full = (tail + 1 == head)
//
// valid data between  [headRead, tail)
// currently filled    [tail, tailWrite)
// currently removed   [head, headRead)

template <typename TValue2, typename TValue, typename TSpec, typename TParallel>
inline bool
tryPopFront(TValue2 & result, ConcurrentQueue<TValue, TSpec> & me, Tag<TParallel> parallelTag)
{
    typedef ConcurrentQueue<TValue, TSpec>              TQueue;
    typedef typename TQueue::TLockTag                   TLockTag;
    typedef typename Host<TQueue>::Type                 TString;
    typedef typename Size<TString>::Type                TSize;
    typedef typename Iterator<TString, Standard>::Type  TIter;

    // try to extract a value
    ScopedReadLock<ReadWriteLock, TLockTag> readLock(me.lock);

    TSize cap = capacity(me.data);
    TSize roundSize = me.roundSize;
    TSize headReadPos;
    TSize newHeadReadPos;
    SpinDelay spinDelay;

    // wait for queue to become filled
    while (true)
    {
        headReadPos = me.headReadPos;
        TSize tailPos = me.tailPos;

        SEQAN_ASSERT_LEQ(headReadPos, tailPos);

        // return if queue is empty
        if (headReadPos == tailPos)
            return false;

        newHeadReadPos = _cyclicInc(headReadPos, cap, roundSize);

        if (atomicCasBool(me.headReadPos, headReadPos, newHeadReadPos, parallelTag))
            break;

        waitFor(spinDelay);
    }

    // extract value and destruct it in the data string
    TIter it = begin(me.data, Standard()) + (headReadPos & (roundSize - 1));
    std::swap(result, *it);
    valueDestruct(it);

    // wait for pending previous reads and synchronize headPos to headReadPos
    spinCas(me.headPos, headReadPos, newHeadReadPos);

    return true;
}

template <typename TValue, typename TSpec>
inline bool
tryPopFront(TValue & result, ConcurrentQueue<TValue, TSpec> & me)
{
    return tryPopFront(result, me, typename DefaultParallelSpec<ConcurrentQueue<TValue, TSpec> >::Type());
}

// ----------------------------------------------------------------------------
// Function waitForWriters()
// ----------------------------------------------------------------------------

/*!
 * @fn ConcurrentQueue#waitForWriters
 * @brief Wait for writers to register.
 *
 * @signature void waitForWriters(queue, writerCount);
 *
 * @param[in] queue       A queue.
 * @param[in] writerCount The minimal required number of registered writers, see @link
 *                        ConcurrentQueue#lockWriting @endlink.
 *
 * If the values are dequeued with @link ConcurrentQueue#popFront2 @endlink,
 * this function is a barrier for all writers to set up completely and should be called before calling @link
 * ConcurrentQueue#appendValue @endlink the first time.
 */

template <typename TValue, typename TSpec>
inline void
waitForWriters(ConcurrentQueue<TValue, TSpec> & me, unsigned writerCount)
{
    SpinDelay spinDelay;
    while (me.writerCount < writerCount)
    {
        waitFor(spinDelay);
    }
    me.virgin = false;
}

/*!
 * @fn ConcurrentQueue#waitForFirstValue
 * @brief Wait for writers to enqueue the first value.
 *
 * @signature void waitForFirstValue(queue);
 *
 * @param[in] queue       A queue.
 *
 * If the values are dequeued with @link ConcurrentQueue#popFront2 @endlink,
 * this function is a barrier for all readers to wait until all writers are set up completely and should be called
 * before calling @link ConcurrentQueue#popFront2 @endlink the first time.
 */

template <typename TValue, typename TSpec>
inline void
waitForFirstValue(ConcurrentQueue<TValue, TSpec> & me)
{
    spinWhileEq(me.virgin, true);
}

// ----------------------------------------------------------------------------
// Function popFront()
// ----------------------------------------------------------------------------

/*!
 * @fn ConcurrentQueue#popFront
 * @headerfile <seqan/parallel.h>
 * @brief Dequeue a value from a queue.
 *
 * @signature bool popFront(result, queue[, parallelTag]);
 *
 *
 * @param[in,out] queue       A queue.
 * @param[out]    result      The dequeued value. If the queue is empty but writers are available the thread spinlocks
 *                            until a value becomes available.
 * @param[in]     parallelTag The concurrency scheme. If multiple threads dequeue values concurrently this tag must be
 *                            @link ParallelismTags#Parallel @endlink. The more efficient @link ParallelismTags#Serial
 *                            @endlink tag can only be used if one thread calls <tt>popFront</tt> at a time.
 *                            Default is @link ParallelismTags#Parallel @endlink.
 * @return        bool        Returns <tt>true</tt> if a value could be dequeued or <tt>false</tt> if no writer is
 *                            available, see @link ConcurrentQueue#waitForWriters @endlink.
 */

// returns if no writer is locked the queue and queue is empty
template <typename TValue, typename TSpec, typename TParallel>
inline bool
popFront(TValue & result, ConcurrentQueue<TValue, TSpec> & me, Tag<TParallel> parallelTag)
{
    SpinDelay spinDelay;

    while (me.writerCount != 0)
    {
        if (tryPopFront(result, me, parallelTag))
            return true;
        waitFor(spinDelay);
    }
    // we have to give it another try if the queue was empty inside the loop
    // but after the check a writer pushes a value and zeroes the writerCount
    return tryPopFront(result, me);
}

template <typename TValue, typename TSpec>
inline bool
popFront(TValue & result, ConcurrentQueue<TValue, TSpec> & me)
{
    return popFront(result, me, typename DefaultParallelSpec<ConcurrentQueue<TValue, TSpec> >::Type());
}

/*!
 * @fn ConcurrentQueue#popFront2
 * @headerfile <seqan/parallel.h>
 * @brief Dequeue a value from a queue.
 *
 * @signature TValue popFront(queue[, parallelTag]);
 *
 *
 * @param[in,out] queue       A queue.
 * @param[in]     parallelTag The concurrency scheme. If multiple threads dequeue values concurrently this tag must be
 *                            @link ParallelismTags#Parallel @endlink. The more efficient @link ParallelismTags#Serial
 *                            @endlink tag can only be used if one thread calls <tt>popFront</tt> at a time.
 *                            Default is @link ParallelismTags#Parallel @endlink.
 * @return        TValue      The dequeued value. If the queue is empty the thread spinlocks until a value becomes
 *                            available.
 */

template <typename TValue, typename TSpec, typename TParallel>
inline TValue SEQAN_FORWARD_RETURN
popFront(ConcurrentQueue<TValue, TSpec> & me, Tag<TParallel> parallelTag)
{
    TValue result;
    SpinDelay spinDelay;

    while (!tryPopFront(result, me, parallelTag))
    {
        waitFor(spinDelay);
    }
    return SEQAN_MOVE(result);
}

template <typename TValue, typename TSpec>
inline TValue SEQAN_FORWARD_RETURN
popFront(ConcurrentQueue<TValue, TSpec> & me)
{
    return popFront(me, typename DefaultParallelSpec<ConcurrentQueue<TValue, TSpec> >::Type());
}

template <typename TValue, typename TSpec, typename TValue2>
inline bool
_queueOverflow(ConcurrentQueue<TValue, TSpec> & me, TValue2 SEQAN_FORWARD_CARG, Insist)
{
    ignoreUnusedVariableWarning(me);
    SEQAN_ASSERT_GT(capacity(me.data), 1u);
    return false;
}

template <typename TValue, typename TSpec, typename TValue2>
inline bool
_queueOverflow(ConcurrentQueue<TValue, TSpec> & me, TValue2 SEQAN_FORWARD_CARG, Limit)
{
    ignoreUnusedVariableWarning(me);
    SEQAN_ASSERT_GT(capacity(me.data), 1u);
    return false;
}

template <typename TValue, typename TSpec, typename TValue2, typename TExpand>
inline bool
_queueOverflow(ConcurrentQueue<TValue, TSpec> & me,
               TValue2 SEQAN_FORWARD_CARG val,
               Tag<TExpand> expandTag)
{
    typedef ConcurrentQueue<TValue, TSpec>                  TQueue;
    typedef typename TQueue::TLockTag                       TLockTag;
    typedef typename Host<TQueue>::Type                     TString;
    typedef typename Size<TString>::Type                    TSize;
    typedef typename Iterator<TString, Standard>::Type      TIter;

    bool queueIsResizable = !IsSameType<TLockTag, Limit>::VALUE;
    ignoreUnusedVariableWarning(queueIsResizable);
    SEQAN_ASSERT(queueIsResizable);

    // try to extend capacity

    ScopedWriteLock<ReadWriteLock, TLockTag> writeLock(me.lock);
    TSize cap = capacity(me.data);
    TSize roundSize = me.roundSize;
    TSize headPos = me.headPos;
    TSize tailPos = me.tailPos;

    SEQAN_ASSERT_EQ(tailPos, me.tailWritePos);
    SEQAN_ASSERT_EQ(headPos, me.headReadPos);

    bool valueWasAppended = false;

    // did we reach the capacity limit (another thread could have done the upgrade already)?
    if (_cyclicInc((TSize)tailPos, cap, roundSize) >= headPos + roundSize)
    {
        if (cap != 0)
        {
            TIter it = begin(me.data, Standard()) + (tailPos & (roundSize - 1));
            valueConstruct(it, SEQAN_FORWARD(TValue2, val));
            tailPos = headPos + roundSize;
            valueWasAppended = true;
        }

        SEQAN_ASSERT_EQ(tailPos, headPos + roundSize);

        // get positions of head/tail in current data sequence
        TSize headIdx = headPos & (roundSize - 1);
        TSize tailIdx = tailPos & (roundSize - 1);

        // increase capacity
        _setLength(me.data, cap);
        reserve(me.data, cap + 1, expandTag);
        TSize delta = capacity(me.data) - cap;
        roundSize = (TSize)1 << (log2(capacity(me.data) - 1) + 1);

        // create a gap of delta many values between tail and head
        _clearSpace(me.data, delta, headIdx, headIdx, expandTag);
        if (cap != 0)
        {
            me.headReadPos = me.headPos = headIdx + delta;
            me.tailWritePos = me.tailPos = tailIdx + roundSize;
        }
        me.roundSize = roundSize;
    }
    return valueWasAppended;
}

// ----------------------------------------------------------------------------
// Function appendValue()
// ----------------------------------------------------------------------------

/*!
 * @fn ConcurrentQueue#appendValue
 * @headerfile <seqan/parallel.h>
 * @brief Enqueue a value to a queue.
 *
 * @signature void appendValue(queue, val[, expandTag[, parallelTag]);
 *
 *
 * @param[in,out] queue       A queue.
 * @param[in]     val         The value to enqueue.
 * @param[in]     expandTag   The overflow strategy. If @link OverflowStrategyTags#Generous @endlink the queue will be
 *                            automatically resized if the capacity is exceeded, otherwise the thread spinlocks until
 *                            the element can be enqueued.
 *                            Default is the @link DefaultOverflowImplicit @endlink result for the <tt>queue</tt> type.
 * @param[in]     parallelTag The concurrency scheme. If multiple threads enqueue values concurrently this tag must be
 *                            @link ParallelismTags#Parallel @endlink. The more efficient @link ParallelismTags#Serial
 *                            @endlink tag can only be used if one thread calls <tt>appendValue</tt> at a time.
 *                            Default is @link ParallelismTags#Parallel @endlink.
 */

template <typename TValue, typename TSpec, typename TValue2, typename TExpand, typename TParallel>
inline void
appendValue(ConcurrentQueue<TValue, TSpec> & me,
            TValue2 SEQAN_FORWARD_CARG val,
            Tag<TExpand> expandTag,
            Tag<TParallel> parallelTag)
{
    typedef ConcurrentQueue<TValue, TSpec>              TQueue;
    typedef typename TQueue::TLockTag                   TLockTag;
    typedef typename Host<TQueue>::Type                 TString;
    typedef typename Size<TString>::Type                TSize;
    typedef typename Iterator<TString, Standard>::Type  TIter;

    SpinDelay spinDelay;

    while (true)
    {
        // try to append the value
        {
            ScopedReadLock<ReadWriteLock, TLockTag> readLock(me.lock);
            TSize cap = capacity(me.data);
            TSize roundSize = me.roundSize;

            while (true)
            {
                TSize tailWritePos = me.tailWritePos;
                TSize newTailWritePos = _cyclicInc(tailWritePos, cap, roundSize);
                TSize headPos = me.headPos;

                SEQAN_ASSERT_LEQ(newTailWritePos, headPos + roundSize + 1);

                // break if we have a wrap around, i.e. queue is full
                if (newTailWritePos >= headPos + roundSize)
                    break;

                if (atomicCasBool(me.tailWritePos, tailWritePos, newTailWritePos, parallelTag))
                {
                    TIter it = begin(me.data, Standard()) + (tailWritePos & (roundSize - 1));
                    valueConstruct(it, SEQAN_FORWARD(TValue2, val));

                    // wait for pending previous writes and synchronize tailPos to tailWritePos
                    spinCas(me.tailPos, tailWritePos, newTailWritePos);

                    return;
                }

                waitFor(spinDelay);
            }
        }

        // if possible extend capacity and return (spin loop otherwise)
        if (_queueOverflow(me, SEQAN_FORWARD(TValue2, val), expandTag))
            return;

        waitFor(spinDelay);
    }
}

template <typename TValue, typename TSpec, typename TValue2, typename TExpand>
inline void
appendValue(ConcurrentQueue<TValue, TSpec> & me,
            TValue2 SEQAN_FORWARD_CARG val,
            Tag<TExpand> expandTag)
{
    appendValue(me, SEQAN_FORWARD(TValue2, val), expandTag, typename DefaultParallelSpec<ConcurrentQueue<TValue, TSpec> >::Type());
}

}  // namespace seqan

#endif  // #ifndef SEQAN_PARALLEL_PARALLEL_QUEUE_H_
