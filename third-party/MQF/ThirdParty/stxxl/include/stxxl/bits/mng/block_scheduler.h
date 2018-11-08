/***************************************************************************
 *  include/stxxl/bits/mng/block_scheduler.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2010-2011 Raoul Steffen <R-Steffen@gmx.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_MNG_BLOCK_SCHEDULER_HEADER
#define STXXL_MNG_BLOCK_SCHEDULER_HEADER

#include <stack>
#include <queue>
#include <limits>

#include <stxxl/bits/mng/block_manager.h>
#include <stxxl/bits/mng/typed_block.h>
#include <stxxl/bits/common/addressable_queues.h>

STXXL_BEGIN_NAMESPACE

//! Virtualization of a block of data.
//! Holds information for allocating and swapping. To use in cooperation with block_scheduler.
//!
//! A swappable_block can be uninitialized, i.e. it holds no data.
//! When access is required, is has to be acquired first, and released afterwards, so it can be swapped in and out as required.
//! If the stored data is no longer needed, it can get uninitialized, freeing both internal and external memory.
//! \tparam ValueType type of contained objects (POD with no references to internal memory).
//! \tparam BlockSize Number of objects in one block.
//!         BlockSize*sizeof(ValueType) must be divisible by 4096.
template <typename ValueType, unsigned BlockSize>
class swappable_block
{
protected:
    static const unsigned_type raw_block_size = BlockSize * sizeof(ValueType);

public:
    typedef typed_block<raw_block_size, ValueType> internal_block_type;
    typedef typename internal_block_type::bid_type external_block_type;

protected:
    external_block_type external_data;      //!external_data.valid if no associated space on disk
    internal_block_type* internal_data;     //NULL if there is no internal memory reserved
    bool dirty;
    int_type reference_count;

    static unsigned_type disk_allocation_offset;

    void get_external_block()
    { block_manager::get_instance()->new_block(striping(), external_data, ++disk_allocation_offset); }

    void free_external_block()
    {
        block_manager::get_instance()->delete_block(external_data);
        external_data = external_block_type(); // make invalid
    }

public:
    //! Create in uninitialized state.
    swappable_block()
        : external_data() /*!valid*/, internal_data(0), dirty(false), reference_count(0) { }

    ~swappable_block() { }

    //! If it has an internal_block. The internal_block implicitly holds valid data.
    bool is_internal() const
    { return (internal_data != NULL); }

    //! If the external_block does not hold valid data.
    bool is_dirty() const
    { return dirty; }

    //! If it has an external_block.
    bool has_external_block() const
    { return external_data.valid(); }

    //! If it has an external_block that holds valid data.
    bool is_external() const
    { return has_external_block() && ! is_dirty(); }

    //! If it is acquired.
    bool is_acquired() const
    { return reference_count > 0; }

    //! If it holds an internal_block but does not need it.
    bool is_evictable() const
    { return ! is_acquired() && is_internal(); }

    //! If it has some valid data (in- or external).
    bool is_initialized() const
    { return is_internal() || is_external(); }

    //! Invalidate external data if true.
    //! \return is_dirty()
    bool make_dirty_if(const bool make_dirty)
    {
        assert(is_acquired());
        return dirty = make_dirty || dirty;
    }

    //! Acquire the block, i.e. add a reference. Has to be internal.
    //! \return A reference to the data-block.
    internal_block_type & acquire()
    {
        assert(is_internal());
        ++reference_count;
        return *internal_data;
    }

    //! Release the block, i.e. subduct a reference. Has to be acquired.
    void release()
    {
        assert(is_acquired());
        --reference_count;
    }

    //! Get a reference to the data-block. Has to be acquired.
    const internal_block_type & get_internal_block() const
    {
        assert(is_acquired());
        return *internal_data;
    }

    //! Get a reference to the data-block. Has to be acquired.
    internal_block_type & get_internal_block()
    {
        assert(is_acquired());
        return *internal_data;
    }

    //! Fill block with default data, is supposed to be overwritten by subclass. Has to be internal.
    void fill_default() { }

    //! Read asyncronusly from external_block to internal_block. Has to be internal and have an external_block.
    //! \return A request pointer to the I/O.
    request_ptr read_async(completion_handler on_cmpl = completion_handler())
    {
        assert(is_internal());
        assert(has_external_block());
        #ifdef RW_VERBOSE
        STXXL_MSG("reading block");
        #endif
        dirty = false;
        return internal_data->read(external_data, on_cmpl);
    }

    //! Read synchronously from external_block to internal_block. Has to be internal and have an external_block.
    void read_sync()
    { read_async()->wait(); }

    //! Write asyncronusly from internal_block to external_block if necessary.
    //! \return A request pointer to the I/O, an invalid request pointer if not necessary.
    request_ptr clean_async(completion_handler on_cmpl = completion_handler())
    {
        if (! is_dirty())
            return request_ptr();
        if (! has_external_block())
            get_external_block();
        #ifdef RW_VERBOSE
        STXXL_MSG("writing block");
        #endif
        dirty = false;
        return internal_data->write(external_data, on_cmpl);
    }

    //! Write synchronously from internal_block to external_block if necessary.
    void clean_sync()
    {
        request_ptr rp = clean_async();
        if (rp.valid())
            rp->wait();
    }

    //! Attach an internal_block, making the block internal. Has to be not internal.
    void attach_internal_block(internal_block_type* iblock)
    {
        assert(! is_internal());
        internal_data = iblock;
    }

    //! Detach the internal_block. Writes to external_block if necessary. Has to be evictable.
    //! \return A pointer to the internal_block.
    internal_block_type * detach_internal_block()
    {
        assert(is_evictable());
        clean_sync();
        internal_block_type* iblock = internal_data;
        internal_data = 0;
        return iblock;
    }

    //! Bring the block in uninitialized state, freeing external and internal memory.
    //! Returns a pointer to the internal_block, NULL if it had none.
    //! \return A pointer to the freed internal_block, NULL if it had none.
    internal_block_type * deinitialize()
    {
        assert(! is_acquired());
        dirty = false;
        // free external_block (so that it becomes invalid and the disk-space can be used again)
        if (has_external_block())
            free_external_block();
        // free internal_block
        internal_block_type* iblock = internal_data;
        internal_data = 0;
        return iblock;
    }

    //! Set the external_block that holds the swappable_block's data. The block gets initialized with it.
    //! \param eblock The external_block holding initial data.
    void initialize(external_block_type eblock)
    {
        assert(! is_initialized());
        external_data = eblock;
    }

    //! Extract the swappable_blocks data in an external_block. The block gets uninitialized.
    //! \return The external_block that holds the swappable_block's data.
    external_block_type extract_external_block()
    {
        assert(! is_internal());
        external_block_type eblock = external_data;
        external_data = external_block_type();
        return eblock;
    }
};

template <typename ValueType, unsigned BlockSize>
unsigned_type swappable_block<ValueType, BlockSize>::disk_allocation_offset = 0;

template <class SwappableBlockType>
class block_scheduler_algorithm;

template <class SwappableBlockType>
class block_scheduler_algorithm_online_lru;

//! Schedules swapping of blocks and provides blocks for temporary storage.
//!
//! In simple mode, it tries to save I/Os through caching only.
//! In simulation mode, it records access patterns into a prediction sequence.
//! The prediction sequence can then be used for prefetching in the (offline) execute mode.
//! This will only work for algorithms with deterministic, data oblivious access patterns.
//! In simulation mode, no I/O is performed; the data provided is accessible but undefined.
//! In execute mode, it does caching, prefetching, and possibly other optimizations.
//! \tparam SwappableBlockType Type of swappable_blocks to manage. Can be some specialized subclass.
template <class SwappableBlockType>
class block_scheduler : private noncopyable
{
protected:
    // tuning-parameter: To acquire blocks, internal memory has to be allocated.
    // This constant limits the number of internal_blocks allocated at once.
    static const int_type max_internal_blocks_alloc_at_once;

    typedef int_type time_type;

public:
    typedef typename SwappableBlockType::internal_block_type internal_block_type;
    typedef typename SwappableBlockType::external_block_type external_block_type;
    typedef typename std::vector<SwappableBlockType>::size_type swappable_block_identifier_type;

    /*/! Mode the block scheduler currently works in
    enum mode
    {
        online,         //serve requests immediately, without any prediction, LRU caching
        simulation,     //record prediction sequence only, do not serve requests, (returned blocks must not be accessed)
        offline_lfd,    //serve requests based on prediction sequence, using longest-forward-distance caching
        offline_lfd_prefetch     //serve requests based on prediction sequence, using longest-forward-distance caching, and prefetching
    };*/

public:
    // -------- prediction_sequence -------
    enum block_scheduler_operation
    {
        op_acquire,
        op_acquire_uninitialized,
        op_release,
        op_release_dirty,
        op_deinitialize,
        op_initialize,
        op_extract_external_block
    };

    struct prediction_sequence_element
    {
        block_scheduler_operation op;
        swappable_block_identifier_type id;
        time_type time;

        prediction_sequence_element(block_scheduler_operation op,
                                    swappable_block_identifier_type id, int_type time)
            : op(op), id(id), time(time) { }
    };

    typedef std::list<prediction_sequence_element> prediction_sequence_type;
    // ---- end prediction_sequence -------

protected:
    template <class SBT>
    friend class block_scheduler_algorithm;

    const int_type max_internal_blocks;
    int_type remaining_internal_blocks;
    //! Stores pointers to arrays of internal_blocks. Used to deallocate them only.
    std::stack<internal_block_type*> internal_blocks_blocks;
    //! holds free internal_blocks with attributes reset.
    std::stack<internal_block_type*> free_internal_blocks;
    //! temporary blocks that will not be needed after algorithm termination.
    mutable std::vector<SwappableBlockType> swappable_blocks;
    //! holds indices of free swappable_blocks with attributes reset.
    std::priority_queue<swappable_block_identifier_type, std::vector<swappable_block_identifier_type>,
                        std::greater<swappable_block_identifier_type> > free_swappable_blocks;
    block_manager* bm;
    block_scheduler_algorithm<SwappableBlockType>* algo;

    //! Get an internal_block from the freelist or a newly allocated one if available.
    //! \return Pointer to the internal_block. NULL if none available.
    internal_block_type * get_free_internal_block()
    {
        if (! free_internal_blocks.empty())
        {
            // => there are internal_blocks in the free-list
            internal_block_type* iblock = free_internal_blocks.top();
            free_internal_blocks.pop();
            return iblock;
        }
        else if (remaining_internal_blocks > 0)
        {
            // => more internal_blocks can be allocated
            int_type num_blocks = std::min(max_internal_blocks_alloc_at_once, remaining_internal_blocks);
            remaining_internal_blocks -= num_blocks;
            internal_block_type* iblocks = new internal_block_type[num_blocks];
            internal_blocks_blocks.push(iblocks);
            for (int_type i = num_blocks - 1; i > 0; --i)
                free_internal_blocks.push(iblocks + i);
            return iblocks;
        }
        else
        {
            // => no internal_block available
            return 0;
        }
    }

    //! Return an internal_block to the freelist.
    void return_free_internal_block(internal_block_type* iblock)
    { free_internal_blocks.push(iblock); }

public:
    //! Create a block_scheduler with empty prediction sequence in simple mode.
    //! \param max_internal_memory Amount of internal memory (in bytes) the scheduler is allowed to use for acquiring, prefetching and caching.
    explicit block_scheduler(const int_type max_internal_memory)
        : max_internal_blocks(div_ceil(max_internal_memory, sizeof(internal_block_type))),
          remaining_internal_blocks(max_internal_blocks),
          bm(block_manager::get_instance()),
          algo(0)
    {
        algo = new block_scheduler_algorithm_online_lru<SwappableBlockType>(*this);
    }

    ~block_scheduler()
    {
        delete algo;
        int_type num_freed_internal_blocks = 0;
        if (free_swappable_blocks.size() != swappable_blocks.size())
        {
            // => not all swappable_blocks are free, at least deinitialize them
            STXXL_ERRMSG("not all swappable_blocks are free, those not acquired will be deinitialized");
            // evictable_blocks would suffice
            for (typename std::vector<SwappableBlockType>::iterator it = swappable_blocks.begin();
                 it != swappable_blocks.end(); ++it)
            {
                if (! it->is_acquired() && it->deinitialize())
                    // count internal_blocks that get freed
                    num_freed_internal_blocks++;
            }
        }
        if (int_type nlost = (max_internal_blocks - remaining_internal_blocks)
                             - (free_internal_blocks.size() + num_freed_internal_blocks)) {
            STXXL_ERRMSG(nlost << " internal_blocks are lost. They will get deallocated.");
        }
        while (! internal_blocks_blocks.empty())
        {
            delete[] internal_blocks_blocks.top();
            internal_blocks_blocks.pop();
        }
    }

    //! Acquire the given block.
    //! Has to be in pairs with release. Pairs may be nested and interleaved.
    //! \return Reference to the block's data.
    //! param sbid Swappable block to acquire.
    internal_block_type & acquire(const swappable_block_identifier_type sbid, const bool uninitialized = false)
    { return algo->acquire(sbid, uninitialized); }

    //! Release the given block.
    //! Has to be in pairs with acquire. Pairs may be nested and interleaved.
    //! \param sbid Swappable block to release.
    //! \param dirty If the data has been changed, invalidating possible data in external storage.
    void release(const swappable_block_identifier_type sbid, const bool dirty)
    { algo->release(sbid, dirty); }

    //! Drop all data in the given block, freeing in- and external memory.
    void deinitialize(const swappable_block_identifier_type sbid)
    { algo->deinitialize(sbid); }

    //! Initialize the swappable_block with the given external_block.
    //!
    //! It will use the the external_block for swapping and take care about
    //! it's deallocation. Has to be uninitialized.
    //! \param sbid identifier to the swappable_block
    //! \param eblock external_block a.k.a. bid
    void initialize(const swappable_block_identifier_type sbid, external_block_type eblock)
    { algo->initialize(sbid, eblock); }

    //! Deinitialize the swappable_block and return it's contents in an external_block.
    //!
    //! \param sbid identifier to the swappable_block
    //! \return external_block a.k.a. bid
    external_block_type extract_external_block(const swappable_block_identifier_type sbid)
    { return algo->extract_external_block(sbid); }

    //! check if the swappable_block is initialized.
    //! \param sbid identifier to the swappable_block
    //! \return if the swappable_block is initialized
    bool is_initialized(const swappable_block_identifier_type sbid) const
    { return algo->is_initialized(sbid); }

    //! Record a timestep in the prediction sequence to seperate consecutive
    //! acquire rsp. release-operations. Has an effect only in simulation mode.
    void explicit_timestep()
    { algo->explicit_timestep(); }

    //! Get a const reference to given block's data. Block has to be already acquired.
    //! \param sbid Swappable block to access.
    internal_block_type & get_internal_block(const swappable_block_identifier_type sbid) const
    { return swappable_blocks[sbid].get_internal_block(); }

    //! Allocate an uninitialized swappable_block.
    //! \return An identifier of the block.
    swappable_block_identifier_type allocate_swappable_block()
    {
        swappable_block_identifier_type sbid;
        if (free_swappable_blocks.empty())
        {
            // create new swappable_block
            sbid = swappable_blocks.size();
            swappable_blocks.resize(sbid + 1);
            algo->swappable_blocks_resize(sbid + 1);
        }
        else
        {
            // take swappable_block from freelist
            sbid = free_swappable_blocks.top();
            free_swappable_blocks.pop();
        }
        return sbid;
    }

    //! Free given no longer used temporary swappable_block.
    //! \param sbid Temporary swappable_block to free.
    void free_swappable_block(const swappable_block_identifier_type sbid)
    {
        deinitialize(sbid);
        free_swappable_blocks.push(sbid);
    }

    //! Returns if simulation mode is on, i.e. if a prediction sequence is being recorded.
    //! \return If simulation mode is on.
    bool is_simulating() const
    { return algo->is_simulating(); }

    //! Switch the used algorithm, e.g. to simulation etc..
    //! \param new_algo Pointer to the new algorithm object. Has to be instantiated to the block scheduler (or the old algorithm object).
    //! \return Pointer to the old algorithm object.
    block_scheduler_algorithm<SwappableBlockType> * switch_algorithm_to(block_scheduler_algorithm<SwappableBlockType>* new_algo)
    {
        assert(&new_algo->bs == this);
        block_scheduler_algorithm<SwappableBlockType>* old_algo = algo;
        algo = new_algo;
        return old_algo;
    }

    //! Return the current algorithm.
    block_scheduler_algorithm<SwappableBlockType> * get_current_algorithm() const
    {
        return algo;
    }

    //! Get the prediction_sequence.
    //! \return reference to the prediction_sequence
    const prediction_sequence_type & get_prediction_sequence() const
    { return algo->get_prediction_sequence(); }

    void flush()
    {
        std::vector<request_ptr> requests;
        while (! algo->evictable_blocks_empty())
        {
            swappable_block_identifier_type sbid = algo->evictable_blocks_pop();
            request_ptr rq = swappable_blocks[sbid].clean_async();
            if (rq.valid())
                requests.push_back(rq);
            return_free_internal_block(swappable_blocks[sbid].detach_internal_block());
        }
        for (typename std::vector<request_ptr>::reverse_iterator it = requests.rbegin();
             it != requests.rend(); ++it)
        {
            (*it)->wait();
        }
    }
};

template <class SwappableBlockType>
const int_type block_scheduler<SwappableBlockType>::max_internal_blocks_alloc_at_once = 128;

//! Interface of a block scheduling algorithm.
template <class SwappableBlockType>
class block_scheduler_algorithm : private noncopyable
{
protected:
    typedef block_scheduler<SwappableBlockType> block_scheduler_type;
    typedef typename block_scheduler_type::internal_block_type internal_block_type;
    typedef typename block_scheduler_type::external_block_type external_block_type;
    typedef typename block_scheduler_type::swappable_block_identifier_type swappable_block_identifier_type;
    typedef typename block_scheduler_type::prediction_sequence_type prediction_sequence_type;
    typedef typename block_scheduler_type::time_type time_type;

public:
    block_scheduler_type& bs;

protected:
    std::vector<SwappableBlockType>& swappable_blocks;
    prediction_sequence_type prediction_sequence;

    block_scheduler_algorithm * get_algorithm_from_block_scheduler()
    { return bs.algo; }

    //! Get an internal_block from the block_scheduler.
    //! \return Pointer to the internal_block. NULL if none available.
    internal_block_type * get_free_internal_block_from_block_scheduler()
    { return bs.get_free_internal_block(); }

    //! Return an internal_block to the block_scheduler.
    void return_free_internal_block_to_block_scheduler(internal_block_type* iblock)
    { bs.return_free_internal_block(iblock); }

public:
    block_scheduler_algorithm(block_scheduler_type& bs)
        : bs(bs),
          swappable_blocks(bs.swappable_blocks)
    { }

    block_scheduler_algorithm(block_scheduler_algorithm* old)
        : bs(old->bs),
          swappable_blocks(bs.swappable_blocks)
    { }

    virtual ~block_scheduler_algorithm() { }

    virtual bool evictable_blocks_empty() = 0;
    virtual swappable_block_identifier_type evictable_blocks_pop() = 0;
    virtual void swappable_blocks_resize(swappable_block_identifier_type /*size*/) { }

    virtual internal_block_type & acquire(const swappable_block_identifier_type sbid, const bool uninitialized = false) = 0;
    virtual void release(swappable_block_identifier_type sbid, const bool dirty) = 0;
    virtual void deinitialize(swappable_block_identifier_type sbid) = 0;
    virtual void initialize(swappable_block_identifier_type sbid, external_block_type eblock) = 0;
    virtual external_block_type extract_external_block(swappable_block_identifier_type sbid) = 0;

    virtual bool is_initialized(const swappable_block_identifier_type sbid) const
    { return swappable_blocks[sbid].is_initialized(); }

    virtual void explicit_timestep() { }
    virtual bool is_simulating() const
    { return false; }
    virtual const prediction_sequence_type & get_prediction_sequence() const
    { return prediction_sequence; }
};

//! Block scheduling algorithm caching via the least recently used policy (online).
template <class SwappableBlockType>
class block_scheduler_algorithm_online_lru : public block_scheduler_algorithm<SwappableBlockType>
{
protected:
    typedef block_scheduler<SwappableBlockType> block_scheduler_type;
    typedef block_scheduler_algorithm<SwappableBlockType> block_scheduler_algorithm_type;
    typedef typename block_scheduler_type::internal_block_type internal_block_type;
    typedef typename block_scheduler_type::external_block_type external_block_type;
    typedef typename block_scheduler_type::swappable_block_identifier_type swappable_block_identifier_type;

    using block_scheduler_algorithm_type::bs;
    using block_scheduler_algorithm_type::swappable_blocks;
    using block_scheduler_algorithm_type::get_algorithm_from_block_scheduler;
    using block_scheduler_algorithm_type::get_free_internal_block_from_block_scheduler;
    using block_scheduler_algorithm_type::return_free_internal_block_to_block_scheduler;

    //! Holds swappable blocks, whose internal block can be freed, i.e. that are internal but unacquired.
    addressable_fifo_queue<swappable_block_identifier_type> evictable_blocks;

    internal_block_type * get_free_internal_block()
    {
        // try to get a free internal_block
        if (internal_block_type* iblock = get_free_internal_block_from_block_scheduler())
            return iblock;
        // evict block
        assert(! evictable_blocks.empty()); // fails it there is not enough memory available
        return swappable_blocks[evictable_blocks.pop()].detach_internal_block();
    }

    void return_free_internal_block(internal_block_type* iblock)
    { return_free_internal_block_to_block_scheduler(iblock); }

    void init()
    {
        if (get_algorithm_from_block_scheduler())
            while (! get_algorithm_from_block_scheduler()->evictable_blocks_empty())
                evictable_blocks.insert(get_algorithm_from_block_scheduler()->evictable_blocks_pop());
    }

public:
    block_scheduler_algorithm_online_lru(block_scheduler_type& bs)
        : block_scheduler_algorithm_type(bs)
    { init(); }

    block_scheduler_algorithm_online_lru(block_scheduler_algorithm_type* old)
        : block_scheduler_algorithm_type(old)
    { init(); }

    virtual ~block_scheduler_algorithm_online_lru()
    {
        if (! evictable_blocks.empty())
            STXXL_ERRMSG("Destructing block_scheduler_algorithm_online that still holds evictable blocks. They get deinitialized.");
        while (! evictable_blocks.empty())
        {
            SwappableBlockType& sblock = swappable_blocks[evictable_blocks.pop()];
            if (internal_block_type* iblock = sblock.deinitialize())
                return_free_internal_block(iblock);
        }
    }

    virtual bool evictable_blocks_empty()
    { return evictable_blocks.empty(); }

    virtual swappable_block_identifier_type evictable_blocks_pop()
    { return evictable_blocks.pop(); }

    virtual internal_block_type & acquire(const swappable_block_identifier_type sbid, const bool uninitialized = false)
    {
        SwappableBlockType& sblock = swappable_blocks[sbid];
        /* acquired => internal -> increase reference count
           internal but not acquired -> remove from evictable_blocks, increase reference count
           not internal => uninitialized or external -> get internal_block, increase reference count
           uninitialized -> fill with default value
           external -> read */
        if (sblock.is_internal())
        {
            if (! sblock.is_acquired())
                // not acquired yet -> remove from evictable_blocks
                evictable_blocks.erase(sbid);
            sblock.acquire();
        }
        else if (sblock.is_initialized())
        {
            // => external but not internal
            //get internal_block
            sblock.attach_internal_block(get_free_internal_block());
            if (! uninitialized)
                //load block synchronously
                sblock.read_sync();
            sblock.acquire();
        }
        else
        {
            // => ! sblock.is_initialized()
            //get internal_block
            sblock.attach_internal_block(get_free_internal_block());
            sblock.acquire();
            //initialize new block
            if (! uninitialized)
                sblock.fill_default();
        }
        return sblock.get_internal_block();
    }

    virtual void release(swappable_block_identifier_type sbid, const bool dirty)
    {
        SwappableBlockType& sblock = swappable_blocks[sbid];
        sblock.make_dirty_if(dirty);
        sblock.release();
        if (! sblock.is_acquired())
        {
            if (sblock.is_dirty() || sblock.is_external())
                // => evictable, put in pq
                evictable_blocks.insert(sbid);
            else
                // => uninitialized, release internal block and put it in freelist
                return_free_internal_block(sblock.detach_internal_block());
        }
    }

    virtual void deinitialize(swappable_block_identifier_type sbid)
    {
        SwappableBlockType& sblock = swappable_blocks[sbid];
        if (sblock.is_evictable())
            evictable_blocks.erase(sbid);
        if (internal_block_type* iblock = sblock.deinitialize())
            return_free_internal_block(iblock);
    }

    virtual void initialize(swappable_block_identifier_type sbid, external_block_type eblock)
    {
        SwappableBlockType& sblock = swappable_blocks[sbid];
        sblock.initialize(eblock);
    }

    virtual external_block_type extract_external_block(swappable_block_identifier_type sbid)
    {
        SwappableBlockType& sblock = swappable_blocks[sbid];
        if (sblock.is_evictable())
            evictable_blocks.erase(sbid);
        if (sblock.is_internal())
            return_free_internal_block(sblock.detach_internal_block());
        return sblock.extract_external_block();
    }
};

//! Pseudo block scheduling algorithm only recording the request sequence.
template <class SwappableBlockType>
class block_scheduler_algorithm_simulation : public block_scheduler_algorithm<SwappableBlockType>
{
protected:
    typedef block_scheduler<SwappableBlockType> block_scheduler_type;
    typedef block_scheduler_algorithm<SwappableBlockType> block_scheduler_algorithm_type;
    typedef typename block_scheduler_type::internal_block_type internal_block_type;
    typedef typename block_scheduler_type::external_block_type external_block_type;
    typedef typename block_scheduler_type::swappable_block_identifier_type swappable_block_identifier_type;
    typedef typename block_scheduler_type::prediction_sequence_element prediction_sequence_element_type;
    typedef typename block_scheduler_algorithm_type::time_type time_type;

    using block_scheduler_algorithm_type::bs;
    using block_scheduler_algorithm_type::prediction_sequence;
    using block_scheduler_algorithm_type::swappable_blocks;
    using block_scheduler_algorithm_type::get_algorithm_from_block_scheduler;
    using block_scheduler_algorithm_type::get_free_internal_block_from_block_scheduler;
    using block_scheduler_algorithm_type::return_free_internal_block_to_block_scheduler;

    //! Holds swappable blocks, whose internal block can be freed, i.e. that are internal but unacquired.
    std::stack<swappable_block_identifier_type> evictable_blocks;
    time_type time_count;
    bool last_op_release;
    std::vector<int_type> reference_counts;
    internal_block_type dummy_block;

    void return_free_internal_block(internal_block_type* iblock)
    { return_free_internal_block_to_block_scheduler(iblock); }

    void init()
    {
        if (get_algorithm_from_block_scheduler())
            while (! get_algorithm_from_block_scheduler()->evictable_blocks_empty())
                evictable_blocks.push(get_algorithm_from_block_scheduler()->evictable_blocks_pop());
        for (swappable_block_identifier_type i = 0; i < reference_counts.size(); ++i)
            reference_counts[i] = swappable_blocks[i].is_initialized();
    }

public:
    block_scheduler_algorithm_simulation(block_scheduler_type& bs)
        : block_scheduler_algorithm_type(bs),
          time_count(0),
          last_op_release(false),
          reference_counts(swappable_blocks.size())
    { init(); }

    block_scheduler_algorithm_simulation(block_scheduler_algorithm_type* old)
        : block_scheduler_algorithm_type(old),
          time_count(0),
          last_op_release(false),
          reference_counts(swappable_blocks.size())
    { init(); }

    virtual ~block_scheduler_algorithm_simulation()
    {
        if (! evictable_blocks.empty())
            STXXL_ERRMSG("Destructing block_scheduler_algorithm_record_prediction_sequence that still holds evictable blocks. They get deinitialized.");
        while (! evictable_blocks.empty())
        {
            SwappableBlockType& sblock = swappable_blocks[evictable_blocks.top()];
            if (internal_block_type* iblock = sblock.deinitialize())
                return_free_internal_block(iblock);
            evictable_blocks.pop();
        }
    }

    virtual bool evictable_blocks_empty()
    { return evictable_blocks.empty(); }

    virtual swappable_block_identifier_type evictable_blocks_pop()
    {
        swappable_block_identifier_type sbid = evictable_blocks.top();
        evictable_blocks.pop();
        return sbid;
    }

    virtual internal_block_type & acquire(const swappable_block_identifier_type sbid, const bool uninitialized = false)
    {
        ++reference_counts[sbid];
        last_op_release = false;
        if (uninitialized)
            prediction_sequence.push_back(
                prediction_sequence_element_type(block_scheduler_type::op_acquire_uninitialized, sbid, time_count)
                );
        else
            prediction_sequence.push_back(
                prediction_sequence_element_type(block_scheduler_type::op_acquire, sbid, time_count)
                );
        return dummy_block;
    }

    virtual void release(swappable_block_identifier_type sbid, const bool dirty)
    {
        --reference_counts[sbid] += dirty;
        time_count += ! last_op_release;
        last_op_release = true;
        if (dirty)
            prediction_sequence.push_back(
                prediction_sequence_element_type(block_scheduler_type::op_release_dirty, sbid, time_count)
                );
        else
            prediction_sequence.push_back(
                prediction_sequence_element_type(block_scheduler_type::op_release, sbid, time_count)
                );
    }

    virtual void deinitialize(swappable_block_identifier_type sbid)
    {
        reference_counts[sbid] = false;
        prediction_sequence.push_back(
            prediction_sequence_element_type(block_scheduler_type::op_deinitialize, sbid, time_count)
            );
    }

    virtual void initialize(swappable_block_identifier_type sbid, external_block_type)
    {
        reference_counts[sbid] = true;
        prediction_sequence.push_back(
            prediction_sequence_element_type(block_scheduler_type::op_initialize, sbid, time_count)
            );
    }

    virtual external_block_type extract_external_block(swappable_block_identifier_type sbid)
    {
        reference_counts[sbid] = false;
        prediction_sequence.push_back(
            prediction_sequence_element_type(block_scheduler_type::op_extract_external_block, sbid, time_count)
            );
        return external_block_type();
    }

    virtual void swappable_blocks_resize(swappable_block_identifier_type size)
    {
        reference_counts.resize(size, 0);
    }

    virtual bool is_initialized(const swappable_block_identifier_type sbid) const
    {
        return reference_counts[sbid] > 0;
    }

    virtual void explicit_timestep()
    { ++time_count; }

    virtual bool is_simulating() const
    { return true; }
};

//! Block scheduling algorithm caching via the longest forward distance policy (offline).
template <class SwappableBlockType>
class block_scheduler_algorithm_offline_lfd : public block_scheduler_algorithm<SwappableBlockType>
{
protected:
    typedef block_scheduler<SwappableBlockType> block_scheduler_type;
    typedef block_scheduler_algorithm<SwappableBlockType> block_scheduler_algorithm_type;
    typedef typename block_scheduler_type::internal_block_type internal_block_type;
    typedef typename block_scheduler_type::external_block_type external_block_type;
    typedef typename block_scheduler_type::swappable_block_identifier_type swappable_block_identifier_type;
    typedef typename block_scheduler_algorithm_type::time_type time_type;
    typedef typename block_scheduler_type::prediction_sequence_type prediction_sequence_type;

    using block_scheduler_algorithm_type::bs;
    using block_scheduler_algorithm_type::swappable_blocks;
    using block_scheduler_algorithm_type::get_algorithm_from_block_scheduler;
    using block_scheduler_algorithm_type::get_free_internal_block_from_block_scheduler;
    using block_scheduler_algorithm_type::return_free_internal_block_to_block_scheduler;

    class priority
    {
        unsigned_type p;

    public:
        priority(const SwappableBlockType& sblock, const std::pair<bool, time_type>& t)
        {
            // p larger => evict earlier
            if (t.first)
            {
                // most significant: next use
                p = unsigned_type(t.second) << 2;
                // less significant: not dirty
                p |= unsigned_type(! sblock.is_dirty()) << 1;
                // less significant: has external_block
                p |= unsigned_type(sblock.has_external_block()) << 0;
            }
            else
            {
                // most significant: next use
                p = std::numeric_limits<unsigned_type>::max() << 2;
                // less significant: next operation: extract > accessed no more > deinitialize
                p |= unsigned_type(t.second) << 0;
            }
        }

        // less => evict earlier
        bool operator < (const priority& right) const
        { return p > right.p; }
    };

    //! Holds swappable blocks, whose internal block can be freed, i.e. that are internal but unacquired.
    addressable_priority_queue<swappable_block_identifier_type, priority> evictable_blocks;
    /*!
     * Stores for the sequence of releases extracted from the prediction_sequence:
     * (true, timestamp of the blocks next acquire) if it is acquired next
     * (false, 0) if it is deinitialized next
     * (false, 1) if it is not accessed any more
     * (false, 2) if it is extracted next
     * (false, 3) if it is initialized next
     */
    std::deque<std::pair<bool, time_type> > next_use;

    internal_block_type * get_free_internal_block()
    {
        // try to get a free internal_block
        if (internal_block_type* iblock = get_free_internal_block_from_block_scheduler())
            return iblock;
        // evict block
        assert(! evictable_blocks.empty()); // fails it there is not enough memory available
        return swappable_blocks[evictable_blocks.pop()].detach_internal_block();
    }

    void return_free_internal_block(internal_block_type* iblock)
    { return_free_internal_block_to_block_scheduler(iblock); }

    void init(block_scheduler_algorithm_type* old_algo)
    {
        std::vector<std::pair<bool, time_type> >
        blocks_next_acquire(swappable_blocks.size(), std::make_pair(false, 1));
        if (old_algo)
        {
            // precomputations for priorities: init next_acquires
            const prediction_sequence_type& ps = old_algo->get_prediction_sequence();
            for (typename prediction_sequence_type::const_reverse_iterator it = ps.rbegin(); it != ps.rend(); ++it)
            {
                switch (it->op)
                {
                case (block_scheduler_type::op_acquire):
                case (block_scheduler_type::op_acquire_uninitialized):
                    blocks_next_acquire[it->id] = std::make_pair(true, it->time);
                    break;
                case (block_scheduler_type::op_release):
                case (block_scheduler_type::op_release_dirty):
                    next_use.push_front(blocks_next_acquire[it->id]);
                    break;
                case (block_scheduler_type::op_deinitialize):
                    blocks_next_acquire[it->id] = std::make_pair(false, 0);
                    break;
                case (block_scheduler_type::op_initialize):
                    blocks_next_acquire[it->id] = std::make_pair(false, 3);
                    break;
                case (block_scheduler_type::op_extract_external_block):
                    blocks_next_acquire[it->id] = std::make_pair(false, 2);
                    break;
                }
            }
        }
        if (get_algorithm_from_block_scheduler())
        {
            while (! get_algorithm_from_block_scheduler()->evictable_blocks_empty())
            {
                // insert already evictable blocks with the right priority
                const swappable_block_identifier_type sbid = get_algorithm_from_block_scheduler()->evictable_blocks_pop();
                evictable_blocks.insert(sbid, priority(swappable_blocks[sbid], blocks_next_acquire[sbid]));
            }
        }
    }

public:
    block_scheduler_algorithm_offline_lfd(block_scheduler_type& bs)
        : block_scheduler_algorithm_type(bs)
    { init(get_algorithm_from_block_scheduler()); }

    // It is possible to keep an old simulation-algorithm object and reuse it's prediction sequence
    block_scheduler_algorithm_offline_lfd(block_scheduler_algorithm_type* old)
        : block_scheduler_algorithm_type(old)
    { init(old); }

    virtual ~block_scheduler_algorithm_offline_lfd()
    {
        if (! evictable_blocks.empty())
            STXXL_ERRMSG("Destructing block_scheduler_algorithm_offline_lfd that still holds evictable blocks. They get deinitialized.");
        while (! evictable_blocks.empty())
        {
            SwappableBlockType& sblock = swappable_blocks[evictable_blocks.pop()];
            if (internal_block_type* iblock = sblock.deinitialize())
                return_free_internal_block(iblock);
        }
    }

    virtual bool evictable_blocks_empty()
    { return evictable_blocks.empty(); }

    virtual swappable_block_identifier_type evictable_blocks_pop()
    { return evictable_blocks.pop(); }

    virtual internal_block_type & acquire(const swappable_block_identifier_type sbid, const bool uninitialized = false)
    {
        SwappableBlockType& sblock = swappable_blocks[sbid];
        /* acquired => internal -> increase reference count
           internal but not acquired -> remove from evictable_blocks, increase reference count
           not intern => uninitialized or external -> get internal_block, increase reference count
           uninitialized -> fill with default value
           external -> read */
        if (sblock.is_internal())
        {
            if (! sblock.is_acquired())
                // not acquired yet -> remove from evictable_blocks
                evictable_blocks.erase(sbid);
            sblock.acquire();
        }
        else if (sblock.is_initialized())
        {
            // => external but not internal
            //get internal_block
            sblock.attach_internal_block(get_free_internal_block());
            if (! uninitialized)
                //load block synchronously
                sblock.read_sync();
            sblock.acquire();
        }
        else
        {
            // => ! sblock.is_initialized()
            //get internal_block
            sblock.attach_internal_block(get_free_internal_block());
            sblock.acquire();
            //initialize new block
            if (! uninitialized)
                sblock.fill_default();
        }
        return sblock.get_internal_block();
    }

    virtual void release(swappable_block_identifier_type sbid, const bool dirty)
    {
        if (next_use.empty())
        {
            STXXL_ERRMSG("block_scheduler_algorithm_offline_lfd got release-request but prediction sequence ended. Switching to block_scheduler_algorithm_online.");
            // switch algorithm
            block_scheduler_algorithm_type* new_algo, * old_algo;
            new_algo = new block_scheduler_algorithm_online_lru<SwappableBlockType>(bs);
            old_algo = bs.switch_algorithm_to(new_algo);
            // redirect call
            new_algo->release(sbid, dirty);
            // delete self
            delete old_algo;
            return;
        }
        SwappableBlockType& sblock = swappable_blocks[sbid];
        sblock.make_dirty_if(dirty);
        sblock.release();
        if (! sblock.is_acquired())
        {
            if (sblock.is_dirty() || sblock.is_external())
                // => evictable, put in pq
                evictable_blocks.insert(sbid, priority(swappable_blocks[sbid], next_use.front()));
            else
                // => uninitialized, release internal block and put it in freelist
                return_free_internal_block(sblock.detach_internal_block());
        }
        next_use.pop_front();
    }

    virtual void deinitialize(swappable_block_identifier_type sbid)
    {
        SwappableBlockType& sblock = swappable_blocks[sbid];
        if (sblock.is_evictable())
            evictable_blocks.erase(sbid);
        if (internal_block_type* iblock = sblock.deinitialize())
            return_free_internal_block(iblock);
    }

    virtual void initialize(swappable_block_identifier_type sbid, external_block_type eblock)
    {
        SwappableBlockType& sblock = swappable_blocks[sbid];
        sblock.initialize(eblock);
    }

    virtual external_block_type extract_external_block(swappable_block_identifier_type sbid)
    {
        SwappableBlockType& sblock = swappable_blocks[sbid];
        if (sblock.is_evictable())
            evictable_blocks.erase(sbid);
        if (sblock.is_internal())
            return_free_internal_block(sblock.detach_internal_block());
        return sblock.extract_external_block();
    }
};

//! Block scheduling algorithm caching via the least recently used policy
//! (offline), and prefetching in addition.
template <class SwappableBlockType>
class block_scheduler_algorithm_offline_lru_prefetching : public block_scheduler_algorithm<SwappableBlockType>
{
protected:
    struct scheduled_block_meta;
    struct write_read_request;

    typedef block_scheduler<SwappableBlockType> block_scheduler_type;
    typedef block_scheduler_algorithm<SwappableBlockType> block_scheduler_algorithm_type;
    typedef typename block_scheduler_type::internal_block_type internal_block_type;
    typedef typename block_scheduler_type::external_block_type external_block_type;
    typedef typename block_scheduler_type::swappable_block_identifier_type swappable_block_identifier_type;
    typedef typename block_scheduler_algorithm_type::time_type time_type;
    typedef typename block_scheduler_type::prediction_sequence_type prediction_sequence_type;
    typedef typename block_scheduler_type::block_scheduler_operation block_scheduler_operation;
    typedef typename std::vector<SwappableBlockType>::iterator swappable_blocks_iterator;

    typedef std::map<swappable_block_identifier_type, scheduled_block_meta> scheduled_blocks_type;
    typedef typename scheduled_blocks_type::iterator scheduled_blocks_iterator;
    typedef typename scheduled_blocks_type::reference scheduled_blocks_reference;
    typedef std::map<swappable_block_identifier_type, write_read_request*> write_scheduled_blocks_type;
    typedef typename write_scheduled_blocks_type::iterator write_scheduled_blocks_iterator;
    typedef typename write_scheduled_blocks_type::reference write_scheduled_blocks_reference;

    using block_scheduler_algorithm_type::bs;
    using block_scheduler_algorithm_type::swappable_blocks;
    using block_scheduler_algorithm_type::get_algorithm_from_block_scheduler;
    using block_scheduler_algorithm_type::get_free_internal_block_from_block_scheduler;
    using block_scheduler_algorithm_type::return_free_internal_block_to_block_scheduler;
    using block_scheduler_algorithm_type::prediction_sequence;

    struct write_read_request
    {
        bool write_done_soon;                          // set by read_after_write, checked by schedule_read()
        bool shall_read;                               // checked by read_after_write, set by schedule_read()
        swappable_blocks_iterator block_to_start_read; // used by read_after_write, set by schedule_read()
        scheduled_blocks_iterator taker;               // read_req set by read_after_write
        request_ptr write_req;                         // completes with read_after_write

        write_read_request()
            : write_done_soon(false), shall_read(false), block_to_start_read(), taker(), write_req(0) { }
    };

    struct read_after_write
    {
        write_read_request* wrr;

        read_after_write(write_read_request* write_read_req)
            : wrr(write_read_req) { }

        void operator () (request*)
        {
            wrr->write_done_soon = true;
            if (wrr->shall_read)
                wrr->taker->second.read_req = wrr->block_to_start_read->read_async();
        }
    };

    struct scheduled_block_meta
    {
        internal_block_type* reserved_iblock;
        std::pair<bool, swappable_block_identifier_type> giver;
        request_ptr read_req;
        std::deque<block_scheduler_operation> operations; // invariant: not empty; front: last scheduled operation, back: upcoming operation

        scheduled_block_meta(block_scheduler_operation op)
            : reserved_iblock(0),
              giver(false, 0),
              read_req(0),
              operations()
        { operations.push_front(op); }
    };

    //! Holds swappable blocks, whose internal block can be freed, i.e. that are internal but unacquired.
    std::set<swappable_block_identifier_type> free_evictable_blocks;
    std::set<swappable_block_identifier_type> scheduled_evictable_blocks;

    //! Holds not internal swappable_blocks, whose next access has already been scheduled.
    scheduled_blocks_type scheduled_blocks;

    //! Holds swappable_blocks, whose internal block has been taken away but the clean did not finish yet.
    write_scheduled_blocks_type write_scheduled_blocks;
    typename prediction_sequence_type::iterator next_op_to_schedule;

    //! Schedule an internal, possibly dirty swappable_block to write.
    //!
    //! The block becomes not dirty. if it was dirty, an entry in write_scheduled_blocks is made referencing the write_read_request.
    //! \param sbid block to write
    //! \return pointer to the write_read_request
    write_read_request * schedule_write(const swappable_block_identifier_type sbid)
    {
        SwappableBlockType& sblock = swappable_blocks[sbid];
        write_read_request* wrr = new write_read_request;
        wrr->write_req = sblock.clean_async(read_after_write(wrr));
        if (wrr->write_req.valid())
        {
            bool t = write_scheduled_blocks.insert(std::make_pair(sbid, wrr)).second;
            STXXL_ASSERT(t);
            return wrr;
        }
        else
        {
            delete wrr;
            return 0;
        }
    }

    //! Try to interrupt a read scheduled in a write_read_request.
    //!
    //! side-effect: possibly erases entry from write_scheduled_blocks, so the iterator writing_block may become invalid
    //! \return if successful
    bool try_interrupt_read(const write_scheduled_blocks_iterator& writing_block)
    {
        // stop read
        writing_block->second->shall_read = false;
        // check if stopped
        if (! writing_block->second->write_done_soon)
            return true;
        // => possibly to late
        scheduled_blocks_reference taker = *writing_block->second->taker;
        // wait
        wait_on_write(writing_block);
        // check if read started
        if (taker.second.read_req.valid())
            // => read started, to late
            return false;
        else
            // => just in time
            return true;
    }

    //! Schedule an internal and external block to read.
    //!
    //! If the giver is still writing, schedule read via its write_read_request.
    void schedule_read(scheduled_blocks_iterator block_to_read)
    {
        // first check if block_to_read is still writing. do not read before write finished
        // wait_on_write(block_to_read->first);

        write_scheduled_blocks_iterator it = write_scheduled_blocks.find(block_to_read->first);
        if (it != write_scheduled_blocks.end())
        {
            scheduled_blocks_iterator other_block_to_read = it->second->taker;
            // check if scheduled to read
            if (it->second->shall_read)
            {
                if (try_interrupt_read(it))
                {
                    // => interrupted, swap internal_blocks
                    std::swap(other_block_to_read->second.giver, block_to_read->second.giver);
                    if (other_block_to_read->second.giver.first)
                    {
                        write_scheduled_blocks_iterator it = write_scheduled_blocks.find(other_block_to_read->second.giver.second);
                        if (it != write_scheduled_blocks.end())
                            it->second->taker = other_block_to_read;
                        else
                            other_block_to_read->second.giver.first = false;
                    }
                    if (block_to_read->second.giver.first)
                    {
                        write_scheduled_blocks_iterator it = write_scheduled_blocks.find(block_to_read->second.giver.second);
                        if (it != write_scheduled_blocks.end())
                            it->second->taker = block_to_read;
                        else
                            block_to_read->second.giver.first = false;
                    }

                    internal_block_type* tmp_iblock = swappable_blocks[block_to_read->first].detach_internal_block();
                    swappable_blocks[block_to_read->first].attach_internal_block(
                        swappable_blocks[other_block_to_read->first].detach_internal_block());
                    swappable_blocks[other_block_to_read->first].attach_internal_block(tmp_iblock);
                    // => this block has its internal_block back, no need to read
                    // reschedule other
                    schedule_read(other_block_to_read);
                    return;
                }
                // else => read already started, but write done -> read this
            }
            else
            {
                // => no read scheduled, swap internal_blocks
                std::swap(other_block_to_read->second.giver, block_to_read->second.giver);
                if (other_block_to_read->second.giver.first)
                {
                    write_scheduled_blocks_iterator it = write_scheduled_blocks.find(other_block_to_read->second.giver.second);
                    if (it != write_scheduled_blocks.end())
                        it->second->taker = other_block_to_read;
                    else
                        other_block_to_read->second.giver.first = false;
                }
                if (block_to_read->second.giver.first)
                {
                    write_scheduled_blocks_iterator it = write_scheduled_blocks.find(block_to_read->second.giver.second);
                    if (it != write_scheduled_blocks.end())
                        it->second->taker = block_to_read;
                    else
                        block_to_read->second.giver.first = false;
                }

                internal_block_type* tmp_iblock = swappable_blocks[block_to_read->first].detach_internal_block();
                swappable_blocks[block_to_read->first].attach_internal_block(
                    other_block_to_read->second.reserved_iblock);
                other_block_to_read->second.reserved_iblock = tmp_iblock;
                // => this block has its internal_block back, no need to read
                return;
            }
        }

        // schedule block_to_read to read
        if (block_to_read->second.giver.first)
        {
            write_scheduled_blocks_iterator writing_block = write_scheduled_blocks.find(block_to_read->second.giver.second);
            if (writing_block != write_scheduled_blocks.end())
            {
                // => there is a write scheduled
                // tell the completion handler that we want a read
                writing_block->second->block_to_start_read = swappable_blocks.begin() + block_to_read->first;
                writing_block->second->taker = block_to_read;
                writing_block->second->shall_read = true;
                // and check if it is not to late
                if (writing_block->second->write_done_soon)
                {
                    // => the completion handler may have missed our wish to read
                    // so wait for it to finish and check
                    wait_on_write(writing_block);
                    block_to_read->second.giver.first = false;
                    if (block_to_read->second.read_req.valid())
                        // read scheduled
                        return;
                }
                else
                    // read scheduled
                    return;
            }
            else
                block_to_read->second.giver.first = false;
        }
        // => read could not be scheduled through the completion handler
        block_to_read->second.read_req = swappable_blocks[block_to_read->first].read_async();
    }

    //! Wait for the write to finish.
    //!
    //! side-effect: erases entry from write_scheduled_blocks
    void wait_on_write(const write_scheduled_blocks_iterator& writing_block)
    {
        writing_block->second->write_req->wait();
        delete writing_block->second;
        write_scheduled_blocks.erase(writing_block);
    }

    //! Wait for the write to finish.
    //!
    //! side-effect: erases entry from write_scheduled_blocks
    void wait_on_write(const swappable_block_identifier_type& writing_block)
    {
        write_scheduled_blocks_iterator it = write_scheduled_blocks.find(writing_block);
        if (it != write_scheduled_blocks.end())
            wait_on_write(it);
    }

    //! Wait for the write of the giver to finish.
    //!
    //! side-effect: erases entry from write_scheduled_blocks
    void wait_on_write(const scheduled_blocks_iterator& schedule_meta)
    {
        if (schedule_meta->second.giver.first)
        {
            wait_on_write(schedule_meta->second.giver.second);
            schedule_meta->second.giver.first = false;
        }
    }

    //! Wait for the read to finish.
    //!
    //! side-effect: erases entry for the write of the giver from write_scheduled_blocks
    void wait_on_read(const scheduled_blocks_iterator& schedule_meta)
    {
        wait_on_write(schedule_meta);
        if (schedule_meta->second.read_req.valid())
        {
            schedule_meta->second.read_req->wait();
            schedule_meta->second.read_req = 0;
        }
    }

    //! Wait for the write of the giver to finish and return reserved internal_block.
    //!
    //! side-effect: erases entry for the write of the giver from write_scheduled_blocks
    internal_block_type * get_ready_block(const scheduled_blocks_iterator& schedule_meta)
    {
        wait_on_write(schedule_meta);
        internal_block_type* r = schedule_meta->second.reserved_iblock;
        schedule_meta->second.reserved_iblock = 0;
        return r;
    }

    bool shall_keep_internal_block(const scheduled_blocks_iterator& schedule_meta, const bool ignore_first = true) const
    {
        // returns true iif there is an acquire or acquire_uninitialized scheduled or there is a deinitialize scheduled and the block is dirty
        for (typename std::deque<block_scheduler_operation>::reverse_iterator
             rit = schedule_meta->second.operations.rbegin() + ignore_first;
             rit != schedule_meta->second.operations.rend(); ++rit)
        {
            switch (*rit)
            {
            case block_scheduler_type::op_acquire:
            case block_scheduler_type::op_acquire_uninitialized:
                return true;
                break;
            case block_scheduler_type::op_release:
            case block_scheduler_type::op_release_dirty:
                break;
            case block_scheduler_type::op_deinitialize:
                if (swappable_blocks[schedule_meta->first].is_dirty()) return true; break;
            case block_scheduler_type::op_initialize:
            case block_scheduler_type::op_extract_external_block:
                break;
            }
        }
        return false;
    }

    // assumes the current operation to be still in operations
    bool shall_be_cleaned(const scheduled_blocks_iterator& schedule_meta) const
    {
        // returns true iif there is an extract_external_block scheduled and no release_dirty, deinitialize or initialize before
        for (typename std::deque<block_scheduler_operation>::reverse_iterator
             rit = schedule_meta->second.operations.rbegin() + 1;
             rit != schedule_meta->second.operations.rend(); ++rit)
        {
            switch (*rit)
            {
            case block_scheduler_type::op_acquire:
            case block_scheduler_type::op_acquire_uninitialized:
            case block_scheduler_type::op_release:
                break;
            case block_scheduler_type::op_release_dirty:
            case block_scheduler_type::op_deinitialize:
            case block_scheduler_type::op_initialize:
                return false;
                break;
            case block_scheduler_type::op_extract_external_block:
                return true;
                break;
            }
        }
        return false;
    }

    bool shall_be_read(const scheduled_blocks_iterator& schedule_meta, const bool ignore_first = true) const
    {
        // returns true iif there is an acquire scheduled next and the block is initialized
        return swappable_blocks[schedule_meta->first].is_initialized()
               && schedule_meta->second.operations.rbegin() + ignore_first != schedule_meta->second.operations.rend()
               && *(schedule_meta->second.operations.rbegin() + ignore_first) == block_scheduler_type::op_acquire;
    }

    void operation_done(scheduled_blocks_iterator& schedule_meta)
    {
        schedule_meta->second.operations.pop_back();
        if (schedule_meta->second.operations.empty())
        {
            assert(! schedule_meta->second.giver.first);
            scheduled_blocks.erase(schedule_meta);
        }
    }

    block_scheduler_algorithm_type * give_up(std::string err_msg = "detected some error in the prediction sequence")
    {
        STXXL_ERRMSG("block_scheduler_algorithm_offline_lru_prefetching: " << err_msg << ". Switching to block_scheduler_algorithm_online.");
        // switch algorithm
        block_scheduler_algorithm_type* new_algo
            = new block_scheduler_algorithm_online_lru<SwappableBlockType>(bs);
        // and delete self
        delete bs.switch_algorithm_to(new_algo);
        return new_algo;
    }

    void return_free_internal_block(internal_block_type* iblock)
    { return_free_internal_block_to_block_scheduler(iblock); }

    void schedule_next_operations()
    {
        while (next_op_to_schedule != prediction_sequence.end())
        {
            // list operation in scheduled_blocks
            std::pair<scheduled_blocks_iterator, bool> ins_res = scheduled_blocks.insert(
                std::make_pair(next_op_to_schedule->id, next_op_to_schedule->op));
            scheduled_blocks_iterator schedule_meta = ins_res.first;
            if (! ins_res.second)
                schedule_meta->second.operations.push_front(next_op_to_schedule->op);
            SwappableBlockType& sblock = swappable_blocks[next_op_to_schedule->id];

            // do appropriate preparations
            if (next_op_to_schedule->op == block_scheduler_type::op_acquire
                || next_op_to_schedule->op == block_scheduler_type::op_acquire_uninitialized)
            {
                if (sblock.is_internal())
                {
                    if (free_evictable_blocks.erase(next_op_to_schedule->id))
                        scheduled_evictable_blocks.insert(next_op_to_schedule->id);
                }
                else
                {
                    if (! schedule_meta->second.reserved_iblock)
                    {
                        // => needs internal_block -> try to get one
                        // -> try to get one from block_scheduler
                        schedule_meta->second.reserved_iblock = get_free_internal_block_from_block_scheduler();
                        if (! schedule_meta->second.reserved_iblock)
                        {
                            // -> try to get one by evicting
                            if (free_evictable_blocks.empty())
                            {
                                // => can not schedule acquire
                                // remove operation from scheduled_blocks
                                if (ins_res.second)
                                    scheduled_blocks.erase(ins_res.first);
                                else
                                    schedule_meta->second.operations.pop_front();
                                // stop scheduling
                                return;
                            }
                            swappable_block_identifier_type giver = pop_begin(free_evictable_blocks);
                            {
                                assert(scheduled_blocks.find(giver) == scheduled_blocks.end() ||
                                       !shall_keep_internal_block(scheduled_blocks.find(giver), false));
                            }
                            write_read_request* wrr = schedule_write(giver);
                            schedule_meta->second.giver.first = (wrr != NULL);
                            schedule_meta->second.giver.second = giver;
                            schedule_meta->second.reserved_iblock = swappable_blocks[giver].detach_internal_block();
                            if (wrr)
                                wrr->taker = schedule_meta;
                        }
                        // read if desired
                        if (shall_be_read(schedule_meta, false))
                        {
                            // => there is no operation scheduled for this block before this acquire and it is initialized
                            // -> start prefetching now
                            sblock.attach_internal_block(schedule_meta->second.reserved_iblock);
                            schedule_meta->second.reserved_iblock = 0;
                            scheduled_evictable_blocks.insert(next_op_to_schedule->id);
                            schedule_read(schedule_meta);
                        }
                    }
                }
            }
            else if (next_op_to_schedule->op == block_scheduler_type::op_deinitialize)
            {
                if (sblock.is_dirty())
                    if (free_evictable_blocks.erase(next_op_to_schedule->id))
                        scheduled_evictable_blocks.insert(next_op_to_schedule->id);
            }

            ++next_op_to_schedule;
        }
        for (typename std::set<swappable_block_identifier_type>::iterator it = free_evictable_blocks.begin();
             it != free_evictable_blocks.end(); ++it)
        {
            if (! write_scheduled_blocks.count(*it))
                schedule_write(*it);
        }
    }

    void init(block_scheduler_algorithm_type* old_algo)
    {
        if (old_algo)
            // copy prediction sequence
            prediction_sequence = old_algo->get_prediction_sequence();
        next_op_to_schedule = prediction_sequence.begin();
        if (get_algorithm_from_block_scheduler())
            while (! get_algorithm_from_block_scheduler()->evictable_blocks_empty())
                free_evictable_blocks.insert(get_algorithm_from_block_scheduler()->evictable_blocks_pop());
        schedule_next_operations();
    }

    void deinit()
    {
        // TODO remove
        if (! scheduled_blocks.empty())
            STXXL_MSG("deinit while scheduled_blocks not empty");
        if (! scheduled_evictable_blocks.empty())
            STXXL_MSG("deinit while scheduled_evictable_blocks not empty");

        // empty scheduled_blocks
        free_evictable_blocks.insert(scheduled_evictable_blocks.begin(), scheduled_evictable_blocks.end());
        //for (typename std::set<swappable_block_identifier_type>::iterator it = scheduled_evictable_blocks.begin();
        //        it != scheduled_evictable_blocks.end(); ++it)
        //    free_evictable_blocks.insert(*it);
        scheduled_evictable_blocks.clear();
        while (! scheduled_blocks.empty())
        {
            scheduled_blocks_iterator it = scheduled_blocks.begin();
            wait_on_read(it);
            if (it->second.reserved_iblock)
                return_free_internal_block(it->second.reserved_iblock);
            scheduled_blocks.erase(it);
        }
        while (! write_scheduled_blocks.empty())
        {
            write_scheduled_blocks_iterator it = write_scheduled_blocks.begin();
            wait_on_write(it);
        }
    }

public:
    block_scheduler_algorithm_offline_lru_prefetching(block_scheduler_type& bs)
        : block_scheduler_algorithm_type(bs)
    { init(get_algorithm_from_block_scheduler()); }

    // It is possible to keep an old simulation-algorithm object and reuse it's prediction sequence
    block_scheduler_algorithm_offline_lru_prefetching(block_scheduler_algorithm_type* old)
        : block_scheduler_algorithm_type(old)
    { init(old); }

    virtual ~block_scheduler_algorithm_offline_lru_prefetching()
    {
        deinit();
        if (! free_evictable_blocks.empty())
            STXXL_ERRMSG("Destructing block_scheduler_algorithm_offline_lru_prefetching that still holds evictable blocks. They get deinitialized.");
        while (! free_evictable_blocks.empty())
        {
            SwappableBlockType& sblock = swappable_blocks[pop_begin(free_evictable_blocks)];
            if (internal_block_type* iblock = sblock.deinitialize())
                return_free_internal_block(iblock);
        }
    }

    virtual bool evictable_blocks_empty()
    {
        deinit();
        return free_evictable_blocks.empty();
    }

    virtual swappable_block_identifier_type evictable_blocks_pop()
    { return pop_begin(free_evictable_blocks); }

    virtual internal_block_type & acquire(const swappable_block_identifier_type sbid, const bool uninitialized = false)
    {
        assert(! prediction_sequence.empty());
        assert(prediction_sequence.front().op ==
               ((uninitialized) ? block_scheduler_type::op_acquire_uninitialized : block_scheduler_type::op_acquire));
        assert(prediction_sequence.front().id == sbid);
        prediction_sequence.pop_front();
        scheduled_blocks_iterator schedule_meta = scheduled_blocks.find(sbid);
        assert(schedule_meta != scheduled_blocks.end());                                                               // acquire not scheduled or out of internal_blocks (i.e. not enough internal memory)
        assert(schedule_meta->second.operations.back() ==
               ((uninitialized) ? block_scheduler_type::op_acquire_uninitialized : block_scheduler_type::op_acquire)); // acquire not scheduled or out of internal_blocks (i.e. not enough internal memory)

        SwappableBlockType& sblock = swappable_blocks[sbid];
        /* acquired => internal -> increase reference count
           internal but not acquired -> remove from scheduled_evictable_blocks, increase reference count
           not internal => uninitialized or external -> get internal_block, increase reference count
           uninitialized -> fill with default value
           external -> read */
        if (sblock.is_internal())
        {
            if (! sblock.is_acquired())
            {
                // not acquired yet -> remove from scheduled_evictable_blocks
                size_t t = scheduled_evictable_blocks.erase(sbid);
                STXXL_ASSERT(t != 0);
                wait_on_read(schedule_meta);
            }
            sblock.acquire();
        }
        else
        {
            assert(uninitialized || ! sblock.is_initialized()); // initialized blocks should be scheduled to read and thus internal
            //get internal_block
            sblock.attach_internal_block(get_ready_block(schedule_meta));
            sblock.acquire();
            //initialize new block
            if (! uninitialized)
                sblock.fill_default();
        }

        operation_done(schedule_meta);
        return sblock.get_internal_block();
    }

    virtual void release(swappable_block_identifier_type sbid, const bool dirty)
    {
        assert(! prediction_sequence.empty());
        assert(prediction_sequence.front().op ==
               ((dirty) ? block_scheduler_type::op_release_dirty : block_scheduler_type::op_release));
        assert(prediction_sequence.front().id == sbid);
        prediction_sequence.pop_front();
        scheduled_blocks_iterator schedule_meta = scheduled_blocks.find(sbid);
        assert(schedule_meta != scheduled_blocks.end());
        assert(schedule_meta->second.operations.back() ==
               ((dirty) ? block_scheduler_type::op_release_dirty : block_scheduler_type::op_release));

        SwappableBlockType& sblock = swappable_blocks[sbid];
        sblock.make_dirty_if(dirty);
        sblock.release();
        if (! sblock.is_acquired())
        {
            if (sblock.is_dirty() || sblock.is_external())
            {
                // => evictable
                if (shall_keep_internal_block(schedule_meta))
                {
                    // => swappable_block shall keep its internal_block
                    scheduled_evictable_blocks.insert(sbid);
                    if (shall_be_cleaned(schedule_meta))
                        schedule_write(sbid);
                }
                else
                {
                    // give block to scheduler
                    free_evictable_blocks.insert(sbid);
                    if (next_op_to_schedule != prediction_sequence.end())
                        schedule_next_operations();
                    else {
                        if (! write_scheduled_blocks.count(sbid))
                            schedule_write(sbid);
                    }
                }
            }
            else
            {
                // => uninitialized
                if (shall_keep_internal_block(schedule_meta))
                    // => swappable_block shall keep its internal_block
                    schedule_meta->second.reserved_iblock = sblock.detach_internal_block();
                else
                {
                    // release internal block and give it to prefetcher
                    return_free_internal_block(sblock.detach_internal_block());
                    if (next_op_to_schedule != prediction_sequence.end())
                        schedule_next_operations();
                }
            }
        }
        operation_done(schedule_meta);
    }

    virtual void deinitialize(swappable_block_identifier_type sbid)
    {
        assert(! prediction_sequence.empty());
        assert(prediction_sequence.front().op == block_scheduler_type::op_deinitialize);
        assert(prediction_sequence.front().id == sbid);
        prediction_sequence.pop_front();
        scheduled_blocks_iterator schedule_meta = scheduled_blocks.find(sbid);
        assert(schedule_meta != scheduled_blocks.end());
        assert(schedule_meta->second.operations.back() == block_scheduler_type::op_deinitialize);

        SwappableBlockType& sblock = swappable_blocks[sbid];
        if (sblock.is_evictable())
        {
            size_t t;
            if (shall_keep_internal_block(schedule_meta, false))
            {
                t = scheduled_evictable_blocks.erase(sbid);
                if (t == 0) {
                    STXXL_ERRMSG("dirty block not scheduled on deinitialize");
                    t = free_evictable_blocks.erase(sbid);
                }
            }
            else
                t = free_evictable_blocks.erase(sbid);
            assert(t != 0);
        }
        if (internal_block_type* iblock = sblock.deinitialize())
        {
            if (shall_keep_internal_block(schedule_meta))
                // => swappable_block shall keep its internal_block
                schedule_meta->second.reserved_iblock = iblock;
            else
            {
                // release internal block and give it to prefetcher
                return_free_internal_block(iblock);
                if (next_op_to_schedule != prediction_sequence.end())
                    schedule_next_operations();
            }
        }
        operation_done(schedule_meta);
    }

    virtual void initialize(swappable_block_identifier_type sbid, external_block_type eblock)
    {
        assert(! prediction_sequence.empty());
        assert(prediction_sequence.front().op == block_scheduler_type::op_initialize);
        assert(prediction_sequence.front().id == sbid);
        prediction_sequence.pop_front();
        scheduled_blocks_iterator schedule_meta = scheduled_blocks.find(sbid);
        assert(schedule_meta != scheduled_blocks.end());
        assert(schedule_meta->second.operations.back() == block_scheduler_type::op_initialize);

        SwappableBlockType& sblock = swappable_blocks[sbid];
        sblock.initialize(eblock);
        if (shall_be_read(schedule_meta))
        {
            sblock.attach_internal_block(schedule_meta->second.reserved_iblock);
            schedule_meta->second.reserved_iblock = 0;
            scheduled_evictable_blocks.insert(sbid);
            schedule_read(schedule_meta);
        }
        operation_done(schedule_meta);
    }

    virtual external_block_type extract_external_block(swappable_block_identifier_type sbid)
    {
        assert(! prediction_sequence.empty());
        assert(prediction_sequence.front().op == block_scheduler_type::op_extract_external_block);
        assert(prediction_sequence.front().id == sbid);
        prediction_sequence.pop_front();
        scheduled_blocks_iterator schedule_meta = scheduled_blocks.find(sbid);
        assert(schedule_meta != scheduled_blocks.end());
        assert(schedule_meta->second.operations.back() == block_scheduler_type::op_extract_external_block);

        SwappableBlockType& sblock = swappable_blocks[sbid];
        wait_on_write(sbid);
        operation_done(schedule_meta);
        return sblock.extract_external_block();
    }
};

STXXL_END_NAMESPACE

#endif // !STXXL_MNG_BLOCK_SCHEDULER_HEADER
