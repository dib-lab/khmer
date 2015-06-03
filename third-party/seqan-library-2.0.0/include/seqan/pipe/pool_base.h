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

#ifndef SEQAN_HEADER_POOL_BASE_H
#define SEQAN_HEADER_POOL_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

/*!
 * @class PoolConfigSize
 * @extends PoolSpec
 * @headerfile <seqan/pipe.h.
 * @brief Configuration of Pool.
 *
 * @signature template <typename TSize[, typename TFile]>
 *            struct PoolConfigSize;
 *
 * @tparam TSize The Pool's size type.
 * @tparam TFile The underlying File type, defaults to <tt>File&lt;&gt;</tt>.
 *
 * @see PoolConfig
 */

    template < typename TSize,
               typename TFile = File<> >                        // default file type
    struct PoolConfigSize {
        typedef TSize SizeType;
        typedef TFile File;
    };

/*!
 * @class PoolConfig
 * @extends PoolSpec
 * @headerfile <seqan/pipe.h>
 * @brief Configuration of Pool.
 *
 * @signature template <typename TFile>
 *            struct PoolConfig;
 *
 * @tparam TFile The underlying File type, defaults to <tt>File&lt;&gt;</tt>.
 *
 * Using this configuration spec., the Pool's size type is <tt>Size&lt;TFile&gt;::Type</tt>.  To use a custom size type
 * PoolConfigSize should be used.
 *
 * @see PoolConfigSize
 */

    template < typename TFile = File<> >                        // default file type
    struct PoolConfig {
        typedef typename Size<TFile>::Type SizeType;
        typedef TFile File;
    };

/*!
 * @class PoolSpec
 * @extends Pool
 * @headerfile <seqan/pipe.h>
 * @brief Stores/Retrieves all elements to/from disk.
 *
 * @signature template <typename TValue, typename TConfig>
 *            class Pool<TValue, PoolSpec<TConfig> >;
 *
 * @tparam TConfig Configuration Spec.  Defines destination function, size type, and file type.
 *                 Types: PoolConfig, PoolConfigSize
 * @tparam TValue  The value type, that is the type of the stream elements.
 *
 * The Pool's input/output type is <tt>TValue</tt> and the size type is determined by the <tt>TConfig</tt>.
 */

    template < typename TConfig = PoolConfig<> >
    struct PoolSpec {
        typedef TConfig Config;
    };

/*!
 * @class Pool
 * @headerfile <seqan/pipe.h>
 * @brief Pools are push- and pop-passive pipeline modules.
 *
 * @signature template <typename TValue[, typename TSpec]>
 *            class Pool;
 *
 * @tparam TSpec  The specializing type. Default: PoolSpec&lt;&gt;, see PoolSpec.
 * @tparam TValue The value type, that is the type of the stream elements.
 *
 * Use Value to get the output type of a given Pipe (returns <tt>Value&lt;TInput&gt;::Type</tt> by default).
 *
 * Use Size to get the size type of a given Pipe (returns <tt>Size&lt;TInput&gt;::Type</tt> by default).
 */

    template < typename TValue,
               typename TSpec = PoolSpec<> >
    struct Pool;


    struct PoolParameters
    {

#ifdef SEQAN_IS_32_BIT
        // in 32bit mode at most 4GB are addressable
        enum { DefaultMemBufferSize     = 384 * 1024,      // low memory config [kB]
               DefaultPageSize          = 32 * 1024,            // [kB]
               DefaultBucketBufferSize  = 64 * 1024,            // [kB]
               DefaultReadAheadBuffers  = 4,
               DefaultWriteBackBuffers  = 4,
               DefaultWriteBackBuckets  = 16,
               DefaultAbsoluteSizes     = true };
#else
        enum { DefaultMemBufferSize     = 8 * 1024*1024,        // max memory config [kB]
               DefaultPageSize          = 1 * 1024*1024,        // [kB]
               DefaultBucketBufferSize  = 2 * 1024*1024,        // [kB]
               DefaultReadAheadBuffers  = 4,
               DefaultWriteBackBuffers  = 4,
               DefaultWriteBackBuckets  = 16,
               DefaultAbsoluteSizes     = true };
#endif

/*
        enum { DefaultMemBufferSize     = 128 * 1024,           // normal memory config [kB]
               DefaultPageSize          = 32 * 1024,            // [kB]
               DefaultBucketBufferSize  = 64 * 1024,            // [kB]
               DefaultReadAheadBuffers  = 2,
               DefaultWriteBackBuffers  = 2,
               DefaultWriteBackBuckets  = 16,
               DefaultAbsoluteSizes     = true };
*/

/*
        enum { DefaultMemBufferSize     = 0*8192,//64 * 1024,    // low memory config [kB]
               DefaultPageSize          = 2 * 1024,             // [kB]
               DefaultBucketBufferSize  = 6 * 1024,             // [kB]
               DefaultReadAheadBuffers  = 4,
               DefaultWriteBackBuffers  = 4,
               DefaultWriteBackBuckets  = 16,
               DefaultAbsoluteSizes     = true };
*/

        size_t  memBufferSize;
        size_t  pageSize;
        size_t  bucketBufferSize;
        size_t  readAheadBuffers;
        size_t  writeBackBuffers;
        size_t  writeBackBuckets;
        bool    absoluteSizes;      // when false, sizes are measured in units of TValue
                                    // when true, sizes are measured in bytes

        PoolParameters():
            memBufferSize((size_t)DefaultMemBufferSize * 1024ul),
            pageSize((size_t)DefaultPageSize * 1024ul),
            bucketBufferSize((size_t)DefaultBucketBufferSize * 1024ul),
            readAheadBuffers(DefaultReadAheadBuffers),
            writeBackBuffers(DefaultWriteBackBuffers),
            writeBackBuckets(DefaultWriteBackBuckets),
            absoluteSizes(DefaultAbsoluteSizes) {}

        template < typename TValue >
        void absolutize(size_t aligning, TValue *)
        {
            if (!absoluteSizes) return;
            memBufferSize = (memBufferSize + sizeof(TValue) - 1) / sizeof(TValue);
            bucketBufferSize = (bucketBufferSize + sizeof(TValue) - 1) / sizeof(TValue);
            pageSize = (pageSize + sizeof(TValue) - 1) / sizeof(TValue);
            pageSize = ((pageSize + aligning - 1) / aligning) * aligning;
        }
    };



    template < typename TBuffer, typename THandler >
    inline TBuffer& processBuffer(TBuffer &h, THandler &) {
        return h;
    }


    //////////////////////////////////////////////////////////////////////////////
    // handler that manages a simple memory buffer
    struct MemorySpec_;
    typedef Tag<MemorySpec_> MemorySpec;

    template <typename TValue, typename TPoolSpec>
    struct BufferHandler<Pool<TValue, TPoolSpec>, MemorySpec>
    {
        typedef Pool<TValue, TPoolSpec> TPool;
        typedef Buffer<TValue>          TBuffer;

        TPool    &pool;
        TBuffer    empty;

        BufferHandler(TPool &_pool):
            pool(_pool) {}

        BufferHandler(TPool &_pool, size_t):
            pool(_pool) {}

        inline TBuffer & first()
        {
            return pool.memBuffer;
        }

        inline TBuffer & next()
        {
            return empty;
        }

        inline void process()
        {
            processBuffer(pool.memBuffer, *this);
        }

        inline void end() {}
        inline void cancel() {}
    };


    //////////////////////////////////////////////////////////////////////////////
    // generic block based asynchronous read handler
    struct ReadFileSpec_;
//IOREV _notio_ not directly relevant to iorev
    typedef Tag<ReadFileSpec_> ReadFileSpec; //IOREV _notio_ not directly relevant to iorev

    template <typename TValue, typename TSpec>
    struct BufferHandler<Pool<TValue, TSpec>, ReadFileSpec >
    {
        typedef Pool<TValue, TSpec>                         TPool;
        typedef typename TPool::File                        TFile;

        typedef Buffer<TValue>                              TBuffer;
        typedef Buffer<TValue, PageFrame<TFile, Dynamic> >  TPageFrame;
        typedef PageChain<TPageFrame>                       TPageChain;

        TPool       &pool;
        TPageChain  chain;
        size_t      pageSize;
        size_t      readPageNo, _pages;
        TBuffer        empty;

        BufferHandler(TPool &_pool):
            pool(_pool),
            chain(_min(_pool.readAheadBuffers, _pool.pages())),
            pageSize(_pool.pageSize) {}

        BufferHandler(TPool &_pool, size_t _requestedBufferSize, size_t _readAheadBuffers = 1):
            pool(_pool),
//            pageSize(alignSize(_min(_pool.size(), _requestedBufferSize), _pool.pageSize)),
            chain(_min(_readAheadBuffers, _pool.pages(pageSize = alignSize(_min(_pool.size(), _requestedBufferSize), _pool.pageSize))))
        {
            #ifdef SEQAN_HEADER_PIPE_DEBUG
                std::cerr << "___BufferHandler___" << std::endl;
                std::cerr << "pagesize: " << pageSize << std::endl;
                std::cerr << "readaheadbuffers: " << chain.maxFrames << std::endl;
                std::cerr << "pages: " << pool.pages(pageSize) << std::endl;
            #endif
        }

        ~BufferHandler()
        {
            end();
        }

        inline TBuffer & first()
        {
            _pages = pool.pages(pageSize);
            if (!_pages) return empty;

            // enqueue reading of the <readAheadBuffers> first blocks
            readPageNo = 0;
            TPageFrame *p = chain.first;
            while (p) {
                p->pageNo = readPageNo++;
                _read(*p);
                p = p->next;
            }

            // retrieve the very first and wait for I/O transfer to complete
            bool waitResult = waitFor(*chain.first);
            if (!waitResult)
                SEQAN_FAIL("%s operation could not be completed: \"%s\"",
                           _pageFrameStatusString(chain.first->status),
                           strerror(errno));

            return processBuffer(*chain.first, *this);
        }

        inline TBuffer & next()
        {
            // step one buffer ahead
            chain.getReadyPage();

            // read ahead
            chain.last->pageNo = readPageNo++;
            _read(*chain.last);

            // retrieve the next buffer in order and wait for I/O transfer to complete
            bool waitResult = waitFor(*chain.first);
            if (!waitResult)
                SEQAN_FAIL("%s operation could not be completed: \"%s\"",
                           _pageFrameStatusString(chain.first->status),
                           strerror(errno));

            return processBuffer(*chain.first, *this);
        }

        inline void end()
        {
            cancel();
        }

        inline void cancel()
        {
            TPageFrame *p = chain.first;
            while (p) {
                seqan::cancel(*p, pool.file);
                freePage(*p, pool.file);
                p = p->next;
            }
        }

        inline void process() {}

    private:
        bool _error() {
//            std::cerr << "Error in BufWriteFileHandler::_read " << pool.file.error() << std::endl;
            return true;
        }

        inline bool _read(TPageFrame &pf)
        {
            if (pf.pageNo < _pages) {
                // alloc if empty
                if (!pf.begin)
                    allocPage(pf, pageSize, pool.file);

                // set buffer size according to read size
//                std::cout << "poolsize="<<pool._size<<" pageno="<<pf.pageNo<<" pagesize="<<pageSize<<" resutl="<<pool.dataSize(pf.pageNo, pageSize)<<std::endl;
                resize(pf, pool.dataSize(pf.pageNo, pageSize));

                // read asynchronously (if possible) from disk
                return readPage(pf, pool.file) || _error();
            } else {
                // free if allocated
                freePage(pf, pool.file);
                return false;
            }
        }
    };


    //////////////////////////////////////////////////////////////////////////////
    // generic block based asynchronous write handler
    struct WriteFileSpec_;
//IOREV _notio_ not directly relevant to iorev
    typedef Tag<WriteFileSpec_> WriteFileSpec; //IOREV _notio_ not directly relevant to iorev

    template < typename TValue, typename TSpec >
    struct BufferHandler<Pool<TValue, TSpec>, WriteFileSpec>
    {
        typedef Pool<TValue, TSpec>                 TPool;
        typedef typename TPool::File                TFile;

        typedef Buffer<TValue>                              TBuffer;
        typedef Buffer<TValue, PageFrame<TFile, Dynamic> >  TPageFrame;
        typedef PageChain<TPageFrame>                       TPageChain;

        TPool       &pool;
        TPageChain  chain;
        size_t      pageSize;
        size_t      writePageNo, _pages;
        TBuffer        empty;

        BufferHandler(TPool &_pool):
            pool(_pool),
            chain(_min(_pool.writeBackBuffers, _pool.pages())),
            pageSize(_pool.pageSize),
            writePageNo(0),
            _pages(0) {}

        BufferHandler(TPool &_pool, size_t _requestedBufferSize, size_t _writeBackBuffers = 1):
            pool(_pool),
            chain(_min(_writeBackBuffers, _pool.pages(pageSize))),
            pageSize(alignSize(_min(_pool.size(), _requestedBufferSize), _pool.pageSize)) {}

        ~BufferHandler() {
            cancel();
        }

    public:

        inline TBuffer & first()
        {
            _pages = pool.pages(pageSize);
            if (!_pages) return empty;

            // get a ready page frame
            TPageFrame & pf = *chain.getReadyPage();

            if (!pf.begin)
                allocPage(pf, pageSize, pool.file);

            writePageNo = 0;
            resize(pf, pool.dataSize(pf.pageNo = writePageNo++, pageSize));
            return pf;
        }

        inline TBuffer & next()
        {
            // write previously provided buffer to disk
            _write(processBuffer(*chain.last, *this));

            // step one buffer ahead
            TPageFrame & pf = *chain.getReadyPage();

            if (!pf.begin)
                allocPage(pf, pageSize, pool.file);

            resize(pf, pool.dataSize(pf.pageNo = writePageNo++, pageSize));
            return pf;
        }

        inline void end()
        {
            // write previously provided buffer to disk
            _write(processBuffer(*chain.last, *this));
            flush();
        }

        inline void flush()
        {
            TPageFrame *p = chain.first;
            while (p)
            {
                // wait for I/O transfer to complete
                bool waitResult = waitFor(*p);
                if (!waitResult)
                    SEQAN_FAIL("%s operation could not be completed: \"%s\"",
                               _pageFrameStatusString(p->status),
                               strerror(errno));

                freePage(*p, pool.file);
                p = p->next;
            }
            seqan::flush(pool.file);
        }

        inline void cancel()
        {
            TPageFrame *p = chain.first;
            while (p) {
                seqan::cancel(*p, pool.file);
                freePage(*p, pool.file);
                p = p->next;
            }
        }

        inline void process() {}

    private:

        bool _error()
        {
//            std::cerr << "Error in BufWriteFileHandler::_write " << pool.file.error() << std::endl;
            return true;
        }

        bool _write(TPageFrame &pf)
        {
            if (pf.pageNo < _pages) {
                // write asynchronously (if possible) to disk
                return writePage(pf, pool.file) || _error();
            } else {
                // free if allocated
                freePage(pf, pool.file);
                return false;
            }
        }
    };

    //////////////////////////////////////////////////////////////////////////////
    // generic buffered multiplex handler
    struct MultiplexSpec_;
    typedef Tag<MultiplexSpec_> MultiplexSpec;

    template < typename TBufferHandler1, typename TBufferHandler2 >
    struct BufferHandler< Bundle2< TBufferHandler1, TBufferHandler2 >, MultiplexSpec >
    {
        typedef typename Value<TBufferHandler1>::Type TBuffer;

        TBufferHandler1 *handler1;
        TBufferHandler2 *handler2;

        template <typename TPool>
        BufferHandler(TPool &_pool) {
            if (_pool.memBuffer.begin || _pool._size == 0) {    // decision to choose handler1 or handler2
                handler1 = new TBufferHandler1(_pool);
                handler2 = NULL;
            } else {
                handler1 = NULL;
                handler2 = new TBufferHandler2(_pool);
            }
        }

        template <typename TPool>
        BufferHandler(TPool &_pool, size_t _requestedBufferSize) {
            if (_pool.memBuffer.begin || _pool._size == 0) {    // decision to choose handler1 or handler2
                handler1 = new TBufferHandler1(_pool, _requestedBufferSize);
                handler2 = NULL;
            } else {
                handler1 = NULL;
                handler2 = new TBufferHandler2(_pool, _requestedBufferSize);
            }
        }

        ~BufferHandler() {
            delete handler1;
            delete handler2;
        }

        inline TBuffer first() {
            if (handler1)    return handler1->first();
            else            return handler2->first();
        }

        inline TBuffer next() {
            if (handler1)    return handler1->next();
            else            return handler2->next();
        }

        inline void end() {
            if (handler1)    handler1->end();
            else            handler2->end();
        }

        inline void process() {
            if (handler1)    handler1->process();
            else            handler2->process();
        }

        inline void cancel() {
            if (handler1)    handler1->cancel();
            else            handler2->cancel();
        }
    };


    template < typename THandler1, typename THandler2 >
    struct Handler< Bundle2< THandler1, THandler2 >, MultiplexSpec >
    {
        typedef typename Value<Handler>::Type TValue;

        THandler1 *handler1;
        THandler2 *handler2;

        template <typename TPool>
        Handler(TPool &_pool)
        {
            if (_pool.memBuffer.begin || _pool._size == 0) {    // decision to choose handler1 or handler2
                handler1 = new THandler1(_pool);
                handler2 = NULL;
            } else {
                handler1 = NULL;
                handler2 = new THandler2(_pool);
            }
        }

        ~Handler() {
            delete handler1;
            delete handler2;
        }

        inline bool begin() {
            if (handler1)    return handler1->begin();
            else            return handler2->begin();
        }

        inline TValue const & front() const {
            if (handler1)    return handler1->front();
            else            return handler2->front();
        }

        inline void pop() {
            if (handler1)    handler1->pop();
            else            handler2->pop();
        }

        inline void pop(TValue &Ref_) {
            if (handler1)    handler1->pop(Ref_);
            else            handler2->pop(Ref_);
        }

        inline void push(TValue const & Val_) {
            if (handler1)    handler1->push(Val_);
            else            handler2->push(Val_);
        }

        inline bool eof() const {
            if (handler1)    return handler1->eof();
            else            return handler2->eof();
        }

        inline void end() {
            if (handler1)    handler1->end();
            else            handler2->end();
        }

        inline void process() {
            if (handler1)    handler1->process();
            else            handler2->process();
        }
    };



    //////////////////////////////////////////////////////////////////////////////
    // character and buffer based handler definitions
    template < typename TValue, typename TSpec >
    struct BufReadHandler< Pool< TValue, TSpec > >
    {
        typedef BufferHandler< Bundle2<
            BufferHandler< Pool< TValue, TSpec >, MemorySpec >,
            BufferHandler< Pool< TValue, TSpec >, ReadFileSpec >
        >, MultiplexSpec > Type;
    };

    template < typename TValue, typename TSpec >
    struct BufWriteHandler< Pool< TValue, TSpec > >
    {
        typedef BufferHandler< Bundle2<
            BufferHandler< Pool< TValue, TSpec >, MemorySpec >,
            BufferHandler< Pool< TValue, TSpec >, WriteFileSpec >
        >, MultiplexSpec > Type;
    };


    template <typename TBufferHandler1, typename TBufferHandler2>
    struct Value< BufferHandler< Bundle2<TBufferHandler1, TBufferHandler2>, MultiplexSpec> > :
        public Value<TBufferHandler1> {};

    template <typename THandler1, typename THandler2>
    struct Value< Handler< Bundle2<THandler1, THandler2>, MultiplexSpec> > :
        public Value<THandler1> {};

    template < typename TPool >
    struct HandlerArgs {
        typedef Nothing Type;
    };


    //////////////////////////////////////////////////////////////////////////////
    // base class of all pool classes like Pool, Mapper, Sorter

    template < typename TValue,
               typename TSpec >
    struct Pool
    {
        typedef typename TSpec::Config      Config;
        typedef typename Config::File       File;
        typedef typename Config::SizeType   SizeType;
        typedef Buffer<TValue>              TBuffer;

        // public handlers to read simultanously buffer- or character-wise
        typedef typename ReadHandler<Pool>::Type    ReadHandler;
        typedef typename WriteHandler<Pool>::Type   WriteHandler;
        typedef typename HandlerArgs<Pool>::Type    HandlerArgs;

        File                file;
        bool                _temporary, _ownFile;
        SizeType            _size;
        size_t              _pages;
        size_t              pageSize;
        size_t              bucketBufferSize;
        size_t              readAheadBuffers;
        size_t              writeBackBuffers;
        size_t              writeBackBuckets;

        TBuffer                memBuffer;
        size_t              memBufferSize;
        HandlerArgs         handlerArgs;

        bool                _partiallyFilled;        // the pool is partially filled (it contains undefined values)
        TValue                undefinedValue;            // value to represent undefined (unwritten) entries

    protected:
        size_t              _lastPageNo;
        size_t              _lastPageSize;

        int                 listeners;

        ReadHandler         *reader;
        WriteHandler        *writer;


    public:

        //////////////////////////////////////////////////////////////////////////////
        // public iterator types

        typedef SizeType                size_type;

        Pool(const PoolParameters &_conf = PoolParameters()) :
            file(NULL), undefinedValue()
        {
            _init(_conf);
            _setSize(0);
        }

        Pool(HandlerArgs const &args, PoolParameters const &_conf = PoolParameters()):
            file(NULL), handlerArgs(args), undefinedValue()
        {
            _init(_conf);
            _setSize(0);
        }

        template < typename TInput, typename TPipeSpec >
        Pool(Pipe<TInput, TPipeSpec> &, const PoolParameters &_conf = PoolParameters()):
            file(NULL), undefinedValue()
        {
            _init(_conf);
            _setSize(0);
        }

        template < typename TInput, typename TPipeSpec >
        Pool(Pipe<TInput, TPipeSpec> &, HandlerArgs const &args, PoolParameters const &_conf = PoolParameters()):
            file(NULL), handlerArgs(args), undefinedValue()
        {
            _init(_conf);
            _setSize(0);
        }

        Pool(File &_file, const PoolParameters &_conf = PoolParameters()):
            file(_file), undefinedValue()
        {
            _init(_conf);
            _ownFile = false;
            _temporary = false;
            memBufferSize = 0;
            _setSize(seqan::length(file) / sizeof(TValue));
        }

        Pool(const char *fileName, const PoolParameters &_conf = PoolParameters()) :
            file(NULL), undefinedValue()
        {
            _init(_conf);
            _temporary = false;
            memBufferSize = 0;
            _ownFile = open(file, fileName);
            if (_ownFile)
                _setSize(seqan::length(file) / sizeof(TValue));
            else
                _setSize(0);
        }

        ~Pool()
        {
            endRead();
            endWrite();
            if (_temporary)
                clear();
            else
                if (_ownFile) close(file);
        }

        inline void clear() {
            resize(0);
        }

        inline size_type size() const {
            return _size;
        }

        // this is not a real resize
        void resize(size_type _newSize) {
            typedef typename Size<File>::Type TFSize SEQAN_UNUSED_TYPEDEF;
            if (_newSize == _size) return;

            _freeHandlers();    // if you forgot to call endRead/endWrite we have no trouble

            if (_temporary && _ownFile) {
                if (_size != 0) {
                    if (memBuffer.begin)
                        freePage(memBuffer, *this);
                    else {
                        close(file);
                        SEQAN_PROSUB(SEQAN_PROIOVOLUME, (_proFloat)((TFSize)_size * (TFSize)sizeof(TValue)));
                    }
                }

                if (_newSize != 0) {
                    if (_newSize <= (size_type)memBufferSize)
                        allocPage(memBuffer, _newSize, *this);
                    else {
                        openTemp(file);
                        SEQAN_PROADD(SEQAN_PROIOVOLUME, (_proFloat)((TFSize)_newSize * (TFSize)sizeof(TValue)));
                    }
                }
            }

            _setSize(_newSize);
        }


        //////////////////////////////////////////////////////////////////////////////
        // auto-disposal interface (deprecated)

        inline void addListener() {
            if (!listeners) return;
            ++listeners;
        }

        inline void delListener() {
            if (!listeners) return;
            if (--listeners == 0) {
                clear();
            }
        }


        //////////////////////////////////////////////////////////////////////////////
        // queue interface

        inline TValue const & front() const
        {
            return reader->front();
        }

        inline TValue const & operator*() const
        {
            return reader->front();
        }

        inline void pop()
        {
            reader->pop();
        }

        inline Pool& operator++()
        {
            reader->pop();
            return *this;
        }

        inline void pop(TValue &Ref_)
        {
            reader->pop(Ref_);
        }

        inline void push(TValue const &Val_)
        {
            writer->push(Val_);
        }

        inline bool eof()
        {
            if (reader) return reader->eof();
            if (writer) return writer->eof();
            return true;
        }


        //////////////////////////////////////////////////////////////////////////////
        // flow control

        bool beginWrite() {
            _freeHandlers(); // if you forgot to call endRead/endWrite we have no trouble
            return (
                (writer = new WriteHandler(*this)) &&
                writer->begin());
        }

        bool endWrite() {
            if (writer) {
                writer->end();
                writer->process();
            }
            delete writer;
            writer = NULL;
            return true;
        }

        bool beginRead() {
            _freeHandlers(); // if you forgot to call endRead/endWrite we have no trouble
            return (
                (reader = new ReadHandler(*this)) &&
                reader->begin());
        }

        bool endRead() {
            if (reader) reader->end();
            delete reader;
            reader = NULL;
            delListener();
            return true;
        }


        inline size_t pages() const {
            return _pages;
        }

        inline size_t pages(size_t pageSize__) const {
            return enclosingBlocks(_size, pageSize__);
        }

        // used by buffer handlers ...
        inline size_t dataSize(size_t pageNo__) const {
            return (pageNo__ != _lastPageNo)? pageSize: _lastPageSize;
        }

        // used by buffer handlers with variable PageSize ...
        inline size_t dataSize(size_t pageNo__, size_t pageSize__) const {
            return (pageNo__ != _size / pageSize__)? pageSize__: _size % pageSize__;
        }

    protected:

        inline void _freeHandlers() {
            if (reader) reader->end();
            if (writer) writer->end();
            delete reader;
            delete writer;
            reader = NULL;
            writer = NULL;
        }

        inline void _setSize(size_type _newSize) {
            _size = _newSize;
            _pages = enclosingBlocks(_size, pageSize);
            _lastPageNo = _size / pageSize;
            _lastPageSize = _size % pageSize;
        }

        void _init(PoolParameters _conf = PoolParameters())
        {
            _conf.absolutize(16*1024/*sectorSize(file)*/, (TValue*)NULL);
            memBufferSize    = _conf.memBufferSize;
            pageSize         = _conf.pageSize;
            bucketBufferSize = _conf.bucketBufferSize;
            readAheadBuffers = _conf.readAheadBuffers;
            writeBackBuffers = _conf.writeBackBuffers;
            writeBackBuckets = _conf.writeBackBuffers;
            _partiallyFilled = true;
            listeners = 0;
            reader = NULL;
            writer = NULL;
            _ownFile = true;
            _temporary = true;
        }

    };

    template < typename TValue, typename TSpec, typename TIteratorSpec>
    struct Iterator< Pool< TValue, TSpec > const, TIteratorSpec> {
        typedef IPipeIterator< Pool< TValue, TSpec > > Type;
    };

    template < typename TValue, typename TSpec, typename TIteratorSpec>
    struct Iterator< Pool< TValue, TSpec >, TIteratorSpec> {
        typedef OPipeIterator< Pool< TValue, TSpec > > Type;
    };

    template < typename TValue, typename TSpec >
    OPipeIterator< Pool< TValue, TSpec > >
    begin(Pool< TValue, TSpec > &pool) {
        return OPipeIterator< Pool< TValue, TSpec > >(pool);
    }

    template < typename TValue, typename TSpec >
    OPipeIterator< Pool< TValue, TSpec > >
    end(Pool< TValue, TSpec > &/*pool*/) {
        return OPipeIterator< Pool< TValue, TSpec > >();
    }



    //////////////////////////////////////////////////////////////////////////////
    // a pool is a queue-like container for a large amount of data
    // in contrast to a queue you can read the whole content more than once
    // but can't pop directly after a push
    // instead, access to a pool looks like the following:
    //
    //  1. resize(<new size>)
    //  2. beginWrite()
    //  3. push(a), push(b), push(c) ...
    //  4. endWrite()
    // (5. do something else)
    //  6. beginRead()
    //  7. pop(a), pop(b), pop(c) ...
    //  8. endRead()
    // (9. clear() - if you want to save memory and refill the pool later)
    //
    // you can repeat steps 2-4 and 6-8 independently as often as you like

    //template < typename TValue, typename TFile = File<> >
    //struct Pool: public Pool< TValue, PoolSpec< PoolConfig< TFile > > >
    //{
    //    typedef TValue                        Type;
    //    typedef TFile                       File;
    //    typedef typename Size<TFile>::Type  SizeType;
    //};


    // seqan namespace traits
    template < typename TValue, typename TSpec >
    struct Value< Pool< TValue, TSpec > > {
        typedef TValue Type;
    };

    template < typename TValue, typename TSpec >
    struct Size< Pool< TValue, TSpec > > {
        typedef typename Pool< TValue, TSpec >::SizeType Type;
    };

    template < typename TValue, typename TSpec >
    struct Position< Pool< TValue, TSpec > > {
        typedef typename Size<Pool< TValue, TSpec > >::Type Type;
    };

    template < typename TValue, typename TSpec >
    struct Difference< Pool< TValue, TSpec > > {
        typedef typename MakeSigned_<typename Size<Pool< TValue, TSpec > >::Type>::Type Type;
    };

/*!
 * @fn Pool#clear
 * @brief Remove all elements from the pool.
 *
 * @signature void clear(pool);
 *
 * @param[in,out] pool Pool to clear.
 */
    template < typename TValue, typename TSpec >
    inline void clear(Pool<TValue, TSpec> &me)
    {
        return me.clear();
    }

    // deprecated
    template < typename TValue, typename TSpec >
    inline typename Size< Pool<TValue, TSpec> >::Type
    size(Pool<TValue, TSpec> const &me)
    {
        return me.size();
    }

/*!
 * @fn Pool#length
 * @headerfile <seqan/pipe.h>
 * @brief Length of the pool.
 *
 * @signature TSize length(pool);
 *
 * @param[in] pool  The Pool to query for its size.
 *
 * @return    TSize The number of elements in the pool.
 */

    template < typename TValue, typename TSpec >
    inline typename Size< Pool<TValue, TSpec> >::Type
    length(Pool<TValue, TSpec> const &me)
    {
        return me.size();
    }

/*!
 * @fn Pool#resize
 * @brief Resize a pool.
 *
 * @signature void resize(pool, len);
 *
 * @param[in,out] pool Pool to resize.
 * @param[in]     len  Length to resize <tt>pool</tt> to.
 */

    template < typename TValue, typename TSpec, typename TSize >
    inline TSize resize(Pool<TValue, TSpec> &me, TSize new_length)
    {
        me.resize(new_length);
        return me.size();
    }

/*!
 * @fn Pool#front
 * @brief Return reference to the first element.
 *
 * @signature TReference front(pool);
 *
 * @param[in] pool The pool to get the first element of.
 *
 * @return TReference A reference to the first element of <tt>seq</tt>.
 */

    template < typename TValue, typename TSpec >
    inline typename Value< Pool<TValue, TSpec> >::Type const & front(Pool<TValue, TSpec> &me) {
        return me.front();
    }
/*!
 * @fn Pool#pop
 * @headerfile <seqan/pipe.h>
 * @brief Pops the first element of the remaining stream.
 *
 * @signature void pop(pool[, ref]);
 *
 * @param[in,out] pool A pop-passive pipeline module.
 * @param[out]    ref    Reference to the result.  Result type is <tt>Value&lt;TObject&gt;::Type</tt> for <tt>object</tt>
 *                       type <tt>TObject</tt>.  Returns the first element of the remaining input stream.
 *
 * In contrast to Pool#front this function also steps one element further.
 *
 * Pool#front or Pool#pop can only be called within a read process surrounded by beginRead and endRead.
 */
    template < typename TValue, typename TSpec >
    inline void pop(Pool<TValue, TSpec> &me) {
        me.pop();
    }

    template < typename TValue, typename TSpec >
    inline void pop(Pool<TValue, TSpec> &me, TValue &Ref_) {
        me.pop(Ref_);
    }

/*!
 * @fn Pool#push
 * @brief Appends an item at the end of an input stream.
 *
 * @signature void push(pool, val);
 *
 * @param[in,out] pool A pool module.
 * @param[in]     val    Item to be pushed.
 *
 * The function <tt>push</tt> can only be called within a write process surrounded by beginWrite and endWrite.
 */

    template < typename TValue, typename TSpec >
    inline void push(Pool<TValue, TSpec> &me, TValue const &Val_) {
        me.push(Val_);
    }

    template < typename TValue, typename TSpec >
    std::ostream& operator<<(std::ostream &out, Pool<TValue, TSpec> &p) {
        beginRead(p);
        while (!eof(p)) {
            out << front(p) << std::endl;
            pop(p);
        }
        endRead(p);
        return out;
    }


    // the pipe interface of pool classes
    //namespace SEQAN_NAMESPACE_PIPELINING
    //{
        //template < typename TValue, typename TSpec >
        //struct Value< Pool< TValue, TSpec > >
        //{
        //    typedef TValue Type;
        //};

        //template < typename TValue, typename TSpec >
        //struct Size< Pool< TValue, TSpec > >
        //{
        //    typedef typename Size< TFile >::Type Type;
        //};

        template < typename TValue, typename TSpec >
        inline bool control(Pool< TValue, TSpec > &me, ControlEof const &) {
            return me.eof();
        }

        template < typename TValue, typename TSpec >
        inline bool control(Pool< TValue, TSpec > &me, ControlEos const &) {
            return me.eof();
        }

        template < typename TValue, typename TSpec >
        inline bool control(Pool< TValue, TSpec > &me, ControlClear const &) {
            me.clear();
            return true;
        }

        template < typename TValue, typename TSpec >
        inline bool control(Pool< TValue, TSpec > &me, ControlBeginRead const &) {
            return me.beginRead();
        }

        template < typename TValue, typename TSpec >
        inline bool control(Pool< TValue, TSpec > &me, ControlEndRead const &) {
            return me.endRead();
        }

/*!
 * @fn Pool#beginWrite
 * @headerfile <seqan/pipe.h>
 * @brief Initiates a write process.
 *
 * @signature bool beginWrite(pool);
 *
 * @param[in,out] pool A pool module.
 *
 * @return bool <tt>true</tt> on success, false on failure.
 *
 * <tt>beginWrite</tt> prepares a pool for succeeding writes.
 *
 * A write process must be terminated with endWrite.  Nested write processes are not allowed.
 *
 * @see Pool#endWrite
*/

        template < typename TValue, typename TSpec >
        inline bool beginWrite(Pool< TValue, TSpec > &me) {
SEQAN_CHECKPOINT
            return me.beginWrite();
        }

/*!
 * @fn Pool#endWrite
 * @brief Terminates a write process.
 *
 * @signature bool endWrite(pool);
 *
 * @param[in,out] pool A push-passive pipeline module.
 *
 * @return bool true on success, false on failure.
 *
 * <tt>endWrite</tt> closes the input stream and frees resources possibly allocated by beginWrite
 *
 * @see Pool#beginWrite
 */

        template < typename TValue, typename TSpec >
        inline bool endWrite(Pool< TValue, TSpec > &me) {
            return me.endWrite();
        }

/*!
 * @fn Pool#atEnd
 * @brief Check whether the @link Pool @endlink object is at end.
 *
 * @signature bool atEnd(pool);
 *
 * @param[in] pool The @link Pool @endlink object to query.
 *
 * @return bool true in case of the pool being at the end, false otherwise.
 */

        template < typename TValue, typename TSpec >
        inline bool eof(Pool< TValue, TSpec > &me) {
SEQAN_CHECKPOINT
            return control(me, ControlEof());
        }

/*!
 * @fn Pool#beginRead
 * @headerfile <seqan/pipe.h>
 * @brief Initiates a read process.
 *
 * @signature bool beginRead(pool);
 *
 * @param[in,out] pool A pool module.
 *
 * @return bool true on success, false on failure.
 *
 * A read process must be terminated with endRead. Nested read processes are not allowed.
 *
 * @see Pool#endRead
 */
        template < typename TValue, typename TSpec >
        inline bool beginRead(Pool< TValue, TSpec > &me) {
SEQAN_CHECKPOINT
            return control(me, ControlBeginRead());
        }
/*!
 * @fn Pool#endRead
 * @headerfile <seqan/pipe.h>
 * @brief Terminates a read process.
 *
 * @signature bool endRead(pool);
 *
 * @param[in,out] pool A pool module.
 *
 * @return bool true on success, false on failure.
 *
 * @see Pool#beginRead
 */
        template < typename TValue, typename TSpec >
        inline bool endRead(Pool< TValue, TSpec > &me) {
SEQAN_CHECKPOINT
            return control(me, ControlEndRead());
        }

    //}

    // pipe/pool -> pool
    template < typename TValue,
               typename TSpec,
               typename TSource >
    inline bool append(Pool<TValue, TSpec> &dest, TSource &src) {
        typename Size<TSource>::Type leftToRead = length(src);
        if (!beginRead(src)) return false;
        while (leftToRead) {
            push(dest, *src);
            ++src;
            --leftToRead;
        }
        endRead(src);
        return true;
    }

    // string -> pool
    template < typename TValue,
               typename TSpec,
               typename TStringSpec >
    inline bool append(Pool<TValue, TSpec> &dest, String<TValue, TStringSpec> &src) {
        typedef typename Iterator< String<TValue, TStringSpec> const, Standard >::Type TIter;
        TIter _cur = begin(src, Standard());
        TIter _end = end(src, Standard());
        while (_cur != _end) {
            push(dest, *_cur);
            ++_cur;
        }
        return true;
    }
/*!
 * @fn Pool#assign
 * @headerfile <seqan/pipe.h>
 * @brief Assigns one object to another object.
 *
 * @signature void assign(target, source);
 *
 * @param[out] target Reference to assign to.
 * @param[in]  source Value to assign.
 *
 * Assign value of source to target.
 */
    template < typename TValue,
               typename TSpec,
               typename TSource >
    inline bool assign(Pool<TValue, TSpec> &dest, TSource &src) {
        typename Size<TSource>::Type _size = length(src);
        resize(dest, _size);
        return beginWrite(dest) && append(dest, src) && endWrite(dest);
    }

    template < typename TValue,
               typename TSpec,
               typename TSource >
    inline bool operator<<(Pool<TValue, TSpec> &dest, TSource &src) {
        return assign(dest, src);
    }



    // pool -> string
    template < typename TValue1,
               typename TStringSpec,
               typename TValue2,
               typename TSpec >
    inline bool assign(String<TValue1, TStringSpec> &dest, Pool<TValue2, TSpec> &src) {
        typedef typename Iterator< String<TValue1, TStringSpec>, Standard >::Type TIter;
        typename Size< String<TValue1, TStringSpec> >::Type _size = length(src);
        resize(dest, _size);
        if (!beginRead(src)) return false;
        TIter _cur = begin(dest, Standard());
        TIter _end = end(dest, Standard());
        while (_cur != _end) {
            *_cur = *src;
            ++_cur;
            ++src;
        }
        endRead(src);
        return true;
    }

    template < typename TValue1,
               typename TStringSpec,
               typename TValue2,
               typename TSpec >
    inline bool operator<<(String<TValue1, TStringSpec> &dest, Pool<TValue2, TSpec> &src) {
        return assign(dest, src);
    }

}

#endif
