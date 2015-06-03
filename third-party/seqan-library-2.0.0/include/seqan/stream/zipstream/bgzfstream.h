/*
zipstream Library License:
--------------------------

The zlib/libpng License Copyright (c) 2003 Jonathan de Halleux.

This software is provided 'as-is', without any express or implied warranty. In no event will the authors be held liable for any damages arising from the use of this software.

Permission is granted to anyone to use this software for any purpose, including commercial applications, and to alter it and redistribute it freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.

2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.

3. This notice may not be removed or altered from any source distribution

Author: Jonathan de Halleux, dehalleux@pelikhan.com, 2003   (original zlib stream)
Author: David Weese, dave.weese@gmail.com, 2014             (extension to parallel block-wise compression in bgzf format)
*/

#ifndef BGZFSTREAM_HPP
#define BGZFSTREAM_HPP


#include <vector>
#include <iostream>
#include <algorithm>
#include <zlib.h>
#include "zutil.h"

namespace seqan {

/// default gzip buffer size,
/// change this to suite your needs
const size_t default_buffer_size = 4096;

const unsigned BGZF_MAX_BLOCK_SIZE = 64 * 1024;
const unsigned BGZF_BLOCK_HEADER_LENGTH = 18;
const unsigned BGZF_BLOCK_FOOTER_LENGTH = 8;
const unsigned ZLIB_BLOCK_OVERHEAD = 5; // 5 bytes block overhead (see 3.2.4. at http://www.gzip.org/zlib/rfc-deflate.html)

// Reduce the maximal input size, such that the compressed data
// always fits in one block even for level Z_NO_COMPRESSION.
const unsigned BGZF_BLOCK_SIZE = BGZF_MAX_BLOCK_SIZE - BGZF_BLOCK_HEADER_LENGTH - BGZF_BLOCK_FOOTER_LENGTH - ZLIB_BLOCK_OVERHEAD;

/// Compression strategy, see bgzf doc.
enum EStrategy
{
    StrategyFiltered = 1,
    StrategyHuffmanOnly = 2,
    DefaultStrategy = 0
};

template<
    typename Elem,
    typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_bgzf_streambuf : public std::basic_streambuf<Elem, Tr>
{
public:
    typedef std::basic_ostream<Elem, Tr>& ostream_reference;
    typedef ElemA char_allocator_type;
    typedef ByteT byte_type;
    typedef ByteAT byte_allocator_type;
    typedef byte_type* byte_buffer_type;
    typedef typename Tr::char_type char_type;
    typedef typename Tr::int_type int_type;

    typedef ConcurrentQueue<size_t, Suspendable<Limit> > TJobQueue;

    struct OutputBuffer
    {
        char    buffer[BGZF_MAX_BLOCK_SIZE];
        size_t  size;
    };

    struct BufferWriter
    {
        ostream_reference ostream;

        BufferWriter(ostream_reference ostream) :
            ostream(ostream)
        {}

        bool operator() (OutputBuffer const & outputBuffer)
        {
            ostream.write(outputBuffer.buffer, outputBuffer.size);
            return ostream.good();
        }
    };

    struct CompressionJob
    {
        typedef std::vector<char_type, char_allocator_type> TBuffer;

        TBuffer         buffer;
        size_t          size;
        OutputBuffer    *outputBuffer;

        CompressionJob() :
            buffer(BGZF_BLOCK_SIZE / sizeof(char_type), 0),
            size(0),
            outputBuffer(NULL)
        {}
    };

    // string of recycable jobs
    size_t                  numThreads;
    size_t                  numJobs;
    String<CompressionJob>  jobs;
    TJobQueue               jobQueue;
    TJobQueue               idleQueue;
    Serializer<
        OutputBuffer,
        BufferWriter>       serializer;

    size_t                  currentJobId;
    bool                    currentJobAvail;


    struct CompressionThread
    {
        basic_bgzf_streambuf            *streamBuf;
        CompressionContext<BgzfFile>    compressionCtx;
        size_t                          threadNum;

        void operator()()
        {
            ScopedReadLock<TJobQueue> readLock(streamBuf->jobQueue);
            ScopedWriteLock<TJobQueue> writeLock(streamBuf->idleQueue);

            // wait for a new job to become available
            bool success = true;
            while (success)
            {
                size_t jobId = -1;
                if (!popFront(jobId, streamBuf->jobQueue))
                    return;

                CompressionJob &job = streamBuf->jobs[jobId];

                // compress block with zlib
                job.outputBuffer->size = _compressBlock(
                    job.outputBuffer->buffer, sizeof(job.outputBuffer->buffer),
                    &job.buffer[0], job.size, compressionCtx);

                success = releaseValue(streamBuf->serializer, job.outputBuffer);
                appendValue(streamBuf->idleQueue, jobId);
            }
        }
    };

    // array of worker threads
    Thread<CompressionThread>   *threads;

    basic_bgzf_streambuf(ostream_reference ostream_,
                         size_t numThreads = 16,
                         size_t jobsPerThread = 8) :
        numThreads(numThreads),
        numJobs(numThreads * jobsPerThread),
        jobQueue(numJobs),
        idleQueue(numJobs),
        serializer(ostream_, numThreads * jobsPerThread)
    {
        resize(jobs, numJobs, Exact());
        currentJobId = 0;

        lockWriting(jobQueue);
        lockReading(idleQueue);
        setReaderWriterCount(jobQueue, numThreads, 1);
        setReaderWriterCount(idleQueue, 1, numThreads);

        for (unsigned i = 0; i < numJobs; ++i)
        {
            bool success = appendValue(idleQueue, i);
            ignoreUnusedVariableWarning(success);
            SEQAN_ASSERT(success);
        }

        threads = new Thread<CompressionThread>[numThreads];
        for (unsigned i = 0; i < numThreads; ++i)
        {
            threads[i].worker.streamBuf = this;
            threads[i].worker.threadNum = i;
            run(threads[i]);
        }

        currentJobAvail = popFront(currentJobId, idleQueue);
        SEQAN_ASSERT(currentJobAvail);

        CompressionJob &job = jobs[currentJobId];
        job.outputBuffer = aquireValue(serializer);
        this->setp(&job.buffer[0], &job.buffer[0] + (job.buffer.size() - 1));
    }

    ~basic_bgzf_streambuf()
    {
        // the buffer is now (after addFooter()) and flush will append the empty EOF marker
        flush(true);

        unlockWriting(jobQueue);
        unlockReading(idleQueue);

        for (unsigned i = 0; i < numThreads; ++i)
            waitFor(threads[i]);
        delete[] threads;
    }

    bool compressBuffer(size_t size)
    {
        // submit current job
        if (currentJobAvail)
        {
            jobs[currentJobId].size = size;
            appendValue(jobQueue, currentJobId);
        }

        // recycle existing idle job
        if (!(currentJobAvail = popFront(currentJobId, idleQueue)))
            return false;

        jobs[currentJobId].outputBuffer = aquireValue(serializer);

        return serializer;
    }

    int_type overflow(int_type c)
    {
        int w = static_cast<int>(this->pptr() - this->pbase());
        if (c != EOF)
        {
            *this->pptr() = c;
            ++w;
        }
        if (compressBuffer(w))
        {
            CompressionJob &job = jobs[currentJobId];
            this->setp(&job.buffer[0], &job.buffer[0] + (job.buffer.size() - 1));
            return c;
        }
        else
        {
            return EOF;
        }
    }

    std::streamsize flush(bool flushEmptyBuffer = false)
    {
        int w = static_cast<int>(this->pptr() - this->pbase());
        if ((w != 0 || flushEmptyBuffer) && compressBuffer(w))
        {
            CompressionJob &job = jobs[currentJobId];
            this->setp(&job.buffer[0], &job.buffer[0] + (job.buffer.size() - 1));
        }
        else
        {
            w = 0;
        }

        // wait for running compressor threads
        waitForMinSize(idleQueue, numJobs - 1);

        serializer.worker.ostream.flush();
        return w;
    }

    int sync()
    {
        if (this->pptr() != this->pbase())
        {
            int c = overflow(EOF);
            if (c == EOF)
                return -1;
        }
        return 0;
    }

    void addFooter()
    {
        // we flush the filled buffer here, so that an empty (EOF) buffer is flushed in the d'tor
        if (this->pptr() != this->pbase())
            overflow(EOF);
    }

    /// returns a reference to the output stream
    ostream_reference get_ostream() const    { return serializer.worker.ostream; };
};

template<
    typename Elem,
    typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_unbgzf_streambuf :
    public std::basic_streambuf<Elem, Tr>
{
public:
    typedef std::basic_istream<Elem, Tr>& istream_reference;
    typedef ElemA char_allocator_type;
    typedef ByteT byte_type;
    typedef ByteAT byte_allocator_type;
    typedef byte_type* byte_buffer_type;
    typedef typename Tr::char_type char_type;
    typedef typename Tr::int_type int_type;
    typedef typename Tr::off_type off_type;
    typedef typename Tr::pos_type pos_type;

    typedef std::vector<char_type, char_allocator_type>     TBuffer;
    typedef ConcurrentQueue<int, Suspendable<Limit> >       TJobQueue;

    static const size_t MAX_PUTBACK = 4;

    struct Serializer
    {
        istream_reference   istream;
        Mutex               lock;
        IOError             *error;
        off_type            fileOfs;

        Serializer(istream_reference istream) :
            istream(istream),
            lock(false),
            error(NULL),
            fileOfs(0u)
        {}

        ~Serializer()
        {
            delete error;
        }
    };

    Serializer serializer;

    struct DecompressionJob
    {
        typedef std::vector<byte_type, byte_allocator_type> TInputBuffer;

        TInputBuffer    inputBuffer;
        TBuffer         buffer;
        off_type        fileOfs;
        int             size;

        CriticalSection cs;
        Condition       readyEvent;
        bool            ready;

        DecompressionJob() :
            inputBuffer(BGZF_MAX_BLOCK_SIZE, 0),
            buffer(MAX_PUTBACK + BGZF_MAX_BLOCK_SIZE / sizeof(char_type), 0),
            fileOfs(),
            size(0),
            readyEvent(cs),
            ready(true)
        {}

        DecompressionJob(DecompressionJob const &other) :
            inputBuffer(other.inputBuffer),
            buffer(other.buffer),
            fileOfs(other.fileOfs),
            size(other.size),
            readyEvent(cs),
            ready(other.ready)
        {}
    };

    // string of recycable jobs
    size_t                      numThreads;
    size_t                      numJobs;
    String<DecompressionJob>    jobs;
    TJobQueue                   runningQueue;
    TJobQueue                   todoQueue;
    int                         currentJobId;

    struct DecompressionThread
    {
        basic_unbgzf_streambuf          *streamBuf;
        CompressionContext<BgzfFile>    compressionCtx;

        void operator()()
        {
            ScopedReadLock<TJobQueue> readLock(streamBuf->todoQueue);
            ScopedWriteLock<TJobQueue> writeLock(streamBuf->runningQueue);

            // wait for a new job to become available
            while (true)
            {
                int jobId = -1;
                if (!popFront(jobId, streamBuf->todoQueue))
                    return;

                DecompressionJob &job = streamBuf->jobs[jobId];
                size_t tailLen = 0;

                // typically the idle queue contains only ready jobs
                // however, if seek() fast forwards running jobs into the todoQueue
                // the caller defers the task of waiting to the decompression threads
                if (!job.ready)
                {
                    ScopedLock<CriticalSection> lock(job.cs);
                    if (!job.ready)
                    {
                        waitFor(job.readyEvent);
                        job.ready = true;
                    }
                }

                {
                    ScopedLock<Mutex> scopedLock(streamBuf->serializer.lock);

                    if (streamBuf->serializer.error != NULL)
                        return;

                    // remember start offset (for tellg later)
                    job.fileOfs = streamBuf->serializer.fileOfs;
                    job.size = -1;

                    // only load if not at EOF
                    if (job.fileOfs != -1)
                    {
                        // read header
                        streamBuf->serializer.istream.read(
                            (char*)&job.inputBuffer[0],
                            BGZF_BLOCK_HEADER_LENGTH);

                        if (!streamBuf->serializer.istream.good())
                        {
                            streamBuf->serializer.fileOfs = -1;
                            if (streamBuf->serializer.istream.eof())
                                goto eofSkip;
                            streamBuf->serializer.error = new IOError("Stream read error.");
                            return;
                        }

                        // check header
                        if (!_bgzfCheckHeader(&job.inputBuffer[0]))
                        {
                            streamBuf->serializer.fileOfs = -1;
                            streamBuf->serializer.error = new IOError("Invalid BGZF block header.");
                            return;
                        }

                        // extract length of compressed data
                        tailLen = _bgzfUnpack16(&job.inputBuffer[0] + 16) + 1u - BGZF_BLOCK_HEADER_LENGTH;

                        // read compressed data and tail
                        streamBuf->serializer.istream.read(
                            (char*)&job.inputBuffer[0] + BGZF_BLOCK_HEADER_LENGTH,
                            tailLen);

                        if (!streamBuf->serializer.istream.good())
                        {
                            streamBuf->serializer.fileOfs = -1;
                            if (streamBuf->serializer.istream.eof())
                                goto eofSkip;
                            streamBuf->serializer.error = new IOError("Stream read error.");
                            return;
                        }

                        streamBuf->serializer.fileOfs += BGZF_BLOCK_HEADER_LENGTH + tailLen;
                        job.ready = false;

                    eofSkip:
                        streamBuf->serializer.istream.clear(
                            streamBuf->serializer.istream.rdstate() & ~std::ios_base::failbit);
                    }

                    if (!appendValue(streamBuf->runningQueue, jobId))
                    {
                        // signal that job is ready
                        {
                            ScopedLock<CriticalSection> lock(job.cs);
                            job.ready = true;
                            signal(job.readyEvent);
                        }
                        return;
                    }
                }

                if (!job.ready)
                {
                    // decompress block
                    job.size = _decompressBlock(
                        &job.buffer[0] + MAX_PUTBACK, capacity(job.buffer),
                        &job.inputBuffer[0], BGZF_BLOCK_HEADER_LENGTH + tailLen, compressionCtx);

                    // signal that job is ready
                    {
                        ScopedLock<CriticalSection> lock(job.cs);
                        job.ready = true;
                        signal(job.readyEvent);
                    }
                }
            }
        }
    };

    // array of worker threads
    Thread<DecompressionThread> *threads;
    TBuffer                     putbackBuffer;

    basic_unbgzf_streambuf(istream_reference istream_,
                           size_t numThreads = 16,
                           size_t jobsPerThread = 8) :
        serializer(istream_),
        numThreads(numThreads),
        numJobs(numThreads * jobsPerThread),
        runningQueue(numJobs),
        todoQueue(numJobs),
        putbackBuffer(MAX_PUTBACK)
    {
        resize(jobs, numJobs, Exact());
        currentJobId = -1;

        lockReading(runningQueue);
        lockWriting(todoQueue);
        setReaderWriterCount(runningQueue, 1, numThreads);
        setReaderWriterCount(todoQueue, numThreads, 1);

        for (unsigned i = 0; i < numJobs; ++i)
        {
            bool success = appendValue(todoQueue, i);
            ignoreUnusedVariableWarning(success);
            SEQAN_ASSERT(success);
        }

        threads = new Thread<DecompressionThread>[numThreads];
        for (unsigned i = 0; i < numThreads; ++i)
        {
            threads[i].worker.streamBuf = this;
            run(threads[i]);
        }
    }

    ~basic_unbgzf_streambuf()
    {
        unlockWriting(todoQueue);
        unlockReading(runningQueue);

        for (unsigned i = 0; i < numThreads; ++i)
            waitFor(threads[i]);
        delete[] threads;
    }

    int_type underflow()
    {
        // no need to use the next buffer?
        if (this->gptr() && this->gptr() < this->egptr())
            return Tr::to_int_type(*this->gptr());

        size_t putback = this->gptr() - this->eback();
        if (putback > MAX_PUTBACK)
            putback = MAX_PUTBACK;

        // save at most MAX_PUTBACK characters from previous page to putback buffer
        if (putback != 0)
            std::copy(
                this->gptr() - putback,
                this->gptr(),
                &putbackBuffer[0]);

        if (currentJobId >= 0)
            appendValue(todoQueue, currentJobId);

        while (true)
        {
            if (!popFront(currentJobId, runningQueue))
            {
                currentJobId = -1;
                SEQAN_ASSERT(serializer.error != NULL);
                if (serializer.error != NULL)
                    throw *serializer.error;
                return EOF;
            }

            DecompressionJob &job = jobs[currentJobId];

            // restore putback buffer
            this->setp(&job.buffer[0], &job.buffer[0] + (job.buffer.size() - 1));
            if (putback != 0)
                std::copy(
                    &putbackBuffer[0],
                    &putbackBuffer[0] + putback,
                    &job.buffer[0] + (MAX_PUTBACK - putback));

            // wait for the end of decompression
            {
                ScopedLock<CriticalSection> lock(job.cs);
                if (!job.ready)
                    waitFor(job.readyEvent);
            }

            size_t size = (job.size != -1)? job.size : 0;

            // reset buffer pointers
            this->setg(
                  &job.buffer[0] + (MAX_PUTBACK - putback),     // beginning of putback area
                  &job.buffer[0] + MAX_PUTBACK,                 // read position
                  &job.buffer[0] + (MAX_PUTBACK + size));       // end of buffer

            if (job.size == -1)
                return EOF;
            else if (job.size > 0)
                return Tr::to_int_type(*this->gptr());      // return next character
        }
    }

    pos_type seekoff(off_type ofs, std::ios_base::seekdir dir, std::ios_base::openmode openMode)
    {
        if ((openMode & (std::ios_base::in | std::ios_base::out)) == std::ios_base::in)
        {
            if (dir == std::ios_base::cur && ofs >= 0)
            {
                // forward delta seek
                while (currentJobId < 0 || this->egptr() - this->gptr() < ofs)
                {
                    ofs -= this->egptr() - this->gptr();
                    if (this->underflow() == EOF)
                        break;
                }

                if (currentJobId >= 0 && ofs <= this->egptr() - this->gptr())
                {
                    DecompressionJob &job = jobs[currentJobId];

                    // reset buffer pointers
                    this->setg(
                          this->eback(),            // beginning of putback area
                          this->gptr() + ofs,       // read position
                          this->egptr());           // end of buffer

                    return pos_type((job.fileOfs << 16) + ((this->gptr() - &job.buffer[MAX_PUTBACK])));
                }

            }
            else if (dir == std::ios_base::beg)
            {
                // random seek
                std::streampos destFileOfs = ofs >> 16;

                // are we in the same block?
                if (currentJobId >= 0 && jobs[currentJobId].fileOfs == (off_type)destFileOfs)
                {
                    DecompressionJob &job = jobs[currentJobId];

                    // reset buffer pointers
                    this->setg(
                          this->eback(),                                        // beginning of putback area
                          &job.buffer[0] + (MAX_PUTBACK + (ofs & 0xffff)),      // read position
                          this->egptr());                                       // end of buffer
                    return ofs;
                }

                // ok, different block
                {
                    ScopedLock<Mutex> scopedLock(serializer.lock);

                    // remove all running jobs and put them in the idle queue unless we
                    // find our seek target

                    if (currentJobId >= 0)
                        appendValue(todoQueue, currentJobId);

                    // Note that if we are here the current job does not represent the sought block.
                    // Hence if the running queue is empty we need to explicitly unset the jobId,
                    // otherwise we would not update the serializers istream pointer to the correct position.
                    if (empty(runningQueue))
                        currentJobId = -1;

                    // empty is thread-safe in serializer.lock
                    while (!empty(runningQueue))
                    {
                        popFront(currentJobId, runningQueue);

                        if (jobs[currentJobId].fileOfs == (off_type)destFileOfs)
                            break;

                        // push back useless job
                        appendValue(todoQueue, currentJobId);
                        currentJobId = -1;
                    }

                    if (currentJobId == -1)
                    {
                        SEQAN_ASSERT(empty(runningQueue));
                        serializer.istream.clear(serializer.istream.rdstate() & ~std::ios_base::eofbit);
                        if (serializer.istream.rdbuf()->pubseekpos(destFileOfs, std::ios_base::in) == destFileOfs)
                            serializer.fileOfs = destFileOfs;
                        else
                            currentJobId = -2;      // temporarily signals a seek error
                    }
                }

                // if our block wasn't in the running queue yet, it should now
                // be the first that falls out after modifying serializer.fileOfs
                if (currentJobId == -1)
                    popFront(currentJobId, runningQueue);
                else if (currentJobId == -2)
                    currentJobId = -1;

                if (currentJobId >= 0)
                {
                    // wait for the end of decompression
                    DecompressionJob &job = jobs[currentJobId];
                    {
                        ScopedLock<CriticalSection> lock(job.cs);
                        if (!job.ready)
                            waitFor(job.readyEvent);
                    }

                    SEQAN_ASSERT_EQ(job.fileOfs, (off_type)destFileOfs);

                    // reset buffer pointers
                    this->setg(
                          &job.buffer[0] + MAX_PUTBACK,                     // no putback area
                          &job.buffer[0] + (MAX_PUTBACK + (ofs & 0xffff)),  // read position
                          &job.buffer[0] + (MAX_PUTBACK + job.size));       // end of buffer
                    return ofs;
                }
            }
        }
        return pos_type(off_type(-1));
    }

    pos_type seekpos(pos_type pos, std::ios_base::openmode openMode)
    {
        return seekoff(off_type(pos), std::ios_base::beg, openMode);
    }

    /// returns the compressed input istream
    istream_reference get_istream()    { return serializer.istream;};
};

/* \brief Base class for zip ostreams

Contains a basic_bgzf_streambuf.
*/
template<
    typename Elem,
    typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_bgzf_ostreambase : virtual public std::basic_ios<Elem,Tr>
{
public:
    typedef std::basic_ostream<Elem, Tr>& ostream_reference;
    typedef basic_bgzf_streambuf<
        Elem,
        Tr,
        ElemA,
        ByteT,
        ByteAT
        > bgzf_streambuf_type;

    basic_bgzf_ostreambase(ostream_reference ostream_)
        : m_buf(ostream_)
    {
        this->init(&m_buf );
    };

    /// returns the underlying zip ostream object
    bgzf_streambuf_type* rdbuf() { return &m_buf; };

    /// returns the bgzf error state
    int get_zerr() const                    {    return m_buf.get_err();};
    /// returns the uncompressed data crc
    long get_crc() const                    {    return m_buf.get_crc();};
    /// returns the compressed data size
    long get_out_size() const                {    return m_buf.get_out_size();};
    /// returns the uncompressed data size
    long get_in_size() const                {    return m_buf.get_in_size();};
private:
    bgzf_streambuf_type m_buf;
};

/* \brief Base class for unzip istreams

Contains a basic_unbgzf_streambuf.
*/
template<
    typename Elem,
    typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_bgzf_istreambase : virtual public std::basic_ios<Elem,Tr>
{
public:
    typedef std::basic_istream<Elem, Tr>& istream_reference;
    typedef basic_unbgzf_streambuf<
        Elem,
        Tr,
        ElemA,
        ByteT,
        ByteAT
        > unbgzf_streambuf_type;

    basic_bgzf_istreambase(istream_reference ostream_)
        : m_buf(ostream_)
    {
        this->init(&m_buf );
    };

    /// returns the underlying unzip istream object
    unbgzf_streambuf_type* rdbuf() { return &m_buf; };

    /// returns the bgzf error state
    int get_zerr() const                    {    return m_buf.get_zerr();};
    /// returns the uncompressed data crc
    long get_crc() const                    {    return m_buf.get_crc();};
    /// returns the uncompressed data size
    long get_out_size() const                {    return m_buf.get_out_size();};
    /// returns the compressed data size
    long get_in_size() const                {    return m_buf.get_in_size();};
private:
    unbgzf_streambuf_type m_buf;
};

/*brief A zipper ostream

This class is a ostream decorator that behaves 'almost' like any other ostream.

At construction, it takes any ostream that shall be used to output of the compressed data.

When finished, you need to call the special method zflush or call the destructor
to flush all the intermidiate streams.

Example:
\code
// creating the target zip string, could be a fstream
ostringstream ostringstream_;
// creating the zip layer
bgzf_ostream zipper(ostringstream_);


// writing data
zipper<<f<<" "<<d<<" "<<ui<<" "<<ul<<" "<<us<<" "<<c<<" "<<dum;
// zip ostream needs special flushing...
zipper.zflush();
\endcode
*/
template<
    typename Elem,
    typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_bgzf_ostream :
    public basic_bgzf_ostreambase<Elem,Tr,ElemA,ByteT,ByteAT>,
    public std::basic_ostream<Elem,Tr>
{
public:
    typedef basic_bgzf_ostreambase<
        Elem,Tr,ElemA,ByteT,ByteAT> bgzf_ostreambase_type;
    typedef std::basic_ostream<Elem,Tr> ostream_type;
    typedef ostream_type& ostream_reference;

    basic_bgzf_ostream(ostream_reference ostream_)
    :
        bgzf_ostreambase_type(ostream_),
        ostream_type(bgzf_ostreambase_type::rdbuf())
    {}

    /// flush inner buffer and zipper buffer
    basic_bgzf_ostream<Elem,Tr>& zflush()
    {
        this->flush(); this->rdbuf()->flush(); return *this;
    };

    ~basic_bgzf_ostream()
    {
        this->rdbuf()->addFooter();
    }

private:
    static void put_long(ostream_reference out_, unsigned long x_);
#ifdef _WIN32
    void _Add_vtordisp1() { } // Required to avoid VC++ warning C4250
    void _Add_vtordisp2() { } // Required to avoid VC++ warning C4250
#endif
};

/* \brief A zipper istream

This class is a istream decorator that behaves 'almost' like any other ostream.

At construction, it takes any istream that shall be used to input of the compressed data.

Simlpe example:
\code
// create a stream on zip string
istringstream istringstream_( ostringstream_.str());
// create unzipper istream
bgzf_istream unzipper( istringstream_);

// read and unzip
unzipper>>f_r>>d_r>>ui_r>>ul_r>>us_r>>c_r>>dum_r;
\endcode
*/
template<
    typename Elem,
    typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_bgzf_istream :
    public basic_bgzf_istreambase<Elem,Tr,ElemA,ByteT,ByteAT>,
    public std::basic_istream<Elem,Tr>
{
public:
    typedef basic_bgzf_istreambase<
        Elem,Tr,ElemA,ByteT,ByteAT> bgzf_istreambase_type;
    typedef std::basic_istream<Elem,Tr> istream_type;
    typedef istream_type& istream_reference;
    typedef char byte_type;

    basic_bgzf_istream(istream_reference istream_)
      :
        bgzf_istreambase_type(istream_),
        istream_type(bgzf_istreambase_type::rdbuf()),
        m_is_gzip(false),
        m_gbgzf_data_size(0)
    {};

    /// returns true if it is a gzip file
    bool is_gzip() const                {    return m_is_gzip;};
    /// return data size check
    bool check_data_size() const        {    return this->get_out_size() == m_gbgzf_data_size;};

    /// return the data size in the file
    long get_gbgzf_data_size() const        {    return m_gbgzf_data_size;};
protected:
    static void read_long(istream_reference in_, unsigned long& x_);

    int check_header();
    bool m_is_gzip;
    unsigned long m_gbgzf_data_size;

#ifdef _WIN32
private:
    void _Add_vtordisp1() { } // Required to avoid VC++ warning C4250
    void _Add_vtordisp2() { } // Required to avoid VC++ warning C4250
#endif
};

/// A typedef for basic_bgzf_ostream<char>
typedef basic_bgzf_ostream<char> bgzf_ostream;
/// A typedef for basic_bgzf_ostream<wchar_t>
typedef basic_bgzf_ostream<wchar_t> bgzf_wostream;
/// A typedef for basic_bgzf_istream<char>
typedef basic_bgzf_istream<char> bgzf_istream;
/// A typedef for basic_bgzf_istream<wchart>
typedef basic_bgzf_istream<wchar_t> bgzf_wistream;

}  // namespace seqan

#include "bgzfstream_impl.h"

#endif
