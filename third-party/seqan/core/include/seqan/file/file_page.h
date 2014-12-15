// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Implements lots of internal data structures for low-level file access used
// by the External String, external sorters and mappers (Pool), FileStream.
// ==========================================================================

#ifndef SEQAN_HEADER_FILE_PAGE_H
#define SEQAN_HEADER_FILE_PAGE_H


/* IOREV
 * _nottested_
 * _nodoc_
 *
 * not tested by any test or app
 * no documentation for the functions
 *
 * contains lots of seemingly important code, however doc is sparse
 * to non-existent. Code needs thorough investigation, to understand
 * how/if this works
 *
 * probably Weese's code according to holtgrew
 *
 *
 */


//////////////////////////////////////////////////////////////////////////////

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TConfig>
struct MMap;


	//////////////////////////////////////////////////////////////////////////////
	// base class for memory buffers

    struct Dynamic;

    template <size_t PAGESIZE>
    struct Fixed;

	template <typename TFile, typename TSpec = Dynamic>
    struct PageFrame;


	template <typename TValue, typename TSpec = Simple>
	struct Buffer
	{
//IOREV _nodoc_
		typedef	typename Size<Buffer>::Type                 TSize;
		typedef	typename Iterator<Buffer, Standard>::Type   TIterator;

		TIterator	begin;      // the beginning of the buffer
        TIterator	end;        // end of valid data
        TSize		pageSize;   // size of allocated memory

        Buffer():
            begin(NULL),
            end(NULL),
            pageSize(0) {}

        Buffer(TIterator _begin, TIterator _end):
            begin(_begin),
            end(_end),
            pageSize(0) {}

        Buffer(TIterator _begin, TSize _size):
            begin(_begin),
            end(_begin + _size),
            pageSize(0) {}

//        template <typename T>
//        Buffer(T &x):
//        {
//            const_cast<int&>(x);
//        }
//
//        Buffer(TSize _pageSize):
//            begin(NULL),
//            end(NULL),
//            pageSize(_pageSize) {}

        inline TValue       & operator[](TSize i)       { return begin[i]; }
        inline TValue const & operator[](TSize i) const { return begin[i]; }
	};


	//////////////////////////////////////////////////////////////////////////////
	// meta-function interface

	template <typename TValue, typename TSpec>
    struct Value<Buffer<TValue, TSpec> >
    {
        typedef TValue Type;
    };

	template <typename TValue, typename TSpec>
	struct Size<Buffer<TValue, TSpec> >
    {
        typedef size_t Type;
    };

	template <typename TValue, typename TSpec>
    struct Iterator<Buffer<TValue, TSpec>, Standard>
    {
        typedef TValue *Type;
    };

	template <typename TValue, typename TSpec>
    struct Iterator<Buffer<TValue, TSpec> const, Standard>
    {
        typedef TValue const *Type;
    };


	//////////////////////////////////////////////////////////////////////////////
	// global interface

	template <typename TValue, typename TSpec>
    inline typename Size<Buffer<TValue, TSpec> >::Type
    capacity(Buffer<TValue, TSpec> const &me)
    {
        return me.pageSize;
    }

	template <typename TValue, typename TFile, size_t PAGESIZE>
    inline typename Size<Buffer<TValue, PageFrame<TFile, Fixed<PAGESIZE> > > >::Type
    capacity(Buffer<TValue, PageFrame<TFile, Fixed<PAGESIZE> > > const &)
    {
        return PAGESIZE;
    }

	template <typename TValue, typename TSpec, typename TSize>
    inline void
    _setCapacity(Buffer<TValue, TSpec> &me, TSize size)
    {
        me.pageSize = size;
    }

	template <typename TValue, typename TSpec>
    inline typename Size<Buffer<TValue, TSpec> >::Type
    length(Buffer<TValue, TSpec> const &me)
    {
        return me.end - me.begin;
    }

	template <typename TValue, typename TFile, size_t PAGESIZE>
    inline typename Size<Buffer<TValue, PageFrame<TFile, Fixed<PAGESIZE> > > >::Type
    length(Buffer<TValue, PageFrame<TFile, Fixed<PAGESIZE> > > const &)
    {
        return PAGESIZE;
    }

	template <typename TValue, typename TSpec, typename TSize>
    inline void
    resize(Buffer<TValue, TSpec> &me, TSize size)
    {
        me.end = me.begin + size;
    }

    template <typename TValue, typename TSpec, typename TSize, typename T>
	inline bool
    allocPage(Buffer<TValue, TSpec> &buffer, TSize size, T const & me)
    {
//IOREV _nodoc_
        allocate(me, buffer.begin, size);
        _setCapacity(buffer, size);
        resize(buffer, size);
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "allocPage: " << ::std::hex << (void*)buffer.begin << ::std::dec << ::std::endl;
		#endif
        return buffer.begin != NULL;
	}

	template <typename TValue, typename TSpec, typename T>
	inline void
    freePage(Buffer<TValue, TSpec> &buffer, T const & me)
    {
//IOREV _nodoc_
		#ifdef SEQAN_VVERBOSE
			if ((void*)buffer.begin)
				::std::cerr << "freePage:  " << ::std::hex << (void*)buffer.begin << ::std::dec << ::std::endl;
		#endif
		deallocate(me, buffer.begin, capacity(buffer));
		buffer.begin = NULL;
        resize(buffer, 0);
        _setCapacity(buffer, 0);
	}

	template <typename TValue, typename TSpec>
	inline typename Iterator<Buffer<TValue, TSpec>, Standard>::Type
    begin(Buffer<TValue, TSpec> &pf, Standard)
    {
		return pf.begin;
	}

	template <typename TValue, typename TSpec>
	inline typename Iterator<Buffer<TValue, TSpec> const, Standard>::Type
    begin(Buffer<TValue, TSpec> const &pf, Standard)
    {
		return pf.begin;
	}

	template <typename TValue, typename TSpec>
	inline typename Iterator<Buffer<TValue, TSpec>, Standard>::Type
    end(Buffer<TValue, TSpec> &pf, Standard)
    {
		return pf.end;
	}

	template <typename TValue, typename TSpec>
	inline typename Iterator<Buffer<TValue, TSpec> const, Standard>::Type
    end(Buffer<TValue, TSpec> const &pf, Standard)
    {
		return pf.end;
	}



    //////////////////////////////////////////////////////////////////////////////
	// a bucket is a structure to represent a small window of a page
    // used by algorithms which need a global view of all pages (merge sort, mapper)

	template <typename TValue>
    struct PageBucket
	{
//IOREV _nodoc_ has some unformatted comments in code, but not enough doc
		typedef	typename Iterator<PageBucket, Standard>::Type TIterator;

        unsigned    pageOfs;                // begin of bucket window with relation to page begin
        TIterator	begin, cur, end;        // begin/end of buckets memory buffer and a pointer
    };

    template <typename TValue>
    struct PageBucketExtended : public PageBucket< TValue > {
//IOREV _nodoc_
		int     	pageNo;		            // related page (needed by merger sort)
    };
    
	template <typename TValue>
    ::std::ostream& operator<<(::std::ostream &out, const PageBucketExtended<TValue> &pb) {
//IOREV _ _nodoc_
        for(TValue *cur = pb.begin; cur != pb.end; cur++)
            out << *cur << " ";
        return out;
    }


	//////////////////////////////////////////////////////////////////////////////
	// meta-function interface

	template <typename TValue>
    struct Value< PageBucket<TValue> >		{ typedef TValue Type; };
//IOREV

    template <typename TValue>
    struct Size< PageBucket<TValue> >		{ typedef size_t Type; };
//IOREV

    template <typename TValue>
    struct Iterator< PageBucket<TValue>, Standard >			{ typedef TValue *Type; };
//IOREV

    template <typename TValue>
    struct Iterator< PageBucket<TValue> const, Standard >	{ typedef TValue const *Type; };
//IOREV

	//////////////////////////////////////////////////////////////////////////////

    template <typename TFilePos, typename TSize>
    struct FilePageRequest
    {
        enum FilePageRequestMode {
            MODE_RDONLY = 1,
            MODE_WRONLY = 2,
            MODE_RDWR   = 3
        };

        TFilePos            filePos;
        TSize               size;
        FilePageRequestMode mode;
    };

	//////////////////////////////////////////////////////////////////////////////

    enum PageFrameStatus
    {
        UNUSED,
        READING,
        PREPROCESSING,
        READY,
        POSTPROCESSING,
        POSTPROCESSED,
        WRITING
    };

	//////////////////////////////////////////////////////////////////////////////
	// page frame of dynamic size

    template <typename TValue, typename TFile>
	struct Buffer<TValue, PageFrame<TFile, Dynamic> > :
        public Buffer<TValue>
	{
//IOREV _nodoc_
		typedef Buffer<TValue>                      TBase;
		typedef TFile                               File;
		typedef typename Size<TFile>::Type          TFileSize;
        typedef typename AsyncRequest<TFile>::Type  AsyncRequest;
        typedef FilePageRequest<TFileSize, size_t>  TFilePageRequest;

        TFilePageRequest pageRequest;

		bool			dirty;		// data needs to be written to disk before freeing
		unsigned		pageNo;		// maps frames to pages (reverse vector mapper)
        AsyncRequest    request;    // request structure of the async io process
		PageFrameStatus status;
        Buffer          *next;      // next buffer in a chained list

        Buffer():
			TBase(),
            dirty(false),
            pageNo(MaxValue<unsigned>::VALUE),
			status(READY),
			next(NULL) {}
    };


	//////////////////////////////////////////////////////////////////////////////
	// page frame of static size

//IOREV

    typedef ::std::list<Position<String<void*> >::Type> PageLRUList;    // least recently usage list
	typedef PageLRUList::iterator PageLRUEntry;

    template < typename TValue,
               typename TFile,
               size_t PAGESIZE_ >
	struct Buffer<TValue, PageFrame<TFile, Fixed<PAGESIZE_> > >
    {
		typedef TFile                                       File;
        typedef typename AsyncRequest<TFile>::Type          AsyncRequest;
		typedef	typename Iterator<Buffer, Standard>::Type   TIterator;

		enum            { PAGESIZE = PAGESIZE_ };
		enum DataStatus	{ ON_DISK = -1, UNINITIALIZED = -2 };
		enum Priority	{ NORMAL_LEVEL = 0, PREFETCH_LEVEL = 1, ITERATOR_LEVEL = 2, PERMANENT_LEVEL = 3 };

        TIterator       begin;
        TIterator       end;
        AsyncRequest    request;    // request structure of the async io process
		PageFrameStatus status;
		DataStatus		dataStatus;
		PageLRUEntry	lruEntry;   // priority based lru
        Priority        priority;
		int     		pageNo;		// maps frames to pages (reverse vector mapper)
		bool			dirty;		// data needs to be written to disk before freeing

		Buffer():
            begin(NULL),
            end(NULL),
			status(READY),
			dataStatus(UNINITIALIZED),
            priority(NORMAL_LEVEL),
            pageNo(-1),
            dirty(false) {}

        template <typename TPos>
		inline TValue &
        operator[] (TPos i)
        {
            return begin[i];
        }

        template <typename TPos>
        inline TValue const &
        operator[] (TPos i) const
        {
            return begin[i];
        }
	};


	//////////////////////////////////////////////////////////////////////////////
	// a memory-mapped page frame
/*
    template <typename TValue, typename TFile, typename TConfig>
	struct PageFrame< TValue, TFile, MMap<TConfig> >: public Buffer<TValue>
    {
		typedef Buffer<TValue> TBase;

        enum Status		{ READY, READING, WRITING };

		unsigned   		pageNo;		// maps frames to pages (reverse vector mapper)
		Status			status;
        PageFrame       *next;      // next buffer in a chained list

        PageFrame():
			TBase(),
            pageNo(MaxValue<unsigned>::VALUE),
			next(NULL) {}
    };
*/

	//////////////////////////////////////////////////////////////////////////////
	// meta-function interface

	template <typename TValue, typename TFile, size_t PAGESIZE>
    struct Iterator<Buffer<TValue, PageFrame<TFile, Fixed<PAGESIZE> > >, Standard >
    {
        typedef VolatilePtr<TValue> Type;
    };

	template <typename TValue, typename TFile, size_t PAGESIZE>
    struct Iterator<Buffer<TValue, PageFrame<TFile, Fixed<PAGESIZE> > > const, Standard >
    {
        typedef VolatilePtr<TValue> const Type;
    };




	//////////////////////////////////////////////////////////////////////////////
	// various page frame methods

	template <typename TValue, typename TFile, typename TSpec>
    inline const char *
    _pageFrameStatusString(Buffer<TValue, PageFrame<TFile, TSpec> > const &pf)
    {
        switch (pf.status)
        {
			case READY:
                return "READY";
			case READING:
                return "READING";
			case WRITING:
                return "WRITING";
            case UNUSED:
                return "UNUSED";
            case PREPROCESSING:
                return "PREPROCESSING";
            case POSTPROCESSING:
                return "POSTPROCESSING";
            case POSTPROCESSED:
                return "POSTPROCESSED";
            default:
                return "UNKNOWN";
        }
    }

	template <typename TValue, typename TFile, typename TSpec>
    ::std::ostream& operator<<(::std::ostream &out, const Buffer<TValue, PageFrame<TFile, TSpec> > &pf)
	{
//IOREV _nodoc_
        out << "PageFrame @ " << pf.pageNo;
        if (pf.dirty)
            out << " DIRTY ";
        else
            out << " CLEAN ";

        out << _pageFrameStatusString(pf);

//        if (pf.dataStatus == pf.ON_DISK)
//            out << " ON_DISK";
//        else
//            out << " UNITIALIZED";
//
//        out << " Prio:" << pf.priority;
        out << " Buffer:" << (void*)pf.begin;

        return out;
	}

	template <typename TValue, typename TFile, typename TSpec, typename T> inline
	void allocPage(Buffer<TValue, PageFrame<TFile, TSpec> > &pf, T const & me) 
	{
//IOREV _nodoc_
		TValue* tmp = NULL;
		allocate(me, tmp, capacity(pf));
		//bzero(tmp, sizeof(TValue) * capacity(pf));
		pf.begin = tmp;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "allocPage: " << ::std::hex << &pf << '\t' << tmp << ::std::dec << ::std::endl;
		#endif
	}

	template <typename TValue, typename TFile, typename TSpec, typename T> inline
	void freePage(Buffer<TValue, PageFrame<TFile, TSpec> > &pf, T const & me) 
	{
//IOREV _nodoc_
		#ifdef SEQAN_VVERBOSE
			if ((void*)pf.begin)
				::std::cerr << "freePage:  " << ::std::hex << &pf << '\t' << (void*)pf.begin << ::std::dec << ::std::endl;
		#endif
        nukeCopies(pf.begin);
		deallocate(me, (TValue*)pf.begin, capacity(pf));
		pf.begin = NULL;
        resize(pf, 0);
	}

	template <typename TValue, typename TFile, typename TSpec> inline
	bool readPage(int pageNo, Buffer<TValue, PageFrame<TFile, TSpec> > &pf, TFile &file)
	{
//IOREV _nodoc_
		typedef typename Position<TFile>::Type TPos;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "readPage:  " << ::std::hex << &pf << '\t' << (void*)pf.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << ::std::endl;
		#endif
		pf.dirty = false;
		pf.status = READING;
//        resize(pf, capacity(pf));

		bool readResult = asyncReadAt(file, (TValue*)pf.begin, length(pf), (TPos)pageNo * (TPos)capacity(pf), pf.request);

        // TODO(weese): Throw an I/O exception
        if (!readResult)
            SEQAN_FAIL("%s operation could not be initiated: \"%s\"", _pageFrameStatusString(pf), strerror(errno));

        return readResult;
	}

	template <typename TValue, typename TFile, typename TSpec> inline
	bool unmapPage(Buffer<TValue, PageFrame<TFile, TSpec> > &pf, TFile &)
	{
        if (pf.begin)
        {
#ifdef PLATFORM_WINDOWS
#else
            munmap(pf.begin, length(pf) * sizeof(TValue));
#endif
			pf.begin = NULL;
        }
        return true;
    }

	template <typename TValue, typename TFile, typename TSpec, typename TSize> inline
	bool mapReadPage(Buffer<TValue, PageFrame<TFile, TSpec> > &pf, TFile &file, TSize size)
	{
		typedef typename Position<TFile>::Type TPos;
        unmapPage(pf, file);
		pf.status = READY;
        SEQAN_ASSERT_GT(size, 0);

#ifdef PLATFORM_WINDOWS
#if 0
		DWORD prot = 0;
		DWORD access = 0;
		if ((me._openMode & OPEN_MASK) == OPEN_RDONLY) 
		{
			prot = PAGE_READONLY;
			access = FILE_MAP_READ;
		} else {
			prot = PAGE_READWRITE;
			access = FILE_MAP_ALL_ACCESS;
		}
        LARGE_INTEGER largeSize;
		largeSize.QuadPart = new_capacity;
		largeSize.QuadPart *= sizeof(TValue);

		me.handle = CreateFileMapping(file.handle, &MMapStringDefaultAttributes, prot, largeSize.HighPart, largeSize.LowPart, NULL);
		if (me.handle == NULL)
		{
            SEQAN_FAIL(
                "CreateFileMapping(%d, %d, 0, 0, 0) failed in mapReadPage: \"%s\"",
                file.handle, access, strerror(errno));
		#if SEQAN_ENABLE_DEBUG
			::std::cerr << "CreateFileMapping failed. (ErrNo=" << GetLastError() << ")" << ::std::endl;
		#endif
			return false;
		}

		pf.begin = (TValue *) MapViewOfFile(pf.handle, access, 0, 0, 0);	
		if (pf.begin == NULL)
		{
            SEQAN_FAIL(
                "MapViewOfFile(%d, %d, 0, 0, 0) failed in mapReadPage: \"%s\"",
                pf.handle, access, strerror(errno));
			return false;
		}
				
		me.data_begin = (TValue *) addr;
#endif
#else
        pf.begin = (TValue*)mmap(NULL, size * sizeof(TValue), PROT_READ | PROT_WRITE, MAP_PRIVATE, file.handle, (TPos)pf.pageNo * (TPos)capacity(pf) * (TPos)sizeof(TValue));

        if (pf.begin == MAP_FAILED)
        {
            pf.begin = pf.end = NULL;
            SEQAN_FAIL(
                "mmap(NULL, %d, PROT_READ | PROT_WRITE, MAP_PRIVATE, %d, %d) failed in mapReadPage: \"%s\"",
                size * sizeof(TValue), file.handle, (TPos)pf.pageNo * (TPos)capacity(pf) * (TPos)sizeof(TValue), strerror(errno));
            return false;
        }
#endif
        pf.end = pf.begin + size;
        return true;
    }

	template <typename TValue, typename TFile, typename TSpec, typename TSize> inline
	bool mapWritePage(Buffer<TValue, PageFrame<TFile, TSpec> > &pf, TFile &file, TSize size) 
	{
		typedef typename Position<TFile>::Type TPos;
        unmapPage(pf, file);
        pf.status = READY;
        SEQAN_ASSERT_GT(size, 0u);
        
#ifdef PLATFORM_WINDOWS
#else
        pf.begin = (TValue*)mmap(NULL, size * sizeof(TValue), PROT_READ | PROT_WRITE, MAP_SHARED, file.handle, (TPos)pf.pageNo * (TPos)capacity(pf) * (TPos)sizeof(TValue));

        if (pf.begin == MAP_FAILED)
        {
            pf.begin = pf.end = NULL;
            SEQAN_FAIL(
                "mmap(NULL, %d, PROT_READ | PROT_WRITE, MAP_SHARED, %d, %d) failed in mapWritePage: \"%s\"",
                size * sizeof(TValue), file.handle, (TPos)pf.pageNo * (TPos)capacity(pf) * (TPos)sizeof(TValue), strerror(errno));
            return false;
        }
#endif
        pf.end = pf.begin + size;
        return true;
    }

	template <typename TValue, typename TFile, typename TSpec> inline
	bool writePage(Buffer<TValue, PageFrame<TFile, TSpec> > &pf, int pageNo, TFile &file) 
	{
//IOREV _nodoc_
		typedef typename Position<TFile>::Type TPos;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "writePage: " << ::std::hex << &pf << '\t' << (void*)pf.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << ::std::endl;
		#endif
		pf.status = WRITING;
//        resize(pf, capacity(pf));

		bool writeResult = asyncWriteAt(file, (TValue*)pf.begin, length(pf), (TPos)pageNo * (TPos)capacity(pf), pf.request);

        // TODO(weese): Throw an I/O exception
        if (!writeResult)
            SEQAN_FAIL("%s operation could not be initiated: \"%s\"", _pageFrameStatusString(pf), strerror(errno));

        return writeResult;
	}

	template < typename TValue, typename TFile, typename TSpec, typename TSize> inline
    bool readLastPage(int pageNo, Buffer<TValue, PageFrame<TFile, TSpec> > &pf, TFile &file, TSize size) 
	{
//IOREV _nodoc_
		typedef typename Position<TFile>::Type TPos;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "readPage:  " << ::std::hex << &pf << '\t' << (void*)pf.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << " size " << size << ::std::endl;
		#endif
		pf.dirty = false;
		pf.status = READY;
//        resize(pf, size);

		bool readResult = readAt(file, (TValue*)pf.begin, size, (TPos)pageNo * (TPos)capacity(pf));

        // TODO(weese): Throw an I/O exception
        if (!readResult)
            SEQAN_FAIL("%s operation could not be initiated: \"%s\"", _pageFrameStatusString(pf), strerror(errno));

        return readResult;
	}

	template <typename TValue, typename TFile, typename TSpec, typename TSize> inline
	bool writeLastPage(Buffer<TValue, PageFrame<TFile, TSpec> > &pf, int pageNo, TFile &file, TSize size) 
	{
//IOREV _nodoc_
		typedef typename Position<TFile>::Type TPos;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "writePage: " << ::std::hex << &pf << '\t' << (void*)pf.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << " size " << size << ::std::endl;
		#endif
		pf.dirty = false;
		pf.status = READY;
//        resize(pf, size);

		bool writeResult = writeAt(file, (TValue*)pf.begin, size, (TPos)pageNo * (TPos)capacity(pf));

        // TODO(weese): Throw an I/O exception
        if (!writeResult)
            SEQAN_FAIL("%s operation could not be initiated: \"%s\"", _pageFrameStatusString(pf), strerror(errno));

        return writeResult;
	}

	template <typename TValue, typename TFile, typename TSpec> inline
	bool waitFor(Buffer<TValue, PageFrame<TFile, TSpec> > &pf) 
	{
//IOREV _nodoc_ equally named functions with different purposes in system_thread.h, system_event.h and file_async.h
		if (pf.status == READY)
            return true;

        bool waitResult = waitFor(pf.request);

        // TODO(weese): Throw an I/O exception
        if (!waitResult)
        {
            //printRequest(pf.request);
            SEQAN_FAIL("%s operation could not be completed: \"%s\"", _pageFrameStatusString(pf), strerror(errno));
        }

        pf.status = READY;
        pf.dirty = false;
        return waitResult;
	}

	template <typename TValue, typename TFile, typename TSpec, typename TTime> inline
	bool waitFor(Buffer<TValue, PageFrame<TFile, TSpec> > &pf, TTime timeOut, bool &inProgress)
	{
//IOREV _nodoc_ equally named functions with different purposes in system_thread.h, system_event.h and file_async.h (none documented)
		if (pf.status == READY)
        {
            inProgress = false;
            return true;
        }

        bool waitResult = waitFor(pf.request, timeOut, inProgress);

        // TODO(weese): Throw an I/O exception
        if (!waitResult)
        {
            //printRequest(pf.request);
            SEQAN_FAIL("%s operation could not be completed: \"%s\"", _pageFrameStatusString(pf), strerror(errno));
        }

        if (!inProgress)
        {
            pf.status = READY;
            pf.dirty = false;
        }
        return waitResult;
	}

	template <typename TValue, typename TFile, typename TSpec> inline
	bool cancel(Buffer<TValue, PageFrame<TFile, TSpec> > &pf, TFile &file) 
	{
//IOREV _nodoc_ equally named functions with different purposes in pipe/pool_*.h, system_thread.h, and file_async.h (only the latter documented)
        bool dummy;
        waitFor(pf, 0, dummy);
		if (pf.status != READY) 
		{
            if (!cancel(file, pf.request)) return false;
            pf.status = READY;
        }
        return true;
	}

	//////////////////////////////////////////////////////////////////////////////
	// page based read/write methods used by Pool classes

    template <typename TValue, typename TFile> inline
	bool readPage(Buffer<TValue, PageFrame<TFile, Dynamic> > &pf, TFile &file) 
	{
//IOREV _nodoc_
        if (length(pf) == capacity(pf))
            return readPage(pf.pageNo, pf, file);
        else
            return readLastPage(pf.pageNo, pf, file, length(pf));
	}

	template <typename TValue, typename TFile> inline
	bool writePage(Buffer<TValue, PageFrame<TFile, Dynamic> > &pf, TFile &file) 
	{
//IOREV _nodoc_
        if (length(pf) == capacity(pf))
            return writePage(pf, pf.pageNo, file);
        else
            return writeLastPage(pf, pf.pageNo, file, length(pf));
	}

	template <typename TValue, typename TFile> inline
	unsigned readBucket(PageBucket<TValue> &b, int pageNo, size_t pageSize, size_t dataSize, TFile &file) 
	{
//IOREV _nodoc_
		typedef typename Position<TFile>::Type TPos;
        size_t readSize = _min(dataSize - b.pageOfs, (size_t)(b.end - b.begin));
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "readBucket:  " << ::std::hex << b.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << " at " << (TPos)pageNo * (TPos)pageSize + b.pageOfs;
			::std::cerr << " size " << readSize << ::std::endl;
		#endif
        if (readSize && readAt(file, b.begin, readSize, (TPos)pageNo * (TPos)pageSize + b.pageOfs)) {
            b.pageOfs += readSize;
            b.cur = b.begin;
            b.end = b.begin + readSize;
            return readSize;
        } else
            return 0;
	}

	template <typename TValue, typename TFile> inline
	bool writeBucket(PageBucket<TValue> &b, int pageNo, size_t pageSize, TFile &file) 
	{
//IOREV _nodoc_
		typedef typename Position<TFile>::Type TPos;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "writeBucket: " << ::std::hex << b.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << " at " << (TPos)pageNo * (TPos)pageSize + b.pageOfs;
			::std::cerr << " size " << b.cur - b.begin << ::std::endl;
		#endif
        if ((b.cur == b.begin) || writeAt(file, b.begin, b.cur - b.begin, (TPos)pageNo * (TPos)pageSize + b.pageOfs)) {
            b.pageOfs += b.cur - b.begin;
            b.cur = b.begin;
            return true;
        } else
            return false;
	}

	template <typename TValue, typename TPageOfs, typename TFile> inline
	bool writeBucket(Buffer<TValue, PageFrame<TFile, Dynamic> > &pf, TPageOfs &pageOfs, TFile &file)
	{
//IOREV _nodoc_
		typedef typename Position<TFile>::Type TPos;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "writeBucket: " << ::std::hex << pf.begin;
			::std::cerr << " from page " << ::std::dec << pf.pageNo << " at " << (TPos)pf.pageNo * (TPos)capacity(pf) + pageOfs;
			::std::cerr << " size " << length(pf) << ::std::endl;
		#endif
        if (pf.end == pf.begin) return true;
        if (asyncWriteAt(file, pf.begin, length(pf), (TPos)pf.pageNo * (TPos)capacity(pf) + pageOfs, pf.request)) {
            pf.status = WRITING;
            pageOfs += length(pf);
            return true;
        } else
            return false;
	}


	//////////////////////////////////////////////////////////////////////////////

    template <typename TPageFrame>
    struct PageChain 
	{
//IOREV _nodoc_
        TPageFrame          *first, *last;
        unsigned            frames, maxFrames;
        
        PageChain(unsigned _maxFrames = 0):
            first(NULL),
            last(NULL),
            frames(0),
            maxFrames(_maxFrames)
        {
            for(unsigned i = 0; i < _maxFrames; ++i)
                pushBack();
        }
        
        ~PageChain()
        {
            while (!empty())
                delete &popFront();
        }
        
        inline TPageFrame & operator[](int k) 
		{
            TPageFrame *p = first;
            for (; k > 0; --k)
                p = p->next;
            return *p;
        }
        
        inline TPageFrame const & operator[](int k) const 
		{
            TPageFrame *p = first;
            for (; k > 0; --k)
                p = p->next;
            return *p;
        }

        inline TPageFrame * getReadyPage() 
		{
            if (empty())
                return &pushBack();

            if (frames < maxFrames)
            {
                bool inProgress;
                bool waitResult = waitFor(*first, 0, inProgress);
                
                // TODO(weese): Throw an I/O exception
                if (!waitResult)
                    SEQAN_FAIL("%s operation could not be completed: \"%s\"", _pageFrameStatusString(*first), strerror(errno));

                if (!inProgress)
                    return &pushBack();
            }

            bool waitResult = waitFor(*first);

            // TODO(weese): Throw an I/O exception
            if (!waitResult)
                SEQAN_FAIL("%s operation could not be completed: \"%s\"", _pageFrameStatusString(*first), strerror(errno));

            return &firstToEnd();
        }

        template <typename TFile>
        inline void cancelAll(TFile &file) 
		{
            TPageFrame *p = first;
            for (; p != NULL; p = p->next)
                cancel(*p, file);
        }

        inline void waitForAll() 
		{
            for (TPageFrame *p = first; p != NULL; p = p->next)
            {
                bool waitResult = waitFor(*p);

                // TODO(weese): Throw an I/O exception
                if (!waitResult)
                    SEQAN_FAIL("%s operation could not be completed: \"%s\"", _pageFrameStatusString(*p), strerror(errno));
            }
        }

    public:

        inline bool
        empty()
        {
            return first == NULL;
        }

        inline unsigned
        size()
        {
            return frames;
        }

        inline TPageFrame &
        front()
        {
            return *first;
        }

        inline TPageFrame &
        back()
        {
            return *last;
        }

        //
        // append a PageFrame at the end of this chain
        //
        inline TPageFrame &
        pushBack(TPageFrame &pageFrame)
		{
            if (last != NULL)
                last->next = &pageFrame;
            else
                first = &pageFrame;
            
            last = &pageFrame;
            pageFrame.next = NULL;
            ++frames;

            return pageFrame;
        }

        //
        // append a *new* PageFrame at the end of this chain
        //
        inline TPageFrame &
        pushBack()
		{
            TPageFrame *pageFrame = new TPageFrame();
            return pushBack(*pageFrame);
        }

        //
        // remove a PageFrame from this chain
        //
        inline TPageFrame &
        erase(TPageFrame & pageFrame)
        {
            TPageFrame *prev = NULL;
            TPageFrame *p;

            for (p = first; p != NULL; p = p->next)
            {
                if (p == &pageFrame)
                {
                    // found it, now detach from neighbors
                    if (prev == NULL)
                        first = p->next;
                    else
                        prev->next = p->next;

                    p->next = NULL;

                    // were we the last?
                    if (last == p)
                        last = prev;

                    --frames;
                    break;
                }
                prev = p;
            }
            // we assert that the page has been found
            SEQAN_ASSERT(p != NULL);
            return pageFrame;
        }

        inline TPageFrame &
        popFront()
        {
            // equivalent to erase(front())
            TPageFrame *p = first;
            if (p != NULL)
            {
                first = p->next;
                p->next = NULL;
                if (first == NULL)
                    last = NULL;
                --frames;
            }
            else
                SEQAN_FAIL("pop() was called on an empty PageChain.");
            return *p;
        }

        //
        // move first page to the end of this chain
        // return the moved page
        //
        inline TPageFrame &
        firstToEnd()
		{
            // equivalent to pushBack(popFront())
            last->next = first;
            last = first;
            first = first->next;
            last->next = NULL;
            return *last;
        }


// call pop(chain.front())

//        inline TPageFrame &
//        popFront()
//		{
//            TPageFrame *p = first;
//            if (p) {
//                first = first->next;
//                if (!first) last = NULL;
//                --frames;
//                delete p;
//            }
//            return p;
//        }
    };

	//////////////////////////////////////////////////////////////////////////////
	// meta-function interface

    template <typename TPageFrame>    
	struct Value< PageChain<TPageFrame> > 
	{
//IOREV
		typedef TPageFrame Type;
	};

    template <typename TPageFrame>    
	struct Size< PageChain<TPageFrame> > 
	{
//IOREV
		typedef unsigned Type;
	};

    template <typename TPageFrame>    
	struct Iterator< PageChain<TPageFrame> > 
	{
//IOREV
		typedef TPageFrame *Type;
	};

    template <typename TPageFrame>    
	struct Iterator< PageChain<TPageFrame> const > 
	{
//IOREV
		typedef TPageFrame const *Type;
	};


	//////////////////////////////////////////////////////////////////////////////
	// page container with lru mechanism
	// the lru strategy uses different priorities
	// the page with the least priority is used first
	// 0..random access pages
	// 1..forward iterator pages
	// 2..quasi permanent pages

    template < typename TPageFrame,
               unsigned FRAMES,
               unsigned PRIORITY_LEVELS = TPageFrame::PERMANENT_LEVEL + 1 >
	struct PageContainer
	{
//IOREV _nodoc_ has some in-code comments
		typedef String<TPageFrame>					TPages;
		typedef typename Position<TPages>::Type		TPos;

		enum { PriorityLevels = PRIORITY_LEVELS };

		TPages			pages;
        PageLRUList		*lruList;

		PageContainer()
		{
			lruList = new PageLRUList[PRIORITY_LEVELS];
			resize(pages, FRAMES, Exact());
            for(TPos i = 0; i < FRAMES; ++i)
			    pages[i].lruEntry = lruList[0].insert(lruList[0].end(), i);
        }

		~PageContainer()
		{
			delete[] lruList;
		}

        inline TPageFrame       & operator[](TPos i)       { return pages[i]; }
        inline TPageFrame const & operator[](TPos i) const { return pages[i]; }

		inline void push_back() 
		{
			TPos last = length(pages);
			resize(pages, last + 1);
			pages[last].lruEntry = lruList[0].insert(lruList[0].end(), last);
		}

		inline void erase(int frameNo) 
		{
			lruList[pages[frameNo].priority].erase(pages[frameNo].lruEntry);
            seqan::erase(pages, frameNo);
		}

        inline void rename(int frameNo) 
		{
            *(pages[frameNo].lruEntry) = frameNo;
        }

		inline void pop_back() 
		{
			lruList[back(pages).priority].erase(back(pages).lruEntry);
			seqan::erase(pages, endPosition(pages) - 1);
		}

        inline void _print()
        {
            std::cout << std::endl;
            for(TPos i = 0; i < FRAMES; ++i)
			    std::cout << pages[i] << std::endl;
        }

		//////////////////////////////////////////////////////////////////////////////
		// lru strategy interface

		inline void upgrade(const TPageFrame &pf) 
		{
			lruList[pf.priority].splice(lruList[pf.priority].begin(), lruList[pf.priority], pf.lruEntry);
			pf.lruEntry = lruList[pf.priority].begin();
		}

		inline void downgrade(const TPageFrame &pf) 
		{
			lruList[pf.priority].splice(lruList[pf.priority].end(), lruList[pf.priority], pf.lruEntry);
			pf.lruEntry = lruList[pf.priority].end();
			--pf.lruEntry;
		}

		inline void upgrade(TPageFrame &pf, int newPriority) 
		{
			lruList[newPriority].splice(lruList[newPriority].begin(), lruList[pf.priority], pf.lruEntry);
			pf.lruEntry = lruList[newPriority].begin();
			pf.priority = static_cast<typename TPageFrame::Priority> (newPriority);
		}

		inline void downgrade(TPageFrame &pf, int newPriority) 
		{
			lruList[newPriority].splice(lruList[newPriority].end(), lruList[pf.priority], pf.lruEntry);
			pf.lruEntry = lruList[newPriority].end();
			--pf.lruEntry;
			pf.priority = static_cast<typename TPageFrame::Priority> (newPriority);
		}

		inline void _dump() 
		{
			for(unsigned i = 0; i < PRIORITY_LEVELS; ++i) {
                ::std::cerr << "|";
                PageLRUList::const_iterator I = lruList[i].end();
                PageLRUList::const_iterator first = lruList[i].begin();
				while (I != first) {
					--I;
                    TPageFrame &pf = pages[*I];
                    ::std::cerr << pf.pageNo;
                    if (pf.dirty) ::std::cerr << "*";
                    else          ::std::cerr << " ";
                    if (pf.status == READY) ::std::cerr << "  ";
                    else                    ::std::cerr << ". ";
				};
            }
            ::std::cerr << ::std::endl;
		}

        // Function is a functor which is called with a PageFrame object,
        // that is dirty or not READY (in an IO transfer)
		template <class Function>
		inline int mru(Function Func_, unsigned maxLevel = PRIORITY_LEVELS - 1) 
		{
			for(unsigned i = 0; i <= maxLevel; ++i) {
                PageLRUList::const_iterator I = lruList[i].end();
                PageLRUList::const_iterator first = lruList[i].begin();
				while (I != first) {
					--I;
					TPageFrame& pf = pages[*I];
					if (pf.status == READY && !pf.dirty)
						return *I;
					else
						if (Func_(pf)) return *I;
				};
            }
			#ifdef SEQAN_VVERBOSE
				::std::cerr << "ALL PAGES DIRTY Or IN USE (try to use const iterators) :-(" << ::std::endl;
			#endif
			return -1;
		}

		inline int mruDirty() 
		{
            for(unsigned i = 0; i < PRIORITY_LEVELS; ++i)
                if (!lruList[i].empty())
                    return lruList[i].back();
			return -1;
		}
	};


	//////////////////////////////////////////////////////////////////////////////
	// meta-function interface

	template <typename TPageFrame, unsigned FRAMES, unsigned PRIORITY_LEVELS>
	struct Host< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> > 
	{
//IOREV
		typedef String<TPageFrame> Type;
	};

	template <typename TPageFrame, unsigned FRAMES, unsigned PRIORITY_LEVELS>
	struct Host< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> const > 
	{
//IOREV
		typedef String<TPageFrame> const Type;
	};

	template <typename TPageFrame, unsigned FRAMES, unsigned PRIORITY_LEVELS>
	struct Value< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> > 
	{
//IOREV
		typedef TPageFrame Type;
	};

	template <typename TPageFrame, unsigned FRAMES, unsigned PRIORITY_LEVELS>
	struct Size< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> >:
		public Size< typename Host< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> >::Type> {};
//IOREV

	template <typename TPageFrame, unsigned FRAMES, unsigned PRIORITY_LEVELS>
	struct Position< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> >:
		public Position< typename Host< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> >::Type> {};
//IOREV

	template <typename TPageFrame, unsigned FRAMES, unsigned PRIORITY_LEVELS>
	struct Iterator< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> >:
		public Iterator< typename Host< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> >::Type> {};
//IOREV

	template <typename TPageFrame, unsigned FRAMES, unsigned PRIORITY_LEVELS>
	struct Iterator< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> const >:
		public Iterator< typename Host< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> const>::Type> {};
//IOREV


	//////////////////////////////////////////////////////////////////////////////
	// global interface

	template <typename TPageFrame, unsigned FRAMES, unsigned PRIORITY_LEVELS, typename TExpand>
	inline void reserve(
		PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> &pageCont,
		unsigned Count_,
		Tag<TExpand> expand)
	{
//IOREV
		reserve(pageCont.pages, Count_, expand);
	}

    template <typename TPageFrame, unsigned FRAMES, unsigned PRIORITY_LEVELS>
	inline void resize(
		PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> &pageCont, 
		unsigned Count_) 
	{
//IOREV
		unsigned Size_ = length(pageCont.pages);
		if (Size_ < Count_) {
			reserve(pageCont.pages, Count_);
            for(unsigned i = Size_; i < Count_; ++i)
                pageCont.push_back();
		} else 
			if (Size_ > Count_)
				for(unsigned i = Count_; i < Size_; ++i)
					pageCont.pop_back();
	}

    template <typename TPageFrame, unsigned FRAMES, unsigned PRIORITY_LEVELS>
	inline typename Size< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> >::Type
	length(PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> const &pageCont) 
	{
//IOREV
		return length(pageCont.pages);
	}

	template <typename TPageFrame, unsigned FRAMES, unsigned PRIORITY_LEVELS>
	inline typename Iterator< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS>, Standard >::Type
	begin(
		PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> &pageCont,
		Standard const) 
	{
//IOREV
		return begin(pageCont.pages, Standard());
	}

	template <typename TPageFrame, unsigned FRAMES, unsigned PRIORITY_LEVELS>
	inline typename Iterator< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> const, Standard >::Type
	begin(
		PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> const &pageCont,
		Standard const) 
	{
//IOREV
		return begin(pageCont.pages, Standard());
	}

	template <typename TPageFrame, unsigned FRAMES, unsigned PRIORITY_LEVELS>
	inline typename Iterator< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS>, Standard >::Type
	end(
		PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> &pageCont,
		Standard const) 
	{
//IOREV
		return end(pageCont.pages, Standard());
	}

	template <typename TPageFrame, unsigned FRAMES, unsigned PRIORITY_LEVELS>
	inline typename Iterator< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> const, Standard >::Type
	end(
		PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> const &pageCont,
		Standard const) 
	{
//IOREV
		return end(pageCont.pages, Standard());
	}




	//////////////////////////////////////////////////////////////////////////////

    template <typename TValue, typename TSize, typename T, class Function>
    inline bool equiDistantDistribution(
        Buffer<TValue> &_clusterBuffer, size_t _bufferSize, T const &me,
        TSize _size, size_t _pageSize,
        Function const &Func_)
    {
//IOREV _nodoc_
//        ::std::cerr << "equiDistantDistribution: size=" << _size << "\tpageSize=" << _pageSize << ::std::endl;
        unsigned _pages         = enclosingBlocks(_size, _pageSize);
        if (!_pages) {
			::std::cerr << "equiDistantDistribution: _pages is null!" << ::std::endl;
            return false;
        }

        if (_bufferSize < _pages) {
			::std::cerr << "equiDistantDistribution: clusterBufferSize is too small -> raised to " << _pages << ::std::endl;
            _bufferSize = _pages;
        }

        size_t lastPageSize = _size % _pageSize;
        unsigned pages      = _pages;

        if ((TSize)_bufferSize > _size)
            _bufferSize = _size;

        allocPage(_clusterBuffer, _bufferSize, me);
        PageBucketExtended<TValue> pb;
        pb.begin = _clusterBuffer.begin;

        size_t clusterSize = _bufferSize / pages;
//        ::std::cerr << "equiDistantDistribution: pages=" << _pages << "\tclusterSize=" << clusterSize << ::std::endl;
        if (lastPageSize > 0 && clusterSize >= lastPageSize) {
            // last page bucket would get more memory than page would need
            // --> exclude from equi-size distribution
            if (--pages) {
                _bufferSize -= lastPageSize;
                clusterSize = _bufferSize / pages;
            }
        }

        if (pages) {
            size_t remainder = _bufferSize % pages;
            for(unsigned i = 0, numerator = 0; i < pages; ++i) {
                pb.end = pb.begin + clusterSize;
                if ((numerator += remainder) >= pages) {    // simple bresenham for distribution
                    numerator -= pages;
                    ++pb.end;
                }
                pb.cur = pb.begin;
                pb.pageOfs = 0;
			    Func_(pb);
                pb.begin = pb.end;
            }
        }

        if (pages < _pages) {
            pb.end = pb.begin + lastPageSize;
            pb.cur = pb.begin;
            pb.pageOfs = 0;
			Func_(pb);
        }

        return true;
    }

    template <typename TValue, typename TSize, typename T, class Function>
    inline unsigned equiDistantAlignedDistribution(
        Buffer<TValue> &_clusterBuffer, size_t aligning, size_t _bufferSize, T const &me,
        TSize _size, size_t _pageSize,
        Function const &Func_)
    {
//IOREV _nodoc_
        unsigned _pages         = enclosingBlocks(_size, _pageSize);
//        ::std::cerr << "equiDistantAlignedDistribution: size=" << _size << "\tpageSize=" << _pageSize << ::std::endl;
        if (!_pages) {
			::std::cerr << "equiDistantAlignedDistribution: _pages is null!" << ::std::endl;
            return 0;
        }

        if (_bufferSize < _pages) {
			::std::cerr << "equiDistantAlignedDistribution: clusterBufferSize is too small -> raised to " << _pages << ::std::endl;
            _bufferSize = _pages;
        }

        size_t lastPageSize = _size % _pageSize;
        unsigned pages      = _pages;

        if ((TSize)_bufferSize > _size)
            _bufferSize = _size;

        size_t clusterSize = _bufferSize / pages;
        size_t aclusterSize = (clusterSize / aligning) * aligning;
        if (clusterSize - aclusterSize > aligning / 2)
            aclusterSize += aligning;

//        ::std::cerr << "equiDistantAlignedDistribution: pages=" << _pages << "\tclusterSize=" << aclusterSize << ::std::endl;

		if (aclusterSize != 0) {

			if (lastPageSize > 0 && aclusterSize > lastPageSize) {
				// last page bucket would get more memory than page would need
				// --> exclude from equi-size distribution
				--pages;
				allocPage(_clusterBuffer, aclusterSize * pages + lastPageSize, me);
			} else
				allocPage(_clusterBuffer, aclusterSize * pages, me);

			PageBucketExtended<TValue> pb;
			pb.begin = _clusterBuffer.begin;

			if (pages) {
				for(unsigned i = 0; i < pages; ++i) {
					pb.end = pb.begin + aclusterSize;
					pb.cur = pb.begin;
					pb.pageOfs = 0;
					Func_(pb);
					pb.begin = pb.end;
				}
			}

			if (pages < _pages) {
				pb.end = pb.begin + lastPageSize;
				pb.cur = pb.begin;
				pb.pageOfs = 0;
				Func_(pb);
			}
		}

        return aclusterSize;
    }

}

#endif
