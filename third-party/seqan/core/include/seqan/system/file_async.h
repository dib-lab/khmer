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

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#ifndef SEQAN_HEADER_FILE_ASYNC_H
#define SEQAN_HEADER_FILE_ASYNC_H

/* IOREV
 *
 * _tested_
 * _nodoc_
 *
 * contains third way of file-io, low level platform specific, i.e. fd for
 * POSIX and handles for WINDOWS
 *
 * not clear in which places this is used exactly
 *
 * SEQAN_DIRECTIO used on UNIX, not sure exactly what it does
 * use of other macros unclear aswell
 *
 * not sure about all the AsyncRequest stuff, as most related calls are
 * deactivated in code, see also file/file_base.h
 *
 */


namespace SEQAN_NAMESPACE_MAIN
{

 
	template <typename TSpec /* = void */>
	struct Async;


#ifdef PLATFORM_WINDOWS

	template <typename TSpec>
	class File<Async<TSpec> >
    {
//IOREV _windows_
    public:

        typedef LONGLONG    FilePtr;
        typedef ULONGLONG   SizeType;
        typedef DWORD       SizeType_;
        typedef HANDLE      Handle;

		Handle              handle, handleAsync;
        bool                noBuffering;

        File():
            handle(INVALID_HANDLE_VALUE) {}

        File(void *): // to be compatible with the FILE*(NULL) constructor
            handle(INVALID_HANDLE_VALUE) {}

        bool open(char const *fileName, int openMode = DefaultOpenMode<File>::VALUE) {
			SEQAN_PROADD(SEQAN_PROOPENFILES, 1);
            noBuffering = (getExtraFlags(openMode | OPEN_ASYNC) & (FILE_FLAG_NO_BUFFERING | FILE_FLAG_OVERLAPPED)) != 0;
            handleAsync = CreateFileA(fileName,
                                getFileAccess(openMode | OPEN_ASYNC),
                                FILE_SHARE_DELETE | FILE_SHARE_READ | FILE_SHARE_WRITE,
                                NULL,
                                getCreationFlags(openMode | OPEN_ASYNC),
                                getExtraFlags(openMode | OPEN_ASYNC),
                                NULL);

            if (handleAsync == INVALID_HANDLE_VALUE) {
				if (!(openMode & OPEN_QUIET))
					std::cerr << "Open failed on file " << fileName << ". (ErrNo=" << GetLastError() << ")" << std::endl;
                return false;
            }
            #ifdef SEQAN_VERBOSE
				if (!(openMode & OPEN_QUIET))
	                std::cerr << "file opened asynchronously " << fileName << " handle " << std::hex << handleAsync << std::dec << std::endl;
            #endif

            if (noBuffering) {
                handle = CreateFileA(fileName,                // in this case io must be sector aligned
                                getFileAccess(openMode),    // so we open a second file, for unaligned access
                                FILE_SHARE_DELETE | FILE_SHARE_READ | FILE_SHARE_WRITE,
                                NULL,
                                OPEN_EXISTING,
                                getExtraFlags(openMode & ~OPEN_ASYNC),
                                NULL);
                if (handle == INVALID_HANDLE_VALUE) {
					if (!(openMode & OPEN_QUIET))
	                	std::cerr << "Open failed on secondary file " << fileName << ". (ErrNo=" << GetLastError() << ")" << std::endl;
                    return false;
                }
	            #ifdef SEQAN_VERBOSE
					if (!(openMode & OPEN_QUIET))
	                	std::cerr << "async file opened  " << fileName << " handle " << std::hex << handle << std::dec << std::endl;
                #endif
            } else
                handle = handleAsync;

            return true;
        }

        bool openTemp(int openMode = DefaultOpenTempMode<File>::VALUE) {
            char szTempName[MAX_PATH];
#ifdef SEQAN_DEFAULT_TMPDIR
            static const char szTempPath[MAX_PATH] = SEQAN_DEFAULT_TMPDIR;
#else
            char szTempPath[MAX_PATH];
            if (!GetTempPathA(MAX_PATH, szTempPath)) {
				if (!(openMode & OPEN_QUIET))
					std::cerr << "Couldn't get a temporary path name. (ErrNo=" << GetLastError() << ")" << std::endl;
                return false;
            }
#endif
            if (!GetTempFileNameA(szTempPath, "GNDX", 0, szTempName)) {
				if (!(openMode & OPEN_QUIET))
					std::cerr << "Couldn't get a temporary file name. (ErrNo=" << GetLastError() << ")" << std::endl;
                return false;
            }
            return open(szTempName, openMode | OPEN_TEMPORARY);
        }

        inline bool close() {
            BOOL result = TRUE;
            #ifdef SEQAN_VERBOSE
                std::cerr << "files closed handles " << std::hex << handleAsync << " and " << handle << std::dec << std::endl;
            #endif
            if (handle != handleAsync)
                result &= CloseHandle(handleAsync);
            result &= CloseHandle(handle);
            handleAsync = INVALID_HANDLE_VALUE;
            handle = INVALID_HANDLE_VALUE;
			SEQAN_PROSUB(SEQAN_PROOPENFILES, 1);
            return result != FALSE;
        }

        inline DWORD read(void *memPtr, SizeType_ count) const {
            SEQAN_PROADD(SEQAN_PROIO, (count + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
            SEQAN_PROTIMESTART(tw);
			DWORD _transferedBytes;
		    ReadFile(handle, memPtr, count, &_transferedBytes, NULL);
            SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
            return _transferedBytes;
        }

        inline DWORD write(void const *memPtr, SizeType_ count) const {
            SEQAN_PROADD(SEQAN_PROIO, (count + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
            SEQAN_PROTIMESTART(tw);
			DWORD _transferedBytes;
		    WriteFile(handle, memPtr, count, &_transferedBytes, NULL);
            SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
            return _transferedBytes;
        }

		inline FilePtr seek(FilePtr _pos, DWORD origin = FILE_BEGIN) {
//          LARGE_INTEGER li = _pos;
//			return SetFilePointer(handleAsync, li.LowPart, &li.HighPart, MoveMethod);
            LARGE_INTEGER new_pos, pos;
            pos.QuadPart = _pos;
            SetFilePointerEx(handle, pos, &new_pos, origin);
//            position = new_pos.QuadPart;
            return new_pos.QuadPart;
		}

		inline FilePtr tell() {
			return seek(0, FILE_CURRENT);
        }

		inline FilePtr size() const {
            LARGE_INTEGER result;
            DWORD dwError, high;
            result.LowPart = GetFileSize(handle, &high);
            result.HighPart = high;
            if (result.LowPart == INVALID_FILE_SIZE && (dwError = GetLastError()) != NO_ERROR) {
				std::cerr << "Couldn't get file size. (ErrNo=" << dwError << ")" << std::endl;
                return 0;
            }
            return result.QuadPart;
        }

        inline bool setEof() const {
            return SetEndOfFile(handle) != FALSE;
        }

		inline static DWORD error() {
			return GetLastError();
		}

        operator bool () const {
            return (handle != INVALID_HANDLE_VALUE) && (handleAsync != INVALID_HANDLE_VALUE);
        }

    protected:

        DWORD getFileAccess(int openMode) {
            switch (openMode & OPEN_MASK) {
                case OPEN_RDONLY:
                    return GENERIC_READ;
                case OPEN_WRONLY:
                    return GENERIC_WRITE;
                case OPEN_RDWR:
                    return GENERIC_READ | GENERIC_WRITE;
				default:
					return 0;
            }
        }

        DWORD getCreationFlags(int openMode) {
            if (openMode & OPEN_CREATE)
                if (openMode & OPEN_APPEND)
                    return OPEN_ALWAYS;
                else
                    return CREATE_ALWAYS;
            else
                return OPEN_EXISTING;
        }

        DWORD getExtraFlags(int openMode) {
            DWORD extra = FILE_ATTRIBUTE_NORMAL | FILE_FLAG_RANDOM_ACCESS;// | FILE_FLAG_WRITE_THROUGH;
            if (openMode & OPEN_ASYNC) {
                extra |= FILE_FLAG_OVERLAPPED;
                #ifdef SEQAN_DIRECTIO
                    extra |= FILE_FLAG_NO_BUFFERING;
                #endif
            }
            if (openMode & OPEN_TEMPORARY)  extra |= FILE_FLAG_DELETE_ON_CLOSE;
            return extra;
        }

    };


    //////////////////////////////////////////////////////////////////////////////
    // (SeqAn adaption)
    //////////////////////////////////////////////////////////////////////////////

    struct aiocb_win32 {
        OVERLAPPED  overlapped;
        Event       xmitDone;
    };

	template <typename TSpec>
    struct AsyncRequest<File<Async<TSpec> > >
    {
        typedef aiocb_win32 Type;
    };
/*
	template <typename TSpec>
    struct aEvent<File<Async<TSpec> > >
    {
        typedef Event Type;
    };


	template <typename TSpec>
    struct aQueue<File<Async<TSpec> > >
    {
        typedef IOQueue Type;
    };

	template <typename TSpec>
    struct aHint<File<Async<TSpec> > >
    {
        typedef typename aQueue<File<Async<TSpec> > >::Type::aHint Type;
    };

	template <typename TSpec>
    struct aCallback<File<Async<TSpec> > >
    {
        typedef typename aQueue<File<Async<TSpec> > >::Type::aCallback Type;
    };*/


	template <typename TSpec>
    inline typename Size<File<Async<TSpec> > >::Type size(File<Async<TSpec> > &me) {
        return me.size();
    }

	template <typename TSpec>
    inline bool setEof(File<Async<TSpec> > &me) {
        return me.setEof();
    }

	template <typename TSpec>
    inline unsigned sectorSize(File<Async<TSpec> > const &) {
        DWORD SpC, nofC, tnoC, aligning;
        if (GetDiskFreeSpace(NULL, &SpC, &aligning, &nofC, &tnoC) == 0)  {
            std::cerr << "Error " << GetLastError() << " while querying cluster size" << std::endl;
            return 4096;
        }
        return aligning;
    }


    template < typename TSpec, typename TValue, typename TSize, typename TPos >
    inline bool asyncReadAt(File<Async<TSpec> > & me, TValue *memPtr, TSize const count, TPos const fileOfs,
        aiocb_win32 &request)
    {
        SEQAN_PROTIMESTART(tw);
        LARGE_INTEGER ofs;
        ofs.QuadPart = fileOfs;
        ofs.QuadPart *= sizeof(TValue);
        request.overlapped.Offset = ofs.LowPart;
        request.overlapped.OffsetHigh = ofs.HighPart;
        if (!request.xmitDone) open(request.xmitDone);
        request.overlapped.hEvent = request.xmitDone.hEvent;
        if (ReadFile(
            me.handleAsync, 
            memPtr, 
            count * sizeof(TValue),
            &ofs.LowPart,
            &request.overlapped) || (me.error() == ERROR_IO_PENDING))
        {
            SEQAN_PROADD(SEQAN_PROIO, (sizeof(TValue) * count + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
            SEQAN_PROADD(SEQAN_PROIWAIT, SEQAN_PROTIMEDIFF(tw));
            return true;
        }
        if (me.error() == ERROR_NO_SYSTEM_RESOURCES) {  // read synchronoulsy instead
            #if SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING
            	std::cerr << "Warning: Falling back to sync. read. :( " << std::endl;
            #endif
			signal(request.xmitDone);
            return readAt(me, memPtr, count, fileOfs);
        }
        return false;
    }
    
    template < typename TSpec, typename TValue, typename TSize, typename TPos >
    inline bool asyncWriteAt(File<Async<TSpec> > & me, TValue const *memPtr, TSize const count, TPos const fileOfs,
        aiocb_win32 &request)
    {
        SEQAN_PROTIMESTART(tw);
        LARGE_INTEGER ofs;
        ofs.QuadPart = fileOfs;
        ofs.QuadPart *= sizeof(TValue);
        request.overlapped.Offset = ofs.LowPart;
        request.overlapped.OffsetHigh = ofs.HighPart;
        if (!request.xmitDone) open(request.xmitDone);
        request.overlapped.hEvent = request.xmitDone.hEvent;
        if (WriteFile(
            me.handleAsync, 
            memPtr, 
            count * sizeof(TValue),
            &ofs.LowPart,
            &request.overlapped) || (me.error() == ERROR_IO_PENDING))
        {
            SEQAN_PROADD(SEQAN_PROIO, (sizeof(TValue) * count + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
            SEQAN_PROADD(SEQAN_PROIWAIT, SEQAN_PROTIMEDIFF(tw));
            return true;
        }
        if (me.error() == ERROR_NO_SYSTEM_RESOURCES) {  // write synchronoulsy instead
            #if SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING
            	std::cerr << "Warning: Falling back to sync. write. :( " << std::endl;
            #endif
			signal(request.xmitDone);
            return writeAt(me, memPtr, count, fileOfs);
        }
        return false;
    }

    //////////////////////////////////////////////////////////////////////
    // queue specific functions

    inline bool waitFor(aiocb_win32 &request) {
//IOREV _doc_ 
        SEQAN_PROTIMESTART(tw);
		bool inProgress;
		bool waitResult = waitFor(request.xmitDone, 60000, inProgress);
		if (inProgress)
            std::cerr << "waitFor timeout" << std::endl;
        SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
        return waitResult;
	}

    template < typename TTime >
    inline bool waitFor(aiocb_win32 &request, TTime timeoutMilliSec, bool &inProgress) {
//IOREV _doc_ 
        SEQAN_PROTIMESTART(tw);
		bool waitResult = waitFor(request.xmitDone, timeoutMilliSec, inProgress);
        SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
        return waitResult;
	}

	template < typename TSize >
	inline int waitForAny(aiocb_win32 const * const contexts[], TSize count, DWORD timeoutMilliSec = Event::Infinite) {
//IOREV _nodoc_ 
        Event::Handle *handles = new Event::Handle[count];
        for(TSize i = 0; i < count; ++i)
            handles[i] = contexts[i]->xmitDone.hEvent;

        SEQAN_PROTIMESTART(tw);
        DWORD result = WaitForMultipleObjects(count, handles, false, timeoutMilliSec);
        SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
		delete[] handles;
        if (/*result >= WAIT_OBJECT_0 && */result < WAIT_OBJECT_0 + count)
    		return result - WAIT_OBJECT_0;
        return -1;
	}

	template <typename TSpec>
    inline bool cancel(File<Async<TSpec> > & me, aiocb_win32 const &request) {
//IOREV _doc_ 
        return CancelIo(me.handleAsync);
    }

	template <typename TSpec>
    inline bool flush(File<Async<TSpec> > & me) {
//IOREV _doc_ 
		if (me.handle != me.handleAsync)	// in case of equality no direct access was done -> no flush needed
        	return FlushFileBuffers(me.handle) != 0;
        else
            return true;
    }

    template < typename TSpec, typename AsyncRequest >
    inline void release(File<Async<TSpec> > & me, AsyncRequest & request) {
//IOREV _nodoc_ 
    }


/*        
    //////////////////////////////////////////////////////////////////////
    // callback based read/write

    template < typename TSpec, typename TValue, typename TSize,
               typename aCallback, typename aHint >
    inline typename AsyncRequest<File<Async<TSpec> > >::Type
    asyncRead(File<Async<TSpec> > & me, TValue *memPtr, TSize const count,
        aCallback* cb, aHint* hint)
    {
        DWORD bsize = (DWORD)(count * sizeof(TValue));
        typename AsyncRequest<File<Async<TSpec> > >::Type request = 
            me.queue->asyncReadAt(
                me.handleAsync,
                me.position,
                memPtr,
                bsize,
                cb,
                hint);
        me.position += bsize;
        return request;
    }
    
    template < typename TSpec, typename TValue, typename TSize,
               typename aCallback, typename aHint >
    inline typename AsyncRequest<File<Async<TSpec> > >::Type
    asyncWrite(File<Async<TSpec> > & me, TValue const *memPtr, TSize const count,
        aCallback* cb, aHint* hint)
    {
        DWORD bsize = (DWORD)(count * sizeof(TValue));
        typename AsyncRequest<File<Async<TSpec> > >::Type request = 
            me.queue->asyncWriteAt(
                memPtr,
                me.handleAsync,
                me.position,
                bsize,
                cb,
                hint);
        me.position += bsize;
        return request;
    }

    template < typename TSpec, typename TValue, typename TSize, typename TPos,
               typename aCallback, typename aHint >
    inline typename AsyncRequest<File<Async<TSpec> > >::Type
    asyncReadAt(File<Async<TSpec> > & me, TValue *memPtr, TSize const count, TPos const fileOfs,
        aCallback* cb, aHint* hint)
    {
        DWORD bsize = (DWORD)(count * sizeof(TValue));
        return me.queue->asyncReadAt(
            me.handleAsync,
            fileOfs * sizeof(TValue),
            memPtr,
            bsize,
            cb,
            hint);
    }
    
    template < typename TSpec, typename TValue, typename TSize, typename TPos,
               typename aCallback, typename aHint >
    inline typename AsyncRequest<File<Async<TSpec> > >::Type
    asyncWriteAt(File<Async<TSpec> > & me, TValue *memPtr, TSize const count, TPos const fileOfs,
        aCallback* cb, aHint* hint)
    {
        DWORD bsize = (DWORD)(count * sizeof(TValue));
        return me.queue->asyncWriteAt(
            memPtr,
            me.handleAsync,
            fileOfs * sizeof(TValue),
            bsize,
            cb,
            hint);
    }


    //////////////////////////////////////////////////////////////////////
    // event based read/write

    template < typename TSpec, typename TValue, typename TSize,
               typename aEvent >
    inline typename AsyncRequest<File<Async<TSpec> > >::Type
    asyncRead(File<Async<TSpec> > & me, TValue *memPtr, TSize const count,
        aEvent &event)
    {
        DWORD bsize = (DWORD)(count * sizeof(TValue));
        typename AsyncRequest<File<Async<TSpec> > >::Type request = 
            me.queue->asyncReadAt(
                me.handleAsync,
                me.position,
                memPtr,
                bsize,
                event);
        me.position += bsize;
        return request;
    }
    
    template < typename TSpec, typename TValue, typename TSize,
               typename aEvent >
    inline typename AsyncRequest<File<Async<TSpec> > >::Type
    asyncWrite(File<Async<TSpec> > & me, TValue *memPtr, TSize const count,
        aEvent &event)
    {
        DWORD bsize = (DWORD)(count * sizeof(TValue));
        typename AsyncRequest<File<Async<TSpec> > >::Type request =  
            me.queue->asyncWriteAt(
                memPtr,
                me.handleAsync,
                me.position,
                bsize,
                event);
        me.position += bsize;
        return request;
    }

    template < typename TSpec, typename TValue, typename TSize, typename TPos,
               typename aEvent >
    inline typename AsyncRequest<File<Async<TSpec> > >::Type
    asyncReadAt(File<Async<TSpec> > & me, TValue *memPtr, TSize const count, TPos const fileOfs,
        aEvent &event)
    {
        DWORD bsize = (DWORD)(count * sizeof(TValue));
        return me.queue->asyncReadAt(
            me.handleAsync,
            fileOfs * sizeof(TValue),
            memPtr,
            bsize,
            event);
    }
    
    template < typename TSpec, TValue, typename TSize, typename TPos,
               typename aEvent >
    inline typename AsyncRequest<File<Async<TSpec> > >::Type
    asyncWriteAt(File<Async<TSpec> > & me, TValue *memPtr, TSize const count, TPos const fileOfs,
        aEvent &event)
    {
        DWORD bsize = (DWORD)(count * sizeof(TValue));
        return me.queue->asyncWriteAt(
            memPtr,
            me.handleAsync,
            fileOfs * sizeof(TValue),
            bsize,
            event);
    }


    //////////////////////////////////////////////////////////////////////
    // queue specific functions

	template <typename TSpec>
    inline void flush(File<Async<TSpec> > & me) {
        me.queue->flush();
    }

    template < typename TSpec, typename AsyncRequest >
    inline void release(File<Async<TSpec> > & me, AsyncRequest & request) {
        me.queue->release(request);
    }
*/

	//////////////////////////////////////////////////////////////////////////////
	// page aligned allocate for direct file io

    struct TagAllocateAligned_;	//< allocate page aligned memory for direct i/o access
    typedef Tag<TagAllocateAligned_> const TagAllocateAligned;

	template <typename T, typename TValue, typename TSize>
	inline void
	allocate(T const &, 
			 TValue * & data,
			 TSize count,
			 TagAllocateAligned const)
	{
//IOREV _doc_ 
		data = (TValue *) VirtualAlloc(NULL, count * sizeof(TValue), MEM_COMMIT, PAGE_READWRITE);
        if (data)
            SEQAN_PROADD(SEQAN_PROMEMORY, count * sizeof(TValue));
        else
			std::cerr << "AlignAllocator: Could not allocate memory of size " << std::hex << count * sizeof(TValue) << std::dec << ". (ErrNo=" << GetLastError() << ")" << std::endl;
	}

	//////////////////////////////////////////////////////////////////////////////
	// page aligned deallocate for direct file io

	template <typename T, typename TValue, typename TSize>
	inline void 
	deallocate( T const &,
				TValue * data, 
				TSize count,
				TagAllocateAligned const)
	{
//IOREV _doc_ 
		if (data) {
			VirtualFree(data, 0, MEM_RELEASE);
			if (count)	// .. to use count if SEQAN_PROFILE is not defined
				SEQAN_PROSUB(SEQAN_PROMEMORY, count * sizeof(TValue));
		}
	}

#else

    
	template <typename TSpec>
    class File<Async<TSpec> > : public File<Sync<TSpec> >
    {
//IOREV _nodoc_ members are not well documented
    public:

        typedef File<Sync<TSpec> >  Base;

        typedef off_t			FilePtr;
		typedef off_t           SizeType;   // type of file size
        typedef size_t          SizeType_;  // type of transfer size (for read or write)
		typedef int				Handle;

        Handle handleAsync;
		using Base::handle;

		File(void * = NULL): 	// to be compatible with the FILE*(NULL) constructor
			handleAsync(-1) {}

        virtual ~File() {}
        
        bool open(char const *fileName, int openMode = DefaultOpenMode<File>::VALUE) {
            handle = ::open(fileName, Base::_getOFlag(openMode & ~OPEN_ASYNC), S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
			if (handle == -1) 
			{
				handleAsync = handle;
				if (!(openMode & OPEN_QUIET))
					std::cerr << "Open failed on file " << fileName << ": \"" << ::strerror(errno) << '"' << std::endl;
				return false;
			}

			if (Base::_getOFlag(openMode | OPEN_ASYNC) & O_DIRECT) 
			{
				handleAsync = ::open(fileName, Base::_getOFlag(openMode | (OPEN_ASYNC & ~OPEN_CREATE)), S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
				if (handleAsync == -1 || errno == EINVAL) {	// fall back to cached access
					#if SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING
						if (!(openMode & OPEN_QUIET))
							std::cerr << "Warning: Direct access openening failed. \"" << ::strerror(errno) << '"' << std::endl;
					#endif
					handleAsync = handle;
				}
				#if SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING
				    else
						if (!(openMode & OPEN_QUIET))
							std::cerr << "Direct access successfully initiated" << std::endl;
				#endif
			} else
				handleAsync = handle;

			if (sizeof(FilePtr) < 8 && !(openMode & OPEN_QUIET))
				// To remove this warning, you have to options:
				// 1. include the following line before including anything in your application
				//    #define _FILE_OFFSET_BITS 64
				// 2. include <seqan/platform.h> or <seqan/sequence.h> before any other include
				std::cerr << "WARNING: FilePtr is not 64bit wide" << std::endl;


			SEQAN_PROADD(SEQAN_PROOPENFILES, 1);
            return true;
        }

		bool close() {
			bool result = true;
			if (handleAsync != handle && handleAsync != -1)
	            result &= (::close(handleAsync) == 0);
            result &= (::close(handle) == 0);
            handleAsync = -1;
            handle = -1;
			SEQAN_PROSUB(SEQAN_PROOPENFILES, 1);
            return result;
        }
    };


    //////////////////////////////////////////////////////////////////////////////
    // (SeqAn adaption)
    //////////////////////////////////////////////////////////////////////////////
/*
	template <typename TSpec>
    struct aQueue<File<Async<TSpec> > >
    {
        typedef void* Type;
    };
*/

	template <typename TSpec>
    struct AsyncRequest<File<Async<TSpec> > >
    {
//IOREV _doc_ 
		typedef aiocb Type;
    };
/*
	template <typename TSpec>
    struct aEvent<File<Async<TSpec> > >
    {
        typedef aiocb Type;
    };
*/

    //////////////////////////////////////////////////////////////////////
    // event based read/write

//    enum { AsyncIOSignal_ = SIGIO };

	inline void printRequest(aiocb &request, const char *_hint)
    {
//IOREV _nodoc_ _notinlined_
		std::cerr << std::hex;
		if (_hint)
			std::cerr << _hint << std::endl;
		std::cerr << "fildes:  " << request.aio_fildes << std::endl;
		std::cerr << "buffer:  " << (unsigned long)request.aio_buf << std::endl;
		std::cerr << "offset:  " << request.aio_offset<< std::endl;
		std::cerr << "nbytes:  " << request.aio_nbytes << std::endl;
		std::cerr << "event:   " << request.aio_sigevent.sigev_notify << std::endl;
		std::cerr << "Raddr:   " << &request << std::endl;
		std::cerr << std::dec;
	}

    inline void printRequest(aiocb &request)
    {
        printRequest(request, NULL);
    }

    template < typename TSpec, typename TValue, typename TSize, typename TPos >
    bool asyncReadAt(File<Async<TSpec> > & me, TValue *memPtr, TSize const count, TPos const fileOfs,
        aiocb &request)
    {
//IOREV _doc_ _notinlined_
        SEQAN_PROTIMESTART(tw);
        memset(&request, 0, sizeof(aiocb));
        request.aio_fildes = me.handleAsync;
        request.aio_buf = memPtr;
        request.aio_offset = fileOfs;
        request.aio_offset *= sizeof(TValue);
        request.aio_nbytes = count * sizeof(TValue);
        request.aio_sigevent.sigev_notify = SIGEV_NONE;
/*      request.aio_sigevent.sigev_notify = SIGEV_SIGNAL;
        request.aio_sigevent.sigev_signo = AsyncIOSignal_;
        request.aio_sigevent.sigev_value.sival_ptr = &request;
		#ifdef SEQAN_VVERBOSE
			printRequest(request, "aio_read():");
		#endif
*/      if (request.aio_nbytes == 0) return true;
		SEQAN_PROADD(SEQAN_PROIO, (request.aio_nbytes + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
		int result = aio_read(&request);
        SEQAN_PROADD(SEQAN_PROIWAIT, SEQAN_PROTIMEDIFF(tw));
        if (result != 0)
		{
			request.aio_nbytes = 0;
			if (errno == EAGAIN) {  // read synchronoulsy instead
#if SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING
            	std::cerr << "Warning: Falling back to sync. read. :( " << std::endl;
#endif
				bool success = readAt(me, memPtr, count, fileOfs);
                if (!success)
                    SEQAN_FAIL(
                        "readAt(%d, %d, %d, %d) failed: \"%s\"",
                        me.handle, (size_t)memPtr, count, fileOfs, strerror(errno));
                return success;
			}
#if SEQAN_ENABLE_DEBUG
			else
				std::cerr << "aio_read failed (asyncReadAt). \"" << ::strerror(errno) << '"' << std::endl;
#endif
        }
		return result == 0;
    }
    
    template < typename TSpec, typename TValue, typename TSize, typename TPos >
    bool asyncWriteAt(File<Async<TSpec> > & me, const TValue *memPtr, TSize const count, TPos const fileOfs,
        aiocb &request)
    {
//IOREV _doc_ _notinlined_
        SEQAN_PROTIMESTART(tw);
        memset(&request, 0, sizeof(aiocb));
        request.aio_fildes = me.handleAsync;
        request.aio_buf = const_cast<TValue*>(memPtr);
        request.aio_offset = fileOfs;
        request.aio_offset *= sizeof(TValue);
        request.aio_nbytes = count * sizeof(TValue);
        request.aio_sigevent.sigev_notify = SIGEV_NONE;
/*      request.aio_sigevent.sigev_notify = SIGEV_SIGNAL;
        request.aio_sigevent.sigev_signo = AsyncIOSignal_;
        request.aio_sigevent.sigev_value.sival_ptr = &request;
		#ifdef SEQAN_VVERBOSE
			printRequest(request, "aio_write():");
		#endif
*/      if (request.aio_nbytes == 0) return true;
		SEQAN_PROADD(SEQAN_PROIO, (request.aio_nbytes + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
		int result = aio_write(&request);
        SEQAN_PROADD(SEQAN_PROIWAIT, SEQAN_PROTIMEDIFF(tw));
        if (result != 0)
		{
			request.aio_nbytes = 0;
            int errorNo = errno;
			if (errorNo == EAGAIN)  // write synchronoulsy instead
            {
#if SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING
            	std::cerr << "Warning: Falling back to sync. write. :( " << std::endl;
#endif
				bool success = writeAt(me, memPtr, count, fileOfs);
                if (!success)
                    SEQAN_FAIL(
                        "writeAt(%d, %d, %d, %d) failed: \"%s\"",
                        me.handle, (size_t)memPtr, count, fileOfs, strerror(errno));
                return success;
			}
#if SEQAN_ENABLE_DEBUG
			else
            {
				std::cerr << "aio_write failed (asyncWriteAt): \"" << ::strerror(errno) << '"' << std::endl;
            }
#endif
        }
        return result == 0;
    }

	template <typename TSpec>
    inline bool flush(File<Async<TSpec> > & me) {
//IOREV _doc_ 
		#if _POSIX_SYNCHRONIZED_IO > 0
			return me.handle == me.handleAsync || fdatasync(me.handle) == 0;
		#else
			return me.handle == me.handleAsync || fsync(me.handle) == 0;
		#endif
    }

    //////////////////////////////////////////////////////////////////////
    // queue specific functions

	inline bool waitFor(aiocb &request)
    {
//IOREV _doc_ 
/*		#ifdef SEQAN_VVERBOSE
			printRequest(request, "aio_suspend():");
		#endif
*/
		if (request.aio_nbytes == 0) return true;
		aiocb * cblist = &request;
        SEQAN_PROTIMESTART(tw);
		int result = aio_suspend(&cblist, 1, NULL);
        ssize_t nbytes = aio_return(&request);
        SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));

//#if SEQAN_ENABLE_DEBUG
        if (result != 0 || nbytes != (ssize_t)request.aio_nbytes)
        {
            int errorNo = aio_error(&request);
            if (errorNo != EINPROGRESS)
            {
                if (errorNo != ECANCELED)
                    errorNo = errno;
				std::cerr << "Asynchronous I/O operation failed (waitFor): \"" << ::strerror(errorNo) << '"' << std::endl;
                printRequest(request);
            }
        }
//#endif

		return (result == 0) && (nbytes == (ssize_t)request.aio_nbytes);
	}

	inline bool waitFor(aiocb &request, long timeoutMilliSec, bool &inProgress)
    {
//IOREV _doc_ 
/*		#ifdef SEQAN_VVERBOSE
			printRequest(request, "aio_suspend_timeout():");
		#endif
*/
		if (request.aio_nbytes == 0)
        {
            inProgress = false;
            return true;
        }

		int result;
		if (timeoutMilliSec != 0)
        {
			aiocb * cblist = &request;
			timespec ts;
			ts.tv_sec = timeoutMilliSec / 1000;
			ts.tv_nsec = (timeoutMilliSec % 1000) * 1000;
			SEQAN_PROTIMESTART(tw);
			result = aio_suspend(&cblist, 1, &ts);
			SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
		}

        result = aio_error(&request);
        inProgress = (result == EINPROGRESS);
        ssize_t nbytes;
        
        if (inProgress)
            result = 0;
        else
        {
            nbytes = aio_return(&request);
            result = (nbytes != (ssize_t)request.aio_nbytes);
        }

        #if SEQAN_ENABLE_DEBUG
			if (result != 0)
            {
                int errorNo = aio_error(&request);
                if (errorNo != EINPROGRESS)
                {
                    if (errorNo != ECANCELED)
                        errorNo = errno;
                    std::cerr << "Asynchronous I/O operation failed (waitFor with timeOut=" << timeoutMilliSec << "ms): \"" << ::strerror(errorNo) << '"' << std::endl;
                    printRequest(request);
                }
			}
		#endif
        return result == 0;
	}

	template < typename TSize >
	inline TSize waitForAny(aiocb const * const contexts[], TSize count) {
//IOREV _nodoc_ 
        SEQAN_PROTIMESTART(tw);
		bool result = aio_suspend(contexts, count, NULL);
        SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
        return result == 0;
	}

	template < typename TSize >
	inline TSize waitForAny(aiocb const * const contexts[], TSize count, long timeoutMilliSec) {
//IOREV _nodoc_ 
        timespec ts;
        ts.tv_sec = timeoutMilliSec / 1000;
        ts.tv_nsec = (timeoutMilliSec % 1000) * 1000;
        SEQAN_PROTIMESTART(tw);
		bool result = aio_suspend(contexts, count, &ts);
        SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
        return result == 0;
	}

	template <typename TSpec>
    inline bool cancel(File<Async<TSpec> > & me, aiocb &request) {
//IOREV _doc_ 
/*		#ifdef SEQAN_VVERBOSE
			printRequest(request, "aio_cancel():");
		#endif
*/      return aio_cancel(me.handleAsync, &request) == 0;
    }

    inline int error(aiocb const & request) {
//IOREV _nodoc_ 
        return aio_error(&request);
    }

    inline int _returnValue(aiocb & request) {
//IOREV _nodoc_ 
        return aio_return(&request);
    }

	template <typename TSpec>
    inline void release(File<Async<TSpec> > & /*me*/, aiocb const & /*request*/) {
//IOREV _nodoc_ 
    }

/*
    typedef void (*sighandler_t)(int);
    static unsigned AsyncIOHandlerRefCount_ = 0;
    static struct sigaction AsyncIOOldSig_;

    inline void AsyncIOHandler_(int sigNo, siginfo_t *info, void *hint) {
        SEQAN_ASSERT_EQ(sigNo, AsyncIOSignal_);
        // TODO: signal respective event
        // currently we don't need async IO handlers because
        // we only wait for single events
    }

    static sighandler_t _addAsyncIOHandler() {
        struct sigaction newSig, oldSig;
        newSig.sa_sigaction = AsyncIOHandler_;
        sigemptyset(&newSig.sa_mask);
        newSig.sa_flags = SA_RESTART + SA_SIGINFO;
        if (sigaction(AsyncIOSignal_, &newSig, &oldSig) < 0)
            return SIG_ERR;
        return oldSig.sa_handler;
    }
*/
    
	//////////////////////////////////////////////////////////////////////////////
	// page aligned allocate for direct file io

    struct TagAllocateAligned_;	//< allocate page aligned memory for direct i/o access
    typedef Tag<TagAllocateAligned_> const TagAllocateAligned;

	template <typename T, typename TValue, typename TSize>
	inline void
	allocate(T const & /*me*/, 
			 TValue * & data,
			 TSize count,
			 TagAllocateAligned const)
	{
//IOREV _doc_ 
		data = (TValue *) ::valloc(count * sizeof(TValue));
#ifdef SEQAN_PROFILE 
        if (data)
			SEQAN_PROADD(SEQAN_PROMEMORY, count * sizeof(TValue));
		else
			std::cerr << "AlignAllocator: Could not allocate memory of size " << std::hex << 
				count * sizeof(TValue) << " with page alignment. (ErrNo=" << std::dec <<
				errno << ")" << std::endl;
#endif
	}

	//////////////////////////////////////////////////////////////////////////////
	// page aligned deallocate for direct file io

	template <typename T, typename TValue, typename TSize>
	inline void 
	deallocate( T const & /*me*/,
				TValue * data,
				TSize
#ifdef SEQAN_PROFILE 
					count
#endif
					,
				TagAllocateAligned const)
	{
//IOREV _doc_ 
#ifdef SEQAN_PROFILE 
        if (data && count)	// .. to use count if SEQAN_PROFILE is not defined
			SEQAN_PROSUB(SEQAN_PROMEMORY, count * sizeof(TValue));
#endif
		::free(data);
	}

    template < typename TSpec, typename TSize >
    inline void resize(File<Async<TSpec> > &me, TSize new_length) {
//IOREV _doc_ 
		if (!me.resize(new_length))
            SEQAN_FAIL(
                "resize(%d, %d) failed: \"%s\"",
                me.handle, new_length, strerror(errno));
    }
	

#endif

    //////////////////////////////////////////////////////////////////////////////
    // global functions

	template <typename TSpec>
    struct Size< File<Async<TSpec> > >
    {
//IOREV
        typedef typename File<Async<TSpec> >::SizeType Type;
    };

	template <typename TSpec>
    struct Position< File<Async<TSpec> > >
    {
//IOREV
        typedef typename File<Async<TSpec> >::FilePtr Type;
    };

	template <typename TSpec>
    struct Difference< File<Async<TSpec> > >
    {
//IOREV
        typedef typename File<Async<TSpec> >::FilePtr Type;
    };



    template < typename TSpec, typename TValue, typename TSize>
	inline void
	allocate( File<Async<TSpec> > const & me,
			  TValue * & data, 
			  TSize count)
	{
//IOREV _doc_ 
		allocate(me, data, count, TagAllocateAligned());
	}

    template <typename TSpec, typename TValue, typename TSize>
	inline void
	deallocate( File<Async<TSpec> > const & me,
				TValue * data, 
				TSize count)
	{
//IOREV _doc_ 
		deallocate(me, data, count, TagAllocateAligned());
	}


}

#endif
