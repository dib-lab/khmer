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

#ifndef SEQAN_HEADER_FILE_SIMPLE_H
#define SEQAN_HEADER_FILE_SIMPLE_H

#include <fcntl.h>          // O_CREAT ..
#include <sys/stat.h>       // 
#include <cstdio>           // tmpnam(..)

#ifdef PLATFORM_WINDOWS
# include <io.h>            // read(..) ..
#else
# include <cstdlib>
# include <cerrno>
# include <unistd.h>
#endif


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
 */



//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{


	template <typename TSpec /* = void */>
	struct Sync;


#ifdef PLATFORM_WINDOWS

    //////////////////////////////////////////////////////////////////////////////
    // Windows rtl file access
	template <typename TSpec>
	class File<Sync<TSpec> >
    {
//IOREV _windows_ _nodoc_
    public:

		typedef __int64			FilePtr;
		typedef __int64         SizeType;   // type of file size
        typedef unsigned int    SizeType_;  // type of transfer size (for read or write)
		typedef int				Handle;

        Handle handle;

        File(void * /*dummy*/ = NULL): // to be compatible with the FILE*(NULL) constructor
            handle(-1) {}

	//File(int posixHandle) : handle(posixHandle) {}

        inline int _getOFlag(int openMode) const 
		{
			int result;
			bool canWrite = false;

			switch (openMode & OPEN_MASK) {
                case OPEN_RDONLY:
                    result = _O_RDONLY;
					break;
                case OPEN_WRONLY:
					canWrite = true;
                    result = _O_WRONLY;
					break;
                case OPEN_RDWR:
				default:
					canWrite = true;
                    result = _O_RDWR;
					break;
			}

			if (openMode & OPEN_CREATE)     result |= _O_CREAT;
			if (canWrite && !(openMode & OPEN_APPEND))	result |= _O_TRUNC;
            if (openMode & OPEN_TEMPORARY)  result |= _O_TEMPORARY;
			return result | _O_BINARY;
        }

        bool open(char const *fileName, int openMode = DefaultOpenMode<File>::VALUE) 
		{
            handle = _open(fileName, _getOFlag(openMode), _S_IREAD | _S_IWRITE);
			if (handle == -1) {
				if (!(openMode & OPEN_QUIET))
					::std::cerr << "Open failed on file " << fileName << ". (" << ::strerror(errno) << ")" << ::std::endl;
				return false;
			}
			SEQAN_PROADD(SEQAN_PROOPENFILES, 1);
            return true;
        }

        bool openTemp(int openMode = DefaultOpenTempMode<File>::VALUE) 
		{
#ifdef SEQAN_DEFAULT_TMPDIR
			char *fileName = _tempnam(SEQAN_DEFAULT_TMPDIR, "SQN");
#else
			char *fileName = _tempnam(NULL, "SQN");
#endif
			if (!fileName) {
				if (!(openMode & OPEN_QUIET))
					::std::cerr << "Cannot create a unique temporary filename" << ::std::endl;
				return false;
			}
            bool result = open(fileName, openMode | OPEN_TEMPORARY);
			free(fileName);
			return result;
        }

        inline bool close() 
		{
            if (_close(handle) != 0)
                return false;
            handle = -1;
			SEQAN_PROSUB(SEQAN_PROOPENFILES, 1);
            return true;
        }

		inline int read(void *buffer, SizeType_ count) const 
		{
            SEQAN_PROADD(SEQAN_PROIO, (count + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
            SEQAN_PROTIMESTART(tw);
		    int result = _read(handle, buffer, count);
            SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
            return result;
		}

		inline int write(void const *buffer, SizeType_ count) const 
		{
            SEQAN_PROADD(SEQAN_PROIO, (count + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
            SEQAN_PROTIMESTART(tw);
		    int result = _write(handle, buffer, count);
            SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
            return result;
		}

		inline FilePtr seek(FilePtr pos, int origin = SEEK_SET) const 
		{
			return _lseeki64(handle, pos, origin);
		}

		inline FilePtr tell() const 
		{
			return _telli64(handle);
		}

		static int error() 
		{
			return errno;
		}

        operator bool () const 
		{
            return handle != -1;
        }
    };

	inline bool fileExists(const char *fileName) 
	{
//IOREV _windows_ _nodoc_
		struct _stat buf;
		return _stat(fileName, &buf) == 0;
	}

	inline bool fileUnlink(const char *fileName) 
	{
//IOREV _windows_ _nodoc_
		return _unlink(fileName) == 0;
	}

#else

    //////////////////////////////////////////////////////////////////////////////
    // Unix file access
	template <typename TSpec>
	class File<Sync<TSpec> >
    {
//IOREV __nodoc_
    public:

		typedef off_t			FilePtr;
		typedef off_t           SizeType;   // type of file size
        typedef size_t          SizeType_;  // type of transfer size (for read or write)
		typedef int				Handle;

        Handle handle;

        File(void * /*dummy*/ = NULL): // to be compatible with the FILE*(NULL) constructor
            handle(-1) {}

	///File(int posixHandle) : handle(posixHandle) {}

        virtual ~File() {}
        
        inline int _getOFlag(int openMode) const {
			int result = O_LARGEFILE;

			switch (openMode & OPEN_MASK) {
                case OPEN_RDONLY:
                    result |= O_RDONLY;
					break;
                case OPEN_WRONLY:
                    result |= O_WRONLY;
					if (!(openMode & OPEN_APPEND))	result |= O_TRUNC;
					break;
                case OPEN_RDWR:
                    result |= O_RDWR;
					if (!(openMode & OPEN_APPEND))	result |= O_TRUNC;
					break;
			}

			if (openMode & OPEN_CREATE)     result |= O_CREAT;
//			if (openMode & OPEN_TEMPORARY)  result |= O_TEMPORARY;
        #ifdef SEQAN_DIRECTIO
    		if (openMode & OPEN_ASYNC)		result |= O_DIRECT;
        #endif
			return result;
        }

        virtual bool open(char const *fileName, int openMode = DefaultOpenMode<File>::VALUE) {
            handle = ::open(fileName, _getOFlag(openMode), S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
			if (handle == -1 && errno == EINVAL) {	// fall back to cached access
	            #ifdef SEQAN_DEBUG_OR_TEST_
					if (!(openMode & OPEN_QUIET))
						::std::cerr << "Warning: Direct access openening failed: " << fileName << "." << ::std::endl;
				#endif			
          	    handle = ::open(fileName, _getOFlag(openMode & ~OPEN_ASYNC), S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
			}
			
			if (handle == -1) {
				if (!(openMode & OPEN_QUIET))
					::std::cerr << "Open failed on file " << fileName << ". (" << ::strerror(errno) << ")" << ::std::endl;
				return false;
			}

			if (sizeof(FilePtr) < 8 && !(openMode & OPEN_QUIET))
				// To remove this warning, you have to options:
				// 1. include the following line before including anything in your application
				//    #define _FILE_OFFSET_BITS 64
				// 2. include <seqan/platform.h> or <seqan/sequence.h> before any other include
				::std::cerr << "WARNING: FilePtr is not 64bit wide" << ::std::endl;

			SEQAN_PROADD(SEQAN_PROOPENFILES, 1);
            return true;
        }

        bool openTemp(int openMode = DefaultOpenTempMode<File>::VALUE) {
            // Construct the pattern for the temporary file.
            //
            // First, try to get the temporary directory from the environment
            // variables TMPDIR, TMP.
            CharString tmpDir;
            if ((getuid() == geteuid()) && (getgid() == getegid())) 
			{
                char * res;
                if ((res = getenv("TMPDIR")) != NULL)
                    tmpDir = res;
                else
                    if ((res = getenv("TMP")) != NULL)
                        tmpDir = res;
            }
            // If this does not work, try to use the constant
            // SEQAN_DEFAULT_TMPDIR, fall back to "/tmp", if this does not
            // work.
#ifdef SEQAN_DEFAULT_TMPDIR
            if (empty(tmpDir))
                tmpDir = SEQAN_DEFAULT_TMPDIR;
#else  // #ifdef SEQAN_DEFAULT_TMPDIR
            if (empty(tmpDir))
                tmpDir = "/tmp";
#endif  // #ifdef SEQAN_DEFAULT_TMPDIR

            // At this point, we have a temporary directory.  Now, we add the
            // file name template to get the full path template.
            append(tmpDir, "/SQNXXXXXX");
            // Open temporary file and unlink it immediately afterwards so the
            // memory is released when the program exits.
            int oldMode = umask(077);  // Create with restrictive permissions.
			if ((handle = ::mkstemp(toCString(tmpDir))) == -1) {
			    umask(oldMode);  // Reset umask mode.
				if (!(openMode & OPEN_QUIET))
					::std::cerr << "Couldn't create temporary file " << tmpDir << ". (" << ::strerror(errno) << ")" << ::std::endl;
				return false;
			}
			if (!(close() && open(toCString(tmpDir), openMode))) 
			{
				umask(oldMode);  // Reset umask mode.
			    return false;
            }
			umask(oldMode);  // Reset umask mode.
            #ifdef SEQAN_DEBUG
				if (::unlink(toCString(tmpDir)) == -1 && !(openMode & OPEN_QUIET))
					::std::cerr << "Couldn't unlink temporary file " << tmpDir << ". (" << ::strerror(errno) << ")" << ::std::endl;
            #else
				::unlink(toCString(tmpDir));
			#endif
			return true;
        }


        virtual bool close() {
            if (::close(this->handle) == -1) return false;
            handle = -1;
			SEQAN_PROSUB(SEQAN_PROOPENFILES, 1);
            return true;
        }

		inline ssize_t read(void *buffer, SizeType_ count) const {
            SEQAN_PROADD(SEQAN_PROIO, (count + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
            SEQAN_PROTIMESTART(tw);
		    ssize_t result = ::read(handle, buffer, count);
            SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
            return result;
		}

		inline ssize_t write(void const *buffer, SizeType_ count) const {
            SEQAN_PROADD(SEQAN_PROIO, (count + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
            SEQAN_PROTIMESTART(tw);
		    ssize_t result = ::write(handle, buffer, count);
            SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
            return result;
		}

		inline FilePtr seek(FilePtr pos, int origin = SEEK_SET) const {
            FilePtr result = ::lseek(handle, pos, origin);
//			#ifdef SEQAN_DEBUG
				if (result < 0)
					::std::cerr << "lseek returned " << result << ". (" << ::strerror(errno) << ")" << ::std::endl;
//			#endif
			return result;
		}

		inline FilePtr tell() const {
            return seek(0, SEEK_CUR);
        }

		inline bool resize(SizeType new_length) const {
			return ftruncate(handle, new_length) == 0;
		}

		static int error() {
            return errno;
		}

        operator bool () const {
            return handle != -1;
        }
    };

	inline bool fileExists(const char *fileName)
	{
//IOREV _nodoc_
		struct stat buf;
		return stat(fileName, &buf) != -1;
	}

	inline bool fileUnlink(const char *fileName)
	{
//IOREV _noddoc_
		return unlink(fileName) == 0;
	}

    template < typename TSpec, typename TSize >
    inline void resize(File<Sync<TSpec> > &me, TSize new_length)
	{
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
    struct Size< File<Sync<TSpec> > >
    {
//IOREV
        typedef typename File<Sync<TSpec> >::SizeType Type;
    };

	template <typename TSpec>
    struct Position< File<Sync<TSpec> > >
    {
//IOREV
        typedef typename File<Sync<TSpec> >::FilePtr Type;
    };

	template <typename TSpec>
    struct Difference< File<Sync<TSpec> > >
    {
//IOREV
        typedef typename File<Sync<TSpec> >::FilePtr Type;
    };

    template < typename TSpec, typename TValue, typename TSize >
    inline bool read(File<Sync<TSpec> > & me, TValue *memPtr, TSize const count) {
//IOREV
		return (int) me.read(memPtr, count * sizeof(TValue)) == (int) (count * sizeof(TValue));
    }
    
    template < typename TSpec, typename TValue, typename TSize >
    inline bool write(File<Sync<TSpec> > & me, TValue const *memPtr, TSize const count) {
//IOREV
		return (int) me.write(memPtr, count * sizeof(TValue)) == (int) (count * sizeof(TValue));
    }

}

#endif
