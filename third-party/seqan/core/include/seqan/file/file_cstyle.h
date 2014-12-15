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
// Adaption of the FILE object from C.
// ==========================================================================

#ifndef SEQAN_HEADER_FILE_CSTYLE_H
#define SEQAN_HEADER_FILE_CSTYLE_H

/* IOREV
 * _tested_
 * _nodoc_
 * 
 * 
 * Tested by tests/file
 * documentation non-existent
 * relation to cstream.h not clear
 * Metafunctions supposedly moved hereto from file_cstyle.h are commented here awell
 */

//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

/*    template <>
    struct Value< FILE* >
    {
	    typedef unsigned char Type;
    };
*/
/*  already defined in "seqan/cstream.h" ...

    template <>
    struct Position< FILE* >
    {
	    typedef long Type;
    };
*/
/*
    template <>
    struct Size< FILE* >
    {
	    typedef size_t Type;
    };
*/
/*

    template <>
    struct Position< FILE* >
    {
	    typedef long Type;
    };
*/
    template <>
    struct Difference< FILE* >
    {
//IOREV shouldnt this be ulong, as the file can be ulong bytes big?
	    typedef long Type;
    };


	inline const char * 
	_getCStyleOpenMode(int openMode) 
	{
//IOREV double check whether this translates FileOpenMode correctly (doesnt look like it)
		switch (openMode & OPEN_MASK) {
            case OPEN_WRONLY:
                if (!(openMode & OPEN_APPEND))
                    if (openMode & OPEN_CREATE)
                        return "w";
                    else
                        return "r+";
                else
                    return "a";
            case OPEN_RDWR:
                if (!(openMode & OPEN_APPEND))
                    if (openMode & OPEN_CREATE)
                        return "w+";
                    else
                        return "r+";
                else
                    return "a+";
			default:
		        return "r";
		}
    }

    inline bool 
	open(FILE* &me, const char *fileName, int openMode) 
	{
//IOREV _duplicate_ of cstream.h's  "open"
		SEQAN_PROADD(SEQAN_PROOPENFILES, 1);
        return (me = fopen(fileName, _getCStyleOpenMode(openMode))) != NULL;
    }

    inline bool 
	open(FILE* &me, const char *fileName) 
	{
//IOREV _duplicate_ of cstream.h's  "open"
		return open(me, fileName, DefaultOpenMode<FILE*>::VALUE);
	}

    inline bool 
	openTemp(FILE* &me) 
	{
//IOREV
		SEQAN_PROSUB(SEQAN_PROOPENFILES, 1);
        return (me = tmpfile()) != NULL;
    }

    inline bool 
	close(FILE* me) 
	{
//IOREV _duplicate_ of cstream.h's  "close"
		SEQAN_PROSUB(SEQAN_PROOPENFILES, 1);
        return fclose(me) == 0;
    }

    inline unsigned 
	sectorSize(FILE* const &) 
	{
//IOREV _duplicate_ _nodoc_ duplicate or identical spec. in file_base.h should'nt this be variable
        return 4096;
    }

    template < typename TPos >
    inline Size<FILE*>::Type 
	seek(FILE* me, TPos const fileOfs, int origin) 
	{
//IOREV _duplicate_ overlaps in function with multiple function in cstream.h
        fseek(me, fileOfs, origin);
		return ftell(me);
    }
    template < typename TPos >
    inline Size<FILE*>::Type 
	seek(FILE* me, TPos const fileOfs) 
	{
//IOREV shouldnt it be SEEK_SET instead of SEEK_BEGIN?
		return seek(me, fileOfs, SEEK_BEGIN);
    }

    inline Size<FILE*>::Type 
	tell(FILE* me) 
	{
//IOREV _duplicate_ overlaps in function with multiple functions in cstream.h
		return ftell(me);
    }

    template < typename TValue, typename TSize >
    inline bool 
	read(FILE* me, TValue *memPtr, TSize const count) 
	{
//IOREV _duplicate_ of cstream.h's  "read"
        SEQAN_PROADD(SEQAN_PROIO, (sizeof(TValue) * count + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
        SEQAN_PROTIMESTART(tw);
        bool result = fread(memPtr, sizeof(TValue), count, me) == (size_t)count;
        SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
        return result;
    }

    template < typename TValue, typename TSize >
    inline bool 
	write(FILE* me, TValue const *memPtr, TSize const count) 
	{
//IOREV
        SEQAN_PROADD(SEQAN_PROIO, (sizeof(TValue) * count + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
        SEQAN_PROTIMESTART(tw);
        bool result = fwrite(memPtr, sizeof(TValue), count, me) == (size_t)count;
        SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
        return result;
    }

    template < typename TValue, typename TSize, typename TPos >
    inline bool 
	readAt(FILE* me, TValue *memPtr, TSize const count, TPos const fileOfs) 
	{
//IOREV
		typedef typename Position<FILE*>::Type pos_t;
		seek(me, (pos_t)fileOfs * (pos_t)sizeof(TValue));
        SEQAN_PROADD(SEQAN_PROIO, (sizeof(TValue) * count + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
        SEQAN_PROTIMESTART(tw);
        bool result = fread(memPtr, sizeof(TValue), count, me) == (size_t)count;
        SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
        return result;
    }
    
    template < typename TValue, typename TSize, typename TPos >
    inline bool 
	writeAt(FILE* me, TValue const *memPtr, TSize const count, TPos const fileOfs) 
	{
//IOREV
		typedef typename Position<FILE*>::Type pos_t;
		seek(me, (pos_t)fileOfs * (pos_t)sizeof(TValue));
        SEQAN_PROADD(SEQAN_PROIO, (sizeof(TValue) * count + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
        SEQAN_PROTIMESTART(tw);
        bool result = fwrite(memPtr, sizeof(TValue), count, me) == (size_t)count;
        SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
        return result;
    }

    inline Size<FILE*>::Type 
	size(FILE* me) 
	{
//IOREV
        Size<FILE*>::Type old_pos = tell(me);
        Size<FILE*>::Type result = 0;
        if (seek(me, 0, SEEK_END) == 0)
            result = tell(me);
        seek(me, old_pos, SEEK_BEGIN);
        return result;
    }

    template < typename TSize >
    inline void 
	resize(FILE* me, TSize new_length) 
	{
//IOREV I think this is not c-standard and only works on POSIX-compliant systems
        Size<FILE*>::Type old_pos = tell(me);
        seek(me, new_length, SEEK_BEGIN);
        seek(me, old_pos, SEEK_BEGIN);
    }

	inline bool 
	flush(FILE*) 
	{
//IOREV _notimplemented_
		return true; 
	}

    template < typename AsyncRequest >
	inline void 
	release(FILE*, AsyncRequest &) 
	{
//IOREV _notimplemented_
	}

    template < typename AsyncRequest >
    inline bool 
	cancel(FILE*, AsyncRequest &) 
	{
//IOREV _notimplemented_
		return true; 
	}

}

#endif
