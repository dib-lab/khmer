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

#ifndef SEQAN_HEADER_CSTREAM_H
#define SEQAN_HEADER_CSTREAM_H
 
#include <cstdio>

/* IOREV
 * _tested_
 * _nodoc_
 * 
 * 
 * Tested by tests/file
 * documentation non-existent
 * relation to file_cstyle.h not clear
 * Metafunctions supposedly moved to file_cstyle.h are commented there too
 */



namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
/**
.Adaption."std::FILE *":
..summary:Standard library C style streams.
*/

//////////////////////////////////////////////////////////////////////////////
// Position is now defined in file/file_cstyle.h
/*
template <>
struct Position<FILE *>
{
	typedef long Type;
};
*/
//////////////////////////////////////////////////////////////////////////////

template <>
struct Value<FILE *>
{
//IOREV
	typedef char Type;
};

//////////////////////////////////////////////////////////////////////////////
/*
template <>
struct Position<FILE *>
{
	typedef ::std::fpos_t Type;
};
*/

//////////////////////////////////////////////////////////////////////////////

template <typename T>
struct IsTellAndSeekStream_;
//IOREV

template <>
struct IsTellAndSeekStream_<FILE *>
{
//IOREV
	typedef True Type;
};

//////////////////////////////////////////////////////////////////////////////

inline bool 
_streamOpen(::std::FILE * & me, String<char> path, bool for_read = true)
{
//IOREV _duplicate_ of file_cstyle's  "open"
SEQAN_CHECKPOINT
	if (for_read)
	{
		me = fopen(toCString(path), "rb");
	}
	else
	{
		me = fopen(toCString(path), "wb");
	}
	return (me != 0);
}


//////////////////////////////////////////////////////////////////////////////

inline void 
_streamClose(::std::FILE * & me)
{
//IOREV _duplicate_ of file_cstyle's  "close"
SEQAN_CHECKPOINT
	if (me)
	{
		fclose(me);
		me = 0;
	}
}

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamEOF.param.stream.type:Adaption."std::FILE *"

inline bool 
_streamEOF(::std::FILE * me)
{
//IOREV
SEQAN_CHECKPOINT
	int c = fgetc(me);
    ungetc(c, me);
	return (c == EOF) || ferror(me);
}

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamRead.param.stream.type:Adaption."std::FILE *"

template <typename TValue>
inline size_t 
_streamRead(TValue * target,
			::std::FILE * source,
			size_t limit)
{
//IOREV _duplicate_ of file_cstyle.h's  "read"
SEQAN_CHECKPOINT
	return ::std::fread(target, sizeof(TValue), limit, source);
}

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamGet.param.stream.type:Adaption."std::FILE *"

inline char 
_streamGet(::std::FILE * source)
{
//IOREV
SEQAN_CHECKPOINT
	return getc(source);
}

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamPut.param.stream.type:Adaption."std::FILE *"

inline void
_streamPut(::std::FILE * target,
		   char character)
{
//IOREV
SEQAN_CHECKPOINT
	putc(character, target);
}


//////////////////////////////////////////////////////////////////////////////

///.Internal._streamPut.param.stream.type:Adaption."std::FILE *"

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamTellG.param.stream.type:Adaption."std::FILE *"

inline Position<FILE *>::Type
_streamTellG(FILE * me)
{
//IOREV _duplicate_ overlaps in function with file_cstyle.h's  "tell"
SEQAN_CHECKPOINT
	return ::std::ftell(me);
}

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamTellP.param.stream.type:Adaption."std::FILE *"

inline Position<FILE *>::Type
_streamTellP(FILE * me)
{
//IOREV _duplicate_ overlaps in function with file_cstyle.h's  "tell"
SEQAN_CHECKPOINT
	return ::std::ftell(me);
}

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamSeekG.param.stream.type:Adaption."std::FILE *"

inline void
_streamSeekG(FILE * me,
			 Position<FILE *>::Type pos)
{
//IOREV _duplicate_ overlaps in function with file_cstyle.h's  "seek"
SEQAN_CHECKPOINT
	::std::fseek(me, pos, SEEK_SET);
}

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamSeekP.param.stream.type:Adaption."std::FILE *"

inline void
_streamSeekP(FILE * me,
			 Position<FILE *>::Type pos)
{
//IOREV _duplicate_ overlaps in function with file_cstyle.h's  "seek"
SEQAN_CHECKPOINT
	::std::fseek(me, pos, SEEK_SET);
}

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamSeek2G.param.stream.type:Adaption."std::FILE *"

inline void
_streamSeek2G(FILE * me,
	 int off)
{
//IOREV _duplicate_ overlaps in function with file_cstyle.h's  "seek"
SEQAN_CHECKPOINT
	::std::fseek(me, off, SEEK_CUR);
}

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamUnget.param.stream.type:Adaption."std::FILE *"

inline void
_streamUnget(::std::FILE * stream)
{
//IOREV _duplicate_ overlaps in function with file_cstyle.h's  "seek"
SEQAN_CHECKPOINT
	_streamSeek2G(stream, -1);
}


//////////////////////////////////////////////////////////////////////////////
// holder<FILE *>

template <typename THolder>
inline void
_holderDeallocate(THolder &, FILE *)
{
//IOREV whats this good for?
}
template <typename THolder>
inline FILE *
_holderAllocatePointer(THolder &, FILE * data)
{
//IOREV whats this good for?
	return data;
}

//////////////////////////////////////////////////////////////////////////////
// Stream operators for FILE *
//////////////////////////////////////////////////////////////////////////////


// ISO C++ operators are only allowed for classes, not for pointers

/*
template <typename TSource>
inline FILE *
operator << (FILE * target, 
			 TSource & source)
{
SEQAN_CHECKPOINT
	write(target, source);
	return target;
}
template <typename TSource>
inline FILE *
operator << (FILE * target, 
			 TSource const & source)
{
SEQAN_CHECKPOINT
	write(target, source);
	return target;
}

//____________________________________________________________________________

template <typename TTarget>
inline FILE *
operator >> (FILE * source, 
			 TTarget & target)
{
SEQAN_CHECKPOINT
	read(source, target);
	return source;
}
template <typename TTarget>
inline FILE *
operator >> (FILE * source, 
			 TTarget const & target)
{
SEQAN_CHECKPOINT
	read(source, target);
	return source;
}
*/

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
