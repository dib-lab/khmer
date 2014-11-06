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

#ifndef SEQAN_HEADER_STREAM_H
#define SEQAN_HEADER_STREAM_H

#include <iosfwd>

/* IOREV
 * _doc_
 * _tested_
 * 
 * well documented
 * tested by tests/file and probably used by a lot of stuff
 * c++ equivalent to  cstream.h
 */


namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
/**
.Adaption."std::iostream":
..summary:Standard library stream classes.
*/

//////////////////////////////////////////////////////////////////////////////
	
template <typename TValue, typename TTraits>
struct Position< ::std::basic_ios<TValue, TTraits> >
{
//IOREV
	typedef typename ::std::basic_ios<TValue, TTraits>::pos_type Type;
};
template <typename TValue, typename TTraits>
struct Position< ::std::basic_streambuf<TValue, TTraits> >
{
//IOREV
	typedef typename ::std::basic_streambuf<TValue, TTraits>::pos_type Type;
};
template <typename TValue, typename TTraits>
struct Position< ::std::basic_istream<TValue, TTraits> >
{
//IOREV
	typedef typename ::std::basic_istream<TValue, TTraits>::pos_type Type;
};
template <typename TValue, typename TTraits>
struct Position< ::std::basic_ostream<TValue, TTraits> >
{
//IOREV
	typedef typename ::std::basic_ostream<TValue, TTraits>::pos_type Type;
};
template <typename TValue, typename TTraits>
struct Position< ::std::basic_iostream<TValue, TTraits> >
{
//IOREV
	typedef typename ::std::basic_iostream<TValue, TTraits>::pos_type Type;
};
template <typename TValue, typename TTraits>
struct Position< ::std::basic_stringbuf<TValue, TTraits> >
{
//IOREV
	typedef typename ::std::basic_stringbuf<TValue, TTraits>::pos_type Type;
};
template <typename TValue, typename TTraits>
struct Position< ::std::basic_istringstream<TValue, TTraits> >
{
//IOREV
	typedef typename ::std::basic_istringstream<TValue, TTraits>::pos_type Type;
};
template <typename TValue, typename TTraits>
struct Position< ::std::basic_ostringstream<TValue, TTraits> >
{
//IOREV
	typedef typename ::std::basic_ostringstream<TValue, TTraits>::pos_type Type;
};
template <typename TValue, typename TTraits>
struct Position< ::std::basic_stringstream<TValue, TTraits> >
{
//IOREV
	typedef typename ::std::basic_stringstream<TValue, TTraits>::pos_type Type;
};
template <typename TValue, typename TTraits>
struct Position< ::std::basic_filebuf<TValue, TTraits> >
{
//IOREV
	typedef typename ::std::basic_filebuf<TValue, TTraits>::pos_type Type;
};
template <typename TValue, typename TTraits>
struct Position< ::std::basic_ifstream<TValue, TTraits> >
{
//IOREV
	typedef typename ::std::basic_ifstream<TValue, TTraits>::pos_type Type;
};
template <typename TValue, typename TTraits>
struct Position< ::std::basic_ofstream<TValue, TTraits> >
{
//IOREV
	typedef typename ::std::basic_ofstream<TValue, TTraits>::pos_type Type;
};
template <typename TValue, typename TTraits>
struct Position< ::std::basic_fstream<TValue, TTraits> >
{
//IOREV
	typedef typename ::std::basic_fstream<TValue, TTraits>::pos_type Type;
};

//////////////////////////////////////////////////////////////////////////////
	
template <typename TValue, typename TTraits>
struct Value< ::std::basic_ios<TValue, TTraits> >
{
//IOREV
	typedef typename ::std::basic_ios<TValue, TTraits>::char_type Type;
};
template <typename TValue, typename TTraits>
struct Value< ::std::basic_streambuf<TValue, TTraits> >
{
//IOREV
	typedef typename ::std::basic_streambuf<TValue, TTraits>::char_type Type;
};
template <typename TValue, typename TTraits>
struct Value< ::std::basic_istream<TValue, TTraits> >
{
//IOREV
	typedef typename ::std::basic_istream<TValue, TTraits>::char_type Type;
};
template <typename TValue, typename TTraits>
struct Value< ::std::basic_ostream<TValue, TTraits> >
{
//IOREV
	typedef typename ::std::basic_ostream<TValue, TTraits>::char_type Type;
};
template <typename TValue, typename TTraits>
struct Value< ::std::basic_iostream<TValue, TTraits> >
{
//IOREV
	typedef typename ::std::basic_iostream<TValue, TTraits>::char_type Type;
};
template <typename TValue, typename TTraits>
struct Value< ::std::basic_stringbuf<TValue, TTraits> >
{
//IOREV
	typedef typename ::std::basic_stringbuf<TValue, TTraits>::char_type Type;
};
template <typename TValue, typename TTraits>
struct Value< ::std::basic_istringstream<TValue, TTraits> >
{
//IOREV
	typedef typename ::std::basic_istringstream<TValue, TTraits>::char_type Type;
};
template <typename TValue, typename TTraits>
struct Value< ::std::basic_ostringstream<TValue, TTraits> >
{
//IOREV
	typedef typename ::std::basic_ostringstream<TValue, TTraits>::char_type Type;
};
template <typename TValue, typename TTraits>
struct Value< ::std::basic_stringstream<TValue, TTraits> >
{
//IOREV
	typedef typename ::std::basic_stringstream<TValue, TTraits>::char_type Type;
};
template <typename TValue, typename TTraits>
struct Value< ::std::basic_filebuf<TValue, TTraits> >
{
//IOREV
	typedef typename ::std::basic_filebuf<TValue, TTraits>::char_type Type;
};
template <typename TValue, typename TTraits>
struct Value< ::std::basic_ifstream<TValue, TTraits> >
{
//IOREV
	typedef typename ::std::basic_ifstream<TValue, TTraits>::char_type Type;
};
template <typename TValue, typename TTraits>
struct Value< ::std::basic_ofstream<TValue, TTraits> >
{
//IOREV
	typedef typename ::std::basic_ofstream<TValue, TTraits>::char_type Type;
};
template <typename TValue, typename TTraits>
struct Value< ::std::basic_fstream<TValue, TTraits> >
{
//IOREV
	typedef typename ::std::basic_fstream<TValue, TTraits>::char_type Type;
};

//////////////////////////////////////////////////////////////////////////////

/**.interal.IsTellAndSeekStream_:
..summary:Determines whether stream supports tell and seek functions.
..cat:Metafunction
*/

template <typename T>
struct IsTellAndSeekStream_
{
//IOREV can't I tell and seek e.g. in an ofstream?
	typedef False Type;
};


template <typename TValue, typename TTraits>
struct IsTellAndSeekStream_< ::std::basic_ifstream<TValue, TTraits> >
{
//IOREV
	typedef True Type;
};
template <typename TValue, typename TTraits>
struct IsTellAndSeekStream_< ::std::basic_fstream<TValue, TTraits> >
{
//IOREV
	typedef True Type;
};

//////////////////////////////////////////////////////////////////////////////

/**
.Internal._streamEOF:
..summary:Test stream for being in eof or error state.
..cat:Streams
..signature:_streamEOF(stream)
..param.stream:A stream object.
...type:Adaption."std::iostream"
..returns:$true$, if stream is at end of file or was set to error state, $false$ otherwise.
*/
template <typename TValue, typename TTraits>
inline bool 
_streamEOF(::std::basic_istream<TValue, TTraits> const & me)
{
//IOREV
SEQAN_CHECKPOINT
	// Andreas missed the fact that eof() of a stream is true after reading the eof character
	// So reading the last character eof() is false, reading beyond eof() is true
	// To fix that we use peek() to get the next character and compare it with the eof character
	typedef ::std::basic_istream<TValue, TTraits> TStream;
	return me.fail() || const_cast<TStream &>(me).peek() == TTraits::eof();
}

//////////////////////////////////////////////////////////////////////////////
 
/**
.Internal._streamRead:
..summary:Read some characters from stream into a buffer.
..cat:Streams
..signature:_streamRead(target, stream, limit)
..param.target:A buffer that is filled.
..param.stream:An input stream.
...type:Adaption."std::iostream"
..param.limit:The maximal number of characters that is read from $stream$.
..returns:The number of characters read from $stream$.
*/
template <typename TValue, typename TTraits>
inline ::std::streamsize 
_streamRead(TValue * target,
			::std::basic_istream<TValue, TTraits> & source,
			::std::streamsize limit)
{
//IOREV
SEQAN_CHECKPOINT
	source.read(target, limit);
	return source.gcount();
}

//////////////////////////////////////////////////////////////////////////////

/**
.Internal._streamGet:
..summary:Read one character from stream.
..cat:Streams
..signature:_streamGet(stream)
..param.stream:An input stream.
...type:Adaption."std::iostream"
..returns:The character read.
*/

template <typename TValue, typename TTraits>
inline TValue 
_streamGet(::std::basic_istream<TValue, TTraits> & source)
{
//IOREV
SEQAN_CHECKPOINT
	return source.get();
}

//////////////////////////////////////////////////////////////////////////////

/**
.Internal._streamPeek:
..summary:Return the next character to be read from stream.
..cat:Streams
..signature:_streamPeek(stream)
..param.stream:An input stream.
...type:Adaption."std::iostream"
..returns:The character to be read.
*/

template <typename TValue, typename TTraits>
inline TValue 
_streamPeek(::std::basic_istream<TValue, TTraits> & source)
{
//IOREV
SEQAN_CHECKPOINT
	return source.peek();
}

//////////////////////////////////////////////////////////////////////////////

/**
.Internal._streamUnget:
..summary:Put the last read character back into stream.
..cat:Streams
..signature:_streamUnget(stream)
..param.stream:An input stream.
...type:Adaption."std::iostream"
*/

template <typename TValue, typename TTraits>
inline void
_streamUnget(::std::basic_istream<TValue, TTraits> & source)
{
//IOREV
SEQAN_CHECKPOINT
	source.unget();
}

//////////////////////////////////////////////////////////////////////////////

/**
.Internal._streamPut:
..summary:Writes one character to stream.
..cat:Streams
..signature:_streamPut(stream, character)
..param.stream:An input stream.
...type:Adaption."std::iostream"
..param.character:A character that is written to $stream$.
*/

template <typename TValue, typename TTraits, typename TChar>
inline void
_streamPut(::std::basic_ostream<TValue, TTraits> & target,
		   TChar character)
{
//IOREV
SEQAN_CHECKPOINT
	target.put(convert<TValue>(character));
}

//////////////////////////////////////////////////////////////////////////////

/**
.Internal._streamTellG:
..cat:Streams
..summary:Gets current position of input stream.
..signature:_streamTellG(stream)
..param.stream:An input stream.
...type:Adaption."std::iostream"
..returns:The current position in $stream$.
*/
template <typename TValue, typename TTraits>
inline typename Position< ::std::basic_istream<TValue, TTraits> >::Type
_streamTellG(::std::basic_istream<TValue, TTraits> & me)
{
//IOREV
SEQAN_CHECKPOINT
	return me.tellg();
}

//////////////////////////////////////////////////////////////////////////////
/**
.Internal._streamTellP:
..cat:Streams
..summary:Gets current position of output stream.
..signature:_streamTellP(stream)
..param.stream:An ouput stream.
...type:Adaption."std::iostream"
..returns:The current position in $stream$.
..see:Internal._streamTellG
*/
template <typename TValue, typename TTraits>
inline typename Position< ::std::basic_ostream<TValue, TTraits> >::Type
_streamTellP(::std::basic_ostream<TValue, TTraits> & me)
{
//IOREV
SEQAN_CHECKPOINT
	return me.tellp();
}

//////////////////////////////////////////////////////////////////////////////

/**
.Internal._streamSeekG:
..summary:Moves input stream to a position.
..cat:Streams
..signature:_streamSeekG(stream, position)
..param.stream:An input stream.
...type:Adaption."std::iostream"
..param.position:A position within the stream.
...remarks:Use @Function._streamTellG@ to get valid stream positions.
..see:Internal._streamTellG
*/
template <typename TValue, typename TTraits>
inline void
_streamSeekG(::std::basic_istream<TValue, TTraits> & me,
	 typename Position< ::std::basic_istream<TValue, TTraits> >::Type pos)
{
//IOREV
SEQAN_CHECKPOINT
	me.clear();
	me.seekg(pos);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Internal._streamSeekP:
..summary:Moves output stream to a position.
..cat:Streams
..signature:_streamSeekP(stream, position)
..param.stream:An output stream.
...type:Adaption."std::iostream"
..param.position:A position within the stream.
...remarks:Use @Function._streamTellP@ to get valid stream positions.
..see:Internal._streamTellP
..see:Internal._streamSeekG
*/
template <typename TValue, typename TTraits>
inline void
_streamSeekP(::std::basic_ostream<TValue, TTraits> & me,
	 typename Position< ::std::basic_ostream<TValue, TTraits> >::Type pos)
{
//IOREV
SEQAN_CHECKPOINT
	me.clear();
	me.seekp(pos);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Internal._streamSeek2G:
..summary:Moves input stream position relative to current position.
..cat:Streams
..signature:_streamSeek2G(stream, offset)
..param.stream:An input stream.
...type:Adaption."std::iostream"
..param.offset:The amout the position is changed.
...remarks:If this value is negative.
..see:Internal._streamSeekG
*/
template <typename TValue, typename TTraits>
inline void
_streamSeek2G(::std::basic_istream<TValue, TTraits> & me,
	 int off)
{
//IOREV
SEQAN_CHECKPOINT
	me.seekg(off, ::std::ios_base::cur);
}


//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
