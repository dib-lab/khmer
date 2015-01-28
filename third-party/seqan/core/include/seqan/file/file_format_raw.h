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

#ifndef SEQAN_HEADER_FILE_RAW_H
#define SEQAN_HEADER_FILE_RAW_H


/* IOREV
 * _tested_
 * _doc_
 *
 *
 * has some documentation, but could be better
 * apperently tested by tests/file/test_file.h, used in some other places
 *
 * unclear why static members are used and functions not overloaded directly
 *
 * raw is defined as default behaviour for read and write, is this desired
 * or should read and write maybe auto-detect fileformat?
 * 
 */

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// File Formats - Raw
//////////////////////////////////////////////////////////////////////////////
/**
.Tag.File Format.tag.Raw:
	The file contains data in a raw format.
..remark:It is supposed that the file contains one single piece of data, 
that is the file cannot store multiple records.
..include:seqan/file.h
*/

struct TagRaw_;
//IOREV
typedef Tag<TagRaw_> Raw; //IOREV



//////////////////////////////////////////////////////////////////////////////
// read
//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TData, typename TTag>
struct ReadRaw_;
//IOREV

//____________________________________________________________________________

template <typename TFile, typename TData>
struct ReadRaw_<TFile, TData, True>
{
//IOREV
	template <typename TSize>
	inline static void
	read_(TFile & file,
		TData & data,
		TSize _limit)
	{
SEQAN_CHECKPOINT
		SEQAN_ASSERT_NOT(_streamEOF(file));

		typename Size<TData>::Type limit = _limit;

		//determine length
		typename Position<TFile>::Type begin_pos = _streamTellG(file);
		typename Size<TData>::Type count;
		typename Position<TData>::Type pos;

		for (count = 0; !_streamEOF(file) && count != limit; ++count)
			_streamGet(file);

		//reserve space
		resize(data, count);

		if (!count) return;

		if (length(data) < count)
			count = length(data);

		//read sequence
		_streamSeekG(file, begin_pos);

		for (pos = 0; pos < count; ++pos)
			assignValue(data, pos, _streamGet(file));
	}

//____________________________________________________________________________

	inline static void
	read_(TFile & file,	
		TData & data)
	{
SEQAN_CHECKPOINT
		typedef typename Size<TData>::Type TSize;
		read(file, data, maxValue<TSize>());
	}

};

//____________________________________________________________________________

template <typename TFile, typename TData>
struct ReadRaw_<TFile, TData, False>
{
//IOREV _notinlined_
	static void
	read_(TFile & file,
		TData & data)
	{
SEQAN_CHECKPOINT

		clear(data);
		if (!_streamEOF(file))
		{
SEQAN_CHECKPOINT
			ChunkCollector_<TData> chunk_collector(data);
			assign(chunk_collector, file);
			append(data, chunk_collector);
		}
	}

//____________________________________________________________________________

	template <typename TSize>
	static void
	read_(TFile & file,
		TData & data,
		TSize limit)
	{
SEQAN_CHECKPOINT

		clear(data);
		if (!_streamEOF(file))
		{
SEQAN_CHECKPOINT
			ChunkCollector_<TData> chunk_collector(data);
			assign(chunk_collector, file, limit);
			append(data, chunk_collector, limit);
		}
	}
};

//////////////////////////////////////////////////////////////////////////////


template <typename TFile, typename TData>
void
read(TFile & file,
	 TData & data,
	 Raw)
{
//IOREV _notinlined_
SEQAN_CHECKPOINT
	ReadRaw_<TFile, TData, typename IsTellAndSeekStream_<TFile>::Type>::read_(file, data);
}

//____________________________________________________________________________

template <typename TFile, typename TData, typename TSize>
void
read(TFile & file,
	 TData & data,
	 TSize limit,
	 Raw)
{
//IOREV _notinlined_
SEQAN_CHECKPOINT
	ReadRaw_<TFile, TData, typename IsTellAndSeekStream_<TFile>::Type>::read_(file, data, limit);
}


//////////////////////////////////////////////////////////////////////////////
// readID
//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TString>
void
readID(TFile & /*file*/,
	   TString & id,
	   Raw)
{
//IOREV _notinlined_
SEQAN_CHECKPOINT
	clear(id);
}

//////////////////////////////////////////////////////////////////////////////
// readMeta
//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TMeta>
void
readMeta(TFile & /*file*/,
		 TMeta & meta,
		 Raw)
{
//IOREV _notinlined_
SEQAN_CHECKPOINT
	clear(meta);
}


//////////////////////////////////////////////////////////////////////////////
// goNext
//////////////////////////////////////////////////////////////////////////////

template <typename TFile>
void
goNext(TFile & file,
	   Raw)
{
//IOREV  _notinlined_ _nodoc_ needs some documented defined behaviour
SEQAN_CHECKPOINT
  (void) file;  // When compiled without assertions.
	SEQAN_ASSERT(!_streamEOF(file));

//??? TODO: set file to eof
}



//////////////////////////////////////////////////////////////////////////////
// write
//////////////////////////////////////////////////////////////////////////////


template <typename TFile, typename TData>
void
write(TFile & file,
	  TData const & data,
	  Raw)
{
//IOREV _notinlined_
SEQAN_CHECKPOINT
	_streamWrite(file, data);
}

//____________________________________________________________________________

template <typename TFile, typename TData, typename TString>
void
write(TFile & file,
	  TData const & data,
	  TString const &,
	  Raw)
{
//IOREV _notinlined_
SEQAN_CHECKPOINT
	_streamWrite(file, data);
}

//____________________________________________________________________________

template <typename TFile, typename TString, typename TData, typename TMeta>
void
write(TFile & file,
	  TData const & data,
	  TString const &,
	  TMeta const &,
	  Raw)
{
//IOREV _notinlined_
SEQAN_CHECKPOINT
	_streamWrite(file, data);
}

//////////////////////////////////////////////////////////////////////////////
// default functions
//////////////////////////////////////////////////////////////////////////////


template <typename TFile, typename TData>
void
read(TFile & file,
	 TData & data)
{
//IOREV _notinlined_
	read(file, data, Raw());
}

template <typename TFile, typename TData, typename TSize>
void
read(TFile & file,
	 TData & data,
	 TSize limit)
{
//IOREV  _notinlined_
SEQAN_CHECKPOINT
	read(file, data, limit, Raw());
}

//____________________________________________________________________________

template <typename TFile, typename TData>
void
write(TFile & file,
	  TData & data)
{
//IOREV _notinlined_
SEQAN_CHECKPOINT
	write(file, data, "", Raw());
}
template <typename TFile, typename TData>
void
write(TFile & file,
	  TData const & data)
{
//IOREV _notinlined_
SEQAN_CHECKPOINT
	write(file, data, "", Raw());
}

//////////////////////////////////////////////////////////////////////////////
} //namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef SEQAN_HEADER_...
