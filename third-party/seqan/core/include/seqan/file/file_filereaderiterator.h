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

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#ifndef SEQAN_HEADER_FILE_FILEREADEITERATOR_H
#define SEQAN_HEADER_FILE_FILEREADEITERATOR_H

/* IOREV
 * _nodoc_
 * 
 * 
 * relation to file_filereader.h not totally clear
 *
 * broken according to holtgrew, but used by many file_format implementations
 * that supposedly work
 *
 * not yet fully understood how this works or is supposed to
 *
 * 
 */



namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

template <typename TFormat, typename TFile = FILE*, typename TSpec = Default>
struct FileReader;
//IOREV

// Manual forward for the skipMeta goBegin() overload.
//template <typename TFile, typename TFormat, typename TFile2, typename TSpec>
//void
//goBegin(Iter<TFile, FileReader<TFormat, TFile2, TSpec> > & it, bool skip_meta);

//////////////////////////////////////////////////////////////////////////////
// FileReader: an iterator that scans through the data of a file
// note: this is not the iterator of the FileReader string (see file_filereader.h)
//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TFormat, typename TFile2, typename TSpec>
class Iter<TFile, FileReader<TFormat, TFile2, TSpec> >
{
//IOREV
public:
	typedef typename Value<TFile>::Type TValue;
	typedef typename Position<TFile>::Type TPosition;

	TPosition data_file_pos;	//position of the last read char relative to data begin (only valid if data_eof == false)
	TFile * data_host;		//the host file
	TValue data_char;		//the last read char
	bool data_eof;			//true if reached end of record
//	TFilePosition data_begin_pos;

	Iter(TFile & file_, bool skip_meta = true):
		data_file_pos(0),
		data_host(& file_),
		data_eof(false)
	{
		data_char = _streamGet(file_);
		goBegin(*this, skip_meta);
	}
	Iter(Iter const & other_):
		data_file_pos(other_.data_file_pos),
		data_host(other_.data_host),
		data_char(other_.data_char),
		data_eof(other_.data_eof)
	{
	}
	~Iter() 
	{
	}

	Iter const &
	operator = (Iter const & other_)
	{
		data_file_pos = other_.data_file_pos;
		data_host = other_.data_host;
		data_char = other_.data_char;
		data_eof = other_.data_eof;
		return *this;
	}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TFormat, typename TFile2, typename TSpec>
struct Value< Iter<TFile, FileReader<TFormat, TFile2, TSpec> > >:
	Value<TFile>
{
//IOREV
};

template <typename TFile, typename TFormat, typename TFile2, typename TSpec>
struct GetValue< Iter<TFile, FileReader<TFormat, TFile2, TSpec> > >
{
//IOREV
	typedef typename Value< Iter<TFile, FileReader<TFormat, TFile2, TSpec> > >::Type Type;
};

template <typename TFile, typename TFormat, typename TFile2, typename TSpec>
struct Reference< Iter<TFile, FileReader<TFormat, TFile2, TSpec> > >
{
//IOREV
	typedef typename Value< Iter<TFile, FileReader<TFormat, TFile2, TSpec> > >::Type & Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TFormat, typename TFile2, typename TSpec>
inline TFile &
host(Iter<TFile, FileReader<TFormat, TFile2, TSpec> > & it)
{
//IOREV
	return *(it.data_host);
}


template <typename TFile, typename TFormat, typename TFile2, typename TSpec>
inline typename Reference<Iter<TFile, FileReader<TFormat, TFile2, TSpec> > >::Type
value(Iter<TFile, FileReader<TFormat, TFile2, TSpec> > & it)
{
//IOREV
	return it.data_char;
}

template <typename TFile, typename TFormat, typename TFile2, typename TSpec>
inline typename GetValue<Iter<TFile, FileReader<TFormat, TFile2, TSpec> > >::Type
getValue(Iter<TFile, FileReader<TFormat, TFile2, TSpec> > & it)
{
//IOREV
	return it.data_char;
}

template <typename TFile, typename TFormat, typename TFile2, typename TSpec>
inline bool
atEnd(Iter<TFile, FileReader<TFormat, TFile2, TSpec> > & it)
{
//IOREV
	return it.data_eof;
}


//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef SEQAN_HEADER_...
