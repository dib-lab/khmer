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

#ifndef SEQAN_HEADER_FILE_FASTA_H
#define SEQAN_HEADER_FILE_FASTA_H

/* IOREV
 * _tested_
 * _nodoc_
 *
 * tested in tests/file/test_file.h
 * tag mentionen in doc, but no further documentation, no link to spec
 * 
 * IMPORTANT: from what I understand: fileReaderIterator does not iterate
 * through records, but through lines of the sequence of one record.
 * goNext() on the iterator goes to beginning of next line and sets
 * data boundaries for next iteration
 *
 * goNext() on the file itself goes to the beginning of the next record.
 */



namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// File Formats - Fasta
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.File Format.tag.Fasta:
	FASTA file format for sequences.
..include:seqan/file.h
*/
struct TagFasta_;
//IOREV
typedef Tag<TagFasta_> Fasta; //IOREV

//////////////////////////////////////////////////////////////////////////////
// Filereader
//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TFile2, typename TSpec>
inline void
goBegin(Iter<TFile, FileReader<Fasta, TFile2, TSpec> > & it, bool skip_meta)
{
//IOREV _postponed_ due to unclear state of filereaderiterator
	if (_streamEOF(host(it)))
	{
		it.data_eof = true;
		return;
	}

	if (skip_meta && (it.data_char == '>'))
	{
		//skip meta line
		_streamSkipLine(host(it), it.data_char);
	}

	//eliminate linebreaks
	while ((it.data_char == '\n') || (it.data_char == '\r'))
	{
		if (_streamEOF(host(it)))
		{
			it.data_eof = true;
			return;
		}
		it.data_char = _streamGet(host(it));
	}

	if (it.data_char == '>')
	{//end of record
		it.data_eof = true;
		_streamUnget(host(it));
		return;
	}

	it.data_file_pos = _streamTellG(host(it)) - 1;
	it.data_eof = _streamEOF(host(it));
}
template <typename TFile, typename TFile2, typename TSpec>
inline void
goBegin(Iter<TFile, FileReader<Fasta, TFile2, TSpec> > & it)
{
//IOREV
	goBegin(it, true);
}


template <typename TFile, typename TFile2, typename TSpec>
inline void
goNext(Iter<TFile, FileReader<Fasta, TFile2, TSpec> > & it)
{
//IOREV goto next line of sequence in current record, set boundaries for
/*
	if (_streamEOF(host(it)))
	{
		it.data_eof = true;
		return;
	}
*/
	it.data_char = _streamGet(host(it));
	++it.data_file_pos;

	if (_streamEOF(host(it)))
	{
		it.data_eof = true;
		return;
	}

	if ((it.data_char == '\n') || (it.data_char == '\r'))
	{//linebreak detected: find begin of next line
		do
		{
			it.data_char = _streamGet(host(it));
			if (_streamEOF(host(it)))
			{
				it.data_eof = true;
				return;
			}
			++it.data_file_pos;
		} while ((it.data_char == '\n') || (it.data_char == '\r'));

		if (it.data_char == '>')
		{//end of record
			_streamUnget(host(it));
			it.data_eof = true;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// FileFormat Interface
//////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////
//count_valid: zaehlt die nicht-Zeilenumbrueche (input/output)
//count_all: zaehlt alle Zeichen incl. Zeilenumbrueche (input/output)
//returns: zuletzt gelesenes Zeichen = das erste hinter dem Zeilenumbruch bzw. eof
//the last read char is not counted!
//count_valid and count_all are not resetted but counted up
template <typename TFile, typename TSize>
inline typename Value<TFile>::Type
_fastaScanLine(TFile & file,
				 TSize & count_valid,
				 TSize & count_all)
{
//IOREV _nodoc_ has doc but in wrong format
SEQAN_CHECKPOINT
	SEQAN_ASSERT_NOT(_streamEOF(file));

	TSize count = 0;
	typename Value<TFile>::Type c;

	do
	{
		c = _streamGet(file);
		if (c == '\n' || c == '\r')
		{
			while (!_streamEOF(file))
			{
				++count_all;
				c = _streamGet(file);
				if (c != '\n' && c != '\r')
					break;
			}
			break;
		}

		++count;
	} while (!_streamEOF(file));
	
	count_valid += count;
	count_all += count;
	return c;
}


/////////////////////////////////////////////////////////////////////////
template <typename TFile, typename TSize>
inline void
_readNCharsFromFile(TFile & file, TSize count)
{
//IOREV shouldn't this check for EOF? Shouldn't it be renamed to sth like skipNChars, since it doesn't save them anywhere?
SEQAN_CHECKPOINT
	for (TSize i = 0; i < count; ++i)
	{
		_streamGet(file);
	}
}


//////////////////////////////////////////////////////////////////////////////
// read
//////////////////////////////////////////////////////////////////////////////

// TODO(holtgrew): The interface of this function does not correlate with documentation where the third parameter is metadata!
template <typename TFile, typename TData, typename TSize>
void
read(TFile & file,
	 TData & data,
	 TSize limit,
	 Fasta)
{
//IOREV _recordreading_ _notinlined_ see comment by holtgrew, signature is wrong for most file formats specializations
SEQAN_CHECKPOINT

	SEQAN_ASSERT_NOT(_streamEOF(file));
	clear(data);

	//determine begin position
	typename Value<TFile>::Type c_first = _streamGet(file);
	SEQAN_ASSERT_NOT(_streamEOF(file));

	typename Position<TFile>::Type begin_pos = _streamTellG(file);
	typename Size<TData>::Type count_valid = 1; //"valid" characters read (without line breaks)
	typename Size<TData>::Type count_all = 1;	//all characters read (with line breaks)

	if (_streamEOF(file))
	{
		return;
	}	

	if (c_first == '>')
	{//there is an id line: skip it
		c_first = _fastaScanLine(file, count_valid, count_all);
	}

	if ((c_first == '>') || _streamEOF(file)) 
	{//another id line = empty entry
		_streamSeekG(file, begin_pos);
		_readNCharsFromFile(file, count_all-1);
		return;
	}

	begin_pos = _streamTellG(file);

	count_valid = 1;
	count_all = 1;
	typename Value<TFile>::Type c;
	bool eof_reached = false;
	//determine length
	while (true)
	{
		c = _fastaScanLine(file, count_valid, count_all);
		if (_streamEOF(file)) 
		{//end of file: stop searching
			eof_reached = true;
			break;
		}
		if (c == '>')
		{//next entry found: stop seaching
			break;
		}
		if ((c != '\n') && (c != '\r'))
		{
			++count_valid; //count c
		}
		++count_all;
	}

	//reserve space
	typename Size<TData>::Type count = count_valid;
	if (count > limit)
		count = limit;
	
	resize(data, count);
	if (count > length(data))
		count = length(data);

	//read sequence
	_streamSeekG(file, begin_pos);

	typename Position<TData>::Type pos = 0;
	c = c_first;
	while (true)
	{
		if ((c != '\n') && (c != '\r'))
		{
			data[pos] = c;
			++pos;
		}
		if (pos >= count) break;

		c = _streamGet(file);
		--count_all;
	}

	//move file ptr to next entry
	_readNCharsFromFile(file, count_all - 1);
	if(eof_reached)
		_streamGet(file);
}

//____________________________________________________________________________

template <typename TFile, typename TData>
void
read(TFile & file,
	 TData & data,
	 Fasta tag)
{
//IOREV _notinlined_
SEQAN_CHECKPOINT
	typedef typename Size<TData>::Type TSize;
	read(file, data, maxValue<TSize>(), tag);
}


//////////////////////////////////////////////////////////////////////////////
// readID
//////////////////////////////////////////////////////////////////////////////
 
//the ID is the complete first line (without the leading '>'-sign)

template <typename TFile, typename TString>
void
readID(TFile & file,
	   TString & seqId,
	   Fasta)
{
//IOREV _nodoc_ _notinlined_ name ambiguous with readId-memfunc
SEQAN_CHECKPOINT
	SEQAN_ASSERT_NOT(_streamEOF(file));

	typename Position<TFile>::Type start_pos = _streamTellG(file);

	typename Value<TFile>::Type c = _streamGet(file);
	if (c != '>')
	{
		clear(seqId);
	}
	else
	{
		typename Size<TString>::Type count_valid = 0;
		typename Size<TString>::Type count_all = 0;
		_fastaScanLine(file, count_valid, count_all);

		if (! count_valid)
		{
			clear(seqId);
		}
		else
		{
			resize(seqId, count_valid);
			if (length(seqId) < count_valid)
			{
				count_valid = length(seqId);
			}

			_streamSeekG(file, start_pos);
			c = _streamGet(file); //pop the '>' character
			for (typename Position<TString>::Type pos = 0; count_valid; --count_valid)
			{
				seqId[pos] = _streamGet(file);
				++pos;
			}
		}
	}
	_streamSeekG(file, start_pos);
}

//short ID, read fasta header up to first whitespace
template <typename TFile, typename TString>
void
readShortID(TFile & file,
	   TString & seqId,
	   Fasta)
{
//IOREV _nodoc_ _notinlined_ better to use _stream* or _parse* calls to get first word
SEQAN_CHECKPOINT
	SEQAN_ASSERT_NOT(_streamEOF(file));

	typename Position<TFile>::Type start_pos = _streamTellG(file);

	typename Value<TFile>::Type c = _streamGet(file);
	if (c != '>')
	{
		clear(seqId);
	}
	else
	{
		typename Size<TString>::Type count_valid = 0;
		typename Size<TString>::Type count_all = 0;
		_fastaScanLine(file, count_valid, count_all);

		if (! count_valid)
		{
			clear(seqId);
		}
		else
		{
			resize(seqId, count_valid);
			if (length(seqId) < count_valid)
			{
				count_valid = length(seqId);
			}

			_streamSeekG(file, start_pos);
			c = _streamGet(file); //pop the '>' character
			for (typename Position<TString>::Type pos = 0; count_valid; --count_valid)
			{
				seqId[pos] = _streamGet(file);
				if(seqId[pos]=='\t' || seqId[pos]=='\b' || seqId[pos]==' ')
				{
					resize(seqId,pos);
					break;
				}
				++pos;
			}
		}
	}
	_streamSeekG(file, start_pos);
}


//////////////////////////////////////////////////////////////////////////////
// readMeta
//////////////////////////////////////////////////////////////////////////////

//Fasta file records have no meta data

template <typename TFile, typename TMeta>
void
readMeta(TFile & file,
		 TMeta & meta,
		 Fasta)
{
//IOREV _nodoc_ _notinlined_ code documentation wrong
SEQAN_CHECKPOINT
	readID(file, meta, Fasta());
}


//////////////////////////////////////////////////////////////////////////////
// goNext
//////////////////////////////////////////////////////////////////////////////

template <typename TFile>
void
goNext(TFile & file,
	   Fasta)
{
//IOREV _notinlined_ goto next record
SEQAN_CHECKPOINT
	SEQAN_ASSERT_NOT(_streamEOF(file));

//	bool found_data = false;
	while (true)
	{
		typename Value<TFile>::Type c = _streamGet(file);

		if (_streamEOF(file)) return;

		if (c == '\n' || c == '\r')
		{
			do {
				c = _streamGet(file);
				if (_streamEOF(file)) return;
			} while (c == '\n' || c == '\r');
/*
// weese: I changed the following lines, as otherwises empty
//        sequences (seq1) will be overread. See:
// >seq1
// >seq2
// ACTGGT
			if (c != '>')
			{
				found_data = true;
			}
			else if (found_data)
			{
				_streamUnget(file);
				return;
			}
*/			if (c == '>')
			{
				_streamUnget(file);
				return;
			}
		}
	}
}



//////////////////////////////////////////////////////////////////////////////
// write
//////////////////////////////////////////////////////////////////////////////


template <typename TFile, typename TString, typename TData>
void
_writeImpl(TFile & file,
			TData & data,
			TString & seqId,
			Fasta)
{
//IOREV _notinlined_
SEQAN_CHECKPOINT
	_streamPut(file, '>');
	_streamWrite(file, seqId);
	_streamPut(file, '\n');

	//typename Iterator<TData, Standard>::Type it = begin(data, Standard());
	//typename Iterator<TData, Standard>::Type it_end = end(data, Standard());
	typename Iterator<TData>::Type it = begin(data);
	typename Iterator<TData>::Type it_end = end(data);

	int i = 0;

	for (; it < it_end; ++it)
	{
		if (i == 60)
		{
			_streamPut(file, '\n');
			i = 0;
		}
		++i;

		_streamPut(file, *it);
	}
	_streamPut(file, '\n');
}

//____________________________________________________________________________

template <typename TFile, typename TString, typename TData, typename TMeta>
void
write(TFile & file,
	  TData & data,
	  Fasta)
{
//IOREV _notinlined_
SEQAN_CHECKPOINT
	_writeImpl(file, data, "", Fasta());
}

//____________________________________________________________________________

template <typename TFile, typename TString, typename TData>
void
write(TFile & file,
	  TData & data,
	  TString & seqId,
	  Fasta)
{
//IOREV _notinlined_
SEQAN_CHECKPOINT
	_writeImpl(file, data, seqId, Fasta());
}


//VisualC++ const array bug workaround
template <typename TFile, typename TString, typename TDataValue>
void
write(TFile & file,
	  TDataValue * data,
	  TString & seqId,
	  Fasta)
{
//IOREV _notinlined_
SEQAN_CHECKPOINT
	_writeImpl(file, data, seqId, Fasta());

}

//____________________________________________________________________________

template <typename TFile, typename TString, typename TData, typename TMeta>
void
write(TFile & file,
	  TData & data,
	  TString & seqId,
	  TMeta &,
	  Fasta)
{
//IOREV _notinlined_
SEQAN_CHECKPOINT
	_writeImpl(file, data, seqId, Fasta());
}


//////////////////////////////////////////////////////////////////////////////
} //namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef SEQAN_HEADER_...
