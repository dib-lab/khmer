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

#ifndef SEQAN_HEADER_FILE_GENBANK_H
#define SEQAN_HEADER_FILE_GENBANK_H

/* IOREV
 * _tested_
 * _nodoc_
 *
 * tested in tests/file/test_file.h
 * tag mentionen in doc, but no further documentation, no link to spec
 *
 * current spec:
 * ftp://ftp.ncbi.nih.gov/genbank/gbrel.txt
 *
 * IMPORTANT: from what I understand: fileReaderIterator does not iterate
 * through records, but through lines of the sequence of one record.
 * goNext() on the iterator goes to beginning of next line and sets
 * data boundaries for next iteration
 *
 * goNext() on the file itself goes to the beginning of the next record.
 *
 * THIS is not intuitive AT ALL.
 */


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// File Formats - Genbank
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.File Format.tag.Genbank:
	Genbank format for sequences from the Genbank database.
..include:seqan/file.h
*/
struct TagGenbank_;
//IOREV
typedef Tag<TagGenbank_> const Genbank; //IOREV


//////////////////////////////////////////////////////////////////////////////
// FileReader Iterator
//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TFile2, typename TSpec>
inline void
goBegin(Iter<TFile, FileReader<Genbank, TFile2, TSpec> > & it, bool skip_meta)
{
//IOREV not sure whether begin_pos is absolute; parsing could be done simpler; no error handling
SEQAN_CHECKPOINT
	String<char> line;

	if (_streamEOF(host(it)))
	{//end of file
		it.data_eof = true;
		return;
	}

	if (skip_meta && (it.data_char != ' '))
	{//skip metadata block
		while (true)
		{
			if (it.data_char == '/')
			{//end of record
				_streamSkipLine(host(it), it.data_char);
				it.data_eof = true;
				return;
			}
			if ((it.data_char == 'O') || (it.data_char == 'o'))
			{
				clear(line);
				_streamAppendLine(host(it), line, it.data_char);
				if ((prefix(line, 6) == "ORIGIN") || (prefix(line, 6) == "origin"))
				{//end of metadata
					break;
				}
			}
			//skip meta line
			_streamSkipLine(host(it), it.data_char);

			if (_streamEOF(host(it)))
			{//end of file
				it.data_eof = true;
				return;
			}
		}
	}

	//find first character
	while (true)
	{
		if (_streamEOF(host(it)))
		{//end of file
			it.data_eof = true;
			return;
		}
		if ((it.data_char != ' ') && ((it.data_char < '0') || (it.data_char > '9')))
		{
			if ((it.data_char != '\n') && (it.data_char != '\r'))
			{//fist char found
				break;
			}

			it.data_char = _streamGet(host(it));
			if (it.data_char == '/')
			{//end of record
				_streamSkipLine(host(it), it.data_char);
				it.data_eof = true;
				return;
			}
		}
		else
		{
			it.data_char = _streamGet(host(it));
		}
	}

	it.data_file_pos = _streamTellG(host(it));
	it.data_file_pos -=1;
	it.data_eof = _streamEOF(host(it));
}


template <typename TFile, typename TFile2, typename TSpec>
inline void
goBegin(Iter<TFile, FileReader<Genbank, TFile2, TSpec> > & it)
{
//IOREV
    SEQAN_CHECKPOINT;
    goBegin(it, true);
}


template <typename TFile, typename TFile2, typename TSpec>
inline void
goNext(Iter<TFile, FileReader<Genbank, TFile2, TSpec> > & it)
{
//IOREV this iterates through the characters of a line in one record
SEQAN_CHECKPOINT
	do
	{
		it.data_char = _streamGet(host(it));
		if (_streamEOF(host(it)))
		{
			it.data_eof = true;
			return;
		}
		it.data_file_pos += 1;

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
				it.data_file_pos += 1;
			} while ((it.data_char == '\n') || (it.data_char == '\r'));

			if (it.data_char == '/')
			{//end of record
				_streamSkipLine(host(it), it.data_char);
				_streamUnget(host(it));
				//it.data_file_pos is invalid now
				it.data_eof = true;
				return;
			}
		}
	} while ((it.data_char == ' ') || ((it.data_char >= '0') && (it.data_char <= '9')));
}

//////////////////////////////////////////////////////////////////////////////
// File Format Access Function
//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TData>
inline void
read(TFile & file,
	 TData & data,
	 Genbank)
{
//IOREV _recordreading_ see above
SEQAN_CHECKPOINT
	Iter<TFile, FileReader<Genbank> > it(file);

	clear(data);
	while (!atEnd(it))
	{
		appendValue(data, getValue(it));
		goNext(it);
	}
}

template <typename TFile, typename TData, typename TSize>
inline void
read(TFile & file,
	 TData & data,
	 TSize limit,
	 Genbank)
{
//IOREV _nodoc_ _recordreading_ see above
SEQAN_CHECKPOINT
	typename Size<TData>::Type siz = length(data);
	Iter<TFile, FileReader<Genbank> > it(file);

	clear(data);
	while (!atEnd(it) && (siz < limit))
	{
		appendValue(data, getValue(it));
		goNext(it);
	}
	while (!atEnd(it))
	{
		goNext(it);
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TMeta>
inline void
readMeta(TFile & file,
		 TMeta & meta,
		 Genbank)
{
//IOREV _ndoc_ _recordreading_ see above
SEQAN_CHECKPOINT
	typedef typename Value<TMeta>::Type TValue;
	String<char> line;

	clear(meta);

	if (_streamEOF(file))
	{
		return;
	}

	TValue c = _streamGet(file);

	while (!_streamEOF(file))
	{
		clear(line);
		_streamAppendLine(file, line, c);

		if (c == '/')
		{//end of record
			_streamUnget(file);
			break;
		}

		append(meta, line);
		appendValue(meta, '\n');

		if ((prefix(line, 6) == "ORIGIN") || (prefix(line, 6) == "origin"))
		{//end of metadata
			_streamUnget(file);
			break;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile>
inline void
goNext(TFile & file,
	   Genbank)
{
//IOREV this seems to iterate to the enxt record
SEQAN_CHECKPOINT
	typedef typename Value<TFile>::Type TValue;

	if (_streamEOF(file))
	{
		return;
	}

	while (!_streamEOF(file))
	{
		TValue c = _streamGet(file);
		if (c == '/')
		{//end of record
			_streamSkipLine(file, c);
			_streamUnget(file);
			return;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
/*
template <typename TFile>
inline void
length(TFile & file,
	   Genbank)
{
SEQAN_CHECKPOINT
}
*/
//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TData>
inline void
write(TFile & file,
	  TData & data,
	  Genbank)
{
//IOREV _recordreading_
SEQAN_CHECKPOINT
	enum
	{
		BLOCK_SIZE = 10,
		BLOCKS_PER_LINE = 6
	};
	char const * NUM_BLOCK_FORMAT = "%9d";

	typedef typename Size<TData>::Type TSize;
	typedef typename Iterator<TData, Standard>::Type TIterator;

	TSize count = 0;
	TIterator it = begin(data, Standard());
	TIterator it_end = end(data, Standard());

	while (it != it_end)
	{
		//write count 
		_streamPutInt(file, count+1, NUM_BLOCK_FORMAT);

		int block_count = 0;
		int char_in_block_count = BLOCK_SIZE;

		//write rest of line
		while (it != it_end)
		{
			if (char_in_block_count == BLOCK_SIZE)
			{//begin new block
				if (block_count >= BLOCKS_PER_LINE)
				{//end of line
					_streamPut(file, '\n');
					break;
				}
				_streamPut(file, ' ');
				char_in_block_count = 0;
				++block_count;
			}

			//write next character
			_streamPut(file, *it);
			++it;
			++count;
			++char_in_block_count;
		}
	}

	write(file, "\n//\n");
}

template <typename TFile, typename TData, typename TMeta>
inline void
write(TFile & file,
	  TData & data,
	  TMeta & meta,
	  Genbank)
{
//IOREV _recordreading_
SEQAN_CHECKPOINT
	write(file, meta);
	write(file, data, Genbank());
}


//////////////////////////////////////////////////////////////////////////////
} //namespace SEQAN_NAMESPACE_MAIN
//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef SEQAN_HEADER_...
