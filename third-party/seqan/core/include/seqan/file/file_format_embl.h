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

#ifndef SEQAN_HEADER_FILE_EMBL_H
#define SEQAN_HEADER_FILE_EMBL_H

/* IOREV
 * _tested_
 * _nodoc_
 *
 * tested in tests/file/test_file.h
 * tag mentioned in doc, but no further documentation, no link to spec
 *
 * should be broken, because based on filereaderiterator which is broken
 * according to holtgrew.
 *
 * IMPORTANT: from what I understand: fileReaderIterator does not iterate
 * through records, but through lines of the sequence of one record.
 * goNext() on the iterator goes to beginning of next line and sets
 * data boundaries for next iteration
 *
 * goNext() on the file itself goes to the beginning of the next record.
 *
 */

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// File Formats - Embl
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.File Format.tag.Embl:EMBL format for sequences.
..include:seqan/file.h
*/
struct TagEmbl_;
//IOREV
typedef Tag<TagEmbl_> const Embl; //IOREV



//////////////////////////////////////////////////////////////////////////////
// FileReader Iterator
//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TFile2, typename TSpec>
inline void
goBegin(Iter<TFile, FileReader<Embl, TFile2, TSpec> > & it, bool skip_meta)
{
//IOREV _bug_ i am not convinced this works. current position is never reset
SEQAN_CHECKPOINT
	if (skip_meta)
	{
		while (true)
		{
			if (_streamEOF(host(it)))
			{
				it.data_eof = true;
				return;
			}
			if (it.data_char == '/')
			{//end of record
				_streamSkipLine(host(it), it.data_char);
				it.data_eof = true;
				return;
			}
			if (it.data_char == ' ')
			{
				break;
			}
			//skip meta line
			_streamSkipLine(host(it), it.data_char);
		}
	}

	//find first character
	while (true)
	{
		if (_streamEOF(host(it)))
		{
			it.data_eof = true;
			return;
		}
		if ((it.data_char != ' ') && ((it.data_char < '0') || (it.data_char > '9')))
		{
			if ((it.data_char != '\n') && (it.data_char != '\r')) break;

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

//	it.data_file_pos = _streamTellG(host(it));
	it.data_file_pos -= 1;
	it.data_eof = _streamEOF(host(it));
}


template <typename TFile, typename TFile2, typename TSpec>
inline void
goBegin(Iter<TFile, FileReader<Embl, TFile2, TSpec> > & it)
{
//IOREV
    SEQAN_CHECKPOINT;
    goBegin(it, true);
}


template <typename TFile, typename TFile2, typename TSpec>
inline void
goNext(Iter<TFile, FileReader<Embl, TFile2, TSpec> > & it)
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
	 Embl)
{
//IOREV _recordreading_
SEQAN_CHECKPOINT
	Iter<TFile, FileReader<Embl> > it(file);

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
	 Embl)
{
//IOREV _recordreading_
SEQAN_CHECKPOINT
	typename Size<TData>::Type siz = length(data);
	Iter<TFile, FileReader<Embl> > it(file);

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
		 Embl)
{
//IOREV _recordreading_ weird! 
SEQAN_CHECKPOINT
	typedef typename Value<TMeta>::Type TValue;

	clear(meta);
	if (_streamEOF(file))
	{
		return;
	}

	TValue c = _streamGet(file);

	while (!_streamEOF(file))
	{
		if (c == ' ')
		{//end of meta data
			_streamUnget(file);
			return;
		}
		if (c == '/')
		{//end of record
			_streamSkipLine(file, c);
			_streamUnget(file);
			return;
		}

		_streamAppendLine(file, meta, c);
		appendValue(meta, '\n');
	}
}

//////////////////////////////////////////////////////////////////////////////



/**
.Function.readLineType:
..cat:Input/Output
..summary:Reads the information belonging to the two-character line code specified.
..signature:readLineType(file,data,key,Embl);
..param.file:The input file or string.
...remarks:This function works on an open file stream or on the string data obtained from calling Function.readMeta
..param.data:The target container that will be filled.
..param.key:The two-character code specifying the file entry to be read, e.g. "AC" for the acession number line or "DE" for the description line. 
..see:Function.readMeta
..see:Function.readFeature
..include:seqan/file.h
*/
template<typename TFile, typename TData, typename TKey>
inline void
readLineType(TFile & file,
			 TData & data,
			 TKey key,
			 Embl)

{
//IOREV _recordreading_ weird! 
SEQAN_CHECKPOINT

	//this function is meant to be used for two letter codes only 
	SEQAN_ASSERT(length(key)==2); ;

	typedef typename Value<TFile>::Type TValue;
	typedef typename Position<TFile>::Type TPosition;

	clear(data);
	if(_streamEOF(file))
		return;
	
	TPosition pos = _streamTellG(file);
	TValue c = _streamGet(file);
	while (!_streamEOF(file))
	{
		if(c == '/')
		{
			_streamSeekG(file,pos);
			return;
		}
		if(c == key[0])
		{
			c = _streamGet(file);
			if(c == key[1])
			{
				for(unsigned int i = 0; i < 4; ++i)
					c = _streamGet(file);
				_streamAppendLine(file, data, c);
				while(!_streamEOF(file) && _streamReadWord(file,c) == key)
				{
					appendValue(data, '\n');
					for(unsigned int i = 0; i < 3; ++i)
						c = _streamGet(file);
					_streamAppendLine(file, data, c);
				}
				_streamSeekG(file,pos);
				return;
			}
		}
		_streamSkipLine(file, c);
	}

	_streamSeekG(file,pos);
	
}



template<typename TData, typename TValue, typename TSpec, typename TKey>
inline void
readLineType(String<TValue,TSpec> & meta,
			 TData & data,
			 TKey key,
			 Embl)

{
//IOREV
	//this function is meant to be used for two letter codes only
	SEQAN_ASSERT(length(key)==2); ;

	typedef typename Iterator<String<TValue,TSpec>,Standard>::Type TIterator;

	clear(data);
	if(empty(meta))
		return;

	TIterator it = begin(meta,Standard());
	TIterator end_it = end(meta,Standard());

	while (it != end_it)
	{
		if(*it == '/')
			return;

		if(*it == key[0])
		{
			++it;
			if(*it == key[1])
			{
				it+=4;
				_stringAppendLine(meta, data, it);
				while(it!=end_it && *it==key[0] && *(++it)==key[1])
				{
					appendValue(data, '\n');
					it+=4;
					_stringAppendLine(meta, data, it);
				}
				return;
			}
		}
		_stringSkipLine(meta, it);
	}

	
}


/**
.Function.readFeature:
..cat:Input/Output
..summary:Finds the first feature specified by 'key' starting from position 'start' in the feature table (the feature table can be
obtained by calling readLineType with the two-character code "FT").
..signature:readFeature(ft_string,start,data,key,Embl);
..param.ft_string:The feature table.
..param.start:Position in feature table where search begins.
..param.data:The target container that will be filled.
..param.key:The key word specifying the feature to be read, e.g. "mRNA" or "CDS".
..return:The position behind the feature if found, 0 otherwise.
..see:Function.readMeta
..see:Function.readFeature
..include:seqan/file.h
*/
//read parts of feature table (those that belong to key)
template<typename TData, typename TKey, typename TString>
inline typename Position<TString>::Type
readFeature(TString & str,
			typename Position<TString>::Type start_pos,
			TData & data,
			TKey key,
			Embl)

{
//IOREV
	typedef typename Iterator<TString,Standard>::Type TIterator;

	clear(data);
	if(empty(str) || start_pos >= length(str))
		return 0;
	TIterator it = iter(str,start_pos,Standard());
	TIterator end_it = end(str,Standard());

	while (it != end_it)
	{
		if(*it == key[0])
		{
			++it;
			bool found = true;
			for(unsigned int i = 1; i < length(key); ++i)
			{
				if(key[i] != *it)
				{
					found = false;
					break;
				}
				++it;
			}
			if(found)
			{
				_stringSkipWhitespace(str,it);
				_stringAppendLine(str, data, it);
				while(it!=end_it && *it == ' ')
				{
					appendValue(data, '\n');
					_stringSkipWhitespace(str,it);
					_stringAppendLine(str, data, it);
				}
				return position(it,str);
			}
		}
		_stringSkipLine(str, it);
	}

	return 0;
	
}
	


//////////////////////////////////////////////////////////////////////////////
template <typename TFile>
inline void
goNext(TFile & file,
	   Embl)
{
//IOREV this iterates to the next record
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
	   Embl)
{
SEQAN_CHECKPOINT
}
*/
//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TData>
inline void
write(TFile & file,
	  TData & data,
	  Embl)
{
//IOREV _recordreading_ weird. see above; for-loops would make code more readable
SEQAN_CHECKPOINT
	enum
	{
		BLOCK_SIZE = 10,
		BLOCKS_PER_LINE = 6,
		FIRST_INDENT = 5
	};
	char const * NUM_BLOCK_FORMAT = "%9d";

	typedef typename Size<TData>::Type TSize;
	typedef typename Iterator<TData, Standard>::Type TIterator;

	TSize count = 0;
	int block_count = 0;
	int char_in_block_count = 0;
	TIterator it = begin(data, Standard());
	TIterator it_end = end(data, Standard());

	while (it != it_end)
	{
		//write indent
		for (int j = 0; j < FIRST_INDENT; ++j) _streamPut(file, ' ');

		//write rest of line
		while (true)
		{
			//write next character
			if (it != it_end)
			{
				_streamPut(file, *it);
				++it;
				++count;
			}
			else
			{//no more chars left: fill up with ' '
				_streamPut(file, ' ');
			}
			++char_in_block_count;

			if (char_in_block_count >= BLOCK_SIZE)
			{//end of block
				char_in_block_count = 0;
				_streamPut(file, ' ');
				++block_count;

				if (block_count >= BLOCKS_PER_LINE)
				{//enough blocks: write end of line
					block_count = 0;
					_streamPutInt(file, count, NUM_BLOCK_FORMAT);
					_streamPut(file, '\n');
					break; //end of line
				}
			}
		}
	}

	write(file, "//\n");
}

template <typename TFile, typename TData, typename TMeta>
inline void
write(TFile & file,
	  TData & data,
	  TMeta & meta,
	  Embl)
{
//IOREV _recordreading_
SEQAN_CHECKPOINT
	write(file, meta);
	write(file, data, Embl());
}

//////////////////////////////////////////////////////////////////////////////

/*leere Vorlage
template <typename TFile, typename TData>
inline void
read(TFile & file,
	 TData & data,
	 Embl)
{
SEQAN_CHECKPOINT
}

template <typename TFile, typename TData, typename TSize>
inline void
read(TFile & file,
	 TData & data,
	 TSize limit,
	 Embl)
{
SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TMeta>
inline void
readMeta(TFile & file,
		 TMeta & meta,
		 Embl)
{
SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile>
inline void
goNext(TFile & file,
	   Embl)
{
SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile>
inline void
length(TFile & file,
	   Embl)
{
SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TData>
inline void
write(TFile & file,
	  TData & data,
	  Embl)
{
SEQAN_CHECKPOINT
}
template <typename TFile, typename TData, typename TMeta>
inline void
write(TFile & file,
	  TData & data,
	  TMeta & meta,
	  Embl)
{
SEQAN_CHECKPOINT
}
*/

//////////////////////////////////////////////////////////////////////////////
} //namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef SEQAN_HEADER_...
