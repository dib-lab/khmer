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

#ifndef SEQAN_HEADER_FILE_FORMAT_H
#define SEQAN_HEADER_FILE_FORMAT_H

/* IOREV
 * _tested_
 * _doc_
 *
 *
 * well documented
 * apperently tested by tests/file/test_file.h
 *
 * 
 * AFAICT this not really used widely, as the FileFormat specialization
 * tags are used diretly e.g.
 *      read(TFile, TData, TMeta, Fasta)
 * and not
 *      read(TFile, TData, TMeta, FileFormat<TFile, TData, TMeta, Fasta>)
 *
 * the global wrappers in this file also break down to this behaviours
 *
 * however it is not evident what the purpose of the FileFormat class is
 * altogether and especially the virtual members
 * (except being return value for guessFileFormat, which could aswell return
 * the tag)
 *
 *------------------
 *
 * contains all sorts of stream-IO helper functions that should go someplace
 * else.
 *
 * For all of these: Should we really look at c passed as argument, how do
 * we know, c is actually current character?
 *
 * For the iterators the current mechanism is error-prone:
 *
 * foobar (TInput & input, TIterator it)
 * {
 *     Iterator<TInput, Standard()> it_end = end(input);
 *     while (it != it_end)
 *          ...
 *
 * }
 *
 * this does not guarantee, that it and it_end are the same type, since
 * it maybe non-standard
 *
 * also some functions set it_end = end(input) - 1; which seems to be related
 * to change of atEnd behaviour http://trac.mi.fu-berlin.de/seqan/wiki/IoRevision
 *
 * IMO all these functions should be changed to not construct a 2nd iterator
 * but use atEnd() directly
 *      
 */


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.File Format:
..cat:Input/Output
..summary:A file format.
..include:seqan/file.h
*/


//////////////////////////////////////////////////////////////////////////////
// Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//Base Class for all FileFormat classes
//////////////////////////////////////////////////////////////////////////////

/**
.Class.FileFormat:
..cat:Input/Output
..summary:Object that stores a file format.
..signature:FileFormat<File, Data [, Format [, Meta] ]>
..see:Tag.File Format
..include:seqan/file.h
*/

template <
	typename TFile, 
	typename TData,
	typename TMeta,
	typename TFormat = void >
struct FileFormat:
	public FileFormat<TFile, TData, TMeta, void>
{
//IOREV
public:
	typedef typename Size<TData>::Type TSize;

	FileFormat() {}
	FileFormat(FileFormat const &) {}
	~FileFormat() {}
	FileFormat const & operator =(FileFormat const &) { return *this; }

	inline void * 
	formatID_() const
	{
SEQAN_CHECKPOINT
		return ClassIdentifier_<TFormat>::getID();
	}

	virtual void
	read_(TFile & file, TData & data) const
	{
SEQAN_CHECKPOINT
		read(file, data, TFormat());
	}
	virtual void
	read_(TFile & file, TData & data, TSize limit) const
	{
SEQAN_CHECKPOINT
		read(file, data, limit, TFormat());
	}

	virtual void
	readMeta_(TFile & file, TMeta & meta) const
	{
SEQAN_CHECKPOINT
		readMeta(file, meta, TFormat());
	}

	virtual void
	goNext_(TFile & file) const
	{
SEQAN_CHECKPOINT
		goNext(file, TFormat());
	}

	virtual TSize
	length_(TFile & file) const
	{
SEQAN_CHECKPOINT
		length(file, TFormat());
	}

	virtual void
	write_(TFile & file, TData & data) const
	{
SEQAN_CHECKPOINT
		write(file, data, TFormat());
	}
	virtual void
	write_(TFile & file, TData & data, TMeta & meta) const
	{
SEQAN_CHECKPOINT
		write(file, data, meta, TFormat());
	}
};

//____________________________________________________________________________

//base class for all file format classes 

template <typename TFile, typename TData, typename TMeta>
struct FileFormat<TFile, TData, TMeta, void>
{
//IOREV
public:
	typedef typename Size<TData>::Type TSize;

	FileFormat() {}
	FileFormat(FileFormat const &) {}
	~FileFormat() {}
	FileFormat const & operator =(FileFormat const &) { return *this; }

	virtual void *
	formatID_() const = 0;

	virtual void
	read_(TFile & file, TData & data) const = 0;
	virtual void
	read_(TFile & file, TData & data, TSize limit) const = 0;

	virtual void
	readMeta_(TFile & file, TMeta & meta) const = 0;

	virtual void
	goNext_(TFile & file) const = 0;

	virtual TSize
	length_(TFile & file) const = 0;

	virtual void
	write_(TFile & file, TData & data) const = 0;
	virtual void
	write_(TFile & file, TData & data, TMeta & meta) const = 0;

};

//////////////////////////////////////////////////////////////////////////////
// Wrapper for functions to virtuals
//////////////////////////////////////////////////////////////////////////////


template <typename TFile, typename TData, typename TMeta, typename TFormat>
inline void *
formatID(FileFormat<TFile, TData, TMeta, TFormat> const & file_format)
{
//IOREV
SEQAN_CHECKPOINT
	return file_format.formatID_();
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.FileFormat#read:
..cat:Input/Output
..summary:Loads a record from file.
..signature:read(file, data [, meta], format)
..signature:read(file, data [, meta], tag)
..param.file:An input file.
..param.data:A container that gets the data read from $file$.
..param.meta:A container that gets meta data from $file$. (optional)
..param.format:A file format object.
...type:Class.FileFormat.File Format object
..param.tag:A file format tag.
...type:Tag.File Format.File Format tag
..remarks:The result of this operation is stored in $data$.
..remarks:The function leaves $file$ at the position for reading the next record.
..see:Function.assign
..include:seqan/file.h
*/
template <typename TFile, typename TData, typename TMeta, typename TFormat>
inline void
read(TFile & file,
	 TData & data,
	 FileFormat<TFile, TData, TMeta, TFormat> const & file_format)
{
//IOREV
SEQAN_CHECKPOINT
	file_format.read_(file, data);
}

template <typename TFile, typename TData, typename TMeta, typename TFormat, typename TSize>
inline void
read(TFile & file,
	 TData & data,
	 TSize limit,
	 FileFormat<TFile, TData, TMeta, TFormat> const & file_format)
{
//IOREV
SEQAN_CHECKPOINT
	file_format.read_(file, data, limit);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.readMeta:
..cat:Input/Output
..summary:Read meta information from file.
..signature:readMeta(file, meta, file_format)
..param.file:A file that contains data in the format specified by $file_format$.
..param.meta:A data structure that is able to store meta informations stored in $file$.
..param.file_format:A file format.
..returns.param.meta:The meta data read from $file$.
...type:Tag.File Format
..include:seqan/file.h
*/

template <typename TFile, typename TData, typename TMeta, typename TFormat>
inline void
readMeta(TFile & file,
		 TMeta & meta,
		 FileFormat<TFile, TData, TMeta, TFormat> const & file_format)
{
//IOREV
SEQAN_CHECKPOINT
	file_format.readMeta_(file, meta);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.goNext:
..cat:Input/Output
..include:seqan/file.h
*/

template <typename TFile, typename TData, typename TMeta, typename TFormat>
inline void
goNext(TFile & file,
	   FileFormat<TFile, TData, TMeta, TFormat> const & file_format)
{
//IOREV
SEQAN_CHECKPOINT
	file_format.goNext_(file);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.length:
..cat:Input/Output
..include:seqan/file.h
*/

template <typename TFile, typename TData, typename TMeta, typename TFormat>
inline void
length(TFile & file,
	   FileFormat<TFile, TData, TMeta, TFormat> const & file_format)
{
//IOREV
SEQAN_CHECKPOINT
	file_format.length_(file);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.FileFormat#write:
..class:Adaption."std::iostream"
..cat:Input/Output
..summary:Writes to stream.
..signature:write(stream, source)
..signature:write(stream, begin, end)
..param.stream: A stream object.
...type:Adaption."std::iostream"
..param.source: Container that is written to $stream$.
..param.begin: Iterator to the first character of the range.
..param.end: Iterator behind the last character of the range.
..remarks:The content of $source$ is written 'as-is' to $stream$.
..include:seqan/file.h
*/

template <typename TFile, typename TData, typename TMeta, typename TFormat>
inline void
write(TFile & file,
	  TData & data,
	  FileFormat<TFile, TData, TMeta, TFormat> const & file_format)
{
//IOREV
SEQAN_CHECKPOINT
	file_format.write_(file, data);
}
template <typename TFile, typename TData, typename TMeta, typename TFormat>
inline void
write(TFile & file,
	  TData & data,
	  TMeta & meta,
	  FileFormat<TFile, TData, TMeta, TFormat> const & file_format)
{
//IOREV
SEQAN_CHECKPOINT
	file_format.write_(file, data, meta);
}




//////////////////////////////////////////////////////////////////////////////
// Comparison of two FileFormat objects
//////////////////////////////////////////////////////////////////////////////

template <typename TFileLeft, typename TDataLeft, typename TMetaLeft, typename TFormatLeft, typename TFileRight, typename TDataRight, typename TMetaRight, typename TFormatRight>
inline bool
operator == (FileFormat<TFileLeft, TDataLeft, TMetaLeft, TFormatLeft> const & left, 
			 FileFormat<TFileRight, TDataRight, TMetaRight, TFormatRight> const & right)
{
//IOREV
SEQAN_CHECKPOINT
	return formatID(left) == formatID(right);
}

template <typename TFile, typename TData, typename TMeta, typename TFormat, typename TFormat2>
inline bool
operator == (FileFormat<TFile, TData, TMeta, TFormat> const & left, 
			 Tag<TFormat2> const)
{
//IOREV
SEQAN_CHECKPOINT
	return formatID(left) == ClassIdentifier_<Tag<TFormat2> const>::getID();
}

template <typename TFile, typename TData, typename TMeta, typename TFormat, typename TFormat2>
inline bool
operator == (Tag<TFormat2> const,
			 FileFormat<TFile, TData, TMeta, TFormat> const & right)
{
//IOREV
SEQAN_CHECKPOINT
	return ClassIdentifier_<Tag<TFormat2> const>::getID() == formatID(right);
}

//____________________________________________________________________________

template <typename TFileLeft, typename TDataLeft, typename TMetaLeft, typename TFormatLeft, typename TFileRight, typename TDataRight, typename TMetaRight, typename TFormatRight>
inline bool
operator != (FileFormat<TFileLeft, TDataLeft, TMetaLeft, TFormatLeft> const & left, 
			 FileFormat<TFileRight, TDataRight, TMetaRight, TFormatRight> const & right)
{
//IOREV
SEQAN_CHECKPOINT
	return formatID(left) != formatID(right);
}

template <typename TFile, typename TData, typename TMeta, typename TFormat, typename TFormat2>
inline bool
operator != (FileFormat<TFile, TData, TMeta, TFormat> const & left, 
			 Tag<TFormat2> const)
{
//IOREV
SEQAN_CHECKPOINT
	return formatID(left) != ClassIdentifier_<Tag<TFormat2> const>::getID();
}

template <typename TFile, typename TData, typename TMeta, typename TFormat, typename TFormat2>
inline bool
operator != (Tag<TFormat2> const,
			 FileFormat<TFile, TData, TMeta, TFormat> const & right)
{
//IOREV
SEQAN_CHECKPOINT
	return ClassIdentifier_<Tag<TFormat2> const>::getID() != formatID(right);
}

//////////////////////////////////////////////////////////////////////////////
// allgemeine Funktionen fuer Streams
//////////////////////////////////////////////////////////////////////////////
//TODO??? Das muss in eine extra Datei


/*
template <typename TStream, typename TIterator>
inline void
write(TStream & target,
	  TIterator begin_,
	  TIterator end_)
{
	while (begin_ != end_)
	{
		_streamPut(target, convert<char>(*begin_));
		++begin_;
	}
}

//____________________________________________________________________________

template <typename TStream, typename TSource>
inline void
write(TStream & target,
	  TSource const & source)
{
	write(target, begin(source), end(source));
}
//TODO???: Spezialisierungen zum blockweise schreiben bei contiguous strings von char
//Anmerkungen: write wird nach dem zweiten Argument (source) spezialisiert!

//____________________________________________________________________________

template <typename TStream, typename TSource>
inline void
write(TStream & target,
	  TSource const & source,
	  typename Size<TSource>::Type limit_)
{
	if (length(source) > limit_)
	{
		write(target, begin(source), begin(source) + limit_);
	}
	else
	{
		write(target, begin(source), end(source));
	}
}

*/
//////////////////////////////////////////////////////////////////////////////

// Helper function for scanning a stream
// c = next character, pass it to the next call of the function

template <typename TFile, typename TString, typename TChar>
inline void
_streamAppendLine(TFile & file,
				   TString & str,
				   TChar & c)
{
//IOREV _nodoc_ _hasCRef_ wrong place
	while (true)
	{
		if (_streamEOF(file)) break;

		if (c == '\r')
		{
			c = _streamGet(file);
			if (c == '\n') 
			{
				c = _streamGet(file);
			}
			break;
		}
		if (c == '\n')
		{
			c = _streamGet(file);
			break;
		}

		appendValue(str, c);

		c = _streamGet(file);
	}
}
//____________________________________________________________________________

template <typename TFile, typename TChar>
inline void
_streamCountLine(TFile & file,
				  TChar & c)

{
//IOREV _nodoc_ _notused_ _hasCRef_ wrong place
	while (true)
	{
		if (_streamEOF(file)) break;

		if (c == '\r')
		{
			c = _streamGet(file);
			if (c == '\n') 
			{
				c = _streamGet(file);
			}
			break;
		}
		if (c == '\n')
		{
			c = _streamGet(file);
			break;
		}

		c = _streamGet(file);
	}
}

//____________________________________________________________________________

template <typename TFile, typename TChar>
inline typename Size<TFile>::Type
_streamSkipLine(TFile & file,
				 TChar & c)

{
//IOREV _nodoc_ _hasCRef_ wrong place
	typename Size<TFile>::Type count = 0;
	while (true)
	{
		if (_streamEOF(file)) break;

		if (c == '\r')
		{
			c = _streamGet(file);
			if (c == '\n') 
			{
				c = _streamGet(file);
			}
			break;
		}
		if (c == '\n')
		{
			c = _streamGet(file);
			break;
		}

		++count;

		c = _streamGet(file);
	}

	return count;
}




////////////////////////////////////////////////////////////////////////////
//new ones

//new ones for streams
template<typename TFile, typename TChar>
inline void 
_streamSkipWhitespace(TFile& file, TChar& c)
{
//IOREV _nodoc_ _hasCRef_ wrong place; according to POSIX \v\f\r\n are also whitespace
	if ((c!=' ') && (c != '\t')) return;
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if ((c!=' ') && (c != '\t')) break;
	}
}

////////////////////////////////////////////////////////////////////////////

template<typename TChar>
inline bool
_streamIsLetter(TChar const c)
{
//IOREV _nodoc_ wrong place; what about return ( ((c >= 'a') && (c <= 'z')) || ((c >= 'A') && (c <= 'Z')) [better performance] replace with _parseIsLetter()
	return ((c == 'a') || (c == 'b') || (c == 'c') || (c == 'd') || (c == 'e') || 
			(c == 'f') || (c == 'g') || (c == 'h') || (c == 'i') || (c == 'j') ||
			(c == 'k') || (c == 'l') || (c == 'm') || (c == 'n') || (c == 'o') || 
			(c == 'p') || (c == 'q') || (c == 'r') || (c == 's') || (c == 't') ||
			(c == 'u') || (c == 'v') || (c == 'w') || (c == 'x') || (c == 'y') || 
			(c == 'z') || (c == 'A') || (c == 'B') || (c == 'C') || (c == 'D') ||
			(c == 'E') || (c == 'F') || (c == 'G') || (c == 'H') || (c == 'I') || 
			(c == 'J') || (c == 'K') || (c == 'L') || (c == 'M') || (c == 'N') ||
			(c == 'O') || (c == 'P') || (c == 'Q') || (c == 'R') || (c == 'S') || 
			(c == 'T') || (c == 'U') || (c == 'V') || (c == 'W') || (c == 'X') ||
			(c == 'Y') || (c == 'Z'));
}


////////////////////////////////////////////////////////////////////////////


template<typename TFile, typename TChar>
inline String<char>
_streamReadWord(TFile & file, TChar& c)
{
//IOREV _nodoc_ _hasCRef_ wrong place
	// Read word
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_streamIsLetter(c)) break;
		append(str, c);
	}
	return str;
}

//////////////////////////////////////////////////////////////////////////////

//new ones for strings

template <typename TString, typename TIter>
inline typename Size<TString>::Type
_stringSkipLine(TString & str,
				 TIter & it)

{
//IOREV _nodoc_ _bug_ wrong place; possible request for *it but it could be it_end in l 672; iterator akwardness, see head of file
	typename Size<TString>::Type count = 0;
	typename Iterator<TString,Standard>::Type end_it = end(str,Standard());
	while (true)
	{
		if (it == end_it) break;

		if (*it == '\r')
		{
			++it;
			if (*it == '\n') 
			{
				++it;
			}
			break;
		}
		if (*it == '\n')
		{
			++it;
			break;
		}

		++count;
		++it;
	}

	return count;
}

/////////////////////////////////////////////////////////////////////////

template <typename TString1, typename TString2, typename TIter>
inline void
_stringAppendLine(TString1 & str,
				   TString2 & a_str,
				   TIter & it)
{
//IOREV _nodoc_ _bug_ wrong place; possible request for *it but it could be it_end in l 708; iterator akwardness, see head of file
	typename Iterator<TString1,Standard>::Type end_it = end(str,Standard());
	while (true)
	{
		if (it == end_it) break;

		if (*it == '\r')
		{
			++it; 
			if (*it == '\n') 
			{
				++it;
			}
			break;
		}
		if (*it == '\n')
		{
			++it;
			break;
		}

		appendValue(a_str, getValue(it));
		++it;
	}
}

////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TIter>
inline void 
_stringSkipWhitespace(TString& str, TIter& it)
{
//IOREV _nodoc_ _bug_ wrong place; according to POSIX \v\f\r\n are also whitespace; iterator akwardness, see head of file
	typename Iterator<TString,Standard>::Type end_it = end(str,Standard())-1;
	while (it != end_it) {
		if ((*it!=' ') && (*it != '\t')) break;
		++it;
	}
}

////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TIter>
inline int
_stringReadNumber(TString & str, TIter& it)
{
//IOREV _nodoc_ wrong place; iterator akwardness, see head of file
	// Read number
	typename Iterator<TString,Standard>::Type end_it = end(str,Standard())-1;
	String<char> numstr(getValue(it));
	while (it != end_it) {
		++it;
		if (!_parseIsDigit(*it)) break;
		append(numstr, getValue(it));
	}
 	return atoi(toCString(numstr));
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef SEQAN_HEADER_...
