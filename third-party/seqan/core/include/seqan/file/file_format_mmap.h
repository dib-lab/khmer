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
// Provides functionality to map a sequence file into memory and split into
// segments (one for each sequence). Implements different sequence file
// formats, e.g. Fastq, Fasta, QSeq, Raw, and auto detection (via the
// TagSelector AutoSeqFormat) based on content or file extension.
// ==========================================================================

#ifndef SEQAN_HEADER_FILE_FORMAT_MMAP_H
#define SEQAN_HEADER_FILE_FORMAT_MMAP_H



/* IOREV
 * _tested_
 *
 *
 * mostly documented
 * contents seems to be used widely in apps
 * 
 * contains all sorts of private _is* and _stream* functions
 *
 * contains many file format specific functions and even two new file format
 * definitions (fastq and qsec) that should definitely go someplace else,
 *
 * contains almost nothing specific to MemoryMapped strings, other than
 * the split()-calls
 */


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// File Formats - Fastq (Fasta extension for quality values)
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.File Format.tag.Fastq:
    FASTQ file format for sequences.
..include:seqan/file.h
*/

    struct TagFastq_;
    //IOREV _doc_ this should ge somewhere else!
    typedef Tag<TagFastq_> Fastq; //IOREV

/**
.Tag.File Format.tag.QSeq:
	QSeq format, used for most of the Illumina read files.
..include:seqan/file.h
*/
	struct QSeq_;
	//IOREV _doc_ this whole part should go someplace different
	typedef Tag<QSeq_> QSeq; //IOREV


    template <typename TFormat, typename T = void>
    struct FileFormatExtensions;
    

    template <typename T>
    struct FileFormatExtensions<Fasta, T>
    {
        static char const * VALUE[6];
    };

    template <typename T>
    char const * FileFormatExtensions<Fasta, T>::VALUE[6] = {
        ".fa",      // default output extension
        ".fasta",
        ".faa",     // FASTA Amino Acid file
        ".ffn",     // FASTA nucleotide coding regions file
        ".fna",     // FASTA Nucleic Acid file
        ".frn" };


    template <typename T>
    struct FileFormatExtensions<Fastq, T>
    {
        static char const * VALUE[2];
    };

    template <typename T>
    char const * FileFormatExtensions<Fastq, T>::VALUE[2] = {
        ".fq",      // default output extension
        ".fastq" };


    template <typename T>
    struct FileFormatExtensions<QSeq, T>
    {
        static char const * VALUE[2];
    };

    template <typename T>
    char const * FileFormatExtensions<QSeq, T>::VALUE[2] = {
        ".txt",     // default output extension
        ".seq" };


    template <typename T>
    struct FileFormatExtensions<Raw, T>
    {
        static char const * VALUE[1];	// default is one extension
    };

    template <typename T>
    char const * FileFormatExtensions<Raw, T>::VALUE[1] = {
        ".qseq" };  // default output extension

/**
.Shortcut.MultiFasta
..summary:A sequence file mapped in memory as a StringSet of concatenated sequence file fragments.
..cat:Input/Output
..signature:MultiSeqFile
..shortcutfor:Spec.ConcatDirect
..shortcutfor:Spec.MMap String
...signature:StringSet<String<char, MMap<> >, Owner<ConcatDirect<> > >
..status:deprecated, will be removed in favour of @Shortcut.MultiSeqFile@
*/

/**
.Shortcut.MultiSeqFile
..summary:A sequence file mapped in memory as a StringSet of concatenated sequence file fragments.
..cat:Input/Output
..signature:MultiSeqFile
..shortcutfor:Spec.ConcatDirect
..shortcutfor:Spec.MMap String
...signature:StringSet<String<char, MMap<> >, Owner<ConcatDirect<> > >
*/

	// define memory mapped stringset
	typedef StringSet<String<char, MMap<> >, Owner<ConcatDirect<> > >	MultiFasta;	//deprecated (use MultiSeqFile instead) //IOREV _doc_ _delcandidate_
	typedef StringSet<String<char, MMap<> >, Owner<ConcatDirect<> > >	MultiSeqFile; //IOREV _doc_


	template <typename TValue>
	inline bool
	_isLineBreak(TValue value)
	{
//IOREV _nodoc_ will return true twice on Windows (ANSI EOL = "\r\n"); unclear whether \r on its own should be considered EOL
		return (value == '\n' || value == '\r');
	}

	template <typename TIterator>
	inline bool
	_seekLineBreak(TIterator &it, TIterator itEnd)
	{
//IOREV _nodoc_
		while (!_isLineBreak(*it))
			if (++it == itEnd) return false;
		return true;
	}

	template <typename TIterator>
	inline bool
	_seekNonLineBreak(TIterator &it, TIterator itEnd)
	{
//IOREV _bug_ _nodoc_ this returns true on [2] of "\n\n\n" (there is no loop)
		if (*it == '\n')
		{
			if (++it == itEnd) return false;
			if (*it == '\r')
				if (++it == itEnd) return false;
		} else
			if (*it == '\r')
			{
				if (++it == itEnd) return false;
				if (*it == '\n')
					if (++it == itEnd) return false;
			}
		return true;
	}

        // Returns true iff value is a whitespace.
	template <typename TValue>
	inline bool
	_isWhiteSpace(TValue value)
	{
//IOREV according to POSIX \v and \f are also whitespace, see also file_format.h:_streamSkipWhiteSpace(), misc_parsing.h: _parseSkipWhitespace() and misc_parsing.h:_parseReadWordUntilWhitespace()
		return (value == ' ' || value == '\t' || value == '\r' || value == '\n');
	}

        // Increment iterator until end of sequence or *it is a whitespace.
	template <typename TIterator>
	inline bool
	_seekWhiteSpace(TIterator &it, TIterator itEnd)
	{
//IOREV _nodoc_
		while (!_isWhiteSpace(*it))
			if (++it == itEnd) return false;
		return true;
	}

	template <typename TIterator>
	inline bool
	_seekTab(TIterator& it, TIterator itEnd)
	{
//IOREV _nodoc_
		for (; it != itEnd; ++it)
			if (*it == '\t') return true;
		return false;
	}


//////////////////////////////////////////////////////////////////////////////
// File Formats - Fasta
//////////////////////////////////////////////////////////////////////////////

/**
.Function.guessFormat:
..summary:Guesses a file format from the contents of a sequence file.
..cat:Input/Output
..signature:guessFormat(text, formatTag)
..param.text:A string storing the contents of a sequence file.
...see:Class.String
..param.formatTag:A file format tag.
...type:Tag.File Format
...type:Class.AutoSeqFormat
..returns:$true$ if the format represented by $formatTag$ was recognized in $text$.
..see:Function.guessFormatFromFilename
..include:seqan/file.h
*/

	// test for Fasta format
	template < typename TSeq >
	inline bool
	guessFormat(
		TSeq const & seq,
		Fasta)
	{
//IOREV _doc_ in competition with file_format_guess.h:guessFileFormat()
		return seq[0] == '>';
	}
	
/**
.Function.guessFormatFromFilename:
..summary:Guesses a file format from a sequence file name.
..cat:Input/Output
..signature:guessFormatFromFilename(fileName, formatTag)
..param.fileName:A filename of a sequence file.
...see:Class.String
..param.formatTag:A file format tag.
...type:Tag.File Format
...type:Class.AutoSeqFormat
..returns:$true$ if the format represented by $formatTag$ was recognized in $fileName$.
..see:Function.guessFormatFromFilename
..include:seqan/file.h
*/

	template <typename TFilename, typename TFormat_>
	inline bool
	guessFormatFromFilename(
		TFilename const & fileName,
		Tag<TFormat_> /*formatTag*/)
	{
		typedef typename Value<TFilename>::Type                                 TValue;
		typedef ModifiedString<TFilename const, ModView<FunctorLowcase<TValue> > >	TLowcase;
		typedef Tag<TFormat_>                                                   TFormat;
		
		TLowcase lowcaseFileName(fileName);
		for (unsigned i = 0; i < sizeof(FileFormatExtensions<TFormat>::VALUE) / sizeof(char*); ++i)
			if (endsWith(lowcaseFileName, FileFormatExtensions<TFormat>::VALUE[i]))
				return true;

		return false;
	}

    template <typename TStringSet, typename TFormat_>
    inline void
    _getFileFormatExtensions(TStringSet &stringSet, Tag<TFormat_> /*formatTag*/)
    {
		typedef Tag<TFormat_> TFormat;
		for (unsigned i = 0; i < sizeof(FileFormatExtensions<TFormat>::VALUE) / sizeof(char*); ++i)
			appendValue(stringSet, FileFormatExtensions<TFormat>::VALUE[i]);
    }

    template <typename TStringSet, typename TTag>
    inline void
    _getFileFormatExtensions(TStringSet &stringSet, TagList<TTag, void> const /*formatTag*/)
    {
        _getFileFormatExtensions(stringSet, TTag());
    }

    template <typename TStringSet, typename TTag, typename TSubList>
    inline void
    _getFileFormatExtensions(TStringSet &stringSet, TagList<TTag, TSubList> const /*formatTag*/)
    {
        _getFileFormatExtensions(stringSet, TTag());
        _getFileFormatExtensions(stringSet, TSubList());
    }

    template <typename TStringSet, typename TTagList>
    inline void
    _getFileFormatExtensions(TStringSet &stringSet, TagSelector<TTagList> const /*formatTag*/)
    {
        _getFileFormatExtensions(stringSet, TTagList());
    }

/**
.Function.split:
..summary:Divides the contents of a sequence file into sequence file fragments separated by a file format specific delimiter.
..cat:Input/Output
..signature:split(stringSet, formatTag)
..param.stringSet:A @Spec.ConcatDirect@ StringSet. The concat member (concatenation string) contains the contents of a sequence file.
...type:Spec.ConcatDirect
...type:Shortcut.MultiSeqFile
..param.formatTag:A file format tag.
...type:Tag.File Format
...type:Class.AutoSeqFormat
..remarks:The @Memvar.ConcatDirect#concat@ member should contain the contents of the sequence file by a prior call of @Function.open@.
..remarks:This function expects a @Spec.ConcatDirect@ StringSet and divides the underlying concatenation string into
sequence fragments separated by a file format specific delimiter.
After calling this function, the StringSet length is the number of sequence fragments and each fragment can be retrieved by @Function.value@ or @Function.getValue@.
..see:Function.guessFormat
..include:seqan/file.h
*/

    // split stringset into single Fasta sequences
    template <typename TValue, typename TSpec, typename TStringSetSpec>
    inline void split(StringSet<String<TValue, TSpec>, TStringSetSpec> &me, Fasta /* tag */)
    {
//IOREV _doc_
        typedef String<TValue, TSpec>                               TString;
        typedef typename Iterator<TString const, Standard>::Type	TIterator;

        clear(me.limits);

        TIterator itBeg = begin(me.concat, Standard());
        TIterator itEnd = end(me.concat, Standard());
        bool newLine = true;
        for (TIterator it = itBeg; it != itEnd; ++it)
        {
            TValue c = *it;
            if (newLine && c == '>')
                appendValue(me.limits, it - itBeg, Generous());
            newLine = _isLineBreak(c);
        }
        if (empty(me.limits))
            appendValue(me.limits, 0);
        appendValue(me.limits, itEnd - itBeg);
    }

/**
.Function.assignSeq:
..summary:Extracts the sequence part of a sequence file fragment.
..cat:Input/Output
..signature:assignSeq(sequence, seqFragment, formatTag)
..param.sequence:The resulting sequence of the fragment.
...type:Class.String
..param.seqFragment:A sequence file fragment.
...type:Class.String
..param.formatTag:A file format tag.
...type:Tag.File Format
...type:Class.AutoSeqFormat
..remarks:After calling @Function.split@ on a @Spec.ConcatDirect@ StringSet to divide a file into fragments, 
this function can be used to extract the sequence of every fragment in the StringSet.
..see:Function.split
..include:seqan/file.h
*/

	template <typename TSeq, typename TFastaSeq>
	inline void
	assignSeq(
		TSeq & dst,
		TFastaSeq const & fasta,
		Fasta)
	{
//IOREV _doc_
		typedef typename Iterator<TFastaSeq const, Standard>::Type	TIterator;
		typedef typename Iterator<TSeq, Standard>::Type				TDstIterator;

		TIterator it = begin(fasta, Standard());
		TIterator itEnd = end(fasta, Standard());

		clear(dst);
		
		// skip Fasta id
		if (it == itEnd) return;
		if (*it == '>')
		{
			if (!_seekLineBreak(it, itEnd)) return;
			if (!_seekNonLineBreak(it, itEnd)) return;
		}

		// copy sequence
		resize(dst, itEnd - it);		
		TDstIterator dit = begin(dst, Standard());
		for (; it != itEnd; ++it)
			if (!_isLineBreak(*it))
			{
				*dit = *it;
				++dit;
			}
		resize(dst, dit - begin(dst, Standard()));
	}

/**
.Function.assignSeqId:
..summary:Extracts the sequence id of a sequence file fragment.
..cat:Input/Output
..signature:assignSeqId(id, seqFragment, formatTag)
..param.id:The resulting sequence id of the fragment (e.g. Fasta Id).
...type:Shortcut.CharString
..param.seqFragment:A sequence file fragment.
...type:Class.String
..param.formatTag:A file format tag.
...type:Tag.File Format
...type:Class.AutoSeqFormat
..remarks:After calling @Function.split@ on a @Spec.ConcatDirect@ StringSet to divide a file into fragments, 
this function can be used to extract the sequence id of every fragment in the StringSet.
..see:Function.split
..include:seqan/file.h
*/

	template <typename TSeq, typename TFastaSeq>
	inline void
	assignSeqId(
		TSeq & dst,
		TFastaSeq const & fasta,
		Fasta)
	{
//IOREV _doc_
		typedef typename Iterator<TFastaSeq const, Standard>::Type	TIterator;

		TIterator itBeg = begin(fasta, Standard());
		TIterator itEnd = end(fasta, Standard());
		TIterator it = itBeg;
		
		clear(dst);
		if (it == itEnd) return;
		if (*it == '>')
		{
			_seekLineBreak(it, itEnd);
			assign(dst, infix(fasta, 1, it - itBeg));
		}
	}

/**
.Function.assignCroppedSeqId:
..summary:Extracts the sequence id up to the first whitespace of a sequence file fragment.
..cat:Input/Output
..signature:assignCroppedSeqId(id, seqFragment, formatTag)
..param.id:The resulting cropped sequence id of the fragment (e.g. Fasta Id).
...note:The resulting id contains no whitespaces.
...type:Shortcut.CharString
..param.seqFragment:A sequence file fragment.
...type:Class.String
..param.formatTag:A file format tag.
...type:Tag.File Format
...type:Class.AutoSeqFormat
..remarks:After calling @Function.split@ on a @Spec.ConcatDirect@ StringSet to divide a file into fragments, 
this function can be used to extract the sequence id up to the first whitespace of every fragment in the StringSet.
..see:Function.split
..include:seqan/file.h
*/

	// Assign sequence id up to first whitespace.
	template <typename TSeq, typename TFastaSeq>
	inline void
	assignCroppedSeqId(
		TSeq & dst,
		TFastaSeq const & fasta,
		Fasta)
	{
//IOREV _doc_
		typedef typename Iterator<TFastaSeq const, Standard>::Type	TIterator;

		TIterator itBeg = begin(fasta, Standard());
		TIterator itEnd = end(fasta, Standard());
		TIterator it = itBeg;
		
		clear(dst);
		if (it == itEnd) return;
		if (*it == '>')
		{
			_seekWhiteSpace(it, itEnd);
			assign(dst, infix(fasta, 1, it - itBeg));
		}
	}

/**
.Function.assignQual:
..summary:Extracts the quality values of a sequence file fragment.
..cat:Input/Output
..signature:assignQual(qualities, seqFragment, formatTag)
..param.qualities:The resulting quality values encoded in ASCII.
...remarks:The quality values are encoded in ASCII and must be manually converted into zero-based values.
...type:Shortcut.CharString
..param.seqFragment:A sequence file fragment.
...type:Class.String
..param.formatTag:A file format tag.
...type:Tag.File Format
...type:Class.AutoSeqFormat
..remarks:After calling @Function.split@ on a @Spec.ConcatDirect@ StringSet to divide a file into fragments, 
this function can be used to extract the sequence quality values of every fragment in the StringSet.
..see:Function.split
..include:seqan/file.h
*/

  	template <typename TSeq, typename TFastaSeq>
	inline void
	assignQual(
		TSeq & dst,
		TFastaSeq const &,
		Fasta)
	{
//IOREV _doc_
		clear(dst);
	}
	
/**
.Function.assignQualId:
..summary:Extracts the quality value id of a sequence file fragment.
..cat:Input/Output
..signature:assignQualId(id, seqFragment, formatTag)
..param.id:The resulting quality value id of a sequence (e.g. Fastq Quality Id).
...type:Shortcut.CharString
..param.seqFragment:A sequence file fragment.
...type:Class.String
..param.formatTag:A file format tag.
...type:Tag.File Format
...type:Class.AutoSeqFormat
..remarks:After calling @Function.split@ on a @Spec.ConcatDirect@ StringSet to divide a file into fragments, 
this function can be used to extract the quality value id of every fragment in the StringSet.
..see:Function.split
..include:seqan/file.h
*/

	template <typename TSeq, typename TFastaSeq>
	inline void
	assignQualId(
		TSeq & dst,
		TFastaSeq const &,
		Fasta)
	{
//IOREV _doc_
		clear(dst);
	}


	// test for Fastq format
	template < typename TSeq >
	inline bool
	guessFormat(
		TSeq const & seq,
		Fastq)
	{
//IOREV _doc_
		return seq[0] == '@';
	}

    // split stringset into single Fasta sequences
    template <typename TValue, typename TSpec, typename TStringSetSpec>
    inline void
    split(StringSet<String<TValue, TSpec>, TStringSetSpec> & me, Fastq /* tag */)
    {
        //IOREV
        typedef String<TValue, TSpec>                               TString;
        typedef typename Iterator<TString const, Standard>::Type	TIterator;

        clear(me.limits);

        TIterator itBeg = begin(me.concat, Standard());
        TIterator itEnd = end(me.concat, Standard());
        bool newLine = true;
        for (TIterator it = itBeg; it != itEnd; ++it)
        {
            if (newLine && *it == '@')
                appendValue(me.limits, it - itBeg, Generous());
            if (newLine && *it == '+')
            {
                // skip qualitity fasta id
                if (!_seekLineBreak(it, itEnd)) break;
                if (!_seekNonLineBreak(it, itEnd)) break;
                // skip qualitity values
                if (!_seekLineBreak(it, itEnd)) break;
            }
            newLine = _isLineBreak(*it);
        }
        if (empty(me.limits))
            appendValue(me.limits, 0);
        appendValue(me.limits, itEnd - itBeg);
    }
    

	template <typename TSeq, typename TFastaSeq>
	inline void
	assignSeq(
		TSeq & dst,
		TFastaSeq const & fasta,
		Fastq)
	{
//IOREV _doc_
		typedef typename Iterator<TFastaSeq const, Standard>::Type	TIterator;
		typedef typename Iterator<TSeq, Standard>::Type				TDstIterator;

		TIterator it = begin(fasta, Standard());
		TIterator itEnd = end(fasta, Standard());

		clear(dst);
		
		// skip Fasta id
		if (it == itEnd) return;
		if (*it == '@')
		{
			if (!_seekLineBreak(it, itEnd)) return;
			if (!_seekNonLineBreak(it, itEnd)) return;
		}

		// copy sequence
		resize(dst, itEnd - it);		
		TDstIterator dit = begin(dst, Standard());
		for (; it != itEnd; ++it) 
		{
			if (_isLineBreak(*it))
			{
				if (!_seekNonLineBreak(it, itEnd)) break;
				if (*it == '+') break;
			}
			*dit = *it;
			++dit;
		}
		resize(dst, dit - begin(dst, Standard()));
	}

	template <typename TSeq, typename TFastaSeq>
	inline void
	assignSeqId(
		TSeq & dst,
		TFastaSeq const & fasta,
		Fastq)
	{
//IOREV _doc_
		typedef typename Iterator<TFastaSeq const, Standard>::Type	TIterator;

		TIterator itBeg = begin(fasta, Standard());
		TIterator itEnd = end(fasta, Standard());
		TIterator it = itBeg;
		
		clear(dst);
		if (it == itEnd) return;
		if (*it == '@')
		{
			_seekLineBreak(it, itEnd);
			assign(dst, infix(fasta, 1, it - itBeg));
		}
	}

        // Assign sequence id up to first whitespace.
	template <typename TSeq, typename TFastaSeq>
	inline void
	assignCroppedSeqId(
		TSeq & dst,
		TFastaSeq const & fasta,
		Fastq)
	{
//IOREV _doc_ probably some more code sharing between these functions possible
		typedef typename Iterator<TFastaSeq const, Standard>::Type	TIterator;

		TIterator itBeg = begin(fasta, Standard());
		TIterator itEnd = end(fasta, Standard());
		TIterator it = itBeg;
		
		clear(dst);
		if (it == itEnd) return;
		if (*it == '@')
		{
			_seekWhiteSpace(it, itEnd);
			assign(dst, infix(fasta, 1, it - itBeg));
		}
	}

	template <typename TSeq, typename TFastaSeq>
	inline void
	assignQual(
		TSeq & dst,
		TFastaSeq const & fasta,
		Fastq)
	{
//IOREV _doc_
		typedef typename Iterator<TFastaSeq const, Standard>::Type	TIterator;
		typedef typename Iterator<TSeq, Standard>::Type				TDstIterator;

		TIterator it = begin(fasta, Standard());
		TIterator itEnd = end(fasta, Standard());

		clear(dst);
		
		if (it == itEnd) return;
		if (*it == '@')
		{
			// seek quality id
			do {
				if (!_seekLineBreak(it, itEnd)) return;
				if (!_seekNonLineBreak(it, itEnd)) return;
			} while (*it != '+');

			// skip quality id
			if (!_seekLineBreak(it, itEnd)) return;
			if (!_seekNonLineBreak(it, itEnd)) return;

			// copy sequence
			resize(dst, itEnd - it);		
			TDstIterator dit = begin(dst, Standard());
			for (; it != itEnd; ++it)
				if (!_isLineBreak(*it))
				{
					*dit = *it;
					++dit;
				}
			resize(dst, dit - begin(dst, Standard()));
		}
	}

	template <typename TSeq, typename TFastaSeq>
	inline void
	assignQualId(
		TSeq & dst,
		TFastaSeq const & fasta,
		Fastq)
	{
//IOREV _doc_
		typedef typename Iterator<TFastaSeq const, Standard>::Type	TIterator;

		TIterator itBeg = begin(fasta, Standard());
		TIterator itEnd = end(fasta, Standard());
		TIterator it1 = itBeg;
		
		clear(dst);
		if (it1 == itEnd) return;
		if (*it1 == '@')
		{
			do {
				if (!_seekLineBreak(it1, itEnd)) return;
				if (!_seekNonLineBreak(it1, itEnd)) return;
			} while (*it1 != '+');
			TIterator it2 = it1;
			_seekLineBreak(it2, itEnd);
			assign(dst, infix(fasta, (it1 - itBeg) + 1, it2 - itBeg));
		}
	}

//////////////////////////////////////////////////////////////////////////////
// File Formats - QSeq (used by Illumina for most of their read files)
//////////////////////////////////////////////////////////////////////////////

	// FIXME The following enum is more or less arbitrary since the information
	// in a QSeq file may differ depending on where they come from. Not sure if
	// this is something that needs to be fixed here, rather than in the Illu-
	// mina pipeline itself.
	//
	// Also, this enum is quite convoluted but I don't feel like spilling a lot
	// of common symbols (e.g. 'X', 'Y') into the SeqAn namespace.
	struct QSeqEntry {
//IOREV _bug_ has FIXME in code comments
		enum {
			MachineName,
			Run,
			Lane,
			Tile,
			X,
			Y,
			Index,
			Read,
			Sequence,
			Quality,
			Filter
		};
	};

	template < typename TString >
	inline bool _isQSeqFile(TString const& filename)
	{
//IOREV _doc_ has some comments in code
        unsigned int const namelen = 19;
        unsigned int const pathlen = length(filename);
        if (pathlen < namelen) return false;
		::std::string str;
		assign(str, suffix(filename, pathlen - namelen));
		::std::istringstream is(str);
		unsigned int num;
		// Format (as regex): /^s_\d_\d_\d{4}_qseq.txt$/
		return is.get() == 's' && is.get() == '_' &&
			isdigit(is.get()) && is.get() == '_' &&
			isdigit(is.get()) && is.get() == '_' &&
			is >> num &&
			getline(is, str) && str == "_qseq.txt" && is.eof();
	}

//	// Needed?
//	template < typename TString >
//	inline TString _getFirstFile(
//		char const* dirname,
//		QSeq)
//	{
//		Directory dir(dirname);
//		
//		for ( ; dir; ++dir)
//			if (_isQSeqFile(*dir))
//				return *dir;
//	}
//
//	template < typename TString >
//	inline String< TString >
//	_getAllFiles(
//		char const* dirname,
//		QSeq)
//	{
//		String< TString > ret;
//		Directory dir(dirname);
//
//		for (; dir; ++dir)
//			if (_isQSeqFile(*dir))
//				appendValue(ret, *dir);
//
//		return ret;
//	}

	// test for QSeq format
	template < typename TSeq >
	inline bool
	guessFormat(
		TSeq const & seq,
		QSeq)
	{
//IOREV _nodoc_ code comments suggest that it is guessing
		typedef typename Iterator<TSeq const>::Type TIter;
		TIter front = begin(seq);
		TIter const back = end(seq);
		if (!_seekTab(front, back)) return false;
		::std::string token_base;
		assign(token_base, seq);
		::std::istringstream is(token_base);
		::std::string mname;
		unsigned int numval;
		// Actual information encoded in qseq file may vary. Take a few guesses:
		//     machine name, run number, lane number, tile number
		is >> mname >> numval >> numval >> numval;
		return is.good();
	}
	
	template < typename TFilename >
	inline bool
	guessFormatFromFilename(
		TFilename const & fname,
		QSeq)
	{
//IOREV _nodoc_ need proper QSeq documentation
		// QSeq files come in a variety of ways throughout the Gerald pipeline.
		// In this simplest case, "sorted.txt" is a file in a fragment genome
		// directory, each corresponding to 10MB worth of DNA.
		static CharString const standalone_name = "sorted.txt";
		typename Size<TFilename>::Type len = length(fname);
        bool const bLen = len >= length(standalone_name);
        bool const bSuff = suffix(fname, len - length(standalone_name)) == standalone_name;
		return (bLen && bSuff) || _isQSeqFile(fname);
	}

	// split stringset into single QSeq sequences
	template <typename TValue, typename TSpec, typename TStringSetSpec>
	inline void
	split(StringSet<String<TValue, TSpec>, TStringSetSpec> & me, QSeq const & /* tag */)
	{
//IOREV
		typedef String<TValue, TSpec>                               TString;
		typedef typename Iterator<TString const, Standard>::Type	TIterator;

		clear(me.limits);

		TIterator const front = begin(me.concat, Standard());
		TIterator const back = end(me.concat, Standard());

		appendValue(me.limits, 0, Generous());
		for (TIterator i = front; i != back; ++i)
			if (_isLineBreak(*i))
				appendValue(me.limits, i - front, Generous());

		if (!_isLineBreak(*(back - 1))) // Ignore final line break.
			appendValue(me.limits, back - front);
	}

	template <typename TSequence, typename TSource>
	void assignQSeqEntry(
		TSequence& destination,
		TSource const& source,
		unsigned int entry
	) {
//IOREV _nodoc_
		typedef typename Iterator<TSource const>::Type TIterator;

		TIterator const front = begin(source, Standard());
		TIterator const back = end(source, Standard());

		TIterator infixStart = front;

		for (unsigned int i = QSeqEntry::MachineName; i < entry; ++i) {
			_seekTab(infixStart, back);
			++infixStart;
		}

		SEQAN_ASSERT(infixStart != back);

		TIterator infixEnd = infixStart + 1;
		_seekTab(infixEnd, back);

		assign(destination, infix(source, infixStart - front, infixEnd - front));
	}

	template <typename TSeq, typename TQSeqSeq>
	inline void
	assignSeq(
		TSeq & dst,
		TQSeqSeq const & fasta,
		QSeq)
	{
//IOREV _nodoc_ why fasta input?
		assignQSeqEntry(dst, fasta, QSeqEntry::Sequence);
	}

	template <typename TSeq, typename TQSeqSeq>
	inline void
	assignSeqId(
		TSeq & dst,
		TQSeqSeq const & fasta,
		QSeq)
	{
//IOREV _nodoc_ why fasta input?
		// For now: just return the whole line.
		typename Position<TQSeqSeq const>::Type front = 0;
		while (_isLineBreak(fasta[front]))
			++front;
		assign(dst, infix(fasta, front, length(fasta)));
	}

        // Assign sequence id up to first whitespace.
	template <typename TSeq, typename TQSeqSeq>
	inline void
	assignCroppedSeqId(
		TSeq & dst,
		TQSeqSeq const & fasta,
		QSeq)
	{
//IOREV _nodoc_ why fasta input?
		// For now: just return the whole line.
		typename Position<TQSeqSeq const>::Type front = 0;
		while (_isWhiteSpace(fasta[front]))
			++front;
		assign(dst, infix(fasta, front, length(fasta)));
	}

	template <typename TSeq, typename TQSeqSeq>
	inline void
	assignQual(
		TSeq & dst,
		TQSeqSeq const & fasta,
		QSeq)
	{
//IOREV _nodoc_ why fasta input
		assignQSeqEntry(dst, fasta, QSeqEntry::Quality);
	}

	template <typename TSeq, typename TQSeqSeq>
	inline void
	assignQualId(
		TSeq & dst,
		TQSeqSeq const & fasta,
		QSeq)
	{
//IOREV _nodoc_ why fasta input
		assignSeqId(dst, fasta, QSeq());
	}

//////////////////////////////////////////////////////////////////////////////
// File Formats - Raw (multiple sequences, separated by line breaks)
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.File Format.tag.Raw:
	Raw file format for sequences.
..include:seqan/file.h
*/
struct TagRaw_;
//IOREV _doc_ whats this doing here?
typedef Tag<TagRaw_> Raw; //IOREV

	// test for Fastq format
	template < typename TSeq >
	inline bool
	guessFormat(
		TSeq const &,
		Raw)
	{
//IOREV _doc_
		return true;
	}
	
    // split stringset into single Fasta sequences
    template <typename TValue, typename TSpec, typename TStringSetSpec>
    inline void split(StringSet<String<TValue, TSpec>, TStringSetSpec> & me, Raw /* tag */)
    {
//IOREV _doc_
        typedef String<TValue, TSpec>                               TString;
        typedef typename Iterator<TString const, Standard>::Type	TIterator;

        clear(me.limits);

        TIterator itBeg = begin(me.concat, Standard());
        TIterator itEnd = end(me.concat, Standard());
        bool newLine = true;
        for (TIterator it = itBeg; it != itEnd; ++it)
        {
            if (newLine)
                appendValue(me.limits, it - itBeg, Generous());
            newLine = _isLineBreak(*it);
        }
        appendValue(me.limits, itEnd - itBeg);
    }

	template <typename TSeq, typename TRawSeq>
	inline void
	assignSeq(
		TSeq & dst,
		TRawSeq const & fasta,
		Raw)
	{
//IOREV _doc_
		typedef typename Iterator<TRawSeq const, Standard>::Type	TIterator;
		typedef typename Iterator<TSeq, Standard>::Type				TDstIterator;

		TIterator it = begin(fasta, Standard());
		TIterator itEnd = end(fasta, Standard());

		clear(dst);
		if (it == itEnd) return;

		// copy sequence
		resize(dst, itEnd - it);		
		TDstIterator dit = begin(dst, Standard());
		for (; it != itEnd; ++it) 
		{
			if (_isLineBreak(*it)) continue;
			*dit = *it;
			++dit;
		}
		resize(dst, dit - begin(dst, Standard()));
	}
	
  	template <typename TSeq, typename TRawSeq>
	inline void
	assignSeqId(
		TSeq & dst,
		TRawSeq const &,
		Raw)
	{
//IOREV _doc_
		clear(dst);
	}
	
  	template <typename TSeq, typename TRawSeq>
	inline void
	assignCroppedSeqId(
		TSeq & dst,
		TRawSeq const &,
		Raw)
	{
//IOREV _doc_
		clear(dst);
	}
	
  	template <typename TSeq, typename TRawSeq>
	inline void
	assignQual(
		TSeq & dst,
		TRawSeq const &,
		Raw)
	{
//IOREV _doc_
		clear(dst);
	}
	
  	template <typename TSeq, typename TRawSeq>
	inline void
	assignQualId(
		TSeq & dst,
		TRawSeq const &,
		Raw)
	{
//IOREV _doc_
		clear(dst);
	}
	

//////////////////////////////////////////////////////////////////////////////
// File Formats - Auto-Format
//////////////////////////////////////////////////////////////////////////////

/**
.Class.AutoSeqFormat
..summary:Auto-detects and stores a file format.
..cat:Input/Output
..general:Class.TagSelector
..signature:AutoSeqFormat
..remarks:Currently, it is defined as $TagSelector<SeqFormats>$, with:
...code:
	typedef
		TagList<Fastq,
		TagList<Fasta,
		TagList<QSeq,
		TagList<Raw> > > > 						SeqFormats;
..include:seqan/file.h
*/

	typedef
		TagList<Fastq,
		TagList<Fasta,
//		TagList<QSeq,   // doesn't work as it uses STL strings and parsers
		TagList<Raw> /* > */ > > SeqFormats;  // if TagSelector is set to -1, the file format is auto-detected

	typedef TagSelector<SeqFormats> AutoSeqFormat;

//____________________________________________________________________________
// guess file format

	template <typename TFileSeq>
	inline bool
	guessFormat(TFileSeq const &, TagSelector<> &)
	{
        // we get here if the file format could not be determined
		return false;
	}
	
	template <typename TFileSeq, typename TTagList>
	inline bool
	guessFormat(TFileSeq const &seq, TagSelector<TTagList> &format)
	{
//IOREV _doc_ as mentionened this recursive method is not intuitive, also inlining doesn't really make sense for recursive functions, does it?
//(weese:) ANSWER: This is not a recursive function as TTagList is different for each called instance.
//                 It will be expanded (inlined) completely.

        typedef typename TTagList::Type TFormatTag;

        if (value(format) == -1 || value(format) == LENGTH<TTagList>::VALUE - 1)
        {
            // if tagId is set to -1 (auto-detect) or the current format (TFormatTag) then test for TFormatTag format
            if (guessFormat(seq, TFormatTag()))
            {
                value(format) = LENGTH<TTagList>::VALUE - 1;
                return true;
            }
        }
		return guessFormat(seq, static_cast<typename TagSelector<TTagList>::Base &>(format));
	}
	
//____________________________________________________________________________
// guess file format from filename

	template <typename TFilename>
	inline bool
	guessFormatFromFilename(TFilename const &, TagSelector<>)
	{
        // we get here if the file format could not be determined
		return false;
	}
	
	template <typename TFilename, typename TTagList>
	inline bool
	guessFormatFromFilename(TFilename const &fname, TagSelector<TTagList> &format)
	{
        typedef typename TTagList::Type TFormatTag;

        if (value(format) == -1 || value(format) == LENGTH<TTagList>::VALUE - 1)
        {
            // if tagId is set to -1 (auto-detect) or the current format (TFormatTag) then test for TFormatTag format
            if (guessFormatFromFilename(fname, TFormatTag()))
            {
                value(format) = LENGTH<TTagList>::VALUE - 1;
                return true;
            }
        }
		return guessFormatFromFilename(fname, static_cast<typename TagSelector<TTagList>::Base &>(format));
	}

//____________________________________________________________________________
// split stringset into single sequences

    template <typename TValue, typename TSpec, typename TStringSetSpec>
    inline void split(StringSet<String<TValue, TSpec>, TStringSetSpec> &,
                      TagSelector<void> const & /* tag */)
    {
//IOREV _doc_
    }

    template <typename TValue, typename TSpec, typename TStringSetSpec, typename TTagList>
    inline void split(StringSet<String<TValue, TSpec>, TStringSetSpec> & me,
                      TagSelector<TTagList> const & format)
    {
//IOREV _doc_
        typedef typename TTagList::Type TFormatTag;

		if (value(format) == LENGTH<TTagList>::VALUE - 1)
            split(me, TFormatTag());
        else
            split(me, static_cast<typename TagSelector<TTagList>::Base const &>(format));
    }
    

//____________________________________________________________________________
// assignSeq

	template <typename TSeq, typename TFileSeq>
	inline void
	assignSeq(
		TSeq &,
		TFileSeq const &,
		TagSelector<> const &)
	{
//IOREV _doc_
	}

	template <typename TSeq, typename TFileSeq, typename TTagList>
	inline void
	assignSeq(
		TSeq & dst,
		TFileSeq const & seq,
		TagSelector<TTagList> const &format)
	{
//IOREV _doc_
        typedef typename TTagList::Type TFormatTag;

		if (value(format) == LENGTH<TTagList>::VALUE - 1)
			assignSeq(dst, seq, TFormatTag());
		else
			assignSeq(dst, seq, static_cast<typename TagSelector<TTagList>::Base const &>(format));
	}

//____________________________________________________________________________
// assignSeqId

	template <typename TSeqId, typename TFileSeq>
	inline void
	assignSeqId(
		TSeqId &,
		TFileSeq const &,
		TagSelector<> const &)
	{
//IOREV _doc_
	}

	template <typename TSeqId, typename TFileSeq, typename TTagList>
	inline void
	assignSeqId(
		TSeqId & dst,
		TFileSeq const & seq,
		TagSelector<TTagList> const &format)
	{
//IOREV _doc_
        typedef typename TTagList::Type TFormatTag;

		if (value(format) == LENGTH<TTagList>::VALUE - 1)
			assignSeqId(dst, seq, TFormatTag());
		else
			assignSeqId(dst, seq, static_cast<typename TagSelector<TTagList>::Base const &>(format));
	}

//____________________________________________________________________________
// assignCroppedSeqId

	template <typename TSeqId, typename TFileSeq>
	inline void
	assignCroppedSeqId(
		TSeqId &,
		TFileSeq const &,
		TagSelector<> const &)
	{
//IOREV _doc_
	}

	template <typename TSeqId, typename TFileSeq, typename TTagList>
	inline void
	assignCroppedSeqId(
		TSeqId & dst,
		TFileSeq const & seq,
		TagSelector<TTagList> const &format)
	{
//IOREV _doc_
        typedef typename TTagList::Type TFormatTag;

		if (value(format) == LENGTH<TTagList>::VALUE - 1)
			assignCroppedSeqId(dst, seq, TFormatTag());
		else
			assignCroppedSeqId(dst, seq, static_cast<typename TagSelector<TTagList>::Base const &>(format));
	}

//____________________________________________________________________________
// assignQual

	template <typename TSeq, typename TFileSeq>
	inline void
	assignQual(
		TSeq &,
		TFileSeq const &,
		TagSelector<> const &)
	{
//IOREV _doc_
	}
	
	template <typename TSeq, typename TFileSeq, typename TTagList>
	inline void
	assignQual(
		TSeq & dst,
		TFileSeq const & seq,
		TagSelector<TTagList> const &format)
	{
//IOREV _doc_
        typedef typename TTagList::Type TFormatTag;

		if (value(format) == LENGTH<TTagList>::VALUE - 1)
			assignQual(dst, seq, TFormatTag());
		else
			assignQual(dst, seq, static_cast<typename TagSelector<TTagList>::Base const &>(format));
	}

//____________________________________________________________________________
// assignQualId

	template <typename TSeq, typename TFileSeq>
	inline void
	assignQualId(
		TSeq &,
		TFileSeq const &,
		TagSelector<> const &)
	{
//IOREV _doc_
	}
	
	template <typename TSeq, typename TFileSeq, typename TTagList>
	inline void
	assignQualId(
		TSeq & dst,
		TFileSeq const & seq,
		TagSelector<TTagList> const &format)
	{
//IOREV _doc_
        typedef typename TTagList::Type TFormatTag;

		if (value(format) == LENGTH<TTagList>::VALUE - 1)
			assignQualId(dst, seq, TFormatTag());
		else
			assignQualId(dst, seq, static_cast<typename TagSelector<TTagList>::Base const &>(format));
	}

//////////////////////////////////////////////////////////////////////////////
// Directory import
//////////////////////////////////////////////////////////////////////////////

/**
.Function.appendSeqs:
..summary:Appends all sequences stored in files of directory to a StringSet.
..cat:Input/Output
..signature:appendSeqs(seqSet, dirName, formatTag)
..param.seqSet:A @Class.StringSet@ of sequences to append to.
...type:Class.StringSet
..param.dirName:A path to a directory or single file.
...type:Class.String
..param.formatTag:A file format tag.
...type:Tag.File Format
...type:Class.AutoSeqFormat
..remarks:This function scans a directory and searches for filenames corresponding to the sequence format store in $formatTag$, 
opens them and append their contained sequences to the $seqSet$.
If $formatTag$ is a @Class.AutoSeqFormat@ object, the file format is set to the first known sequence format guessed from a file name.
..see:Function.assignSeq
..include:seqan/file.h
*/

	template <typename TSeqSet, typename TFilename, typename TSeqFormat>
	inline void
	appendSeqs(
		TSeqSet &seqSet,
		TFilename &dirname,
		TSeqFormat format)
	{
//IOREV _doc_
		typedef typename Value<TSeqSet>::Type TSeq;
		
		Directory		dir(dirname);
		TSeq			seq;
		MultiSeqFile	multiSeqFile;
		
		if (!atEnd(dir))
		{
			// dirname is path of a directory
			CharString fname = dirname;
#ifdef PLATFORM_WINDOWS
			appendValue(fname, '\\');
#else
			appendValue(fname, '/');
#endif
			size_t len = length(fname);

			for (; !atEnd(dir); goNext(dir))
			{
                replace(fname, length(fname) - len, length(fname), value(dir));
				if (guessFormatFromFilename(fname, format))
				{
					if (!open(multiSeqFile.concat, toCString(fname), OPEN_RDONLY)) continue;

					split(multiSeqFile, format);
					unsigned seqCount = length(multiSeqFile);
					
					reserve(seqSet, length(seqSet) + seqCount, Generous());					
					for(unsigned i = 0; i < seqCount; ++i)
					{
						assignSeq(seq, multiSeqFile[i], format);
						appendValue(seqSet, seq, Generous());
					}
					close(multiSeqFile.concat);
				}
			}
		} 
		else
		{
			// dirname is path of a file
			if (guessFormatFromFilename(dirname, format))
			{
				if (!open(multiSeqFile.concat, toCString(dirname), OPEN_RDONLY)) return;

				split(multiSeqFile, format);
				unsigned seqCount = length(multiSeqFile);
				
				reserve(seqSet, length(seqSet) + seqCount, Generous());					
				for(unsigned i = 0; i < seqCount; ++i)
				{
					assignSeq(seq, multiSeqFile[i], format);
					appendValue(seqSet, seq, Generous());
				}
				close(multiSeqFile.concat);
			}
		}
	}
	
}

#endif
