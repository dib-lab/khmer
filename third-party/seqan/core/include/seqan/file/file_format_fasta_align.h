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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Support for writing and reading FASTA alignment files.
// ==========================================================================

#ifndef SEQAN_FILE_FILE_FORMAT_FASTA_ALIGN_H_
#define SEQAN_FILE_FILE_FORMAT_FASTA_ALIGN_H_

/* IOREV
 * _tested_
 * _nodoc_
 *
 * tested in tests/file/test_file.h
 * tag mentionen in doc, but no further documentation, no link to spec
 *
 *
 * bacthreading and -writing
 *
 * does not treat IDs as Meta-Data, while format_fasta does
 * 
 */


namespace seqan {

// ===========================================================================
// Forward Declarations
// ===========================================================================

//forward declarations
template <typename T>
struct Row;
//IOREV

template <typename T>
struct Rows;
//IOREV

// ===========================================================================
// Tags, Enums, Classes, Specializations
// ===========================================================================

/**
.Tag.File Format.tag.Fasta alignment:
	FASTA alignment file format for sequences.
..include:seqan/file.h
*/
struct FastaAlign_;
//IOREV
typedef Tag<FastaAlign_> FastaAlign; //IOREV

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

template <typename TFile, typename TSize>
void _fastaAlignScanLine(TFile & file, TSize & count) {
//IOREV _notinlined_ no EOF check in beginning; non-default EOL-handling

	SEQAN_CHECKPOINT;
	SEQAN_ASSERT_NOT(_streamEOF(file));

	while (true) {
		typename Value<TFile>::Type c = _streamGet(file);

		if (_streamEOF(file)) return;
		if (c == '\n') return;

		if ((c != '\r') && (c!='-'))
			++count;
	}
}

//////////////////////////////////////////////////////////////////////////////
// read
//////////////////////////////////////////////////////////////////////////////
template <typename TFile, typename TSource, typename TSpec>
void read(TFile & file, Align<TSource, TSpec> & align, FastaAlign const &) {
//IOREV _notinlined_ _batchreading_ not sure whether begin_pos is absolute
    SEQAN_CHECKPOINT;

	SEQAN_ASSERT_NOT(_streamEOF(file));
	
	typedef typename Value<TSource>::Type TSourceValue;
	typedef typename Size<TSourceValue>::Type TSize;
	typedef typename Position<TFile>::Type TFilePos;
	typedef Triple<TFilePos, TFilePos, TSize> TTriple;
	TSize limit = maxValue<TSize>();

	//Determine begin position, end position and length of each sequence
	String<TTriple> beg_end_length;
	
	TFilePos begin_pos;
	TFilePos end_pos;
	typename Value<TFile>::Type c;
	TSize count;

	while (!_streamEOF(file)) {
		begin_pos = _streamTellG(file);
		count = 0;
		SEQAN_ASSERT_NOT(_streamEOF(file));

	
		c = _streamGet(file);
		
		// Skip id
		if (c == '>') {
			_fastaAlignScanLine(file, count);
			begin_pos = _streamTellG(file);
			count = 0;
		} else {  //If no id first letter belongs to sequence
			count = 1;
		}

		// Count letters
		while (true) {
			_fastaAlignScanLine(file, count);

			typename Value<TFile>::Type c = _streamGet(file);
			if (c == '>') {
				_streamSeek2G(file, -1);
				end_pos = _streamTellG(file);
				break;
			}
			if (_streamEOF(file)) {
				end_pos = _streamTellG(file);
				break;
			}
			if ((c != '\n') && (c != '\r') && (c!='-'))	{
				++count;
			}
		}
		if (count > limit) {
			count = limit;
		}

		appendValue(beg_end_length, TTriple(begin_pos, end_pos, count));
	}

	// Resize alignment data structure
	TSize numRows=length(beg_end_length);
	resize(rows(align), numRows);	//rows

	//typedef Align<TSource, TSpec> TAlign;

	for(TSize i=0;i<numRows;++i) {
		TSize begin = beg_end_length[i].i1;
//		TSize end = beg_end_length[i].i2;
		count = beg_end_length[i].i3;

		//Reserve space
        TSource buffer;
        resize(buffer, count);
        assignSource(row(align, i), buffer);
		// clear(row(align,i));
		// createSource(row(align,i));
		// resize(source(row(align,i)),count);
		// if (length(source(row(align,i))) < count) {
		//	count = length(source(row(align,i)));
		// }
		// setClippedEndPosition(row(align,i),count);

		//Read sequence
		_streamSeekG(file, begin);

		typename Position<TSource>::Type pos, viewPos;
		for (pos = 0, viewPos = 0; pos < count; ++viewPos) {
			c = _streamGet(file);
			if ((c != '\n') && (c != '\r') && (c != '-'))	{
				source(row(align,i))[pos] = c;
				++pos;
			}
			if (c=='-') {
				insertGap(row(align,i), viewPos);
			}
		}
	}

	_streamSeekG(file, 0);
}

//////////////////////////////////////////////////////////////////////////////
// readIDs
//////////////////////////////////////////////////////////////////////////////
 
template <typename TFile, typename TStringContainer>
void readIDs(TFile& file, TStringContainer& ids, FastaAlign) {
//IOREV _notinlined_ _batchreading_ _nodoc_
	
	SEQAN_CHECKPOINT;
	
    SEQAN_ASSERT_NOT(_streamEOF(file));

	typedef typename Value<TStringContainer>::Type TString;
	typename Position<TFile>::Type start_pos;
	typename Value<TFile>::Type c;


	TString seqId;
	while(true) {
		c = _streamGet(file);
		while ((!_streamEOF(file)) && (c != '>')) c = _streamGet(file);
		if (!_streamEOF(file)) {
			start_pos = _streamTellG(file);
			typename Size<TString>::Type count = 0;
			_fastaAlignScanLine(file, count);
			if (! count) clear(seqId);
			else {
				resize(seqId, count);
				if (length(seqId) < count)	{
					count = length(seqId);
				}
				_streamSeekG(file, start_pos);
				for (typename Position<TString>::Type pos = 0; pos<count; ++pos) {
					seqId[pos] = _streamGet(file);
				}
			}
			appendValue(ids, seqId);
		} else {
			break;
		}
	}
	_streamSeekG(file, 0);
}

//////////////////////////////////////////////////////////////////////////////
// readMeta
//////////////////////////////////////////////////////////////////////////////

//Fasta file records have no meta data

template <typename TFile, typename TMeta>
void readMeta(TFile & /*file*/, TMeta & meta, FastaAlign) {
//IOREV _notinlined_ _bug_ ids is meta-data, at least it is handled this way in format_fasta
	SEQAN_CHECKPOINT
	clear(meta);
}


//////////////////////////////////////////////////////////////////////////////
// goNext
//////////////////////////////////////////////////////////////////////////////
template <typename TFile>
void goNext(TFile & file, FastaAlign) {
//IOREV
	SEQAN_CHECKPOINT;
	(void) file; // When compiled without assertions.
	SEQAN_ASSERT_NOT(_streamEOF(file));
	
	return;
}


//////////////////////////////////////////////////////////////////////////////
// write
//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TStringContainer, typename TSource, typename TSpec>
void _writeImpl(TFile & file, Align<TSource, TSpec> const & align, TStringContainer const & ids, FastaAlign const &) {
//IOREV _notinlined_ _batchreading_
	SEQAN_CHECKPOINT

	typedef Align<TSource, TSpec> const TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Position<typename Rows<TAlign>::Type>::Type TRowsPosition;
	TRowsPosition row_count = length(rows(align));

	for(TRowsPosition i=0;i<row_count;++i) {
		TRow & row_ = row(align, i);

		typedef typename Iterator<typename Row<TAlign>::Type const, Standard>::Type TIter;
		TIter begin_ = iter(row_, beginPosition(cols(align)));
		TIter end_ = iter(row_, endPosition(cols(align)));

		_streamPut(file, '>');
		_streamWrite(file, getValue(ids,i));
		_streamPut(file, '\n');

		int chars=0;
		while(begin_ != end_) {
			if (chars == 60) {
				_streamPut(file, '\n');
				chars = 0;
			}
			if (isGap(begin_)) _streamPut(file, gapValue<char>());
			else _streamPut(file, getValue(begin_));
			chars++;
			++begin_;
		}
		_streamPut(file, '\n');
	}
}

//____________________________________________________________________________

template <typename TFile, typename TSource, typename TSpec>
void write(TFile & file, Align<TSource, TSpec> const & align, FastaAlign const & ) {
//IOREV _notinlined_ _batchreading_
	SEQAN_CHECKPOINT
	_writeImpl(file, align, String<String<char> >(), FastaAlign());
}

//____________________________________________________________________________

template <typename TFile, typename TStringContainer, typename TSource, typename TSpec>
void write(TFile & file, Align<TSource, TSpec> const & align, TStringContainer const & ids, FastaAlign const & ) {
//IOREV _notinlined_ _batchreading_
	SEQAN_CHECKPOINT
	_writeImpl(file, align, ids, FastaAlign());
}


//VisualC++ const array bug workaround
// TODO(holtgrew): Superflous?!
template <typename TFile, typename TStringContainer, typename TSource, typename TSpec>
void write(TFile & file, Align<TSource, TSpec> const * align, TStringContainer const & ids, FastaAlign const & ) {
//IOREV  _notinlined_ _batchreading_ _windows_
	SEQAN_CHECKPOINT
	_writeImpl(file, align, ids, FastaAlign());
}

//____________________________________________________________________________

template <typename TFile, typename TStringContainer, typename TSource, typename TSpec, typename TMeta>
void write(TFile & file, Align<TSource, TSpec> const & align, TStringContainer const & ids, TMeta &, FastaAlign const & ) {
//IOREV _notinlined_ _batchreading_
	SEQAN_CHECKPOINT;
	_writeImpl(file, align, ids, FastaAlign());
}

}  // namespace seqan

#endif   // #ifndef SEQAN_FILE_FILE_FORMAT_FASTA_ALIGN_H_
