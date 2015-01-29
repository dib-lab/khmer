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

#ifndef SEQAN_HEADER_FILE_CGVIZ_H
#define SEQAN_HEADER_FILE_CGVIZ_H

/* IOREV
 * _tested_
 * _nodoc_
 *
 * tested in tests/file/test_file.h
 * tag mentionen in doc, but no further documentation, no link to spec
 * 
 */


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// File Formats - CGViz
//////////////////////////////////////////////////////////////////////////////


/**
.Tag.File Format.tag.CGViz:
	CGViz file format for sequences. Only output.
..include:seqan/file.h
*/
struct TagCGViz_;
//IOREV
typedef Tag<TagCGViz_> const CGViz; //IOREV

/////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////
// goNext
//////////////////////////////////////////////////////////////////////////////

template <typename TFile>
void goNext(TFile & file, CGViz) {
//IOREV _notinlined_ purpose not clear to me
	SEQAN_CHECKPOINT;
    (void) file; // When compiled without assertions.
	SEQAN_ASSERT_NOT(streamEof(file));
	
	return;
}



//////////////////////////////////////////////////////////////////////////////
// write
//////////////////////////////////////////////////////////////////////////////

// Returns status code.

template <typename TFile, typename TStringContainer, typename TSource, typename TSpec>
int _writeImpl(TFile& target, Align<TSource, TSpec>& align, TStringContainer& ids, CGViz)
{
    typedef Align<TSource, TSpec> const TAlign;
    typedef typename Position<typename Rows<TAlign>::Type>::Type TRowsPosition;
    typedef typename Position<TAlign>::Type TPosition;
    TRowsPosition row_count = length(rows(align));
    if (row_count < 2)
        return 1;

    unsigned int pair=1;
    unsigned int count=0;
    for (TRowsPosition i=0;i<row_count-1;++i)
    {
        for (TRowsPosition j=i+1;j<row_count;++j)
        {
            // Print header
            streamPut(target, "{DATA dat");
            streamPut(target, pair);
            streamPut(target, '\n');
            streamPut(target, "[__GLOBAL__] dimension=2:\n");

            TPosition begin_ = beginPosition(cols(align));
            TPosition end_ = endPosition(cols(align));

            bool match = false;
            while(begin_ < end_)
            {
                if ((row(align, i)[begin_]==row(align, j)[begin_]) && (row(align, i)[begin_]!='-'))
                {
                    if (!match)
                    {
                        match = true;
                        streamPut(target, toSourcePosition(row(align,i),begin_+1));
                        streamPut(target, ' ');
                        streamPut(target, toSourcePosition(row(align,j),begin_+1));
                        streamPut(target, ' ');
                    }
                }
                if ((row(align, i)[begin_]!=row(align, j)[begin_]) || (row(align, i)[begin_]=='-') ||
                    (row(align, j)[begin_]=='-'))
                {
                    if (match)
                    {
                        streamPut(target, toSourcePosition(row(align,i),begin_));
                        streamPut(target, ' ');
                        streamPut(target, toSourcePosition(row(align,j),begin_));
                        streamPut(target, '\n');
                        match = false;
                    }
                }
                begin_++;
            }
            if (match)
            {
                streamPut(target, toSourcePosition(row(align,i),begin_));
                streamPut(target, ' ');
                streamPut(target, toSourcePosition(row(align,j),begin_));
                streamPut(target, '\n');
                match = false;
            }
            streamPut(target, '}');
            streamPut(target, '\n');

            // Write footer
            streamPut(target, "{GLYPH Glyph");
            streamPut(target, pair);
            streamPut(target, '\n');
            streamPut(target, "drawerName=Lines\n");
            streamPut(target, "lineWidth=3\n");
            streamPut(target, '}');
            streamPut(target, '\n');
            streamPut(target, "{PANE Pane");
            streamPut(target, pair);
            streamPut(target, '\n');
            streamPut(target, "uLabel=");
            streamPut(target, getValue(ids,i));
            streamPut(target, '\n');
            streamPut(target, "uStop=");
            streamPut(target, length(source(row(align,i))));
            streamPut(target, '\n');
            streamPut(target, "vLabel=");
            streamPut(target, getValue(ids,j));
            streamPut(target, '\n');
            streamPut(target, "vStop=");
            streamPut(target, length(source(row(align,j))));
            streamPut(target, '\n');
            streamPut(target, '}');
            streamPut(target, '\n');
            streamPut(target, "{WINDOW Window");
            streamPut(target, pair);
            streamPut(target, '\n');
            streamPut(target, '}');
            streamPut(target, '\n');
            streamPut(target, "{FEEDER Feeder<");
            streamPut(target, pair);
            streamPut(target, '>');
            streamPut(target, ' ');
            streamPut(target, count);
            streamPut(target, ' ');
            streamPut(target, count+1);
            streamPut(target, '\n');
            streamPut(target, '}');
            streamPut(target, '\n');
            ++count;
            streamPut(target, "{THREADER Threader<");
            streamPut(target, pair);
            streamPut(target, '>');
            streamPut(target, ' ');
            streamPut(target, count);
            streamPut(target, ' ');
            streamPut(target, count+1);
            streamPut(target, '\n');
            streamPut(target, '}');
            streamPut(target, '\n');
            ++count;
            streamPut(target, "{ANCHOR Anchor<");
            streamPut(target, pair);
            streamPut(target, '>');
            streamPut(target, ' ');
            streamPut(target, count);
            streamPut(target, ' ');
            streamPut(target, count+1);
            streamPut(target, '\n');
            streamPut(target, '}');
            streamPut(target, '\n');
            count+=2;
            ++pair;
        }
    }

    return streamError(target);
}


//____________________________________________________________________________

template <typename TFile, typename TSource, typename TSpec>
int write(TFile & file, Align<TSource, TSpec>& align, CGViz)
{
	return _writeImpl(file, align, String<String<char> >(), CGViz());
}

//____________________________________________________________________________

template <typename TFile, typename TStringContainer, typename TSource, typename TSpec>
int write(TFile & file, Align<TSource, TSpec> & align, TStringContainer& ids, CGViz)
{
    return _writeImpl(file, align, ids, CGViz());
}


// TODO(holtgrew): Is this still necessary? I believe this is a VS2005 issue.
//VisualC++ const array bug workaround
template <typename TFile, typename TStringContainer, typename TSource, typename TSpec>
int write(TFile & file, Align<TSource, TSpec>* align, TStringContainer & ids, CGViz)
{
	return _writeImpl(file, align, ids, CGViz());
}

//____________________________________________________________________________

template <typename TFile, typename TStringContainer, typename TSource, typename TSpec, typename TMeta>
int write(TFile & file, Align<TSource, TSpec> & align, TStringContainer& ids, TMeta &, CGViz)
{
	return _writeImpl(file, align, ids, CGViz());
}

}  // namespace seqan

#endif //#ifndef SEQAN_HEADER_...
