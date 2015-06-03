// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: Tobias Rausch <rausch@embl.de>
// ==========================================================================
// Input and Output of an alignment graph, Tree Reading im Newick Format
// ==========================================================================

#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_IO_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_IO_H

#include <seqan/stream.h>

namespace SEQAN_NAMESPACE_MAIN
{

// ============================================================================
// Forwards
// ============================================================================

struct FastaAlign_;
typedef Tag<FastaAlign_> FastaAlign;

// ============================================================================
// Tags
// ============================================================================

//////////////////////////////////////////////////////////////////////////////

/*!
 * @tag AlignmentGraphFormatTags#TCoffeLib
 * @brief T-Coffee library format to read and write an @link AlignmentGraph @endlink.
 *
 * @signature typedef Tag<TCoffeeLib_> const TCoffeeLib;
 *
 * @tag AlignmentGraphFormatTags#BlastLib
 * @brief A BLAST library format to read and write an @link AlignmentGraph @endlink.
 *
 * @signature typedef Tag<BlastLib_> const BlastLib;
 *
 * @tag AlignmentGraphFormatTags#MummerLib
 * @brief A MUMMER library format to read and write an @link AlignmentGraph @endlink.
 *
 * @signature typedef Tag<MummerLib_> const MummerLib;
 *
 * @tag AlignmentGraphFormatTags#NewickFormat
 * @brief NewickFormat to write a guide tree.
 *
 * @signature typedef Tag<NewickFormat_> const NewickFormat;
 */

struct TCoffeeLib_;
typedef Tag<TCoffeeLib_> const TCoffeeLib;

//////////////////////////////////////////////////////////////////////////////

struct BlastLib_;
typedef Tag<BlastLib_> const BlastLib;

//////////////////////////////////////////////////////////////////////////////

struct MummerLib_;
typedef Tag<MummerLib_> const MummerLib;


//////////////////////////////////////////////////////////////////////////////

struct NewickFormat_;
typedef Tag<NewickFormat_> const NewickFormat;

typedef OrFunctor<EqualsChar<';'>,
                  OrFunctor<EqualsChar<')'>,
                  OrFunctor<EqualsChar<','>, EqualsChar<':'> > > > IsNewickNodeStop_;

typedef OrFunctor<EqualsChar<'['>,
                  OrFunctor<EqualsChar<']'>,
                  OrFunctor<EqualsChar<'('>, IsWhitespace> > > IsForbiddenNewickIdentifierChar_;

/////////////////////////////////////////////////////////////////////////////
// T-Coffee Library Reading / Writing
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////



template<typename TFile, typename TFragment, typename TSpec, typename TScoreValue, typename TSpec2, typename TNames>
void
read(TFile & file,
     String<TFragment, TSpec>& matches,
     String<TScoreValue, TSpec2>& scores,
     TNames const& names,
     TCoffeeLib)
{
//IOREV _nodoc_ _notinlined_ specialization not documented
    typedef typename Size<TFragment>::Type TSize;
    //typedef typename Value<TFile>::Type TValue;
    typedef typename Value<TNames>::Type TName;

    typename DirectionIterator<TFile, Input>::Type reader = directionIterator(file, Input());
    CharString buffer;

    // Map the names to slots
    typedef std::map<TName, TSize> TNameToPosition;
    TNameToPosition namePosMap;
    for(TSize i = 0;i<length(names);++i) namePosMap.insert(std::make_pair(names[i], i));

    // Remember the correct spots
    if (atEnd(reader))
        return;
    skipLine(reader);
    clear(buffer);
    readUntil(buffer, reader, NotFunctor<IsDigit>());
    TSize nseq = lexicalCast<TSize>(buffer);
    skipLine(reader);

    // Read sequences
    typedef String<TSize> TMapping;
    TName myName;
    TMapping posMap;
    resize(posMap, nseq);
    for (TSize i=0; i<nseq; ++i)
    {
        clear(myName);
        readUntil(myName, reader, IsWhitespace());
        value(posMap, i) = namePosMap.find(myName)->second;
        skipLine(reader);
    }

    bool seq1ToN = false;
    if (atEnd(reader))
        return;

    typedef std::pair<std::pair<TSize, TSize>, TScoreValue> TResiduePair;
    typedef std::set<TResiduePair> TResiduePairSet;
    String<TResiduePairSet> resPair;
    resize(resPair, nseq * nseq);
    TSize seq1 = 0;
    TSize seq2 = 0;
    bool firstPass = true;
    while (!atEnd(reader))
    {
        skipUntil(reader, NotFunctor<IsWhitespace>());
        if (atEnd(reader))
            break;
        if (value(reader) == '#') {
            goNext(reader);
            skipUntil(reader, NotFunctor<IsWhitespace>());
            clear(buffer);
            readUntil(buffer, reader, IsWhitespace());
            seq1 = lexicalCast<TSize>(buffer);
            skipUntil(reader, NotFunctor<IsWhitespace>());
            clear(buffer);
            readUntil(buffer, reader, IsWhitespace());
            seq2 = lexicalCast<TSize>(buffer);
            if (firstPass) {
                firstPass = false;
                if ((seq1 != 0) && (seq2 != 0)) seq1ToN = true;
            }
            if (seq1ToN) {
                --seq1;
                --seq2;
            }
            seq1 = value(posMap, seq1);
            seq2 = value(posMap, seq2);
        } else if (value(reader) == '!') {
            skipUntil(reader, NotFunctor<IsWhitespace>());
        } else {
            clear(buffer);
            readUntil(buffer, reader, IsWhitespace());
            TSize res1 = lexicalCast<TSize>(buffer);
            skipUntil(reader, NotFunctor<IsWhitespace>());
            clear(buffer);
            readUntil(buffer, reader, IsWhitespace());
            TSize res2 = lexicalCast<TSize>(buffer);
            skipUntil(reader, NotFunctor<IsWhitespace>());
            clear(buffer);
            readUntil(buffer, reader, IsWhitespace());
            TScoreValue weight = lexicalCast<TScoreValue>(buffer);
            skipLine(reader);

            if (seq1 < seq2) {
                TSize index = seq1 * nseq + seq2;
                resPair[index].insert(std::make_pair(std::make_pair(--res1,--res2), weight));
            } else {
                TSize index = seq2 * nseq + seq1;
                resPair[index].insert(std::make_pair(std::make_pair(--res2,--res1), weight));
            }
        }
    }
    for(TSize i = 0; i<length(resPair); ++i) {
        if (resPair[i].empty()) continue;
        TSize seq1 = i / nseq;
        TSize seq2 = i % nseq;
        //std::cout << "#" << seq1 << ',' << seq2 << std::endl;
        typename TResiduePairSet::const_iterator pos = resPair[i].begin();
        typename TResiduePairSet::const_iterator posEnd = resPair[i].end();
        TSize startMatch1 = pos->first.first;
        TSize startMatch2 = pos->first.second;
        TScoreValue carg = pos->second;
        TSize len = 1;
        //std::cout << pos->first.first << ',' << pos->first.second << ',' << pos->second << std::endl;
        ++pos;
        while(pos != posEnd) {
            if ((startMatch1 + len == pos->first.first) &&
                (startMatch2 + len == pos->first.second) &&
                (carg / (TScoreValue) len == pos->second)) {
                    carg += pos->second;
                    ++len;
            } else {
                appendValue(matches, TFragment(seq1, startMatch1, seq2, startMatch2, len));
                appendValue(scores, carg);
                startMatch1 = pos->first.first;
                startMatch2 = pos->first.second;
                carg = pos->second;
                len = 1;
            }
            //std::cout << pos->first.first << ',' << pos->first.second << ',' << pos->second << std::endl;
            ++pos;
        }
        appendValue(matches, TFragment(seq1, startMatch1, seq2, startMatch2, len));
        appendValue(scores, carg);
    }
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TString, typename TSpec, typename TNames>
void
read(TFile & file,
     StringSet<TString, TSpec>& oriStr,
     TNames& names,
     TCoffeeLib)
{
//IOREV _nodoc_ _notinlined_ specialization not documented
    typedef typename Size<TNames>::Type TSize;
    //typedef typename Position<TFile>::Type TPosition;
    //typedef typename Value<TFile>::Type TValue;

    CharString buffer;
    typename DirectionIterator<TFile, Input>::Type reader = directionIterator(file, Input());

    if (atEnd(reader))
        return;

    // Ignore first line
    skipLine(reader);

    // Read number of sequences
    clear(buffer);
    readUntil(buffer, reader, IsWhitespace());
    TSize nSeq = lexicalCast<TSize>(buffer);
    resize(oriStr, nSeq);
    skipLine(reader);

    // Read sequences
    for(TSize i = 0; i < nSeq; ++i)
    {
        clear(buffer);
        readUntil(buffer, reader, IsWhitespace());
        appendValue(names, buffer);
        skipUntil(reader, NotFunctor<IsWhitespace>());
        skipUntil(reader, IsWhitespace());
        // Note, not used.
        skipUntil(reader, NotFunctor<IsWhitespace>());
        readUntil(oriStr[i], reader, IsWhitespace());
        skipLine(reader);
    }
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TStringSet, typename TCargo, typename TSpec, typename TNames>
void write(TFile & file,
           Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
           TNames& names,
           TCoffeeLib)
{
//IOREV _nodoc_ _notinlined_ specialization not documented
    SEQAN_CHECKPOINT
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    //typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef typename Value<TStringSet>::Type TString;
    typedef typename Size<TStringSet>::Type TSize;

    typedef std::pair<std::pair<TSize, TSize>, TCargo> TResiduePair;
    typedef std::set<TResiduePair> TResiduePairSet;
    TSize nseq = length(stringSet(g));
    String<TResiduePairSet> resPair;
    resize(resPair, nseq * nseq);

    typedef typename Iterator<TGraph, EdgeIterator>::Type TIter;
    TIter it(g);
    for(;!atEnd(it);++it) {
        TVertexDescriptor sV = sourceVertex(it);
        TVertexDescriptor tV = targetVertex(it);
        TSize fragLen = fragmentLength(g,sV);
        TSize fragPos1 = fragmentBegin(g,sV);
        TSize fragPos2 = fragmentBegin(g,tV);
        TSize seq1 = sequenceId(g,sV);
        TSize seq2 = sequenceId(g,tV);
        TCargo my_carg =  getCargo(*it);
        if (my_carg <= 0) my_carg = 1;
        else my_carg = (TCargo) ((double) my_carg / (double) fragLen);
        for(TSize i = 0; i<fragLen; ++i) {
            resPair[seq1 * nseq + seq2].insert(std::make_pair(std::make_pair(fragPos1 + i, fragPos2 + i), my_carg));
        }
    }

    write(file, "! TC_LIB_FORMAT_01\n");
    TSize len = length(getStringSet(g));
    appendNumber(file, (int)len);
    writeValue(file, '\n');
    for(TSize i=0;i<len;++i) {
        write(file, names[i]);
        writeValue(file, ' ');
        TString str = value(getStringSet(g), i);
        appendNumber(file, length(str));
        writeValue(file, ' ');
        write(file, str);
        writeValue(file, '\n');
    }

    for(TSize i = 0; i<length(resPair); ++i) {
        if (resPair[i].empty()) continue;
        TSize seq1 = i / nseq;
        TSize seq2 = i % nseq;
        writeValue(file, '#');
        appendNumber(file, seq1 + 1);
        writeValue(file, ' ');
        appendNumber(file, seq2 + 1);
        writeValue(file, '\n');
        typename TResiduePairSet::const_iterator pos = resPair[i].begin();
        typename TResiduePairSet::const_iterator posEnd = resPair[i].end();
        while(pos != posEnd) {
            appendNumber(file, pos->first.first + 1);
            writeValue(file, ' ');
            appendNumber(file, pos->first.second + 1);
            writeValue(file, ' ');
            appendNumber(file, pos->second);
            writeValue(file, '\n');
            ++pos;
        }
    }
    write(file, "! SEQ_1_TO_N");
    writeValue(file, '\n');
}


/////////////////////////////////////////////////////////////////////////////
// FastaAlign Reading
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TString, typename TSpec, typename TNames>
void
read(TFile & file,
     StringSet<TString, TSpec>& oriStr,
     TNames& names,
     FastaAlign)
{
//IOREV _nodoc_ _notinlined_ specialization not documented
    SEQAN_CHECKPOINT
    //typedef typename Size<TNames>::Type TSize;
    //typedef TSize TWord;
    //typedef typename Value<TFile>::Type TValue;

    CharString buffer;
    typename DirectionIterator<TFile, Input>::Type reader = directionIterator(file, Input());
    if (atEnd(reader))
        return;

    // Read sequences
    TString seq;
    while(!atEnd(reader))
    {
        skipUntil(reader, NotFunctor<IsWhitespace>());
        if (atEnd(reader))
            break;
        if (value(reader) == '>')
        {
            if (!empty(seq))
            {
                appendValue(oriStr, seq);
                clear(seq);
            }
            goNext(reader);
            clear(buffer);
            readUntil(buffer, reader, IsWhitespace());
            appendValue(names, buffer);
            skipLine(reader);
        }
        else if ((value(reader) == '-') || (value(reader) == '.') || (value(reader) == '\n') ||
                 (value(reader) == '\r'))
        {
            goNext(reader);
        }
        else
        {
            appendValue(seq, value(reader));
            goNext(reader);
        }
    }
    if (!empty(seq))
        appendValue(oriStr, seq);
}


template<typename TValue, typename TSpec2, typename TFragment, typename TSpec, typename TScores, typename TSize>
void
_collectSegmentMatches(String<TValue, TSpec2> const& mat,
                       String<TFragment, TSpec>& matches,
                       TScores& scores,
                       TSize nseq)
{
    SEQAN_CHECKPOINT
    TSize len = length(mat) / nseq;
    TValue gapChar = gapValue<TValue>();

    // Create the anchor graph
    //typedef String<TFragment, TSpec> TFragmentString;
    //typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
    typedef std::pair<TSize, TSize> TResiduePair;
    typedef std::set<TResiduePair> TResiduePairSet;
    String<TResiduePairSet> resPair;
    resize(resPair, nseq * nseq);
    for(TSize seq1 = 0; seq1 < nseq - 1; ++seq1) {
        for(TSize seq2 = seq1 + 1; seq2 < nseq; ++seq2) {
            TSize index = seq1 * nseq + seq2;
            TSize offset1 = 0;
            TSize offset2 = 0;
            for(TSize col = 0; col<len; ++col) {
                if (value(mat, seq1 * len + col) != gapChar) {
                    if (value(mat, seq2 * len + col) != gapChar) {
                        resPair[index].insert(std::make_pair(offset1, offset2));
                        ++offset1;
                        ++offset2;
                    } else ++offset1;
                } else if (value(mat, seq2 * len + col) != gapChar) ++offset2;
            }
        }
    }
    for(TSize i = 0; i<length(resPair); ++i) {
        if (resPair[i].empty()) continue;
        TSize seq1 = i / nseq;
        TSize seq2 = i % nseq;
        typename TResiduePairSet::const_iterator pos = resPair[i].begin();
        typename TResiduePairSet::const_iterator posEnd = resPair[i].end();
        TSize startMatch1 = pos->first;
        TSize startMatch2 = pos->second;
        TSize len = 1;
        ++pos;
        while(pos != posEnd) {
            if ((startMatch1 + len == pos->first) && (startMatch2 + len == pos->second)) ++len;
            else {
                appendValue(matches, TFragment(seq1, startMatch1, seq2, startMatch2, len));
                appendValue(scores, len);
                startMatch1 = pos->first;
                startMatch2 = pos->second;
                len = 1;
            }
            ++pos;
        }
        appendValue(matches, TFragment(seq1, startMatch1, seq2, startMatch2, len));
        appendValue(scores, len);
    }
}


//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TFragment, typename TSpec, typename TScoreValue, typename TSpec2, typename TNames>
void
read(TFile & file,
     String<TFragment, TSpec>& matches,
     String<TScoreValue, TSpec2>& scores,
     TNames const& origNames,
     FastaAlign)
{
//IOREV _nodoc_ _notinlined_ specialization not documented
    SEQAN_CHECKPOINT
    typedef typename Size<TNames>::Type TSize;
    typedef typename Value<TFile>::Type TValue;
    typedef typename Value<TNames>::Type TName;

    CharString buffer;
    typename DirectionIterator<TFile, Input>::Type reader = directionIterator(file, Input());
    if (atEnd(reader))
        return;

    // Read sequences
    String<TValue> mat;
    TNames names;
    TName nextSeq;
    while (!atEnd(reader))
    {
        skipUntil(reader, NotFunctor<IsWhitespace>());
        if (atEnd(reader))
            break;
        if (value(reader) == '>')
        {
            goNext(reader);
            clear(nextSeq);
            readUntil(nextSeq, reader, IsWhitespace());
            appendValue(names, nextSeq);
            skipLine(reader);
        }
        else if ((value(reader) == '\n') || (value(reader) == '\r'))
        {
            goNext(reader);
        }
        else
        {
            appendValue(mat, value(reader));
            goNext(reader);
        }
    }
    // Reorder rows according to names order
    String<TValue> finalMat = mat;
    TSize nseq = length(names);
    TSize len = length(mat) / nseq;
    for(TSize i = 0; i<nseq; ++i) {
        if (value(names, i) == value(origNames, i)) continue;
        else {
            for(TSize j = 0; j<length(origNames); ++j) {
                if (value(names, i) != value(origNames, j)) continue;
                // Copy the whole row
                replace(finalMat, j * len, j * len + len, infix(mat, i * len, i*len + len));
                break;
            }
        }
    }
    clear(mat);

    // Collect the segment matches
    _collectSegmentMatches(finalMat, matches, scores, nseq);
}

/////////////////////////////////////////////////////////////////////////////
// BLAST Library Reading / Writing
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TSizeSpec, typename TSpec1, typename TSpec2, typename TSize>
inline void
_appendFragment(String<Fragment<TSizeSpec, ExactReversableFragment<TSpec1> >, TSpec2>& matches,
                  TSize seq1Id,
                  TSize beg1,
                  TSize seq2Id,
                  TSize beg2,
                  TSize len,
                  bool reversed)
{
    SEQAN_CHECKPOINT
    typedef Fragment<TSizeSpec, ExactReversableFragment<TSpec1> > TFragment;
    appendValue(matches, TFragment(seq1Id, beg1, seq2Id, beg2, len, reversed));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSizeSpec, typename TSpec1, typename TSpec2, typename TSize>
inline void
_appendFragment(String<Fragment<TSizeSpec, ExactFragment<TSpec1> >, TSpec2>& matches,
                  TSize seq1Id,
                  TSize beg1,
                  TSize seq2Id,
                  TSize beg2,
                  TSize len,
                  bool reversed)
{
    SEQAN_CHECKPOINT
    typedef Fragment<TSizeSpec, ExactFragment<TSpec1> > TFragment;
    if (!reversed) appendValue(matches, TFragment(seq1Id, beg1, seq2Id, beg2, len));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TFragment, typename TSpec1, typename TScoreValue, typename TSpec2, typename TNames>
inline void
read(TFile & file,
     String<TFragment, TSpec1>& matches,
     String<TScoreValue, TSpec2>& scores,
     TNames& names,
     BlastLib)
{
//IOREV _nodoc_ specialization not documented
    SEQAN_CHECKPOINT
    typedef typename Size<TNames>::Type TSize;
    //typedef typename Value<TFile>::Type TValue;
    typedef typename Value<TNames>::Type TName;

    // Map the names to slots
    typedef std::map<TName, TSize> TNameToPosition;
    TNameToPosition namePosMap;
    for(TSize i = 0;i<length(names);++i) namePosMap.insert(std::make_pair(names[i], i));

    CharString buffer;
    typename DirectionIterator<TFile, Input>::Type reader = directionIterator(file, Input());

    TName seq1;
    TName seq2;
    while (!atEnd(reader))
    {
        clear(seq1);
        clear(seq2);
        skipUntil(reader, NotFunctor<IsWhitespace>());
        readUntil(seq1, reader, IsWhitespace());

        skipUntil(reader, NotFunctor<IsWhitespace>());
        readUntil(seq2, reader, IsWhitespace());

        skipUntil(reader, NotFunctor<IsWhitespace>());
        if (seq1 == seq2)
        {
            skipLine(reader);
            continue;
        }
        TSize seq1Id = namePosMap[seq1];
        TSize seq2Id = namePosMap[seq2];
        skipUntil(reader, NotFunctor<IsWhitespace>());

        // We skip the data from the file, buffer is only rarely used.
        for (unsigned i = 0; i < 4; ++i)
        {
            skipUntil(reader, IsWhitespace());
            skipUntil(reader, NotFunctor<IsWhitespace>());
        }

        clear(buffer);
        readUntil(buffer, reader, IsWhitespace());
        TSize beg1 = lexicalCast<TSize>(buffer);
        skipUntil(reader, NotFunctor<IsWhitespace>());

        clear(buffer);
        readUntil(buffer, reader, IsWhitespace());
        TSize end1 = lexicalCast<TSize>(buffer);
        skipUntil(reader, NotFunctor<IsWhitespace>());

        TSize len = end1 - beg1 + 1;
        clear(buffer);
        readUntil(buffer, reader, IsWhitespace());
        TSize beg2 = lexicalCast<TSize>(buffer);

        skipUntil(reader, NotFunctor<IsWhitespace>());
        clear(buffer);
        readUntil(buffer, reader, IsWhitespace());
        TSize end2 = lexicalCast<TSize>(buffer);

        skipUntil(reader, IsWhitespace());
        skipUntil(reader, NotFunctor<IsWhitespace>());

        clear(buffer);
        readUntil(buffer, reader, IsWhitespace());
        TScoreValue rawScore = static_cast<TScoreValue>(lexicalCast<double>(buffer));

        bool reversed = false;
        if (beg1 > end1) { TSize tmp = beg1; beg1 = end1; end1 = tmp; reversed = !reversed; }
        if (beg2 > end2) { TSize tmp = beg2; beg2 = end2; end2 = tmp; reversed = !reversed; }
        //// Debug code
        //std::cout << seq1Id << ',' << beg1 << ',' << seq2Id << ',' << beg2 << ',' << len << std::endl;
        //std::cout << infix(strSet[seq1Id], beg1, beg1+len) << std::endl;
        //std::cout << infix(strSet[seq2Id], beg2, beg2+len) << std::endl;

        _appendFragment(matches, seq1Id, --beg1, seq2Id, --beg2, len, reversed);
        appendValue(scores, rawScore);
        skipLine(reader);
    }
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TStringSet, typename TCargo, typename TSpec, typename TNames, typename TEdgeMap>
void write(TFile & file,
           Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
           TNames& names,
           TEdgeMap& edgeMap,
           BlastLib)
{
//IOREV _nodoc_ _notinlined_ specialization not documented
    SEQAN_CHECKPOINT
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    //typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    //typedef typename Value<TStringSet>::Type TString;
    typedef typename Size<TStringSet>::Type TSize;

    TStringSet& str = stringSet(g);

    typedef typename Iterator<TGraph, EdgeIterator>::Type TIter;
    TIter it(g);
    for(;!atEnd(it);++it) {
        TVertexDescriptor sV = sourceVertex(it);
        TVertexDescriptor tV = targetVertex(it);
        if (sequenceId(g,sV) < sequenceId(g,tV)) {
            TVertexDescriptor tmp = sV;
            sV = tV;
            tV = tmp;
        }
        TSize fragLen = fragmentLength(g,sV);
        TSize fragPos1 = fragmentBegin(g,sV);
        TSize fragPos2 = fragmentBegin(g,tV);
        TSize seq1 = idToPosition(str, sequenceId(g,sV));
        TSize seq2 = idToPosition(str, sequenceId(g,tV));
        TCargo my_carg =  getCargo(*it);
        write(file, names[seq1]);
        writeValue(file, '\t');
        write(file, names[seq2]);
        write(file, "\t0\t");
        appendNumber(file, fragLen);
        write(file, "\t0\t0\t");
        appendNumber(file, fragPos1+1);
        writeValue(file, '\t');
        appendNumber(file, fragPos1 + fragLen);
        writeValue(file, '\t');
        if (!property(edgeMap, *it)) {
            appendNumber(file, fragPos2+1);
            writeValue(file, '\t');
            appendNumber(file, fragPos2 + fragLen);
            writeValue(file, '\t');
        } else {
            appendNumber(file, fragPos2 + fragLen);
            writeValue(file, '\t');
            appendNumber(file, fragPos2+1);
            writeValue(file, '\t');
        }
        write(file, "0\t");
        write(file, my_carg);
        writeValue(file, '\n');
    }
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TStringSet, typename TCargo, typename TSpec, typename TNames>
void write(TFile & file,
           Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
           TNames& names,
           BlastLib)
{
//IOREV _nodoc_ _notinlined_ specialization not documented
    SEQAN_CHECKPOINT
    String<bool> edgeMap;
    resize(edgeMap, getIdUpperBound(_getEdgeIdManager(g)), false);
    write(file, g, names, edgeMap, BlastLib());
}



/////////////////////////////////////////////////////////////////////////////
// MUMMER Format Reading
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TPos, typename TSpec2, typename TSpec1, typename TScores, typename TId, typename TSize>
inline void
_appendNewMatch(String<Fragment<TPos, ExactReversableFragment<TSpec2> >, TSpec1>& matches,
                 TScores& scores,
                 TId seq1Id,
                 TId seq2Id,
                 TSize beg1,
                 TSize beg2,
                 TSize len,
                 bool)
{
    typedef Fragment<TPos, ExactReversableFragment<TSpec2> > TFragment;
    appendValue(matches, TFragment(seq1Id, beg1, seq2Id, beg2, len));
    appendValue(scores, len * len);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TPos, typename TSpec2, typename TSpec1, typename TScores, typename TId, typename TSize>
inline void
_appendNewMatch(String<Fragment<TPos, ExactFragment<TSpec2> >, TSpec1>& matches,
                 TScores& scores,
                 TId seq1Id,
                 TId seq2Id,
                 TSize beg1,
                 TSize beg2,
                 TSize len,
                 bool reversed)
{
    typedef Fragment<TPos, ExactFragment<TSpec2> > TFragment;
    if (!reversed) {
        appendValue(matches, TFragment(seq1Id, beg1, seq2Id, beg2, len));
        appendValue(scores, len * len);
    }
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TFragment, typename TSpec1, typename TScoreValue, typename TSpec2, typename TStringSet, typename TNames>
inline void
read(TFile & file,
     String<TFragment, TSpec1>& matches,
     String<TScoreValue, TSpec2>& scores,
     TStringSet const& strSet,
     TNames const& names,
     MummerLib)
{
//IOREV _nodoc_ specialization not documented
    SEQAN_CHECKPOINT
    typedef typename Size<TNames>::Type TSize;
    //typedef typename Value<TFile>::Type TValue;
    typedef typename Value<TNames>::Type TName;

    CharString buffer;
    typename DirectionIterator<TFile, Input>::Type reader = directionIterator(file, Input());

    // Map the names to slots
    typedef std::map<TName, TSize> TNameToPosition;
    TNameToPosition namePosMap;
    for(TSize i = 0;i<length(names);++i) namePosMap.insert(std::make_pair(value(names, i), i));

    // Read the Mummer file
    if (atEnd(reader))
        return ;
    TName seq1;
    TName seq2;
    TSize seq1Id = 0;
    TSize seq2Id = 0;
    bool reversed = false;
    while (!atEnd(reader))
    {
        if (value(reader) == '>')
        {
            goNext(reader);
            skipUntil(reader, NotFunctor<IsWhitespace>());
            clear(seq1);
            readUntil(seq1, reader, IsWhitespace());
            seq1Id = namePosMap[seq1];
            skipUntil(reader, NotFunctor<IsWhitespace>());
            if (value(reader) == 'R')
            {
                 reversed = true;
                 skipLine(reader);
            }
            else
            {
                reversed = false;
            }
        }
        else
        {
            skipUntil(reader, NotFunctor<IsWhitespace>());
            if (atEnd(reader))
                break;
            clear(seq2);
            readUntil(seq2, reader, IsWhitespace());
            seq2Id = namePosMap[seq2];
            skipUntil(reader, NotFunctor<IsWhitespace>());
            clear(buffer);
            readUntil(buffer, reader, IsWhitespace());
            TSize beg2 = lexicalCast<TSize>(buffer);
            skipUntil(reader, NotFunctor<IsWhitespace>());
            clear(buffer);
            readUntil(buffer, reader, IsWhitespace());
            TSize beg1 = lexicalCast<TSize>(buffer);
            skipUntil(reader, NotFunctor<IsWhitespace>());
            clear(buffer);
            readUntil(buffer, reader, IsWhitespace());
            TSize len = lexicalCast<TSize>(buffer);
            skipLine(reader);
            if (seq1Id == seq2Id)
                continue;

            if (!reversed)
                _appendNewMatch(matches, scores, seq1Id, seq2Id, --beg1, --beg2, len, reversed);
            else
                _appendNewMatch(matches, scores, seq1Id, seq2Id, (length(value(strSet, seq1Id)) - (--beg1 + len)), --beg2, len, reversed);
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
// Newick Format
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TCargo, typename TSpec, typename TNames>
void
read(TFile & file,
     Graph<Tree<TCargo, TSpec> >& guideTree,
     TNames& names,
     NewickFormat)
{
//IOREV _nodoc_ _notinlined_ specialization not documented
    typedef Graph<Tree<TCargo, TSpec> > TGuideTree;
    typedef typename VertexDescriptor<TGuideTree>::Type TVertexDescriptor;
    //typedef typename EdgeDescriptor<TGuideTree>::Type TEdgeDescriptor;
    //typedef typename Size<TGuideTree>::Type TSize;
    typedef typename Id<TGuideTree>::Type TId;
    //typedef typename Position<TFile>::Type TPosition;
    //typedef typename Value<TFile>::Type TValue;
    typedef typename Value<TNames>::Type TName;
    TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();


    if (length(names) < 3)
    {
        TVertexDescriptor v1 = addVertex(guideTree);
        TVertexDescriptor v2 = addVertex(guideTree);
        TVertexDescriptor internalVertex = addVertex(guideTree);
        addEdge(guideTree, internalVertex, v1, (TCargo) 1 * SEQAN_DISTANCE_UNITY);
        addEdge(guideTree, internalVertex, v2, (TCargo) 1 * SEQAN_DISTANCE_UNITY);
        assignRoot(guideTree, internalVertex);
        return;
    }

    typedef std::map<TName, TId> TNameToId;
    TNameToId nameToId;
    for(TId i=0; i<length(names);++i) {
        addVertex(guideTree);    // Create the sequence vertices
        nameToId.insert(std::make_pair(names[i], i));
    }

    CharString buffer;
    typename DirectionIterator<TFile, Input>::Type reader = directionIterator(file, Input());
    if (atEnd(reader))
        return;

    TVertexDescriptor lastVertex = nilVertex;
    TVertexDescriptor lastChild = nilVertex;
    while (!atEnd(reader))
    {
        if (value(reader) == '(')
        {
            if (lastVertex == nilVertex) {
                lastVertex = addVertex(guideTree);
                assignRoot(guideTree, lastVertex);
            } else {
                TVertexDescriptor ch = addChild(guideTree, lastVertex);
                lastVertex = ch;
            }
            goNext(reader);
            skipUntil(reader, NotFunctor<IsWhitespace>());
        }
        else if (value(reader) == ')')
        {
            if (!isRoot(guideTree, lastVertex)) {
                lastChild = lastVertex;
                lastVertex = parentVertex(guideTree, lastVertex);
            } else {
                lastChild = lastVertex;
                lastVertex = nilVertex;
            }
            goNext(reader);
            skipUntil(reader, NotFunctor<IsWhitespace>());
        } else if (value(reader) == ',') {
            goNext(reader);
            skipUntil(reader, NotFunctor<IsWhitespace>());
        } else if (value(reader) == ':') {
            goNext(reader);
            clear(buffer);
            readUntil(buffer, reader, IsNewickNodeStop_());
            double tmpDouble = lexicalCast<double>(buffer);
            cargo(findEdge(guideTree, lastVertex, lastChild)) = (TCargo) (tmpDouble * SEQAN_DISTANCE_UNITY);
        } else if (value(reader) == ';') {
            goNext(reader);
            skipUntil(reader, NotFunctor<IsWhitespace>());
        } else {
            TName tmp;
            readUntil(tmp, reader, IsNewickNodeStop_(),
                                   AssertFunctor<NotFunctor<IsForbiddenNewickIdentifierChar_>, ParseError>());
            //std::cout << tmp << std::endl;
            if (lastVertex == nilVertex) {
                // Tree is rooted at a leaf
                // Create artificial root node
                lastVertex = length(names);
                assignRoot(guideTree, addVertex(guideTree));
                addEdge(guideTree, getRoot(guideTree), lastVertex);
                addEdge(guideTree, getRoot(guideTree), nameToId[tmp]);
            } else {
                addEdge(guideTree, lastVertex, nameToId[tmp]);
            }
            lastChild = nameToId[tmp];
        }
    }

    // Root the tree if necessary
    if (outDegree(guideTree, lastChild) > 2) {
        TVertexDescriptor myRoot = addVertex(guideTree);
        assignRoot(guideTree, myRoot);
        typedef typename Iterator<TGuideTree, OutEdgeIterator>::Type TOutEdgeIterator;
        TOutEdgeIterator it(guideTree, lastChild);
        goNext(it); goNext(it);
        TVertexDescriptor tV = targetVertex(it);
        TCargo c = cargo(*it);
        removeEdge(guideTree, lastChild, tV);
        addEdge(guideTree, myRoot, tV, (TCargo) (c / 2));
        addEdge(guideTree, myRoot, lastChild, (TCargo) (c / 2));
    }

    //std::fstream strm1; // Alignment graph as dot
    //strm1.open("tree23.dot", std::ios_base::out | std::ios_base::trunc);
    //writeRecords(strm1,guideTree,DotDrawing());
    //strm1.close();

    //std::cout << guideTree << std::endl;
}


//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TNames, typename TNewickString, typename TVertexDescriptor>
void
_buildNewickString(Graph<Tree<TCargo, TSpec> >& guideTree,
                   TNames& names,
                   TNewickString& str,
                   TVertexDescriptor v,
                   bool collapseRoot)
{
    typedef Graph<Tree<TCargo, TSpec> > TGuideTree;
    typedef typename EdgeDescriptor<TGuideTree>::Type TEdgeDescriptor;
    typedef typename Iterator<TGuideTree, AdjacencyIterator>::Type TAdjIter;
    if (isLeaf(guideTree, v)) {
        append(str, names[v], Generous());
    } else {
        if ((collapseRoot) && (isRoot(guideTree, v))) {
            TAdjIter adjIterRoot(guideTree, v);
            TVertexDescriptor v1 = *adjIterRoot;
            goNext(adjIterRoot);
            TVertexDescriptor v2 = *adjIterRoot;
            if (isLeaf(guideTree, v2)) {
                TVertexDescriptor tmp = v1;
                v1 = v2;
                v2 = tmp;
            }
            String<char> subStr;
            appendValue(subStr, '(', Generous());
            _buildNewickString(guideTree, names, subStr, v1, collapseRoot);
            TEdgeDescriptor e1 = findEdge(guideTree, v, v1);
            TEdgeDescriptor e2 = findEdge(guideTree, v, v2);
            appendValue(subStr, ':', Generous());
            std::ostringstream weight;
            weight << (cargo(e1) + cargo(e2));
            append(subStr, weight.str().c_str(), Generous());
            TAdjIter adjIter(guideTree, v2);
            for(;!atEnd(adjIter); goNext(adjIter)) {
                appendValue(subStr, ',', Generous());
                _buildNewickString(guideTree, names, subStr, *adjIter, collapseRoot);
                TEdgeDescriptor e = findEdge(guideTree, v2, *adjIter);
                appendValue(subStr, ':', Generous());
                std::ostringstream weight;
                weight << cargo(e);
                append(subStr, weight.str().c_str(), Generous());
            }
            appendValue(subStr, ')', Generous());
            append(str, subStr, Generous());
        } else {
            TAdjIter adjIter(guideTree, v);
            String<char> subStr;
            appendValue(subStr, '(', Generous());
            _buildNewickString(guideTree, names, subStr, *adjIter, collapseRoot);
            TEdgeDescriptor e = findEdge(guideTree, v, *adjIter);
            appendValue(subStr, ':', Generous());
            std::ostringstream weight;
            weight << cargo(e);
            append(subStr, weight.str().c_str(), Generous());
            goNext(adjIter);
            for(;!atEnd(adjIter); goNext(adjIter)) {
                appendValue(subStr, ',', Generous());
                _buildNewickString(guideTree, names, subStr, *adjIter, collapseRoot);
                TEdgeDescriptor e = findEdge(guideTree, v, *adjIter);
                appendValue(subStr, ':', Generous());
                std::ostringstream weight;
                weight << cargo(e);
                append(subStr, weight.str().c_str(), Generous());
            }
            appendValue(subStr, ')', Generous());
            append(str, subStr, Generous());
        }
    }
}


//////////////////////////////////////////////////////////////////////////////

template<typename TTarget, typename TCargo, typename TSpec, typename TNames>
inline void
writeRecords(TTarget & target,
             Graph<Tree<TCargo, TSpec> >& guideTree,
             TNames& names,
             bool collapseRoot,
             NewickFormat)
{
//IOREV _nodoc_ specialization not documented
    //typedef Graph<Tree<TCargo, TSpec> > TGuideTree;
    //typedef typename Size<TGuideTree>::Type TSize;
    //typedef typename VertexDescriptor<TGuideTree>::Type TVertexDescriptor;
    typename DirectionIterator<TTarget, Output>::Type iter = directionIterator(target, Output());

    String<char> myNewickString;
    _buildNewickString(guideTree, names, myNewickString, getRoot(guideTree), collapseRoot);
    write(iter, myNewickString);
    write(iter, ";\n");
}

//////////////////////////////////////////////////////////////////////////////

template<typename TTarget, typename TCargo, typename TSpec, typename TNames>
inline void
writeRecords(TTarget & target,
             Graph<Tree<TCargo, TSpec> >& guideTree,
             TNames& names,
             NewickFormat)
{
//IOREV _nodoc_ specialization not documented
    write(target, guideTree, names, false, NewickFormat());
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
