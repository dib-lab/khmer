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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_INDEX_SA_QSORT_H
#define SEQAN_HEADER_INDEX_SA_QSORT_H

namespace SEQAN_NAMESPACE_MAIN
{

    struct SAQSort {};

    // compare two suffices of a given text
    template < typename TSAValue, typename TText >
    struct SuffixLess_ :
        public std::binary_function < TSAValue, TSAValue, bool >
    {
        typedef typename Iterator<TText const, Standard>::Type TIter;
        TIter _begin, _end;

        SuffixLess_(TText const &text):
            _begin(begin(text, Standard())),
            _end(end(text, Standard())) {}

        // skip the first <offset> characters
        template <typename TSize>
        SuffixLess_(TText const &text, TSize offset):
            _begin(begin(text, Standard()) + offset),
            _end(end(text, Standard())) {}

        inline bool operator() (TSAValue const a, TSAValue const b) const
        {
            if (a == b) return false;
            TIter itA = _begin + a;
            TIter itB = _begin + b;
            if (a <= b) {
                for(; itB != _end; ++itB, ++itA) {
                    if (ordLess(*itA, *itB)) return true;
                    if (ordLess(*itB, *itA)) return false;
                }
                return false;
            } else {
                for(; itA != _end; ++itA, ++itB) {
                    if (ordLess(*itA, *itB)) return true;
                    if (ordLess(*itB, *itA)) return false;
                }
                return true;
            }
        }
    };

    // compare two suffices of a given text
    template < typename TSAValue, typename TString, typename TSetSpec >
    struct SuffixLess_<TSAValue, StringSet<TString, TSetSpec> const > :
        public std::binary_function < TSAValue, TSAValue, bool >
    {
        typedef StringSet<TString, TSetSpec> const TText;

        typename Size<TString>::Type _offset;
        TText &_text;

        SuffixLess_(TText &text):
            _text(text) {}

        // skip the first <offset> characters
        template <typename TSize>
        SuffixLess_(TText &text, TSize offset):
            _text(text),
            _offset(offset) {}

        inline bool operator() (TSAValue const a, TSAValue const b) const
        {
            typedef typename Iterator<TString const, Standard>::Type TIter;
            if (a == b) return false;
            TIter itA = begin(getValue(_text, getSeqNo(a)), Standard()) + getSeqOffset(a);
            TIter itB = begin(getValue(_text, getSeqNo(b)), Standard()) + getSeqOffset(b);
            TIter itAEnd = end(getValue(_text, getSeqNo(a)), Standard());
            TIter itBEnd = end(getValue(_text, getSeqNo(b)), Standard());
            if (itAEnd - itA < itBEnd - itB) {
                for(; itA != itAEnd; ++itA, ++itB) {
                    if (ordLess(*itA, *itB)) return true;
                    if (ordLess(*itB, *itA)) return false;
                }
                return true;
            } else {
                for(; itB != itBEnd; ++itB, ++itA) {
                    if (ordLess(*itA, *itB)) return true;
                    if (ordLess(*itB, *itA)) return false;
                }
                if (itA != itAEnd)
                    return false;
                return getSeqNo(a) > getSeqNo(b);
            }
        }
    };

    // compare two suffices of a given text and skip the first <lcp> characters
    template < typename TSAValue, typename TText >
    struct SuffixLessOffset_: SuffixLess_<TSAValue, TText>
    {
        // skip the first <offset> characters
        template <typename TSize>
        SuffixLessOffset_(TText &text, TSize offset):
            SuffixLess_<TSAValue, TText> (text, offset) {}
    };



    template < typename TSA,
            typename TText,
            typename TSize>
    void _sortBucketQuickSort(
        TSA &sa,
        TText &text,
        TSize lcp)
    {
    SEQAN_CHECKPOINT
        // sort bucket with quicksort
        std::sort(
            begin(sa, Standard()),
            end(sa, Standard()),
            SuffixLessOffset_<typename Value<TSA>::Type, TText>(text, lcp));
    }

    template < typename TSA,
               typename TText >
    inline void createSuffixArray(
        TSA &SA,
        TText const &s,
        SAQSort const &)
    {
    SEQAN_CHECKPOINT
        typedef typename Size<TSA>::Type TSize;
        typedef typename Iterator<TSA, Standard>::Type TIter;

        // 1. Fill suffix array with a permutation (the identity)
        TIter it = begin(SA, Standard());
        TIter itEnd = end(SA, Standard());
        for(TSize i = 0; it != itEnd; ++it, ++i)
            *it = i;

        // 2. Sort suffix array with quicksort
        std::sort(
            begin(SA, Standard()),
            end(SA, Standard()),
            SuffixLess_<typename Value<TSA>::Type, TText const>(s));
    }

    template < typename TSA,
               typename TString,
               typename TSSetSpec >
    inline void createSuffixArray(
        TSA &SA,
        StringSet< TString, TSSetSpec > const &s,
        SAQSort const &)
    {
    SEQAN_CHECKPOINT
        typedef StringSet< TString, TSSetSpec > TText;
        typedef typename Size<TSA>::Type TSize;
        typedef typename Iterator<TSA, Standard>::Type TIter;

        // 1. Fill suffix array with a permutation (the identity)
        TIter it = begin(SA, Standard());
        for(unsigned j = 0; j < length(s); ++j)
        {
            TSize len = length(s[j]);
            for(TSize i = 0; i < len; ++i, ++it)
                *it = Pair<unsigned, TSize>(j, i);
        }

        // 2. Sort suffix array with quicksort
        std::sort(
            begin(SA, Standard()),
            end(SA, Standard()),
            SuffixLess_<typename Value<TSA>::Type, TText const>(s));
    }

    //////////////////////////////////////////////////////////////////////////////
    // suffix quicksort pipe
    template < typename TInput >
    struct Pipe< TInput, SAQSort >
    {
        typedef typename Value<TInput>::Type    TValue;
        typedef typename SAValue<TInput>::Type    TSAValue;

        typedef String<TValue, Alloc<> >        TText;
        typedef String<TSAValue, Alloc<> >        TSA;
        typedef Pipe<TSA, Source<> >            TSource;

        TSA        sa;
        TSource    in;

        Pipe(TInput &_textIn):
            in(sa)
        {
            TText text;
            text << _textIn;

            resize(sa, length(_textIn), Exact());
            createSuffixArray(sa, text, SAQSort());
        }

        inline typename Value<TSource>::Type const & operator*() {
            return *in;
        }

        inline Pipe& operator++() {
            ++in;
            return *this;
        }
    };

}

#endif
