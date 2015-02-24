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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_FIND_SIMPLE_H
#define SEQAN_HEADER_FIND_SIMPLE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Simple Finder
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class SimplePattern
 * @extends Pattern
 * @headerfile <seqan/find.h>
 * @brief A brute force online searching algorithm.
 *
 * @signature template <typename TNeedle>
 *            class Pattern<TNeedle, Simple>;
 *
 * @tparam TNeedle The needle type.  Types: @link ContainerConcept @endlink.
 *
 * This specialization should only be used if no other is applicable.
 */

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
class Pattern<TNeedle, Simple> {
//____________________________________________________________________________
public:
    Holder<TNeedle> data_host;

//____________________________________________________________________________

    Pattern() {}

    Pattern(Pattern const & other):
        data_host(other.data_host)
    {
    }

    template <typename TNeedle2>
    Pattern(TNeedle2 const & ndl)
    {
        setHost(*this, ndl);
    }

    ~Pattern(){}

    Pattern const &
    operator = (Pattern const & other)
    {
        data_host = other.data_host;
        return *this;
    }

//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TNeedle2>
void setHost (Pattern<TNeedle, Simple> & me,
              TNeedle2 & needle)
{
    setValue(me.data_host, needle);
}
template <typename TNeedle, typename TNeedle2>
void setHost (Pattern<TNeedle, Simple> & me,
              TNeedle2 const & needle)
{
    setValue(me.data_host, needle);
}

//____________________________________________________________________________


template <typename TFinder, typename TNeedle>
inline bool find(TFinder & finder,
                 Pattern<TNeedle, Simple> & me)
{
    typedef typename Haystack<TFinder>::Type THaystack;
    typedef typename Parameter_<THaystack const>::Type TParamHaystack;
    typedef typename Iterator<THaystack const, Standard>::Type THaystackIterator;

    if (empty(finder))
    {
        _setFinderLength(finder, length(needle(me)));
        _finderSetNonEmpty(finder);
    }
    else ++finder;

    TParamHaystack hstk = haystack(finder);
    TNeedle const & ndl = needle(me);

    THaystackIterator res = std::search(begin(hstk, Standard())+position(finder), end(hstk, Standard()), begin(ndl, Standard()), end(ndl, Standard()));

    if (res == end(hstk, Standard())) return false;

    _setFinderEnd(finder, (res - begin(hstk, Standard())) + length(ndl));
    setPosition(finder, beginPosition(finder));
    return true;

/*
    TSize n = length(hstk);
    TSize m = length(ndl);
    while (position(finder)+m <= n)
    {
        if (ndl == infix(hstk, position(finder), position(finder)+m))
        {
            _setFinderEnd(finder, position(finder)+m);
            return true;
        }
        ++finder;
    }
    return false;
*/
}

//____________________________________________________________________________

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
