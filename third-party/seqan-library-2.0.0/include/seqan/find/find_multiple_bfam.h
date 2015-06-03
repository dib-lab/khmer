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

#ifndef SEQAN_HEADER_FIND_MULTIPLE_BFAM_H
#define SEQAN_HEADER_FIND_MULTIPLE_BFAM_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// MultiBfam
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class MultiBfamPattern
 * @extends Pattern
 * @headerfile <seqan/find.h>
 * @brief Multi-Pattern Backward Factor Automaton Matching.
 *
 * @signature template <typename TNeedle, typename TAutomaton>
 *            class Pattern<TNeedle, MultiBfam<TAutomaton> >;
 *
 * @tparam TAutomaton A tag that specifies the used automaton. Default: @link OracleMultiBfamPattern @endlink.
 * @tparam TNeedle    The needle type. Types: String
 *
 * @see BfamPattern
 */

/*!
 * @class OracleMultiBfamPattern
 * @extends MultiBfamPattern
 * @headerfile <seqan/find.h>
 * @brief Multi-Pattern Backward Factor Automaton Matching using an oracle automaton.
 *
 * @signature template <typename TNeedle>
 *            class Pattern<TNeedle, MultiBfam<Oracle> >;
 *
 * @tparam TNeedle The needle type. Types: String
 *
 * @see OracleBfamPattern
 */

//struct Oracle //defined in find_bom.h

template <typename TSpec = Oracle>
struct MultiBfam; //multiple backward factor automaton searching

typedef MultiBfam<Oracle> SBomAlgo; //deprecated

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec>
class Pattern<TNeedle, MultiBfam<TSpec> >
{
//____________________________________________________________________________
public:
    typedef typename Value<TNeedle>::Type TKeyword;
    typedef typename Position<TNeedle>::Type TNeedlePosition;
    typedef typename Size<TKeyword>::Type TSize;
    typedef typename Value<TKeyword>::Type TValue;
    typedef Graph<Automaton<TValue, void, WithoutEdgeId> > TGraph;

    //searching data: these members are initialized in _patternInit or during search
    TNeedlePosition * position;        //pointer to last found position
    TNeedlePosition * position_end; //end of list in verify
            //note: if to_verify_begin == to_verify_end then searching in Haystack must go on

    //preprocessed data: these members are initialized in setHost
    Holder<TNeedle> needle;
    TGraph automaton; //automaton of the reverse lmin-prefixes of the keywords
    String<String<TNeedlePosition> > terminals; //map of terminal states in automaton

    TSize lmin;    //min length of keyword

//____________________________________________________________________________

    Pattern():
        lmin(0)
    {
    }

    template <typename TNeedle2>
    Pattern(TNeedle2 const & ndl)
    {
        SEQAN_CHECKPOINT
        setHost(*this, ndl);
    }

    ~Pattern()
    {
    }
//____________________________________________________________________________

private:
    Pattern(Pattern const& other);
    Pattern const & operator=(Pattern const & other);

//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////


template <typename TNeedle, typename TStrs>
inline void
_buildAutomatonMultiBfam(Pattern<TNeedle, MultiBfam<Oracle> > & me,
                          TStrs const & strs)
{
    createSetOracle(me.automaton, me.terminals, strs);
}

template <typename TNeedle, typename TStrs>
inline void
_buildAutomatonMultiBfam(Pattern<TNeedle, MultiBfam<Trie> > & me,
                          TStrs const & strs)
{
    createSetSuffixTrie(me.automaton, me.terminals, strs);
}

//____________________________________________________________________________

template <typename TNeedle, typename TAutomaton, typename TNeedle2>
void _setHostMultiBfam(Pattern<TNeedle, MultiBfam<TAutomaton> > & me,
                        TNeedle2 const & needle_)
{
    typedef typename Value<TNeedle>::Type TKeyword;
    typedef typename Value<TKeyword>::Type TValue;
    typedef typename Size<TKeyword>::Type TSize;

    //me.needle
    setValue(me.needle, needle_);
    TNeedle & ndl = needle(me);

    //determine lmin
    TSize len = length(ndl);
    if (len == 0)
    {
        me.lmin = 0;
        return;
    }

    me.lmin = length(ndl[0]);
    for (TSize i = 1; i < len; ++i)
    {
        TSize len = length(ndl[i]);
        if (len < me.lmin)
        {
            me.lmin = len;
        }
    }

    if (me.lmin == 0) return;

    //collect reverse prefixes for automaton
    String<String<TValue> > strs;
    resize(strs, len);
    for (unsigned int i = 0; i < len; ++i)
    {
        strs[i] = prefix(ndl[i], me.lmin);
        reverse(strs[i]);
    }

    //build automaton
    _buildAutomatonMultiBfam(me, strs);
}

template <typename TNeedle, typename TAutomaton, typename TNeedle2>
void setHost (Pattern<TNeedle, MultiBfam<TAutomaton> > & me,
              TNeedle2 const & needle)
{
    _setHostMultiBfam(me, needle);
}

template <typename TNeedle, typename TAutomaton, typename TNeedle2>
inline void
setHost(Pattern<TNeedle, MultiBfam<TAutomaton> > & me,
        TNeedle2 & needle)
{
    _setHostMultiBfam(me, needle);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TAutomaton>
inline typename Host<Pattern<TNeedle, MultiBfam<TAutomaton> > >::Type &
host(Pattern<TNeedle, MultiBfam<TAutomaton> > & me)
{
SEQAN_CHECKPOINT
    return value(me.needle);
}

template <typename TNeedle, typename TAutomaton>
inline typename Host<Pattern<TNeedle, MultiBfam<TAutomaton> > const>::Type &
host(Pattern<TNeedle, MultiBfam<TAutomaton> > const & me)
{
SEQAN_CHECKPOINT
    return value(me.needle);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TAutomaton>
inline typename Size<TNeedle>::Type
position(Pattern<TNeedle, MultiBfam<TAutomaton> > & me)
{
    return *(me.position);
}

//////////////////////////////////////////////////////////////////////////////

//called when search begins
template <typename TNeedle, typename TAutomaton>
inline void _patternInit (Pattern<TNeedle, MultiBfam<TAutomaton> > & me)
{
SEQAN_CHECKPOINT
    me.position = 0;
    me.position_end = 0;
}

//////////////////////////////////////////////////////////////////////////////

//test whether the me.lmin prefix of the *(me.position)-th pattern matches haystack beginning at *it
//the Haystack is long enough that *(it+me.lmin) is a valid value.
//default implementation: it is assumed that if me.automaton parses a me.lmin-length string S, then S is a pattern
template <typename TNeedle, typename TAutomaton, typename THaystackIterator>
inline bool
_startVerifyMultiBfam(Pattern<TNeedle, MultiBfam<TAutomaton> > &,
                       THaystackIterator)
{
    return true;
}

//specialization for oracles: must test explicitely, since set-oracles may also parse me.lmin-length strings that are no patterns
template <typename TNeedle, typename THaystackIterator>
inline bool
_startVerifyMultiBfam(Pattern<TNeedle, MultiBfam<Oracle> > & me,
                       THaystackIterator tit)
{
    typedef typename Value<TNeedle>::Type TKeyword;
    typedef typename Iterator<TKeyword, Standard>::Type TKeywordIterator;
    TKeyword & keyword = needle(me)[*(me.position)];
    TKeywordIterator kit = begin(keyword);
    TKeywordIterator kit_end = kit + me.lmin;
    while (kit != kit_end)
    {
        if (*kit != *tit)
        {
            return false;
        }
        ++kit;
        ++tit;
    }

    return true;
}

//____________________________________________________________________________

template <typename TFinder, typename TAutomaton, typename TNeedle>
inline bool find(TFinder & finder,
                 Pattern<TNeedle, MultiBfam<TAutomaton> > & me)
{
SEQAN_CHECKPOINT
    typedef typename Haystack<TFinder>::Type THaystack;
    typedef typename Iterator<THaystack, Standard>::Type THaystackIterator;
    typedef typename Value<TNeedle>::Type TKeyword;
    typedef typename Value<TKeyword>::Type TValue;
    typedef typename Size<TKeyword>::Type TSize;
    typedef typename Iterator<TKeyword, Standard>::Type TKeywordIterator;
    typedef Graph<Automaton<TValue, void, WithoutEdgeId> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

    if (me.lmin == 0) return false;

    //some variables
    THaystackIterator haystack_end = end(haystack(finder));
    THaystackIterator it1;
    THaystackIterator it1_end = haystack_end - me.lmin + 1;
    THaystackIterator it2;
    TGraph & automaton = me.automaton;
    TVertexDescriptor root = getRoot(automaton);
    TVertexDescriptor nil_ = getNil<TVertexDescriptor>();
    TVertexDescriptor current;
    TKeywordIterator kit;
    TKeywordIterator kit_end;
    THaystackIterator tit;
    TKeyword * p_keyword;
    TSize len;

    if (empty(finder))
    {
//START
        _patternInit(me);
        _finderSetNonEmpty(finder);
        it1 = hostIterator(finder);
    }
    else
    {
//RESUME
        it1 = hostIterator(finder);
        goto VERIFY_NEXT;
    }

//SEARCH
    while (it1 < it1_end)
    {
        it2 = it1 + me.lmin; //me.lmin > 0 => it1 != it2
        current = root;
        while (true)
        {
            --it2;
            current = getSuccessor(automaton, current, *it2);
            if (current == nil_)
            {
//SKIP
                it1 = it2 + 1;
                break;
            }
            if (it1 == it2)
            {
//VERIFY
                me.position = begin(property(me.terminals, current), Standard());
                me.position_end = end(property(me.terminals, current), Standard());
                SEQAN_ASSERT_NEQ(me.position, me.position_end);

                if (_startVerifyMultiBfam(me, it1)) //this returns true if the lmin-length prefixe matches
                {
                    while (me.position != me.position_end)
                    {
                        p_keyword = & needle(me)[*me.position];
                        len = length(*p_keyword);
                        if ((it1 + len) <= haystack_end)
                        {
                            //compare rest
                            kit = begin(*p_keyword) + me.lmin;
                            kit_end = end(*p_keyword);
                            tit = it1 + me.lmin;
                            while (true)
                            {
                                if (kit == kit_end)
                                {
//MATCH FOUND
                                    setPosition(finder, it1 - begin(haystack(finder), Standard()));
                                    _setFinderLength(finder, length(needle(*p_keyword)));
                                    _setFinderEnd(finder, position(finder) + length(finder));
                                    return true;
                                }
                                if (*kit != *tit) break;
                                ++kit;
                                ++tit;
                            }
                        }
VERIFY_NEXT:
                        ++me.position;
                    }
                }
                ++it1;
                break;
            }

        }
    }
    return false;
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
