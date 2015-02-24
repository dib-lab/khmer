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

#ifndef SEQAN_HEADER_FIND_HORSPOOL_H
#define SEQAN_HEADER_FIND_HORSPOOL_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Horspool
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class HorspoolPattern
 * @extends Pattern
 * @headerfile <seqan/find.h>
 *
 * @brief Exact string matching using Horspool's algorithm (1980).
 *
 * @signature template <typename TNeedle>
 *            class Pattern<TNeedle, Horspool>;
 *
 * @tparam TNeedle The needle type. Types: String
 */

struct Horspool_;
typedef Tag<Horspool_> Horspool;

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
class Pattern<TNeedle, Horspool>
{
//____________________________________________________________________________

public:
    typedef typename Size<TNeedle>::Type TSize;

    Holder<TNeedle>        data_host;
    String<TSize>        data_map;

//____________________________________________________________________________

public:
    Pattern() {}

    template <typename TNeedle2>
    Pattern(TNeedle2 const & ndl)
    {
        setHost(*this, ndl);
    }
//____________________________________________________________________________
};


template <typename TNeedle, typename TNeedle2>
void
setHost(Pattern<TNeedle, Horspool> & me, TNeedle2 const & ndl)
{
    typedef typename Value<TNeedle>::Type TValue;
    typedef typename Size<TNeedle>::Type TSize;

    SEQAN_ASSERT_NOT(empty(ndl));

    TSize value_size = ValueSize<TValue>::VALUE;

    //make room for map
    resize(me.data_map, value_size);

    //fill map
    typename Value<String<TSize> >::Type jump_width = length(ndl); //das ist so umstaendlich wegen VC++ 2003
    arrayFill(begin(me.data_map, Standard()), begin(me.data_map, Standard()) + value_size, jump_width);

    typename Iterator<TNeedle2 const, Standard>::Type it;
    it = begin(ndl, Standard());
    while (jump_width > 1)
    {
        --jump_width;
        unsigned int pos_ = *it; //conversion value type to unsigned int
        me.data_map[pos_] = jump_width;
        ++it;
    }

    me.data_host = ndl;
}

template <typename TNeedle, typename TNeedle2>
void
setHost(Pattern<TNeedle, Horspool> & horsp, TNeedle2 & ndl)
{
    setHost(horsp, reinterpret_cast<TNeedle2 const &>(ndl));
}

//____________________________________________________________________________


template <typename TNeedle>
inline void _patternInit (Pattern<TNeedle, Horspool> &) {}


//____________________________________________________________________________

template <typename TFinder, typename TNeedle2>
bool
_findHorspool(TFinder & finder,
              Pattern<TNeedle2, Horspool> & me,
              bool find_first)
{
SEQAN_CHECKPOINT
    typedef typename Haystack<TFinder>::Type THaystack;
    typedef typename Parameter_<THaystack>::Type TParamHaystack;

    TParamHaystack hayst = haystack(finder);

    typedef Pattern<TNeedle2, Horspool> TPattern;
    typedef typename Needle<TPattern>::Type TNeedle;
    TNeedle & ndl = needle(me);

    typedef typename Size<TNeedle>::Type TNeedleSize;
    TNeedleSize ndl_size = length(ndl);

    typedef typename Iterator<THaystack, Standard>::Type THaystackIterator;
    THaystackIterator haystack_end = end(hayst, Standard());
    THaystackIterator it = begin(hayst, Standard());
    it += position(finder) + ndl_size - 1; //it points to the last character
    THaystackIterator it_next = it;

    typedef typename Iterator<TNeedle, Standard>::Type TNeedleIterator;
    TNeedleIterator nit; //needle iterator
    TNeedleIterator nit_begin = begin(ndl, Standard());
    TNeedleIterator nit_end = end(ndl, Standard()) - 1; //here the verification begins

    unsigned int char_i;

    if (find_first)
    {
        goto VALIDATE;
    }

MOVE_FURTHER:
    //move to next position
    char_i = *it; //conversion to unsigned integer
    it_next = it + me.data_map[char_i];
    if (it_next >= haystack_end)
    {//found nothing
        return false;
    }

    it = it_next;

VALIDATE:
    //validate current position
    for (nit = nit_end; nit >= nit_begin; --nit)
    {
        if (*nit != *it_next)
        {//invalid!
            goto MOVE_FURTHER;
        }
        --it_next;
    }

    //valid! return hit
    _setFinderEnd(finder, it - begin(hayst, Standard()) + 1);
    setPosition(finder, beginPosition(finder));
    return true;
}

//____________________________________________________________________________
// Sentinel variant (not used at the moment)
//TODO: if not enough space at the end of the haystack: call non-sentinel search
/*
template <typename TFinder, typename TNeedle2>
bool
find_horspool_sentinel(TFinder & finder,
                       Pattern<TNeedle2, Horspool> & me,
                       bool find_first)
{
SEQAN_CHECKPOINT
    typedef typename Haystack<TFinder>::Type THaystack;
    THaystack & hayst = haystack(finder);


    typedef Pattern<TNeedle2, Horspool> TPattern;
    typedef typename Needle<TPattern>::Type TNeedle;
    TNeedle & ndl = needle(me);

    //implant sentinel
    typename Size<THaystack>::Type old_haystack_size = length(hayst);
    if (find_first)
    {
        typedef typename Position<TFinder>::Type TFinderPosition;
        TFinderPosition finder_pos = position(finder);

        append(hayst, ndl, Exact());
        if (length(hayst) != old_haystack_size + length(ndl))
        {//not enough place in haystack
//TODO!!!
printf("error!");
return false;
        }
        setPosition(finder, finder_pos);
    }
    else
    {
        _setLength(hayst, old_haystack_size + length(ndl));
    }

    typedef typename Size<TNeedle>::Type TNeedleSize;
    TNeedleSize ndl_size = length(ndl);

    typedef typename Iterator<THaystack, Standard>::Type THaystackIterator;
    THaystackIterator it = begin(hayst, Standard());
    THaystackIterator haystack_end = it + old_haystack_size;
    it += position(finder) + ndl_size - 1; //it points to the last character

    typedef typename Iterator<TNeedle, Standard>::Type TNeedleIterator;
    TNeedleIterator nit; //needle iterator
    TNeedleIterator nit_begin = begin(ndl, Standard());
    TNeedleIterator nit_end = end(ndl, Standard()) - 1; //here the verification begins

    typedef typename Value<TNeedle>::Type TNeedleValue;
    TNeedleValue char_needle_last = *nit_end;
    TNeedleValue char_haystack_last;

    char_haystack_last = *it;

    if (find_first)
    {
        goto VALIDATE;
    }

    //main loop
MOVE_FURTHER:
    it += me.data_map[_ord(char_haystack_last)];
    char_haystack_last = *it;
    if (char_haystack_last != char_needle_last) goto MOVE_FURTHER;


    if (it >= haystack_end)
    {//found nothing
        resize(hayst, old_haystack_size);
        return false;
    }

VALIDATE:
    //validate current position
    THaystackIterator it_back = it;
    for (nit = nit_end; nit >= nit_begin; --nit)
    {
        if (*nit != *it_back)
        {//invalid!
            goto MOVE_FURTHER;
        }
        --it_back;
    }

    //valid! return hit
    setPosition(finder, it - begin(hayst, Standard()) - ndl_size + 1);
    resize(hayst, old_haystack_size);
    return true;
}
*/
//____________________________________________________________________________
//spec for file reader haystacks

template <typename TFormat, typename TFile, typename TSpec>
struct FileReader;
//IOREV

template <typename TValue, typename TFormat, typename TFile, typename FileReaderTSpec, typename TFinderSpec, typename TNeedle2>
bool
_findHorspool(Finder<String<TValue, FileReader<TFormat, TFile, FileReaderTSpec> >, TFinderSpec > & finder,
              Pattern<TNeedle2, Horspool> & me,
              bool find_first)
{
SEQAN_CHECKPOINT
    typedef Finder<String<TValue, FileReader<TFormat, TFile, FileReaderTSpec> >, TFinderSpec > TFinder;
    typedef typename Haystack<TFinder>::Type THaystack;
    THaystack & hayst = haystack(finder);

    typedef Pattern<TNeedle2, Horspool> TPattern;
    typedef typename Needle<TPattern>::Type TNeedle;
    TNeedle & ndl = needle(me);

    typedef typename Size<TNeedle>::Type TNeedleSize;
    TNeedleSize ndl_size = length(ndl);

    typedef typename Iterator<THaystack, Standard>::Type THaystackIterator;
    THaystackIterator it(hayst, position(finder) + ndl_size - 1); //it points to the last character

    typedef typename Iterator<TNeedle, Standard>::Type TNeedleIterator;
    TNeedleIterator nit; //needle iterator
    TNeedleIterator nit_begin = begin(ndl, Standard());
    TNeedleIterator nit_end = end(ndl, Standard()) - 1; //here the verification begins

    unsigned int char_i;

    if (find_first)
    {
        goto VALIDATE;
    }

MOVE_FURTHER:
    //move to next position
    char_i = *it; //conversion to unsigned integer
    it += me.data_map[char_i];
    if (atEnd(it))
    {//found nothing
        return false;
    }

VALIDATE:
    //validate current position
    for (nit = nit_end; nit >= nit_begin; --nit)
    {
        if (*nit != *it)
        {//invalid!
            it += (nit_end - nit);
            goto MOVE_FURTHER;
        }
        --it;
    }

    //valid! return hit
    setPosition(finder, it - begin(hayst, Standard()) + 1);
    return true;

}

//____________________________________________________________________________
/* groepl variante
(Beruht vermutlich auf einem Missverstehen von Navarro/Raffinot Seite 26:
Mit "the main loop can be 'unrolled'" ist dort naemlich
"the INNER loop can be 'unrolled'" gemeint.)

template <typename TFinder, typename TNeedle2>
bool
_findHorspool(TFinder & finder,
    Pattern<TNeedle2, Horspool> & me,
    bool find_first)
{
SEQAN_CHECKPOINT
    typedef typename Haystack<TFinder>::Type THaystack;
    THaystack & hayst = haystack(finder);

    typedef Pattern<TNeedle2, Horspool> TPattern;
    typedef typename Needle<TPattern>::Type TNeedle;
    TNeedle & ndl = needle(me);

    typename Size<TNeedle>::Type ndl_size = length(ndl);

    typedef typename Iterator<THaystack, Standard>::Type THaystackIterator;
    THaystackIterator haystack_end = end(hayst, Standard());
    THaystackIterator it = begin(hayst, Standard());
    it += position(finder) + ndl_size - 1; //it points to the last character
    THaystackIterator it2;

    typedef typename Iterator<TNeedle, Standard>::Type TNeedleIterator;
    TNeedleIterator nit; //needle iterator
    TNeedleIterator nit_begin = begin(ndl, Standard());
    TNeedleIterator nit_end = end(ndl, Standard()) - 2; //here the verification begins

    typedef typename Value<TNeedle>::Type TNeedleValue;
    TNeedleValue last_needle_char = value(ndl, ndl_size - 1);

    unsigned int char_i;

    if (!find_first)
    {
        ++it;
    }

MOVE_FURTHER:
    //scan for the last character
    while (true)
    {
        if (it >= haystack_end)
        {//found nothing
            return false;
        }
        if (*it == last_needle_char)
        {
            break;
        }
        ++it;
    }

VALIDATE:
    it2 = it;
    //validate current position
    for (nit = nit_end; nit >= nit_begin; --nit)
    {
        --it2;
        if (*nit != *it2)
        {//invalid! skip!
            char_i = *it; //conversion to unsigned integer
            it += me.data_map[char_i];
            goto MOVE_FURTHER;
        }
    }

    //valid! return hit
    setPosition(finder, it - begin(hayst, Standard()) - ndl_size + 1);
    return true;
}
*/

template <typename TFinder, typename TNeedle2>
bool
find(TFinder & finder, Pattern<TNeedle2, Horspool> & me)
{
SEQAN_CHECKPOINT
    bool find_first = empty(finder);
    if (find_first)
    {
        _patternInit(me);
        _setFinderLength(finder, length(needle(me)));
        _finderSetNonEmpty(finder);
    }

    SEQAN_ASSERT_GT(length(needle(me)), 0u);

    return _findHorspool(finder, me, find_first);
}

//////////////////////////////////////////////////////////////////////////////
// Host
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
struct Host< Pattern<TNeedle, Horspool> >
{
    typedef TNeedle Type;
};

template <typename TNeedle>
struct Host< Pattern<TNeedle, Horspool> const>
{
    typedef TNeedle const Type;
};


//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
