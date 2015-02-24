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
// Definition of the class Pattern and supporting functions and
// metafunctions.
// ==========================================================================

#ifndef SEQAN_HEADER_FIND_PATTERN_BASE_H
#define SEQAN_HEADER_FIND_PATTERN_BASE_H


//////////////////////////////////////////////////////////////////////////////

namespace seqan {

//////////////////////////////////////////////////////////////////////////////

/*!
 * @class Pattern
 * @headerfile <seqan/find.h>
 * @brief Holds the needle and preprocessing data (depends on algorithm).
 *
 * @signature template <typename TNeedle[, typename TSpec]>
 *            class Pattern;
 *
 * @tparam TNeedle The needle type.  Types: @link TextConcept @endlink.
 * @tparam TSpec   A tag that specifies the online algorithm to use for the search.  Defaults to the result of
 *                 @link DefaultPattern @endlink.
 *
 * If <tt>Needle</tt> is a StringSet then <tt>position(pattern)</tt> returns a @link Pair @endlink with the index of the currently
 * matching needle and the position in the needle.
 */

template < typename TNeedle, typename TSpec = typename DefaultPattern<TNeedle>::Type >
class Pattern;

//default implementation
template < typename TNeedle >
class Pattern<TNeedle, void>
{
public:
    typedef typename Position<TNeedle>::Type TNeedlePosition;

    Holder<TNeedle> data_host;
    TNeedlePosition data_begin_position;
    TNeedlePosition data_end_position;

    Pattern() : data_begin_position(), data_end_position()
    {}

    template <typename TNeedle_>
    Pattern(TNeedle_ & ndl) : data_host(ndl), data_begin_position(), data_end_position()
    {}

    template <typename TNeedle_>
    Pattern(TNeedle_ const & ndl) : data_host(ndl), data_begin_position(), data_end_position()
    {}

};
//////////////////////////////////////////////////////////////////////////////

/*!
 * @mfn Pattern#Container
 * @brief Returns the needle type of the pattern.
 *
 * @signature Container<TPattern>::Type;
 *
 * @tparam TPattern The pattern to query for its needle type.
 *
 * @return Type The needle type.
 */

template <typename TNeedle, typename TSpec>
struct Container< Pattern<TNeedle, TSpec> > {
    typedef TNeedle Type;
};

template <typename TNeedle, typename TSpec>
struct Container< Pattern<TNeedle, TSpec> const > {
    typedef TNeedle const Type;
};

/*!
 * @mfn Pattern#Host
 * @brief Returns the host type of the pattern.
 *
 * @signature Host<TPattern>::Type;
 *
 * @tparam TPattern The pattern to query for its host type.
 *
 * @return Type The host type.
 */

template <typename TNeedle, typename TSpec>
struct Host< Pattern<TNeedle, TSpec> > {
    typedef TNeedle Type;
};

template <typename TNeedle, typename TSpec>
struct Host< Pattern<TNeedle, TSpec> const > {
    typedef TNeedle const Type;
};

/*!
 * @mfn Pattern#Value
 * @brief Returns the value type of the underlying pattern.
 *
 * @signature Value<TPattern>::Type;
 *
 * @tparam TPattern The Pattern to query.
 *
 * @return Type The value type.
 */

template <typename TPattern, typename TSpec>
struct Value< Pattern<TPattern, TSpec> > {
    typedef typename Value<TPattern>::Type Type;
};

/*!
 * @mfn Pattern#Position
 * @brief Returns the position type of the underlying pattern.
 *
 * @signature Position<TPattern>::Type;
 *
 * @tparam TPattern The Pattern to query.
 *
 * @return Type The position type.
 */

template <typename TPattern, typename TSpec>
struct Position< Pattern<TPattern, TSpec> > {
    typedef typename Position<TPattern>::Type Type;
};

/*!
 * @mfn Pattern#Difference
 * @brief Returns the difference type of the underlying pattern.
 *
 * @signature Difference<TPattern>::Type;
 *
 * @tparam TPattern The Pattern to query.
 *
 * @return Type The difference type.
 */

template <typename TPattern, typename TSpec>
struct Difference< Pattern<TPattern, TSpec> > {
    typedef typename Difference<TPattern>::Type Type;
};

/*!
 * @mfn Pattern#Size
 * @brief Returns the size type of the underlying pattern.
 *
 * @signature Size<TPattern>::Type;
 *
 * @tparam TPattern The Pattern to query.
 *
 * @return Type The size type.
 */

template <typename TPattern, typename TSpec>
struct Size< Pattern<TPattern, TSpec> > {
    typedef typename Size<TPattern>::Type Type;
};


//////////////////////////////////////////////////////////////////////////////

/*!
 * @mfn Pattern#ScoringScheme
 * @brief Returns the scoring scheme type of an approximate search algorithm.
 *
 * @signature ScoringScheme<TPattern>::Type;
 *
 * @tparam TPattern The Pattern to query for its scoring scheme type.  Default: EditDistanceScore.
 */

template <typename TNeedle>
struct ScoringScheme
{
    typedef EditDistanceScore Type;
};
template <typename TNeedle>
struct ScoringScheme<TNeedle const>:
    ScoringScheme<TNeedle>
{
};

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec>
inline Holder<TNeedle> &
_dataHost(Pattern<TNeedle, TSpec> & me)
{
    return me.data_host;
}
template <typename TNeedle, typename TSpec>
inline Holder<TNeedle> &
_dataHost(Pattern<TNeedle, TSpec> const & me)
{
    return const_cast<Holder<TNeedle> &>(me.data_host);
}

//host access: see basic_host.h


//???TODO: Diese Funktion entfernen! (sobald setHost bei anderen pattern nicht mehr eine Art "assignHost" ist)
template <typename TNeedle, typename TSpec, typename TNeedle2>
inline void
setHost(Pattern<TNeedle, TSpec> & me,
        TNeedle2 const & ndl)
{
     me.data_host = ndl; //assign => Pattern haelt eine Kopie => doof!
}
template <typename TNeedle, typename TSpec, typename TNeedle2>
inline void
setHost(Pattern<TNeedle, TSpec> & me,
        TNeedle2 & ndl)
{
     me.data_host = ndl; //assign => Pattern haelt eine Kopie => doof!
}
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec>
inline typename Position<Pattern<TNeedle, TSpec> >::Type &
beginPosition(Pattern<TNeedle, TSpec> & me)
{
    return me.data_begin_position;
}
template <typename TNeedle, typename TSpec>
inline typename Position<Pattern<TNeedle, TSpec> const >::Type &
beginPosition(Pattern<TNeedle, TSpec> const & me)
{
    return me.data_begin_position;
}


template <typename TNeedle, typename TSpec, typename TPosition>
inline void
setBeginPosition(Pattern<TNeedle, TSpec> & me,
                 TPosition _pos)
{
    me.data_begin_position = _pos;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec>
inline typename Position<Pattern<TNeedle, TSpec> >::Type &
endPosition(Pattern<TNeedle, TSpec> & me)
{
    return me.data_end_position;
}
template <typename TNeedle, typename TSpec>
inline typename Position<Pattern<TNeedle, TSpec> const >::Type &
endPosition(Pattern<TNeedle, TSpec> const & me)
{
    return me.data_end_position;
}

template <typename TNeedle, typename TSpec, typename TPosition>
inline void
setEndPosition(Pattern<TNeedle, TSpec> & me,
               TPosition _pos)
{
    me.data_end_position = _pos;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec>
inline typename Infix<TNeedle>::Type
segment(Pattern<TNeedle, TSpec> & me)
{
    typedef typename Infix<TNeedle>::Type TInfix;
    return TInfix(host(me), me.data_begin_position, me.data_end_position);
}
template <typename TNeedle, typename TSpec>
inline typename Infix<TNeedle>::Type
segment(Pattern<TNeedle, TSpec> const & me)
{
    typedef typename Infix<TNeedle>::Type TInfix;
    return TInfix(host(me), me.data_begin_position, me.data_end_position);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Pattern#host
 * @brief Query a Pattern for its host.
 *
 * @signature THost host(pattern);
 *
 * @param[in] pattern The Pattern to query for its host.
 *
 * @return THost Reference to the host.
 */

template <typename TNeedle, typename TSpec>
inline typename Host<Pattern<TNeedle, TSpec> >::Type &
host(Pattern<TNeedle, TSpec> & me)
{
SEQAN_CHECKPOINT
    return value(me.data_host);
}

template <typename TNeedle, typename TSpec>
inline typename Host<Pattern<TNeedle, TSpec> const>::Type &
host(Pattern<TNeedle, TSpec> const & me)
{
SEQAN_CHECKPOINT
    return value(me.data_host);
}


//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Pattern#needle
 * @brief Returns the needle of a @link Pattern @endlink object (not implemented for some online-algorithms).
 *
 * @signature TNeedle needle(pattern);
 *
 * @param[in] pattern The Pattern to query for its needle.
 *
 * @return TNeedle Reference of the needle object.
 *
 * TNeedle is the result of the Needle metafunction of TPattern.  This is an alias to the function @link Pattern#host @endlink.
 */


template < typename TObject >
inline typename Needle<TObject>::Type &
needle(TObject &obj)
{
    return obj;
}

template < typename TObject >
inline typename Needle<TObject const>::Type &
needle(TObject const &obj)
{
    return obj;
}


/*!
 * @fn Pattern#position
 * @brief Return the position of the last match in the pattern.
 *
 * @signature TPosition position(pattern);
 *
 * @param[in] pattern The Pattern to query for its position.
 *
 * @return TPosition The position of the last match in the pattern.
 */

template < typename TNeedle, typename TSpec >
inline typename Needle< Pattern<TNeedle, TSpec> >::Type &
needle(Pattern<TNeedle, TSpec> & obj)
{
    return host(obj);
}

template < typename TNeedle, typename TSpec >
inline typename Needle< Pattern<TNeedle, TSpec> const>::Type &
needle(Pattern<TNeedle, TSpec> const & obj)
{
    return host(obj);
}

/*!
 * @fn Pattern#setNeedle
 * @brief Sets the needle of a Pattern object and optionall induces preprocessing.
 *
 * @signature void setNeedle(pattern, needle);
 *
 * @param[in,out] pattern The pattern to set the needle for.
 * @param[in]     needle  The needle to set.
 */

template < typename TNeedle, typename TSpec >
inline void
setNeedle(Pattern<TNeedle, TSpec> &obj, TNeedle const &ndl) {
    setHost(obj, ndl);
}


//____________________________________________________________________________

/*!
 * @fn Pattern#scoringScheme
 * @brief The scoring scheme used for finding or aligning.
 *
 * @signature TScoringScheme scoringScheme(pattern);
 *
 * @param[in] pattern The Pattern to query for its scoring scheme.
 *
 * @return TScoringScheme The scoring scheme of the pattern.
 */

template <typename TNeedle, typename TSpec>
inline typename ScoringScheme<Pattern<TNeedle, TSpec> >::Type
scoringScheme(Pattern<TNeedle, TSpec> &)
{
SEQAN_CHECKPOINT
    return typename ScoringScheme<Pattern<TNeedle, TSpec> >::Type();
}
template <typename TNeedle, typename TSpec>
inline typename ScoringScheme<Pattern<TNeedle, TSpec> const>::Type
scoringScheme(Pattern<TNeedle, TSpec> const &)
{
SEQAN_CHECKPOINT
    return typename ScoringScheme<Pattern<TNeedle, TSpec> const>::Type();
}

//____________________________________________________________________________

/*!
 * @fn Pattern#setScoringScheme
 * @brief Sets the scoring scheme used for finding or aligning.
 *
 * @signature void setScoringScheme(pattern, score);
 *
 * @param[in,out] pattern The pattern to set the scoring scheme for.
 * @param[in]     score   The scoring scheme to set.
 */

template <typename TNeedle, typename TSpec, typename TScore2>
inline void
setScoringScheme(Pattern<TNeedle, TSpec> & /*me*/,
                 TScore2 & /*score*/)
{
//dummy implementation for compatibility reasons
}
//////////////////////////////////////////////////////////////////////////////

}  // namespace seqan

#endif //#ifndef SEQAN_HEADER_...
