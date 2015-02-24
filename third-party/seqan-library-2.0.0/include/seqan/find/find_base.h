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
// Definition of the class Finder and supporting tags and metafunctions.
// ==========================================================================

#ifndef SEQAN_HEADER_FIND_BASE_H
#define SEQAN_HEADER_FIND_BASE_H

// TODO(holtgrew): Should Finder implement the RootedIteratorConcept?

//////////////////////////////////////////////////////////////////////////////

namespace seqan {

//////////////////////////////////////////////////////////////////////////////
// Tags

// TODO(holtgrew): Define as Tag<FindInfix_>?

/*!
 * @defgroup ApproximateFinderSearchTypeTags Approximate Finder Search Type Tags
 * @brief Tags for specifying the @link Finder @endlink search type (prefix or infix).
 *
 * There are two interesting kinds of search for approximate string search algorithms: (1) search for a matching
 * substring anywhere in the text, infix search, and (2) search for a matching prefix, prefix search.  The tags
 * in this group can be used to select this search for approximate string search algorithms.
 */

/*!
 * @tag ApproximateFinderSearchTypeTags#FindInfix
 * @headerfile <seqan/find.h>
 *
 * @brief Find needle as a substring of the haystack (infix search).
 *
 * @signature struct FindInfix;
 *
 * @see ApproximateFinderSearchTypeTags#FindPrefix
 */

struct FindInfix;

/*!
 * @tag ApproximateFinderSearchTypeTags#FindPrefix
 * @headerfile <seqan/find.h>
 *
 * @brief Find needle as a substring of the haystack (prefix search).
 *
 * @signature struct FindPrefix;
 *
 * @see ApproximateFinderSearchTypeTags#FindInfix
 */

struct FindPrefix;

//////////////////////////////////////////////////////////////////////////////

/*!
 * @mfn DefaultFinder
 * @headerfile <seqan/find.h>
 * @brief Default @link Finder @endlink specialization type.
 *
 * @signature DefaultFinder<THaystack>::Type;
 *
 * @tparam THaystack The given haystack type.
 *
 * @return Type The Finder specialization.  By default, this is <tt>void</tt> and <tt>FinderMlr</tt> is an @link Index
 *              @endlink.
 */

template < typename TObject >
struct DefaultFinder
{
    typedef void Type;
};

/*!
 * @mfn DefaultPattern
 * @headerfile <seqan/find.h>
 * @brief Default @link Pattern @endlink specialization type.
 *
 * @signature DefaultPattern<TNeedle>::Type;
 *
 * @tparam TNeedle The given needle type.
 *
 * @return Type Is <tt>void</tt> by default
 */

template < typename TObject >
struct DefaultPattern
{
    typedef void Type;
};

/*!
 * @mfn Finder#Haystack
 * @headerfile <seqan/file.h>
 * @brief Returns the haystack type of a @link Finder @endlink type.
 *
 * @signature Haystack<TFinder>::Type;
 *
 * @tparam TFinder The finder to query.
 *
 * @return Type The haystack type of <tt>TFinder</tt>, i.e. <tt>THaystack</tt> for <tt>Finder&lt;THaystack,
 *              TSpec&gt;</tt>.  This is an alias to function <tt>host()</tt> of the pattern function.
 */

template <typename TFinder>
struct Haystack
{
    typedef typename Container<TFinder>::Type Type;
};

/*!
 * @mfn Pattern#Needle
 * @headerfile <seqan/find.h>
 * @brief Returns the needle type of a @link Pattern @endlink type.
 *
 * @signature Needle<TPattern>::Type;
 *
 * @tparam TPattern The pattern type to query.
 *
 * @return Type The needle type of <tt>TPattern</tt., i.e. <tt>TNeedle</tt> for <tt>Pattern&lt;TNeedle, TSpec&gt;</tt>.
 */

template <typename TPattern>
struct Needle
{
    typedef typename Host<TPattern>::Type Type;
};

template <typename THost, typename TSpec>
struct Needle<Segment<THost, TSpec> >
{
    typedef Segment<THost, TSpec> Type;
};

template <typename THost, typename TSpec>
struct Needle<Segment<THost, TSpec> const>
{
    typedef Segment<THost, TSpec> const Type;
};


//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Finder#find
 * @headerfile <seqan/find.h>
 * @brief Search for a @link Pattern @endlink in a @link Finder @endlink object.
 *
 * @signature bool find(finder, pattern[, k]);
 *
 * @param[in,out] finder  The @link Finder @endlink object to search through.
 * @param[in]     pattern The @link Pattern @endlink to search for.  For index finders, pattern can also be a text.
 *                        Types: @link Pattern @endlink, @link TextConcept @endlink.
 * @param[in]     k       Desired minimal score (for approximate matching).  <tt>k</tt> is a number <tt>&lt;= 0</tt>.
 *                        Differences are deletions, insertions, and substitutions.
 *
 * @return bool <tt>true</tt> if an occurence was found and <tt>false</tt> if not.
 *
 * Repeated calls of this function iterate through all occurences of <tt>pattern</tt>.
 *
 * @section Examples
 *
 * The following example shows how one can search online for a pattern in a haystack.  Note that it is neccessary to
 * reset the finder befor searching for another pattern.
 *
 * @include demos/find/finder_online.cpp
 *
 * The output is as follows.
 *
 * @include demos/find/finder_online.cpp.stdout
 *
 * In contrast to the example above the code below shows how one can use a Finder with an index as base.  Again, note
 * that it is neccessary to reset the finder befor searching for another pattern.
 *
 * @include demos/find/finder_index.cpp
 *
 * The output is as follows.
 *
 * @include demos/find/finder_index.cpp.stdout
 */


/*!
 * @class Finder
 *
 * @headerfile <seqan/find.h>
 *
 * @brief Holds the haystack and a current search context.
 *
 * @signature template <typename THaystack[, typename TSpec]>
 *            class Finder;
 *
 * @tparam TSpec The index-algorithm to search with (Optional).Leave empty for
 *               online pattern matching (see @link Pattern @endlink).If
 *               <tt>THaystack</tt> is an @link Index @endlink, then
 *               <tt>TSpec</tt> specifies the index search algorithm. Types:
 *               Pigeonhole, Swift, Backtracking Default: The result of @link
 *               DefaultFinder @endlink
 * @tparam THaystack The haystack type. Types: String, Index
 *
 * <tt>position(finder)</tt> returns the position of the current hit in the haystack.  If <tt>THaystack</tt> is a set of
 * strings or an index of a set of strings, then <tt>position(finder)</tt> returns a @link Pair @endlink <tt>(hayNo,
 * pos)</tt>, in which <tt>hayNo</tt> is the haystack index and <tt>pos</tt> the local position of the hit.
 *
 * To reset the finder object and use it on another text or different text position, use <tt>clear(finder)</tt> Note
 * that <tt>clear(finder)</tt> doesn't move the text iterator. To start the search from the beginning or somewhere else
 * in the text, use @link Finder#goBegin @endlink or @link Finder#setPosition @endlink.
 *
 * @section Examples
 *
 * The following example shows how to restart a search from the beginning of a
 * text.
 *
 * @code{.cpp}
 * CharString hstck = "I spy with my little eye something that is yellow";
 * Finder<CharString> finder(hstck);
 *
 * Pattern<CharString, Horspool> p1("y");
 * findAll(finder, p1);
 *
 * goBegin(finder);    // move Finder to the beginning of the text
 * clear(finder);      // reset Finder
 *
 * Pattern<CharString, Horspool> p2("t");
 * findAll(finder, p2);
 * @endcode
 * Demo: Demo.Index Finder StringSet
 *
 * Demo: Demo.Index Finder
 */

template < typename THaystack, typename TSpec = typename DefaultFinder<THaystack>::Type >
class Finder
{
public:
    typedef typename Iterator<THaystack, Rooted>::Type TIterator;
    typedef typename Position<THaystack>::Type TPosition;
    typedef typename Size<THaystack>::Type TSize;

    TIterator data_iterator;
    TPosition data_endPos; //note: we need this since iterator could point to begin or end (depending on pattern type)
    TSize data_length;
    bool _needReinit;                    // if true, the Pattern needs to be reinitialized
    bool _beginFind_called;                    // if false, then findBegin was not yet called for this match position (see findBegin default implementation)

/*!
 * @fn Finder::Finder
 * @brief Constructor
 *
 * @signature Finder::Finder();
 * @signature Finder::Finder(other);
 * @signature Finder::Finder(haystack);
 * @signature Finder::Finder(iter);
 *
 * @param[in] other    Other Finder of the same type (copy constructor).
 * @param[in] haystack The haystack to work on, of type <tt>THaystack</tt>.
 * @param[in] iter     The iter to work on on, either const or non-const.
 */

    Finder()
        : data_endPos(0)
        , data_length(0)
        , _needReinit(true)
        , _beginFind_called(false)
    {}

    Finder(THaystack & haystack)
        : data_iterator(begin(haystack, Rooted()))
        , data_endPos(0)
        , data_length(0)
        , _needReinit(true)
        , _beginFind_called(false)
    {}

    Finder(TIterator &iter)
        : data_iterator(iter)
        , data_endPos(0)
        , data_length(0)
        , _needReinit(true)
        , _beginFind_called(false)
    {}

    Finder(TIterator const &iter)
        : data_iterator(iter)
        , data_endPos(0)
        , data_length(0)
        , _needReinit(true)
        , _beginFind_called(false)
    {}

    Finder(Finder const &orig)
        : data_iterator(orig.data_iterator)
        , data_endPos(orig.data_endPos)
        , data_length(orig.data_length)
        , _needReinit(orig._needReinit)
        , _beginFind_called(orig._beginFind_called)
    {}

    ~Finder() {}

//____________________________________________________________________________

    Finder const &
    operator = (Finder const & other)
    {
        data_iterator = other.data_iterator;
        data_endPos = other.data_endPos;
        data_length = other.data_length;
        _needReinit = other._needReinit;
        _beginFind_called = other._beginFind_called;
        return *this;
    }

//____________________________________________________________________________

    inline typename Reference<TIterator>::Type
    operator* ()
    {
SEQAN_CHECKPOINT
        return value(hostIterator(*this));
    }

    inline typename Reference<TIterator const>::Type
    operator* () const
    {
SEQAN_CHECKPOINT
        return value(hostIterator(*this));
    }

//____________________________________________________________________________

    operator TIterator () const
    {
SEQAN_CHECKPOINT
        return data_iterator;
    }

//____________________________________________________________________________

};

//////////////////////////////////////////////////////////////////////////////

template <typename T>
inline void
_setFinderEnd(T &) {}

template <typename T, typename TPosition>
inline void
_setFinderEnd(T &, TPosition) {}


template <typename THaystack, typename TSpec>
inline void
_setFinderEnd(Finder<THaystack, TSpec> & me)
{//shortcut: move end position to iterator position +1
    me._beginFind_called = false;
    me.data_endPos = position(me)+1;
}
template <typename THaystack, typename TSpec, typename TPosition>
inline void
_setFinderEnd(Finder<THaystack, TSpec> & me,
              TPosition end_pos)
{
    me._beginFind_called = false;
    me.data_endPos = end_pos;
}

//____________________________________________________________________________

template <typename T, typename TSize>
inline void
_setFinderLength(T &, TSize) {}

template <typename THaystack, typename TSpec, typename TSize>
inline void
_setFinderLength(Finder<THaystack, TSpec> & me,
                 TSize _length)
{
    me.data_length = _length;
}

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn Finder#beginPosition
 * @brief Return begin position of match.
 *
 * @signature TPosition beginPosition(finder);
 *
 * @param[in] finder The Finder to query.
 *
 * @return TPosition The begin position of the finder.  TPosition is the position type of THaystack.
 */

template <typename THaystack, typename TSpec>
inline typename Position<THaystack>::Type
beginPosition(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
    return me.data_endPos - me.data_length;
}

template <typename THaystack, typename TSpec>
inline typename Position<THaystack const>::Type
beginPosition(Finder<THaystack, TSpec> const & me)
{
SEQAN_CHECKPOINT
    return me.data_endPos - me.data_length;
}

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn Finder#begin
 * @brief Return begin iterator of the match in the haystack.
 *
 * @signature TIter begin(finder[, tag]);
 *
 * @param[in] finder The Finder to query.
 * @param[in] tag    The tag to select the iterator type.
 *
 * @return TIter The iterator to the begin of the match in the haystack.  TIter is the same type as returned by
 *               <tt>begin(haystack[, tag])</tt> where <tt>haystack</tt> is the haystack.
 */

template <typename THaystack, typename TSpec, typename TTag>
inline typename Iterator<THaystack, Tag<TTag> const>::Type
begin(Finder<THaystack, TSpec> & me,
      Tag<TTag> const tag)
{
SEQAN_CHECKPOINT
    return iter(haystack(me), beginPosition(me), tag);
}

template <typename THaystack, typename TSpec, typename TTag>
inline typename Iterator<THaystack const, Tag<TTag> const>::Type
begin(Finder<THaystack, TSpec> const & me,
      Tag<TTag> const tag)
{
SEQAN_CHECKPOINT
    return iter(haystack(me), beginPosition(me), tag);
}

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn Finder#endPosition
 * @brief Return end position of match.
 *
 * @signature TPosition endPosition(finder);
 *
 * @param[in] finder The Finder to query.
 *
 * @return TPosition The end position of the finder.  TPosition is the position type of THaystack.
 */

template <typename THaystack, typename TSpec>
inline typename Position<THaystack>::Type
endPosition(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
    return me.data_endPos;
}

template <typename THaystack, typename TSpec>
inline typename Position<THaystack const>::Type
endPosition(Finder<THaystack, TSpec> const & me)
{
SEQAN_CHECKPOINT
    return me.data_endPos;
}

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn Finder#end
 * @brief Return end iterator of the match in the haystack.
 *
 * @signature TIter end(finder[, tag]);
 *
 * @param[in] finder The Finder to query.
 * @param[in] tag    The tag to select the iterator type.
 *
 * @return TIter The iterator to the end of the match in the haystack.  TIter is the same type as returned by
 *               <tt>end(haystack[, tag])</tt> where <tt>haystack</tt> is the haystack.
 */

template <typename THaystack, typename TSpec, typename TTag>
inline typename Iterator<THaystack, Tag<TTag> const>::Type
end(Finder<THaystack, TSpec> & me,
    Tag<TTag> const tag)
{
SEQAN_CHECKPOINT
    return iter(haystack(me), endPosition(me), tag);
}

template <typename THaystack, typename TSpec, typename TTag>
inline typename Iterator<THaystack const, Tag<TTag> const>::Type
end(Finder<THaystack, TSpec> const & me,
    Tag<TTag> const tag)
{
SEQAN_CHECKPOINT
    return iter(haystack(me), endPosition(me), tag);
}

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn Finder#length
 * @brief Return the length of the match.
 *
 * @signature TSize length(finder);
 *
 * @param[in] finder The finder to query for its match length.
 */

template <typename THaystack, typename TSpec>
inline typename Size<THaystack>::Type
length(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
    return me.data_length;
}
template <typename THaystack, typename TSpec>
inline typename Size<THaystack const>::Type
length(Finder<THaystack, TSpec> const & me)
{
SEQAN_CHECKPOINT
    return me.data_length;
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Finder#infix
 * @brief Returns the segment of the last found match in the haystack.
 *
 * @signature TInfix infix(finder);
 *
 * @param[in] finder The Finder to query.
 *
 * @return TInfix The @link SegmentableConcept#Infix @endlink of the match in the haystack.
 *
 * This function works only correct if the begin position of the match was already found, see @link Finder#findBegin @endlink
 *
 * For finders or patterns of filtering algorithms (e.g. @Spec.Swift@) the returned infix is a potential match.
 */

template <typename THaystack, typename TSpec>
inline typename Infix<THaystack>::Type
infix(Finder<THaystack, TSpec> & me)
{
    return infix(haystack(me), beginPosition(me), endPosition(me));
}

template <typename THaystack, typename TSpec>
inline typename Infix<THaystack const>::Type
infix(Finder<THaystack, TSpec> const & me)
{
    return infix(haystack(me), beginPosition(me), endPosition(me));
}

//////////////////////////////////////////////////////////////////////////////

template <typename THaystack, typename TSpec>
inline typename Parameter_<THaystack>::Type
host(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
    return container(hostIterator(me));
}

template <typename THaystack, typename TSpec>
inline typename Parameter_<THaystack>::Type
host(Finder<THaystack, TSpec> const & me)
{
SEQAN_CHECKPOINT
    return container(hostIterator(me));
}

template <typename THaystack, typename TSpec>
inline typename Parameter_<THaystack>::Type
container(Finder<THaystack, TSpec> & me)
{
    return container(hostIterator(me));
}

template <typename THaystack, typename TSpec>
inline typename Parameter_<THaystack>::Type
container(Finder<THaystack, TSpec> const & me)
{
SEQAN_CHECKPOINT
    return container(hostIterator(me));
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline void
setHost(Finder<THaystack, TSpec> & me,
        typename Parameter_<THaystack>::Type container_)
{
SEQAN_CHECKPOINT
    setContainer(hostIterator(me), container_);
    goBegin(me);
}

template <typename THaystack, typename TSpec>
inline void
setContainer(Finder<THaystack, TSpec> & me,
             typename Parameter_<THaystack>::Type container_)
{
SEQAN_CHECKPOINT
    setContainer(hostIterator(me), container_);
    goBegin(me);
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline typename Iterator<THaystack, Rooted>::Type &
hostIterator(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
    return me.data_iterator;
}

template <typename THaystack, typename TSpec>
inline typename Iterator<THaystack, Rooted>::Type const &
hostIterator(Finder<THaystack, TSpec> const & me)
{
SEQAN_CHECKPOINT
    return me.data_iterator;
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline bool
empty(Finder<THaystack, TSpec> const & me)
{
SEQAN_CHECKPOINT
    return me._needReinit;
}

/*!
 * @fn Finder#clear
 * @headerfile <seqan/find.h>
 * @brief Clear the Finder.
 *
 * @signature void clear(finder);
 *
 * @param[in,out] finder The Finder to clear.
 */

template <typename THaystack, typename TSpec>
inline void
clear(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
    me._needReinit = true;
}

//____________________________________________________________________________

template <typename T>
inline void
_finderSetNonEmpty(T & me)
{
SEQAN_CHECKPOINT
    goBegin(me);
}


template <typename THaystack, typename TSpec>
inline void
_finderSetNonEmpty(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
    me._needReinit = false;
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline bool
atBegin(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
    return (!empty(me) && atBegin(hostIterator(me)));
}

template <typename THaystack, typename TSpec>
inline bool
atEnd(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
    return (!empty(me) && atEnd(hostIterator(me)));
}

//____________________________________________________________________________

/*!
 * @fn Finder#goBegin
 * @brief Go to the beginning of the text.
 *
 * @signature void goBegin(finder);
 *
 * @param[in] finder The finder to reset to the beginning of the text.
 */

template <typename THaystack, typename TSpec>
inline void
goBegin(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
    //_finderSetNonEmpty(me);
    goBegin(hostIterator(me));
}

/*!
 * @fn Finder#goEnd
 * @brief Go to the end of the text.
 *
 * @signature void goEnd(finder);
 *
 * @param[in] finder The finder to reset to the end of the text.
 */

template <typename THaystack, typename TSpec>
inline void
goEnd(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
    //_finderSetNonEmpty(me);
    goEnd(hostIterator(me));
}

//____________________________________________________________________________

/*!
 * @fn Finder#position
 * @brief Return current position of the finder in the haystack.
 *
 * @signature TPosition position(finder);
 *
 * @param[in] finder The Finder to query.
 *
 * @return TPosition The current position.  TPosition is the position type of the haystack.
 */

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, TSpec> >::Type
position(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
    if (empty(me)) return 0;
    return position(hostIterator(me));
}

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, TSpec> >::Type
position(Finder<THaystack, TSpec> const & me)
{
SEQAN_CHECKPOINT
    if (empty(me)) return 0;
    return position(hostIterator(me));
}

//____________________________________________________________________________

/*!
 * @fn Finder#setPosition
 * @brief Sets the position of a finder.
 *
 * @signature void setPosition(finder, pos);
 *
 * @param[in,out] finder The Findre to set the position for.
 * @param[in]     pos    The position to set the finder to.
 */

template <typename THaystack, typename TSpec, typename TPosition>
inline void
setPosition(Finder<THaystack, TSpec> & me, TPosition pos_)
{
SEQAN_CHECKPOINT
    setPosition(hostIterator(me), pos_);
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline Finder<THaystack, TSpec> &
operator--(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
    --hostIterator(me);
    return me;
}

template <typename THaystack, typename TSpec>
inline Finder<THaystack, TSpec> &
operator++(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
/*            if (beforeBegin()) {
        goBegin(hostIterator(me));
    } else*/
        ++hostIterator(me);
    return me;
}

//////////////////////////////////////////////////////////////////////////////
// operator +
//////////////////////////////////////////////////////////////////////////////

template <typename THaystack, typename TSpec, typename TIntegral>
inline Finder<THaystack, TSpec> const
operator + (Finder<THaystack, TSpec> const & left, TIntegral right)
{
SEQAN_CHECKPOINT
    return Finder<THaystack, TSpec>(hostIterator(left) + right);
}

//////////////////////////////////////////////////////////////////////////////
// operator +=
//////////////////////////////////////////////////////////////////////////////

template <typename THaystack, typename TSpec, typename TIntegral>
inline Finder<THaystack, TSpec> &
operator += (Finder<THaystack, TSpec> & left,
                TIntegral right)
{
SEQAN_CHECKPOINT
    hostIterator(left) += right;
    return left;
}

//////////////////////////////////////////////////////////////////////////////
// operator -
//////////////////////////////////////////////////////////////////////////////

template <typename THaystack, typename TSpec, typename TIntegral>
inline Finder<THaystack, TSpec> const
operator - (Finder<THaystack, TSpec> const & left, TIntegral right)
{
SEQAN_CHECKPOINT
    return Finder<THaystack, TSpec>(hostIterator(left) - right);
}

template <typename THaystack, typename TSpec, typename TIntegral>
inline typename Difference<Finder<THaystack, TSpec> const>::Type
operator - (Finder<THaystack, TSpec> const & left, Finder<THaystack, TSpec> const & right)
{
SEQAN_CHECKPOINT
    return hostIterator(left) - hostIterator(right);
}

//////////////////////////////////////////////////////////////////////////////
// operator -=
//////////////////////////////////////////////////////////////////////////////

template <typename THaystack, typename TSpec, typename TIntegral>
inline Finder<THaystack, TSpec> &
operator -= (Finder<THaystack, TSpec> & left,
                TIntegral right)
{
SEQAN_CHECKPOINT
    hostIterator(left) -= right;
    return left;
}

//____________________________________________________________________________

/*!
 * @fn Finder#setHaystack
 * @brief Sets the haystack of a @link Finder @endlink object.
 *
 * @signature void setHaystack(finder, haystack);
 *
 * @param[in,out] finder   The finder to set the haystack for.
 * @param[in]     haystack The haystack to set.
 */

template < typename THaystack, typename TSpec >
inline void
setHaystack(Finder<THaystack, TSpec> &obj, THaystack const &hstk) {
    setHost(obj, hstk);
}

template < typename THaystack, typename TSpec >
inline void
setHaystack(Finder<THaystack, TSpec> &obj, THaystack &hstk)
{
    setHost(obj, hstk);
}


/*!
 * @fn Finder#haystack
 * @brief Returns the haystack of a Finder.
 *
 * @signature THaystack haystack(finder);
 *
 * @param[in] finder The Finder to query for its haystack.
 *
 * @return THaystack The result type can be retrieved using @link Finder#Haystack @endlink.
 */

template < typename TObject >
inline typename Parameter_<typename Haystack<TObject>::Type>::Type
haystack(TObject &obj) {
    return container(obj);
}

template < typename TObject >
inline typename Parameter_<typename Haystack<TObject const>::Type>::Type
haystack(TObject const &obj) {
    return container(obj);
}



//////////////////////////////////////////////////////////////////////////////

template <typename THaystack, typename TSpec>
struct Container< Finder<THaystack, TSpec> > {
    typedef THaystack Type;
};

template <typename THaystack, typename TSpec>
struct Container< Finder<THaystack, TSpec> const> {
    typedef THaystack const Type;
};


template <typename THaystack, typename TSpec>
struct Host< Finder<THaystack, TSpec> > {
    typedef THaystack Type;
};

template <typename THaystack, typename TSpec>
struct Host< Finder<THaystack, TSpec> const> {
    typedef THaystack const Type;
};


template <typename THaystack, typename TSpec>
struct Value< Finder<THaystack, TSpec> > {
    typedef typename Value<THaystack>::Type Type;
};

template <typename THaystack, typename TSpec>
struct Reference< Finder<THaystack, TSpec> >
{
    typedef typename Reference<THaystack>::Type Type;
};

template <typename THaystack, typename TSpec>
struct Reference< Finder<THaystack, TSpec> const>
{
    typedef typename Reference<THaystack const>::Type Type;
};


template <typename THaystack, typename TSpec>
struct Position< Finder<THaystack, TSpec> >:
    Position<THaystack> {};


template <typename THaystack, typename TSpec>
struct Difference< Finder<THaystack, TSpec> > {
    typedef typename Difference<THaystack>::Type Type;
};

template <typename THaystack, typename TSpec>
struct Size< Finder<THaystack, TSpec> > {
    typedef typename Size<THaystack>::Type Type;
};


template <typename THaystack, typename TSpec, typename TIteratorSpec>
struct Iterator< Finder<THaystack, TSpec>, TIteratorSpec >
{
    typedef typename Iterator<THaystack>::Type Type;
};

template <typename THaystack, typename TSpec, typename TIteratorSpec>
struct Iterator< Finder<THaystack, TSpec> const, TIteratorSpec >
{
    typedef typename Iterator<THaystack>::Type Type;
};


// TODO(holtgrew): Document DefaultGetIterator at main location, first.
// .Metafunction.DefaultGetIterator.param.T.type:Class.Finder
template <typename THaystack, typename TSpec>
struct DefaultGetIteratorSpec< Finder<THaystack, TSpec> >:
    DefaultGetIteratorSpec< THaystack >
{
};

template <typename THaystack, typename TSpec>
struct DefaultGetIteratorSpec< Finder<THaystack, TSpec> const>:
    DefaultGetIteratorSpec< THaystack const>
{
};

//////////////////////////////////////////////////////////////////////////////

}  // namespace seqan

#endif //#ifndef SEQAN_HEADER_...
