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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Declarations related to and implementation of the Segment class.
// ==========================================================================

#ifndef SEQAN_HEADER_SEGMENT_BASE_H
#define SEQAN_HEADER_SEGMENT_BASE_H


namespace SEQAN_NAMESPACE_MAIN
{

/*!
 * @concept SegmentableConcept
 * @headerfile <seqan/sequence.h>
 * @brief A concept for sequences that can be used as the host of a @link SegmentConcept segment @endlink.
 *
 * @signature concept Segmentable;
 * @brief Returns prefix type in a infix fashion.
 */

/*!
 * @mfn SegmentableConcept#Prefix
 * @brief Return prefix type in a flattening fashion.
 *
 * @signature Prefix<TSeq>::Type
 *
 * @tparam TSeq The segmentable sequence type to get infix type for.
 * @return Type The prefix type.
 *
 * The prefix type of a prefix is a suffix, the prefix of any other segment type is an infix.
 */

/*!
 * @fn SegmentableConcept#prefix
 * @brief Returns the prefix of a Segmentable type.
 *
 * @signature TPrefix prefix(s, endPos);
 *
 * @param s      Segmentable sequence to return the prefix for (type <tt>TSeq</tt>).
 * @param endPos End position must be convertible to <tt>Position<TSeq>::Type</tt>.
 *
 * @return TPrefix The prefix of length <tt>endPos</tt>.  Type as returned by @link Segmentable#Prefix @endlink for TSeq.
 */

/*!
 * @mfn SegmentableConcept#Infix
 * @brief Returns infix type in a flattening fashion.
 *
 * @signature Suffix<TSeq>::Type
 *
 * @tparam TSeq The segmentable sequence type to get infix type for.
 * @return Type The infix type.
 *
 * The infix any segment is an infix.
 */

/*!
 * @fn SegmentableConcept#infixWithLength
 * @brief Returns the infix of a Segmentable type.
 *
 * @signature TPrefix infixWithLength(s, beginPos, len);
 *
 * @param s        Segmentable sequence to return the infix for (type <tt>TSeq</tt>).
 * @param beginPos Begin position must be convertible to <tt>Position<TSeq>::Type</tt>.
 * @param len      Length of the prefix, must be convertible to <tt>Size<TSeq>::Type</tt>.
 *
 * Equivalent to <tt>infix(s, beginPos, beginPos + len)</tt>.
 */

/*!
 * @fn SegmentableConcept#infix
 * @brief Returns the infix of a Segmentable type.
 *
 * @signature TPrefix infix(s, beginPos, endPos);
 *
 * @param s        Segmentable sequence to return the infix for (type <tt>TSeq</tt>).
 * @param beginPos Begin position must be convertible to <tt>Position<TSeq>::Type</tt>.
 * @param endPos   End position must be convertible to <tt>Position<TSeq>::Type</tt>.
 */

/*!
 * @mfn SegmentableConcept#Suffix
 * @brief Returns suffix type in a flattening fashion.
 *
 * @signature Suffix<TSeq>::Type
 *
 * @tparam TSeq The segmentable sequence type to get suffix type for.
 * @return Type The suffix type.
 * 
 * The suffix type of a suffix is a suffix, the suffix of any other segment type is an infix.
 */

/*!
 * @fn SegmentableConcept#suffix
 * @brief Returns the suffix of a Segmentable type.
 *
 * @signature TPrefix suffix(s, beginPos);
 *
 * @param s        The segmentable type to get the suffix of.
 * @param beginPos Begin position must be convertible to <tt>Position<TSeq>::Type</tt>.
 *
 * @return TSuffix The suffix type as returned by @link Segmentable#Suffix @endlink.
 */

//////////////////////////////////////////////////////////////////////////////
// Segment
//////////////////////////////////////////////////////////////////////////////

// TODO(holtgrew): We need to document how to get reference/pointer type of segments.

/*!
 * @class Segment
 * @implements RandomAccessContainerConcept
 * @implements SegmentableConcept
 * @headerfile <seqan/sequence.h>
 * @brief A contiguous part of a sequence.
 *
 * @signature template <typename THost, typename TSpec>
 *            class Segment;
 *
 * @tparam THost The underlying @link SequenceConcept sequence@ type.
 * @tparam TSpec The tag to use for selecting the Segment specialization.
 *
 * Segments are lightweight representations of an underlying sequence (host).  Only a pointer to the host and begin
 * and/or end position have to be stored.
 *
 * Segments support element access (reading and writing) as well as random access iteration.
 *
 * @snippet core/demos/sequence/segment.cpp basic operations
 *
 * You can get the type of the infix/prefix/suffix of a sequence using @link Sequence#Infix @endlink,
 * @link Sequence#Prefix @endlink, and @link Sequence#Suffix @endlink.  These metafunctions will
 * "flatten" the type such that using these metafunctions, the infix of an infix is an infix and not
 * an Infix Segment with an Infix Segment as its host.  Instead, it will again be an Infix Segment
 * of the host of the inner type.
 *
 * A suffix of a suffix remains a suffix, a prefix of a prefix remains a prefix.  Any other combination leads to
 * the resulting type being a infix segment.
 *
 * @snippet core/demos/sequence/segment.cpp metafunction examples
 *
 * If you explicitely need a segment and keep the underlying sequence as it is, explicitely use the <tt>Segment</tt>
 * template class.
 *
 * @snippet core/demos/sequence/segment.cpp explicit segment
 */

/**
.Class.Segment:
..cat:Sequences
..summary:A contiguous part of a sequence.
..signature:Segment<THost, TSpec>
..param.THost:Type of the whole sequence.
...metafunction:Metafunction.Host
...text:Instances of $Segment<THost, TSpec>$ are subsequences of $THost$ objects.
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...default:@Spec.InfixSegment@.
..example.file:demos/sequence/infix.cpp
..example.text:The output is as follows:
..example.output:
Infix: CGCG
..include:seqan/sequence.h
*/

struct InfixSegment {};

template <typename THost, typename TSpec = InfixSegment>
class Segment
{
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Host.param.T.type:Class.Segment
///.Metafunction.Host.class:Class.Segment

template <typename THost, typename TSpec>
struct Host<Segment<THost, TSpec> >
{
    typedef THost Type;
};

template <typename THost, typename TSpec>
struct Host<Segment<THost, TSpec> const >
{
    typedef THost Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec>
struct Pointer_<Segment<THost, TSpec> >
{
    typedef Segment<THost, TSpec> Type;
};

template <typename THost, typename TSpec>
struct Pointer_<Segment<THost, TSpec> const >
{
    typedef Segment<THost, TSpec> Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec>
struct Parameter_<Segment<THost, TSpec> >
{
    typedef Segment<THost, TSpec> Type;
};

template <typename THost, typename TSpec>
struct Parameter_<Segment<THost, TSpec> const >
{
    typedef Segment<THost, TSpec> Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Spec.param.T.type:Class.Segment
///.Metafunction.Spec.class:Class.Segment

template <typename THost, typename TSpec>
struct Spec<Segment<THost, TSpec> >
{
    typedef TSpec Type;
};
template <typename THost, typename TSpec>
struct Spec<Segment<THost, TSpec> const>
{
    typedef TSpec Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Class.Segment
///.Metafunction.Value.class:Class.Segment

template <typename THost, typename TSpec>
struct Value<Segment<THost, TSpec> > :
    Value<THost> {};

template <typename THost, typename TSpec>
struct Value<Segment<THost, TSpec> const > :
    Value<THost const> {};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.GetValue.param.T.type:Class.Segment
///.Metafunction.GetValue.class:Class.Segment

template <typename THost, typename TSpec>
struct GetValue<Segment<THost, TSpec> > :
    GetValue<THost> {};

template <typename THost, typename TSpec>
struct GetValue<Segment<THost, TSpec> const > :
    GetValue<THost const> {};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Reference.param.T.type:Class.Segment
///.Metafunction.Reference.class:Class.Segment

template <typename THost, typename TSpec>
struct Reference<Segment<THost, TSpec> > :
    Reference<THost> {};

template <typename THost, typename TSpec>
struct Reference<Segment<THost, TSpec> const > :
    Reference<THost> {};

// TODO(weese): It's a philosophical question whether const Views like Segments are allowed to modify the host or not
//              It would be more restrictive but consistent in generic algorithms if Segments would behave like Strings (immutable host)
//              On the other hand segments are prone to unintentionally be given as a const to a function (by-const-reference, not by copy)

//template <typename THost, typename TSpec>
//struct Reference<Segment<THost, TSpec> const > :
//    Reference<THost const> {};

//////////////////////////////////////////////////////////////////////////////

// TODO(holtgrew): Should the iterators of const segments be iterators with the constness of the host.

///.Metafunction.Iterator.param.T.type:Class.Segment
///.Metafunction.Iterator.class:Class.Segment

template <typename THost, typename TSpec>
struct Iterator<Segment<THost, TSpec>, Rooted>
{
    typedef Segment<THost, TSpec> TSequence_;
    typedef typename Iterator<THost, Standard>::Type TIterator_;
    typedef Iter<TSequence_, AdaptorIterator<TIterator_> > Type;
};
template <typename THost, typename TSpec>
struct Iterator<Segment<THost, TSpec> const, Rooted>
{
    typedef Segment<THost, TSpec> TSequence_;
    typedef typename Iterator<THost, Standard>::Type TIterator_;
    typedef Iter<TSequence_, AdaptorIterator<TIterator_> > Type;
};

template <typename THost, typename TSpec>
struct Iterator<Segment<THost, TSpec>, Standard>:
    Iterator<THost, Standard> {};

template <typename THost, typename TSpec>
struct Iterator<Segment<THost, TSpec> const, Standard>:
    Iterator<THost , Standard> {};



//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Size.param.T.type:Class.Segment
///.Metafunction.Size.class:Class.Segment

template <typename THost, typename TSpec>
struct Size<Segment<THost, TSpec> > :
    Size<THost> {};

template <typename THost, typename TSpec>
struct Size<Segment<THost, TSpec> const > :
    Size<THost> {};

template <typename THost, typename TSpec>
struct Position<Segment<THost, TSpec> > :
    Position<THost> {};

template <typename THost, typename TSpec>
struct Position<Segment<THost, TSpec> const > :
    Position<THost> {};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.DefaultOverflowImplicit.param.T.type:Class.Segment
///.Metafunction.DefaultOverflowImplicit.class:Class.Segment

template <typename THost, typename TSpec>
struct DefaultOverflowImplicit<Segment<THost, TSpec > > :
    DefaultOverflowImplicit<THost> {};

template <typename THost, typename TSpec>
struct DefaultOverflowImplicit<Segment<THost, TSpec > const > :
    DefaultOverflowImplicit<THost> {};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.DefaultOverflowExplicit.param.T.type:Class.Segment
///.Metafunction.DefaultOverflowExplicit.class:Class.Segment

template <typename THost, typename TSpec>
struct DefaultOverflowExplicit<Segment<THost, TSpec > > :
    DefaultOverflowExplicit<THost> {};

template <typename THost, typename TSpec>
struct DefaultOverflowExplicit<Segment<THost, TSpec > const > :
    DefaultOverflowExplicit<THost> {};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.IsContiguous.param.T.type:Class.Segment
///.Metafunction.IsContiguous.class:Class.Segment

template <typename THost, typename TSpec>
struct IsContiguous< Segment<THost, TSpec> > :
    IsContiguous<THost> {};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.IsSequence.param.T.type:Class.Segment
///.Metafunction.IsSequence.class:Class.Segment

template <typename THost, typename TSpec>
struct IsSequence< Segment<THost, TSpec> > :
    True {};

//////////////////////////////////////////////////////////////////////////////

///.Function.atBegin.param.iterator.type:Class.Segment
///.Function.atBegin.class:Class.Segment
///.Function.atEnd.param.iterator.type:Class.Segment
///.Function.atEnd.class:Class.Segment
///.Function.goBegin.param.iterator.type:Class.Segment
///.Function.goBegin.class:Class.Segment
///.Function.goEnd.param.iterator.type:Class.Segment
///.Function.goEnd.class:Class.Segment
///.Function.goNext.param.iterator.type:Class.Segment
///.Function.goNext.class:Class.Segment
///.Function.goPrevious.param.iterator.type:Class.Segment
///.Function.goPrevious.class:Class.Segment
///.Function.value.param.container.type:Class.Segment
///.Function.value.class:Class.Segment

///.Function.shareResources.param.sequence1, sequence2.type:Class.Segment
///.Function.shareResources.class:Class.Segment

//////////////////////////////////////////////////////////////////////////////
// functions for all Segment classes
//////////////////////////////////////////////////////////////////////////////

///.Function.getObjectId.param.object.type:Class.Segment
///.Function.getObjectId.class:Class.Segment

template <typename THost, typename TSpec>
inline void const *
getObjectId(Segment<THost, TSpec> const & me)
{
SEQAN_CHECKPOINT
    return getObjectId(host(me));
}

// --------------------------------------------------------------------------
// Function _toPointer()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
inline typename Pointer_<Segment<THost, TSpec> >::Type
_toPointer(Segment<THost, TSpec> & me)
{
    return me;
}

template <typename THost, typename TSpec>
inline typename Pointer_<Segment<THost, TSpec> const >::Type
_toPointer(Segment<THost, TSpec> const & me)
{
    return me;
}

// --------------------------------------------------------------------------
// Function _fromPointer()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
Segment<THost, TSpec> _fromPointer(Segment<THost, TSpec> & me)
{
    return me;
}

template <typename THost, typename TSpec>
Segment<THost, TSpec> _fromPointer(Segment<THost, TSpec> const & me)
{
    return me;
}


//////////////////////////////////////////////////////////////////////////////

///.Function.empty.param.object.type:Class.Segment

template <typename THost, typename TSpec>
inline bool
empty(Segment<THost, TSpec> const & me)
{
    SEQAN_CHECKPOINT;
    return (beginPosition(me) == endPosition(me));
}

//////////////////////////////////////////////////////////////////////////////

///.Function.length.param.object.type:Class.Segment
///.Function.length.class:Class.Segment

template <typename THost, typename TSpec>
inline typename Size<Segment<THost, TSpec> const>::Type
length(Segment<THost, TSpec> const & me)
{
SEQAN_CHECKPOINT
    return endPosition(me) - beginPosition(me);
}

//////////////////////////////////////////////////////////////////////////////

///.Function.capacity.param.object.type:Class.Segment
///.Function.capacity.class:Class.Segment

template <typename THost, typename TSpec>
inline typename Size< Segment<THost, TSpec> const>::Type
capacity(Segment<THost, TSpec> const & me)
{
SEQAN_CHECKPOINT
    return capacity(host(me)) + length(me) - length(host(me));
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec>
inline bool
hasNoHost(Segment<THost, TSpec> const & target)
{
    return !_toPointer(host(target));
}


//////////////////////////////////////////////////////////////////////////////

// operationToSet testet, ob statt einer assign/append-Funktion
// eine set Funktion aufgerufen werden soll. Das ist nur dann der Fall,
// wenn target keinen Host hat und source einen kompatiblen Host
// anbietet.
// returns true:  set Funktion verwendet,
//                wurde bereits von test_operation_2_set gemacht
// returns false: keine set Funktion verwenden
//                muss noch assign/append gemacht werden.

template <typename TSameSpec, typename TTargetInfix>
struct SegmentSetImpl_
{
    template <typename TTarget, typename TSource>
    static inline bool operationToSet(TTarget &, TSource &)
    {
        return false;
    }
};

template <typename TTargetInfix>
struct SegmentSetImpl_<True, TTargetInfix>
{
    template <typename TTarget, typename TSource>
    static inline bool operationToSet(TTarget & target, TSource & source)
    {
        set(target, source);
        return true;
    }
};

template <>
struct SegmentSetImpl_<False, True>
{
    template <typename TTarget, typename TSource>
    static inline bool operationToSet(TTarget & target, TSource & source)
    {
        set(target, host(source), beginPosition(source), endPosition(source));
        return true;
    }
};

//////////////////////////////////////////////////////////////////////////////
// assign
//////////////////////////////////////////////////////////////////////////////

/**
.Function.assign:
..class:Class.Segment
..remarks:If $target$ is a @Class.Segment@ object, then
$limit$ denotes the maximal length of @Function.host.$host(target)$@ after the operation.
..param.target.type:Class.Segment
..param.source.type:Class.Segment
*/

// TODO(holtgrew): We'd rather only have one version.

template<typename THost, typename TSpec>
inline void
assign(Segment<THost, TSpec> & target,
       Segment<THost, TSpec> const & source)
{
    typedef Segment<THost, TSpec> TSegment;
    assign(target, source, typename DefaultOverflowImplicit<TSegment>::Type());
}

template<typename THost, typename TSpec>
inline void
assign(Segment<THost, TSpec> & target,
       Segment<THost, TSpec> & source)
{
    typedef Segment<THost, TSpec> TSegment;
    assign(target, source, typename DefaultOverflowImplicit<TSegment>::Type());
}

//////////////////////////////////////////////////////////////////////////////

template <typename TExpand>
struct AssignSegment_
{
    template <typename THost, typename TSpec, typename TSource>
    static inline void
    assign_(
        Segment<THost const, TSpec> & target,
        TSource & source)
    {
        set(target, source);
    }

    template <typename THost, typename TSpec, typename TSource>
    static inline void
    assign_(
        Segment<THost const, TSpec> & target,
        TSource & source,
        typename Size< Segment<THost const, TSpec> >::Type /*limit*/)
    {
        set(target, source);
    }

    template <typename THost, typename TSpec, typename TSource>
    static inline void
    assign_(
        Segment<THost, TSpec> & target,
        TSource & source)
    {
        set(target, source);
    }

    template <typename THost, typename TSpec, typename TSource>
    static inline void
    assign_(
        Segment<THost, TSpec> & target,
        TSource & source,
        typename Size< Segment<THost, TSpec> >::Type limit)
    {
        (void)limit;
        set(target, source);
    }

    template <typename THost, typename TSpec, typename TSource>
    static inline void
    assign_(
        Segment<THost, TSpec> const & target,
        TSource & source)
    {
        set(target, source);
    }

    template <typename THost, typename TSpec, typename TSource>
    static inline void
    assign_(
        Segment<THost, TSpec> const & target,
        TSource & source,
        typename Size< Segment<THost, TSpec> >::Type limit)
    {
        (void)limit;
        set(target, source);
    }
};

//____________________________________________________________________________

template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void
assign(Segment<THost, TSpec> & target,
       TSource & source,
       Tag<TExpand>)
{
SEQAN_CHECKPOINT
    AssignSegment_<Tag<TExpand> >::assign_(target, source);
}
template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void
assign(Segment<THost, TSpec> & target,
       TSource const & source,
       Tag<TExpand>)
{
SEQAN_CHECKPOINT
    AssignSegment_<Tag<TExpand> >::assign_(target, source);
}

template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void
assign(Segment<THost, TSpec> & target,
       TSource & source,
       typename Size< Segment<THost, TSpec> >::Type limit,
       Tag<TExpand>)
{
SEQAN_CHECKPOINT
    AssignSegment_<Tag<TExpand> >::assign_(target, source, limit);
}
template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void
assign(Segment<THost, TSpec> & target,
       TSource const & source,
       typename Size< Segment<THost, TSpec> >::Type limit,
       Tag<TExpand>)
{
SEQAN_CHECKPOINT
    AssignSegment_<Tag<TExpand> >::assign_(target, source, limit);
}

//(for temporary targets)

template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void
assign(Segment<THost, TSpec> const & target,
       TSource & source,
       Tag<TExpand>)
{
SEQAN_CHECKPOINT
    AssignSegment_<Tag<TExpand> >::assign_(target, source);
}
template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void
assign(Segment<THost, TSpec> const & target,
       TSource const & source,
       Tag<TExpand>)
{
SEQAN_CHECKPOINT
    AssignSegment_<Tag<TExpand> >::assign_(target, source);
}

template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void
assign(Segment<THost, TSpec> const & target,
       TSource & source,
       typename Size< Segment<THost, TSpec> >::Type limit,
       Tag<TExpand>)
{
SEQAN_CHECKPOINT
    AssignSegment_<Tag<TExpand> >::assign_(target, source, limit);
}
template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void
assign(Segment<THost, TSpec> const & target,
       TSource const & source,
       typename Size< Segment<THost, TSpec> >::Type limit,
       Tag<TExpand>)
{
SEQAN_CHECKPOINT
    AssignSegment_<Tag<TExpand> >::assign_(target, source, limit);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename TLeftHost, typename TLeftSpec, typename TRight >
inline bool
operator == (Segment<TLeftHost, TLeftSpec> const & left,
            TRight const & right)
{
SEQAN_CHECKPOINT
    typename Comparator<Segment<TLeftHost, TLeftSpec> >::Type _lex(left, right);
    return isEqual(_lex);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TLeftHost, typename TLeftSpec, typename TRight >
inline bool
operator != (Segment<TLeftHost, TLeftSpec> const & left,
            TRight const & right)
{
SEQAN_CHECKPOINT
    typename Comparator<Segment<TLeftHost, TLeftSpec> >::Type _lex(left, right);
    return isNotEqual(_lex);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TLeftHost, typename TLeftSpec, typename TRight>
inline bool
operator < (Segment<TLeftHost, TLeftSpec> const & left,
            TRight const & right)
{
SEQAN_CHECKPOINT
    return isLess(left, right, typename DefaultPrefixOrder<Segment<TLeftHost, TLeftSpec> >::Type());
}

//////////////////////////////////////////////////////////////////////////////

template <typename TLeftHost, typename TLeftSpec, typename TRight>
inline bool
operator <= (Segment<TLeftHost, TLeftSpec> const & left,
             TRight const & right)
{
SEQAN_CHECKPOINT
    return isLessOrEqual(left, right, typename DefaultPrefixOrder<Segment<TLeftHost, TLeftSpec> >::Type());
}
//////////////////////////////////////////////////////////////////////////////

template <typename TLeftHost, typename TLeftSpec, typename TRight>
inline bool
operator > (Segment<TLeftHost, TLeftSpec> const & left,
            TRight const & right)
{
SEQAN_CHECKPOINT
    return isGreater(left, right, typename DefaultPrefixOrder<Segment<TLeftHost, TLeftSpec> >::Type());
}

//////////////////////////////////////////////////////////////////////////////

template <typename TLeftHost, typename TLeftSpec, typename TRight>
inline bool
operator >= (Segment<TLeftHost, TLeftSpec> const & left,
        TRight const & right)
{
SEQAN_CHECKPOINT
    return isGreaterOrEqual(left, right, typename DefaultPrefixOrder<Segment<TLeftHost, TLeftSpec> >::Type());
}


//////////////////////////////////////////////////////////////////////////////
// stream operators
//////////////////////////////////////////////////////////////////////////////

template <typename TStream, typename THost, typename TSpec>
inline TStream &
operator << (TStream & target,
             Segment<THost, TSpec> const & source)
{
SEQAN_CHECKPOINT
    write(target, source);
    return target;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TStream, typename THost, typename TSpec>
inline TStream &
operator >> (TStream & source,
             Segment<THost, TSpec> & target)
{
SEQAN_CHECKPOINT
    read(source, target);
    return source;
}
template <typename TStream, typename THost, typename TSpec>
inline TStream &
operator >> (TStream & source,
             Segment<THost, TSpec> const & target)
{
SEQAN_CHECKPOINT
    read(source, target);
    return source;
}

//////////////////////////////////////////////////////////////////////////////

// this function doesn't do anything as we are not allowed to change the host (only its elements)
// it is, however, implemented for algorithms that get a sequence to work on
// and need to make sure that it has a certain length

template <typename THost, typename TSpec, typename TSize, typename TExpand>
inline typename Size< Segment<THost, TSpec> >::Type
resize(
    Segment<THost, TSpec> & me,
    TSize new_length,
    Tag<TExpand>)
{
    ignoreUnusedVariableWarning(new_length);

    SEQAN_ASSERT_EQ(new_length, length(me));
    return length(me);
}

//////////////////////////////////////////////////////////////////////////////

// TODO(singer): moveValue still works. Should make the compiler throw an error.

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
