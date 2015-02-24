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

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_FRAGMENT_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_FRAGMENT_H_

namespace seqan
{

//////////////////////////////////////////////////////////////////////////////
// Fragment Specs
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class ExactFragment
 * @extends Fragment
 * @headerfile <seqan/align.h>
 * @brief A type for ungapped, pairwise segment matches.
 *
 * @signature template <[typename TSize[, typename TSpec]]>
 *            class Fragment<TSize, ExactFragment<TSpec> >;
 *
 * @tparam TSize The size type of the underlying sequence.  Default: <tt>Size&lt;CharString&gt;Type</tt>.
 * @tparam TSpec Specializing type.  Default: <tt>ExactFragment&lt;&gt;</tt>.
 */

template<typename TSpec = Default>
struct ExactFragment;


/*!
 * @class ExactReversableFragment
 * @extends Fragment
 * @headerfile <seqan/align.h>
 * @brief A type for ungapped, pairwise segment matches that maybe in reverse orientation.
 *
 * Compared to the @link ExactFragment @endlink specialzing type of @link Fragment @endlink, a @link
 * ExactReversableFragment @endlink stores an additional bool value to indicate whether a match is in reverse
 * orientation or not.
 *
 * @signature template <[typename TSize[, typename TSpec]]>
 *            class Fragment<TSize, ExactReversableFragment<TSpec> >;
 *
 * @tparam TSize The size type of the underlying sequence.  Default: <tt>Size&lt;CharString&gt;Type</tt>.
 * @tparam TSpec Specializing type.  Default: <tt>ExactFragment&lt;&gt;</tt>.
 */

template<typename TSpec = Default>
struct ExactReversableFragment;


//////////////////////////////////////////////////////////////////////////////
// Default Fragment is the exact one
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class Fragment
 * @headerfile <seqan/align.h>
 * @brief A type for pairwise segment matches.
 *
 * @signature template <[typename TSize[, typename TSpec]]>
 *            class Fragment;
 *
 * @tparam TSize The size type of the underlying sequence.  Default: <tt>Size&lt;CharString&gt;Type</tt>.
 * @tparam TSpec Specializing type.  Default: <tt>ExactFragment&lt;&gt;</tt>.
 *
 * @section Examples
 *
 * @code{.cpp}
 * // Construct fragment.
 * unsigned seqId1 = 0, beg1 = 0, seqId2 = 32, beg2 = 42, len = 33;
 * Fragment<> fragment(seqId1, beg1, seqId2, beg2, len);
 *
 * // Update fragment's properties.
 * fragmentBegin(fragment, 0) = 10;
 * fragmentBegin(fragment, 1) = 10;
 * sequenceId(fragment, 0) = 33;
 * sequenceId(fragment, 1) = 44;
 * fragmentLength(fragment) += 42;
 * @endcode
 */


template<typename TSize = typename Size<String<char> >::Type, typename TSpec = ExactFragment<> >
class Fragment;


//////////////////////////////////////////////////////////////////////////////
// Size Metafunction
//////////////////////////////////////////////////////////////////////////////

template<typename TSize, typename TSpec>
struct Size<Fragment<TSize, TSpec> > {
    typedef TSize Type;
};


template<typename TSize, typename TSpec>
struct Size<Fragment<TSize, TSpec> const> {
    typedef TSize Type;
};


//////////////////////////////////////////////////////////////////////////////
// Exact Fragment
//////////////////////////////////////////////////////////////////////////////


template<typename TSize, typename TSpec>
class Fragment<TSize, ExactFragment<TSpec> > {
public:
    typedef typename Id<Fragment>::Type TId;

    TId seqId1;
    TSize begin1;
    TId seqId2;
    TSize begin2;
    TSize len;

/*!
 * @fn ExactFragment::Fragment
 * @brief Constructor.
 *
 * @signature Fragment::Fragment();
 * @signature Fragment::Fragment(seqID1, beg1, seqID2, beg2, l);
 *
 * @param[in] seqID1  ID of the first sequence.  Type: <tt>TId</tt>.
 * @param[in] beg1    Begin position of segment match in first sequence.  Type: <tt>TSize</tt>.
 * @param[in] seqID2  ID of the second sequence.  Type: <tt>TId</tt>.
 * @param[in] beg2    Begin position of segment match in second sequence.  Type: <tt>TSize</tt>.
 * @param[in] l       The length of the segment match.  Type: <tt>TSize</tt>.
 */

    Fragment() : seqId1(0), begin1(0), seqId2(0), begin2(0), len(0) {}

    Fragment(TId sqId1, TSize beg1, TId sqId2, TSize beg2, TSize l) :
            seqId1(sqId1), begin1(beg1), seqId2(sqId2), begin2(beg2), len(l)
    {}

};

template<typename TSize, typename TSpec>
inline bool
operator==(Fragment<TSize, ExactFragment<TSpec> > const & left,
           Fragment<TSize, ExactFragment<TSpec> > const & right)
{
    return (left.seqId1 == right.seqId1 &&
            left.begin1 == right.begin1 &&
            left.seqId2 == right.seqId2 &&
            left.begin2 == right.begin2 &&
            left.len == right.len);
}

template<typename TSize, typename TSpec>
inline bool
operator<(Fragment<TSize, ExactFragment<TSpec> > const & left,
          Fragment<TSize, ExactFragment<TSpec> > const & right)
{
    if (left.seqId1 < right.seqId1)
        return true;
    if (left.seqId1 > right.seqId1)
        return false;
    if (left.begin1 < right.begin1)
        return true;
    if (left.begin1 > right.begin1)
        return false;
    if (left.seqId2 < right.seqId2)
        return true;
    if (left.seqId2 > right.seqId2)
        return false;
    if (left.begin2 < right.begin2)
        return true;
    if (left.begin2 > right.begin2)
        return false;
    if (left.len < right.len)
        return true;
    // if (left.len > right.len)
    //     return false;
    return false;
}

//////////////////////////////////////////////////////////////////////////////
// Exact Fragment that is a forward or reverse match
//////////////////////////////////////////////////////////////////////////////


template<typename TSize, typename TSpec>
class Fragment<TSize, ExactReversableFragment<TSpec> > {
public:
    typedef typename Id<Fragment>::Type TId_;

    TId_ seqId1;
    TSize begin1;
    TId_ seqId2;
    TSize begin2;
    TSize len;
    bool reversed;

/*!
 * @fn ExactReversableFragment::Fragment
 * @brief Constructor.
 *
 * @signature Fragment::Fragment();
 * @signature Fragment::Fragment(seqID1, beg1, seqID2, beg2, l[, reversed]);
 *
 * @param[in] seqID1   ID of the first sequence.  Type: <tt>TId</tt>.
 * @param[in] beg1     Begin position of segment match in first sequence.  Type: <tt>TSize</tt>.
 * @param[in] seqID2   ID of the second sequence.  Type: <tt>TId</tt>.
 * @param[in] beg2     Begin position of segment match in second sequence.  Type: <tt>TSize</tt>.
 * @param[in] l        The length of the segment match.  Type: <tt>TSize</tt>.
 * @param[in] reversed A bool; <tt>true</tt> if the segments match in reverse orientation, <tt>false</tt> otherwise.
 */


    Fragment() : seqId1(0), begin1(0), seqId2(0), begin2(0), len(0), reversed(false) {}

    Fragment(TId_ sqId1, TSize beg1, TId_ sqId2, TSize beg2, TSize l) :
            seqId1(sqId1), begin1(beg1), seqId2(sqId2), begin2(beg2), len(l), reversed(false)
    {}

    Fragment(TId_ sqId1, TSize beg1, TId_ sqId2, TSize beg2, TSize l, bool rev) :
            seqId1(sqId1), begin1(beg1), seqId2(sqId2), begin2(beg2), len(l), reversed(rev)
    {}
};

template<typename TSize, typename TSpec>
inline bool
operator==(Fragment<TSize, ExactReversableFragment<TSpec> > const & left,
           Fragment<TSize, ExactReversableFragment<TSpec> > const & right)
{
    return (left.seqId1 == right.seqId1 &&
            left.begin1 == right.begin1 &&
            left.seqId2 == right.seqId2 &&
            left.begin2 == right.begin2 &&
            left.len == right.len &&
            left.reversed == right.reversed);
}

template<typename TSize, typename TSpec>
inline bool
operator<(Fragment<TSize, ExactReversableFragment<TSpec> > const & left,
          Fragment<TSize, ExactReversableFragment<TSpec> > const & right)
{
    if (left.seqId1 < right.seqId1)
        return true;
    if (left.seqId1 > right.seqId1)
        return false;
    if (left.begin1 < right.begin1)
        return true;
    if (left.begin1 > right.begin1)
        return false;
    if (left.seqId2 < right.seqId2)
        return true;
    if (left.seqId2 > right.seqId2)
        return false;
    if (left.begin2 < right.begin2)
        return true;
    if (left.begin2 > right.begin2)
        return false;
    if (left.len < right.len)
        return true;
    if (left.len > right.len)
        return false;
    if (left.reversed < right.reversed)
        return true;
    // if (left.reversed > right.reversed)
    //     return false;
    return false;
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Fragment#label
 * @brief Access to the Fragment's label.
 *
 * @signature TInfix label(frag, stringSet, seqID);
 *
 * @param[in] frag      The Fragment to query.
 * @param[in] stringSet The @link StringSet @endlink with the sequences.
 * @param[in] seqID     The id of the sequence for which the label should be retrieved.
 */

template<typename TSize, typename TSpec, typename TStringSet, typename TVal>
inline typename Infix<typename Value<TStringSet>::Type>::Type
label(Fragment<TSize, TSpec> const& f,
      TStringSet& str,
      TVal const seqId)
{
    SEQAN_CHECKPOINT
    typedef typename Id<Fragment<TSize, TSpec> >::Type TId;
    return ((TId) seqId == (f.seqId1)) ? infix(getValueById(str, (TId) seqId), f.begin1, f.begin1 + f.len) : infix(getValueById(str, (TId) seqId), f.begin2, f.begin2 + f.len);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Fragment#sequenceId
 * @brief Access to the sequence ID of a fragment.
 *
 * @signature TId sequenceId(frag, seqNum);
 *
 * @param[in] frag   A Fragment.
 * @param[in] seqNum The sequence number fo rwhich the id should be retrieved.  Note that @link Fragment @endlink
                     stores information about exactly two sequences which can be accessed with seqNum 0 or 1 but whose
                     ids may differ from their seqNum.
 *
 * @return TId Reference to the sequence fragment id member.
 */

template<typename TSize, typename TSpec, typename TVal>
inline typename Id<Fragment<TSize, TSpec> >::Type &
sequenceId(Fragment<TSize, TSpec> const& f,
           TVal const seqId)
{
    SEQAN_CHECKPOINT
    typedef typename Id<Fragment<TSize, TSpec> >::Type TId;
    return ((TId) seqId == 0) ? const_cast<TId &>(f.seqId1) : const_cast<TId &>(f.seqId2);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Fragment#fragmentBegin
 * @brief Return fragment begin.
 *
 * @signature TSize fragmentBegin(frag, seqId);
 *
 * @param[in] frag  The Fragment to query.
 * @param[in] seqId The id of the sequence to get the begin for.
 *
 * @return TSize Reference to the fragment begin position member.
 */

template<typename TSize, typename TSpec, typename TVal>
inline TSize&
fragmentBegin(Fragment<TSize, TSpec> const& f,
              TVal const seqId)
{
    SEQAN_CHECKPOINT
    typedef typename Id<Fragment<TSize, TSpec> >::Type TId;
    return ((TId) seqId == f.seqId1) ? const_cast<TSize&>(f.begin1) : const_cast<TSize&>(f.begin2);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSize, typename TSpec, typename TVal>
inline TSize&
fragmentLength(Fragment<TSize, TSpec> const& f,
               TVal const)
{
    SEQAN_CHECKPOINT
    return const_cast<TSize&>(f.len);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Fragment#fragmentLength
 * @brief Return fragment begin.
 *
 * @signature TSize fragmentLength(frag);
 *
 * @param[in] frag The Fragment to query for its length.
 *
 * @return TSize Reference to the Fragment's length.
 */

template<typename TSize, typename TSpec>
inline TSize&
fragmentLength(Fragment<TSize, TSpec> const& f)
{
    SEQAN_CHECKPOINT
    return const_cast<TSize&>(f.len);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSize, typename TSpec, typename TId1, typename TPosition1, typename TId2, typename TPosition2>
inline void
getProjectedPosition(Fragment<TSize, ExactFragment<TSpec> > const& f,
                     TId1 const seqId,
                     TPosition1 const pos,
                     TId2& seqId2,
                     TPosition2& pos2)
{
    SEQAN_CHECKPOINT
    typedef typename Id<Fragment<TSize, TSpec> >::Type TId;

    if ((TId) seqId == f.seqId1) {
        SEQAN_ASSERT((TPosition1)f.begin1<=pos);
        SEQAN_ASSERT(pos - f.begin1 < f.len)    ;
        pos2 = f.begin2 + (pos - f.begin1);
        seqId2 = f.seqId2;
        return;
    } else {
        SEQAN_ASSERT((TPosition1)f.begin2<=pos);
        SEQAN_ASSERT(pos - f.begin2 < f.len);
        pos2 = f.begin1 + (pos - f.begin2);
        seqId2 = f.seqId1;
        return;
    }
}


//////////////////////////////////////////////////////////////////////////////

template<typename TSize, typename TSpec, typename TValue, typename TId1, typename TPosition1, typename TId2, typename TPosition2>
inline void
getProjectedPosition(Fragment<TSize, ExactFragment<TSpec> > const& f,
                    TValue seg_num,
                     TId1 const seqId,
                     TPosition1 const pos,
                     TId2& seqId2,
                     TPosition2& pos2)
{
    (void) seqId;  // When compiled without assertions.
    SEQAN_ASSERT((seg_num == 0 && seqId == f.seqId1) || (seg_num == 1 && seqId == f.seqId2));

    if (seg_num == 0) {
        SEQAN_ASSERT((TPosition1)f.begin1<=pos);
        SEQAN_ASSERT(pos - f.begin1 < f.len)    ;
        pos2 = f.begin2 + (pos - f.begin1);
        seqId2 = f.seqId2;
        return;
    } else {
        SEQAN_ASSERT((TPosition1)f.begin2<=pos);
        SEQAN_ASSERT(pos - f.begin2 < f.len);
        pos2 = f.begin1 + (pos - f.begin2);
        seqId2 = f.seqId1;
        return;
    }
}



//////////////////////////////////////////////////////////////////////////////

template<typename TSize, typename TSpec, typename TId1, typename TPosition1, typename TId2, typename TPosition2>
inline void
getProjectedPosition(Fragment<TSize, ExactReversableFragment<TSpec> > const& f,
                     TId1 const seqId,
                     TPosition1 const pos,
                     TId2& seqId2,
                     TPosition2& pos2)
{
    SEQAN_CHECKPOINT
    typedef typename Id<Fragment<TSize, TSpec> >::Type TId;

    if ((TId) seqId == f.seqId1) {
        SEQAN_ASSERT((TPosition1)f.begin1<=pos);
        SEQAN_ASSERT(pos - f.begin1 < f.len)    ;
        if (f.reversed) pos2 = (f.begin2 + f.len - 1) - (pos - f.begin1);
        else pos2 = f.begin2 + (pos - f.begin1);
        seqId2 = f.seqId2;
        return;
    } else {
        SEQAN_ASSERT((TPosition1)f.begin2<=pos);
        SEQAN_ASSERT(pos - f.begin2 < f.len);
        if (f.reversed) pos2 = (f.begin1 + f.len - 1) - (pos - f.begin2);
        else pos2 = f.begin1 + (pos - f.begin2);
        seqId2 = f.seqId1;
        return;
    }
}


/////////////////////////////////////////////////////////////

template<typename TSize, typename TSpec, typename TValue, typename TId1, typename TPosition1, typename TId2, typename TPosition2>
inline void
getProjectedPosition(Fragment<TSize, ExactReversableFragment<TSpec> > const& f,
                     TValue seg_num,
                     TId1 const seqId,
                     TPosition1 const pos,
                     TId2& seqId2,
                     TPosition2& pos2)
{
    SEQAN_CHECKPOINT
    (void) seqId;  // When compiled without assertions.
    SEQAN_ASSERT((seg_num == 0 && seqId==f.seqId1) || (seg_num == 1 && seqId==f.seqId2));

    if (seg_num == 0) {
        SEQAN_ASSERT((TPosition1)f.begin1<=pos);
        SEQAN_ASSERT(pos - f.begin1 < f.len)    ;
        if (f.reversed) pos2 = (f.begin2 + f.len - 1) - (pos - f.begin1);
        else pos2 = f.begin2 + (pos - f.begin1);
        seqId2 = f.seqId2;
        return;
    } else {
        SEQAN_ASSERT((TPosition1)f.begin2<=pos);
        SEQAN_ASSERT(pos - f.begin2 < f.len);
        if (f.reversed) pos2 = (f.begin1 + f.len - 1) - (pos - f.begin2);
        else pos2 = f.begin1 + (pos - f.begin2);
        seqId2 = f.seqId1;
        return;
    }
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn ExactReversableFragment#isReversed
 * @brief Return whether a segment match is in reverse orientation.
 *
 * @signature bool isReversed(frag);
 *
 * @param[in] frag The Fragment to query for reverseness.
 *
 * @return bool <tt>true</tt> if the fragment is reversed and <tt>false</tt> otherwise.
 */

template<typename TSize, typename TSpec>
inline bool
isReversed(Fragment<TSize, ExactReversableFragment<TSpec> > const& f)
{
    SEQAN_CHECKPOINT
    return f.reversed;
}

// Compare lexicographically as tuple.

template<typename TSize, typename TSpec>
inline bool operator>(Fragment<TSize, ExactFragment<TSpec> > const & lhs,
                      Fragment<TSize, ExactFragment<TSpec> > const & rhs)
{
    if (lhs.seqId1 > rhs.seqId1)
        return true;
    if (lhs.seqId1 == rhs.seqId1 && lhs.begin1 > rhs.begin1)
        return true;
    if (lhs.seqId1 == rhs.seqId1 && lhs.begin1 == rhs.begin1 && lhs.seqId2 > rhs.seqId2)
        return true;
    if (lhs.seqId1 == rhs.seqId1 && lhs.begin1 == rhs.begin1 && lhs.seqId2 == rhs.seqId2 && lhs.begin2 > rhs.begin2)
        return true;
    if (lhs.seqId1 == rhs.seqId1 && lhs.begin1 == rhs.begin1 && lhs.seqId2 == rhs.seqId2 && lhs.begin2 == rhs.begin2 && lhs.len > rhs.len)
        return true;
    return false;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_FRAGMENT_H_
