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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#ifndef INCLUDE_SEQAN_JOURNALED_SET_JOURNALED_SET_IMPL_H_
#define INCLUDE_SEQAN_JOURNALED_SET_JOURNALED_SET_IMPL_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class StringSet                                        [Owner<JournaledSet>]
// ----------------------------------------------------------------------------

template <typename TString>
class StringSet<TString, Owner<JournaledSet> >
{
public:
    typedef typename Position<StringSet>::Type TPosition_;
    typedef typename Host<TString>::Type THostString_;
    typedef String<TString>  TStrings_;
    typedef typename StringSetLimits<StringSet>::Type TLimits_;

    Holder<THostString_> _globalRefHolder;
    TStrings_   strings;
    TLimits_    limits;
    bool        limitsValid;

    StringSet() : limitsValid(true)
    {
        appendValue(limits,0);
    }

    ~StringSet() {}

    template <typename TPosition>
    inline typename Reference<StringSet>::Type
    operator[](TPosition const pos)
    {
        return value(*this, pos);
    }

    template <typename TPosition>
    inline typename Reference<StringSet const>::Type
    operator[](TPosition const pos) const
    {
        return value(*this, pos);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Host
// ----------------------------------------------------------------------------

template <typename TString>
struct Host<StringSet<TString, Owner<JournaledSet> > >
{
    typedef typename Host<TString>::Type Type;
};

template <typename TString>
struct Host<StringSet<TString, Owner<JournaledSet> > const >
{
    typedef typename Host<TString>::Type const Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TString, typename TPosition>
inline typename Reference<StringSet<TString, Owner<JournaledSet> > >::Type
value(StringSet<TString, Owner<JournaledSet> > & me,
      TPosition const pos)
{
    return me.strings[pos];
}

template <typename TString, typename TPosition>
inline typename Reference<StringSet<TString, Owner<JournaledSet> > const >::Type
value(StringSet<TString, Owner<JournaledSet> > const & me,
      TPosition const pos)
{
    return me.strings[pos];
}

// ----------------------------------------------------------------------------
// Function appendValue()
// ----------------------------------------------------------------------------

template <typename TString, typename TString2, typename TExpand>
void
appendValue(StringSet<TString, Owner<JournaledSet> > & journalSet,
            TString2 const & newElement,
            Tag<TExpand> tag)
{
    resize(journalSet, length(journalSet) + 1, tag);
    assignValue(journalSet, length(journalSet) - 1 , newElement);
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TString>
inline void
clear(StringSet<TString, Owner<JournaledSet> > & journalSet)
{
    clear(journalSet._globalRefHolder);
    clear(journalSet.strings);
    resize(journalSet.limits, 1, Exact());
    journalSet.limitsValid = true;
}

// ----------------------------------------------------------------------------
// Function erase()
// ----------------------------------------------------------------------------

template <typename TString, typename TPos>
inline typename Size<StringSet<TString, Owner<JournaledSet> > >::Type
erase(StringSet<TString, Owner<JournaledSet> > & journalSet, TPos pos)
{
    erase(journalSet.strings, pos);
    journalSet.limitsValid = false;
    return length(journalSet);
}

template <typename TString, typename TPos, typename TPosEnd>
inline typename Size<StringSet<TString, Owner<JournaledSet> > >::Type
erase(StringSet<TString, Owner<JournaledSet> > & journalSet, TPos pos, TPosEnd posEnd)
{
    erase(journalSet.strings, pos, posEnd);
    journalSet.limitsValid = false;
    return length(journalSet);
}

// ----------------------------------------------------------------------------
// Function resize()
// ----------------------------------------------------------------------------

template <typename TString, typename TSize, typename TValue, typename TExpandTag>
inline typename Size<StringSet<TString, Owner<JournaledSet> > >::Type
resize(StringSet<TString, Owner<JournaledSet> > & journalSet,
       TSize const & newSize,
       TValue const & fillValue,
       Tag<TExpandTag> const & expansionTag)
{
    resize(journalSet.strings, newSize, fillValue, expansionTag);
    resize(journalSet.limits, newSize + 1);
    journalSet.limitsValid = true;
    return length(journalSet);
}

// --------------------------------------------------------------------------
// Function assignValue()
// --------------------------------------------------------------------------

// No journaled strings as values.
template <typename TString, typename TPos,  typename TString2>
inline void assignValue(
    StringSet<TString, Owner<JournaledSet> > & journalSet,
    TPos pos,
    TString2 const & newElement)
{
    SEQAN_ASSERT_GEQ(pos, static_cast<TPos>(0));
    SEQAN_ASSERT_LT(pos, static_cast<TPos>(length(journalSet)));

    assign(journalSet[pos], newElement);
}
// ----------------------------------------------------------------------------
// Function host()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledSet#host
 * @brief Returns the global reference sequence of a @link JournaledSet @endlink.
 *
 * @signature THost host(stringSet);
 *
 * @param[in] stringSet The JournaledStringSet that stores the sequences. Types: @link  JournaledSet  @endlink
 *
 * @return THost Reference to the host.
 */

template <typename TString>
inline typename Host<StringSet<TString, Owner<JournaledSet> > >::Type const &
host(StringSet<TString, Owner<JournaledSet> > const & journalSet)
{
    return value(journalSet._globalRefHolder);
}


template <typename TString>
inline typename Host<StringSet<TString, Owner<JournaledSet> > >::Type &
host(StringSet<TString, Owner<JournaledSet> > & journalSet)
{
    return value(journalSet._globalRefHolder);
}

// ----------------------------------------------------------------------------
// Function setHost()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledSet#setHost
 * @brief Sets the global reference of a @link JournaledSet @endlink.
 *
 * @signature void setHost(stringSet, ref);
 *
 * @param[in,out] stringSet The string set that stores the sequences. Types: @link JournaledSet @endlink
 * @param[in]     ref       The new reference sequence of the @link JournaledSet  @endlink.
 *
 * @section Remarks
 *
 * Uses an @link Holder @endlink to store a reference to the new global reference sequence instead of copying it.
 */

template <typename TString, typename THost>
inline void
setHost(StringSet<TString, Owner<JournaledSet> > & journalSet,
        THost & newGlobalRef)
{
    setValue(journalSet._globalRefHolder, newGlobalRef);
}

// ----------------------------------------------------------------------------
// Function createHost()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledSet#createHost
 * @brief Creates the global reference of a @link JournaledSet @endlink.
 *
 * @signature void createHost(stringSet, ref);
 *
 * @param[in,out] stringSet The JournaledStringSet that stores the sequences.
 * @param[in]     ref       The new reference sequence of the JournaledSet.  Stores a copy of the passed global
 *                          reference sequence.
 */

template <typename TString>
inline void
createHost(StringSet<TString, Owner<JournaledSet> > & journalSet,
                   typename Host<TString>::Type const & newGlobalRef)
{
    create(journalSet._globalRefHolder, newGlobalRef);
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_SET_JOURNALED_SET_IMPL_H_
