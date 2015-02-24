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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Journal entries implementation using a sorted string of journal elements.
// ==========================================================================

#ifndef SEQAN_SEQUENCE_JOURNALED_JOURNAL_ENTRIES_SORTED_ARRAY_H_
#define SEQAN_SEQUENCE_JOURNALED_JOURNAL_ENTRIES_SORTED_ARRAY_H_

namespace seqan {

// ============================================================================
// Tags, Classes
// ============================================================================

// Tag: SortedArray.
struct SortedArray {};


template <typename TNode, typename TTreeSpec>
class JournalEntries;


template <typename TCargo_>
class JournalEntries<TCargo_, SortedArray>
{
public:
    typedef TCargo_ TCargo;
    typedef typename Size<TCargo>::Type TSize;

    String<TCargo> _journalNodes;
    TSize _originalStringLength;

    JournalEntries() : _originalStringLength(0) {}
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TCargo>
struct Iterator<JournalEntries<TCargo, SortedArray>, Standard>
{
    typedef typename Iterator<String<TCargo>, Standard>::Type Type;
};

template <typename TCargo>
struct Iterator<JournalEntries<TCargo, SortedArray> const, Standard>
{
    typedef typename Iterator<String<TCargo> const, Standard>::Type Type;
};

template <typename TCargo>
struct Value<JournalEntries<TCargo, SortedArray> >
{
    typedef TCargo Type;
};

template <typename TCargo>
struct Value<JournalEntries<TCargo, SortedArray> const>
{
    typedef TCargo const Type;
};

template <typename TCargo>
struct Reference<JournalEntries<TCargo, SortedArray> >
{
    typedef TCargo & Type;
};

template <typename TCargo>
struct Reference<JournalEntries<TCargo, SortedArray> const>
{
    typedef TCargo const & Type;
};

// ----------------------------------------------------------------------------
// Metafunction
// ----------------------------------------------------------------------------

template <typename TJournalEntries>
struct GetCargoString_{};

template <typename TCargo>
struct GetCargoString_<JournalEntries<TCargo, SortedArray> >
{
    typedef String<TCargo> Type;
};

template <typename TCargo>
struct GetCargoString_<JournalEntries<TCargo, SortedArray> const>
{
    typedef String<TCargo> const Type;
};

// ============================================================================
// Functions
// ============================================================================

template <typename TStream, typename TNode>
inline
TStream &
operator<<(TStream & stream, JournalEntries<TNode, SortedArray> const & tree)
{
    stream << "JournalEntries(";
    for (unsigned i = 0; i < length(tree._journalNodes); ++i) {
        if (i > 0) stream << ", ";
        stream << tree._journalNodes[i];
    }
    stream << ")";
    return stream;
}

template <typename TCargo>
bool _checkSortedArrayTree(JournalEntries<TCargo, SortedArray> const & tree)
{
    // std::cerr << "TREE\n" << tree._journalNodes << "\nEND OF TREE\n";
    typedef typename Iterator<String<TCargo> const, Standard>::Type TIterator;
    if (length(tree._journalNodes) == 0)
        return true;
    if (tree._journalNodes[0].virtualPosition != 0)
        return false;
    if (tree._journalNodes[0].length == 0)
        return false;

    if (tree._journalNodes[0].segmentSource == SOURCE_ORIGINAL)
    {
        if (tree._journalNodes[0].physicalPosition != tree._journalNodes[0].physicalOriginPosition)
            return false;
    }
    else
    {
        if (tree._journalNodes[0].physicalOriginPosition != 0)
            return false;
    }

    for (TIterator it = begin(tree._journalNodes, Standard()) + 1, itend = end(tree._journalNodes, Standard()); it != itend; ++it)
    {
        if (it->length == 0)
            return false;
        if ((it - 1)->virtualPosition >= it->virtualPosition)
            return false;
        if ((it - 1)->virtualPosition + (it - 1)->length != it->virtualPosition)
            return false;

        if (it->segmentSource == SOURCE_ORIGINAL)
        {
            if (it->physicalPosition != it->physicalOriginPosition)
                return false;
        }
        else
        {
            if ((it - 1)->physicalOriginPosition != it->physicalOriginPosition)
                return false;
        }
    }
    return true;
}

template <typename TCargo>
typename Iterator<JournalEntries<TCargo, SortedArray>, Standard>::Type
begin(JournalEntries<TCargo, SortedArray> & journalTree, Standard const &)
{
    SEQAN_CHECKPOINT;
    return begin(journalTree._journalNodes, Standard());
}


template <typename TCargo>
typename Iterator<JournalEntries<TCargo, SortedArray> const, Standard>::Type
begin(JournalEntries<TCargo, SortedArray> const & journalTree, Standard const &)
{
    SEQAN_CHECKPOINT;
    return begin(journalTree._journalNodes, Standard());
}


template <typename TCargo>
typename Iterator<JournalEntries<TCargo, SortedArray>, Standard>::Type
end(JournalEntries<TCargo, SortedArray> & journalTree, Standard const &)
{
    SEQAN_CHECKPOINT;
    return end(journalTree._journalNodes, Standard());
}


template <typename TCargo>
typename Iterator<JournalEntries<TCargo, SortedArray> const, Standard>::Type
end(JournalEntries<TCargo, SortedArray> const & journalTree, Standard const &)
{
    SEQAN_CHECKPOINT;
    return end(journalTree._journalNodes, Standard());
}


template <typename TCargo>
inline
void reinit(JournalEntries<TCargo, SortedArray> & tree,
            typename Size<TCargo>::Type originalStringLength)
{
    SEQAN_CHECKPOINT;
    clear(tree._journalNodes);
    appendValue(tree._journalNodes, TCargo(SOURCE_ORIGINAL, 0, 0, 0, originalStringLength));
    tree._originalStringLength = originalStringLength;
}


template <typename TCargo, typename TPos>
inline
typename Iterator<JournalEntries<TCargo, SortedArray> const, Standard>::Type
findInJournalEntries(JournalEntries<TCargo, SortedArray> const & journalEntries,
                     TPos pos)
{
    typedef typename Size<TCargo>::Type TSize;
    typedef typename Position<TCargo>::Type TPos_;
    typedef typename Iterator<String<TCargo> const, Standard>::Type TIterator;
    typedef JournalEntryLtByVirtualPos<TPos_, TSize> TCmp;

    if (pos >= back(journalEntries._journalNodes).virtualPosition)
        return end(journalEntries._journalNodes, Standard()) - 1;

    TCargo refCargo;
    refCargo.virtualPosition = pos;
    TIterator iter = std::upper_bound(begin(journalEntries._journalNodes, Standard()),
                                      end(journalEntries._journalNodes, Standard()),
                                      refCargo,
                                      TCmp());
    SEQAN_ASSERT(iter != begin(journalEntries._journalNodes, Standard()));
    --iter;

    return iter;
}

template <typename TCargo, typename TPos>
inline
typename Iterator<JournalEntries<TCargo, SortedArray>, Standard>::Type
findInJournalEntries(JournalEntries<TCargo, SortedArray> & journalEntries,
                     TPos pos)
{
    typedef typename Size<TCargo>::Type TSize;
    typedef typename Position<TCargo>::Type TPos_;
    typedef typename Iterator<String<TCargo>, Standard>::Type TIterator;
    typedef JournalEntryLtByVirtualPos<TPos_, TSize> TCmp;

    if (pos >= static_cast<TPos>(back(journalEntries._journalNodes).virtualPosition))
        return end(journalEntries._journalNodes, Standard()) - 1;

    TCargo refCargo;
    refCargo.virtualPosition = pos;
    TIterator iter = std::upper_bound(begin(journalEntries._journalNodes, Standard()),
                                      end(journalEntries._journalNodes, Standard()),
                                      refCargo,
                                      TCmp());
    SEQAN_ASSERT(iter != begin(journalEntries._journalNodes, Standard()));
    --iter;

    return iter;
}


template <typename TCargo, typename TPos>
inline
TCargo const &
findJournalEntry(JournalEntries<TCargo, SortedArray> const & journalEntries,
                 TPos pos)
{
    return *findInJournalEntries(journalEntries, pos);
}

template <typename TCargo, typename TPos>
inline
TCargo &
findJournalEntry(JournalEntries<TCargo, SortedArray> & journalEntries,
                 TPos pos)
{
    return *findInJournalEntries(journalEntries, pos);
}


template <typename TCargo>
inline
void recordInsertion(JournalEntries<TCargo, SortedArray> & tree,
                     typename Position<TCargo>::Type virtualPosition,
                     typename Position<TCargo>::Type physicalBeginPos,
                     typename Size<TCargo>::Type len)
{
    typedef typename Position<TCargo>::Type TPos;
    typedef typename Iterator<String<TCargo>, Standard>::Type TIterator;

    //std::cerr << __FILE__ << ":" << __LINE__ << " -- INSERT(" << virtualPosition << ", " << physicalBeginPos << ", " << len << ")" << std::endl;
    //std::cerr << __FILE__ << ":" << __LINE__ << " -- " << tree << std::endl;

    // Handle special case that the entry list is empty.
    if (empty(tree._journalNodes))
    {
        SEQAN_ASSERT_EQ(virtualPosition, 0u);
        if (len == 0)
            return;
        appendValue(tree._journalNodes, TCargo(SOURCE_PATCH, physicalBeginPos, virtualPosition, 0, len));
        return;
    }

    // Find position in sorted array of nodes to insert in.
    TIterator iter = findInJournalEntries(tree, virtualPosition);
    // TODO(holtgrew): Maybe move and update entries right of pos at the same time?

    // Create new journal entries.
    if (iter->virtualPosition + iter->length > virtualPosition)
    {
        TPos pos = iter - begin(tree._journalNodes, Standard());
        TPos shiftRightOf = pos;
        // Found node that contains virtualPos.
        SEQAN_ASSERT_LEQ(iter->virtualPosition, virtualPosition);
        if (iter->virtualPosition == virtualPosition)
        {
            // Simple case:  Insert left of iter.
            TPos physicalOriginPosition = 0;
            if (iter != begin(tree._journalNodes, Standard()))
                physicalOriginPosition = (iter - 1)->physicalOriginPosition;
            insertValue(tree._journalNodes, pos, TCargo(SOURCE_PATCH, physicalBeginPos, virtualPosition, physicalOriginPosition, len));
            shiftRightOf += 1;
        }
        else
        {
            // Harder case:  Split current and insert new node.
            String<TCargo, Array<3> > buffer;
            resize(buffer, 3);
            TPos offset = virtualPosition - iter->virtualPosition;
            buffer[0] = TCargo(iter->segmentSource, iter->physicalPosition, iter->virtualPosition, iter->physicalOriginPosition, offset);
            TPos physicalOriginPos1 = iter->physicalOriginPosition;
            buffer[1] = TCargo(SOURCE_PATCH, physicalBeginPos, virtualPosition, physicalOriginPos1, len);
            if (iter->segmentSource == SOURCE_ORIGINAL)
                physicalOriginPos1 += offset;
            buffer[2] = TCargo(iter->segmentSource, iter->physicalPosition + offset, virtualPosition + len, physicalOriginPos1, iter->length - offset);
            // Insert new journal entries.
            replace(tree._journalNodes, pos, pos + 1, buffer);
            shiftRightOf += 3;
        }
        // Update journal entries right of pos.
        for (TIterator it = begin(tree._journalNodes, Standard()) + shiftRightOf, itend = end(tree._journalNodes, Standard()); it != itend; ++it)
        {
            it->virtualPosition += len;
            if (it != begin(tree._journalNodes, Standard()) && it->segmentSource == SOURCE_PATCH)
                it->physicalOriginPosition = (it - 1)->physicalOriginPosition;
        }
    }
    else
    {
        // Insert at end.
        SEQAN_ASSERT_EQ(virtualPosition, iter->virtualPosition + iter->length);
        TPos physicalOriginPosition = back(tree._journalNodes).physicalOriginPosition;
        appendValue(tree._journalNodes, TCargo(SOURCE_PATCH, physicalBeginPos, virtualPosition, physicalOriginPosition, len));
    }
    //std::cerr << __FILE__ << ":" << __LINE__ << " -- " << tree << std::endl;

    SEQAN_ASSERT(_checkSortedArrayTree(tree));
}

template <typename TCargo>
inline
void recordErase(JournalEntries<TCargo, SortedArray> & tree,
                 typename Position<TCargo>::Type pos,
                 typename Position<TCargo>::Type posEnd)
{
    typedef typename Size<TCargo>::Type TSize;
    typedef typename Position<TCargo>::Type TPos;
    typedef typename Iterator<String<TCargo>, Standard>::Type TIter;
//    std::cerr << __FILE__ << ":" << __LINE__ << " -- ERASE(" << pos << ", " << posEnd << ")" << std::endl;
//    std::cerr << __FILE__ << ":" << __LINE__ << " -- " << tree << std::endl;

    // Handle special case of removing all of the singleton existing entry.
    if (length(tree._journalNodes) == 1 && pos == 0 && (TPos)front(tree._journalNodes).length == posEnd)
    {
        clear(tree._journalNodes);
        return;
    }
    // Handle case of an empty journal.
    if (empty(tree._journalNodes))
    {
        SEQAN_ASSERT_EQ(pos, 0u);
        SEQAN_ASSERT_EQ(posEnd, 0u);
        return;
    }

    // Find node.
    TIter it = findInJournalEntries(tree, pos);

    // We will shift the virtual positions of all entries right of and
    // including beginShiftPos by delta positions to the left.
    TPos delta = 0;
    TPos beginShiftPos = 0;

    TPos itPos = it - begin(tree._journalNodes, Standard());
    if (it->virtualPosition == pos && (TPos)it->length == posEnd - pos)
    {
        // Remove the whole entry.
        erase(tree._journalNodes, itPos);
        delta = posEnd - pos;
        beginShiftPos = itPos;
    }
    else if (it->virtualPosition == pos && (TPos)it->length > posEnd - pos)
    {
        // Remove a prefix of the entry.
        SEQAN_ASSERT_LT(pos, it->virtualPosition + it->length);
        delta = posEnd - pos;
        it->physicalPosition += delta;
        if (it->segmentSource == SOURCE_ORIGINAL)
            it->physicalOriginPosition += delta;
        it->length -= delta;
        beginShiftPos = itPos + 1;
    }
    else if (it->virtualPosition < pos && it->virtualPosition + it->length == posEnd)
    {
        // Remove a suffix of the entry.
        SEQAN_ASSERT_GT(pos, it->virtualPosition);
        delta = posEnd - pos;
        it->length -= delta;
        beginShiftPos = itPos + 1;
    }
    else if (it->virtualPosition < pos && it->virtualPosition + it->length > posEnd)
    {
        // Remove a true infix of the entry.
        TSize prefixLength = pos - it->virtualPosition;
        TSize suffixLength = it->length - prefixLength - (posEnd - pos);
        TSize removedInfixLength = posEnd - pos;
        // Update the left part, this must be done before the iterator is possibly invalidate because of copied memory.
        it->length -= removedInfixLength + suffixLength;
        // Insert a new entry for the right part.
        TPos physicalHostPos = it->physicalPosition;
        if (it->segmentSource == SOURCE_ORIGINAL)
            physicalHostPos += prefixLength + removedInfixLength;
        TCargo tmpEntry(it->segmentSource, it->physicalPosition + prefixLength + removedInfixLength, it->virtualPosition + prefixLength, physicalHostPos, suffixLength);
        insertValue(tree._journalNodes, itPos + 1, tmpEntry);
        // Set shift position and delta.
        delta = removedInfixLength;
        beginShiftPos = itPos + 2;
    }
    else
    {
        // Remove more than one entry.
        TPos rmBeginPos = itPos;
        TPos rmEndPos = itPos;
        if (it->virtualPosition != pos)
        {
            // Do not remove all of first.
            delta += it->length - (pos - it->virtualPosition);
            rmBeginPos += 1;
            rmEndPos += 1;
            it->length = (pos - it->virtualPosition);
        }
        else
        {
            // Remove all of first.
            delta = it->length;
            rmEndPos += 1;
        }
        it += 1;
        while (posEnd > it->virtualPosition + it->length)
        {
            rmEndPos += 1;
            delta += it->length;
            it += 1;
        }
        if (it->virtualPosition + it->length == posEnd)
        {
            // Remove all of last.
            rmEndPos += 1;
            delta += it->length;
            beginShiftPos = rmBeginPos;
        }
        else
        {
            // Do not remove all of last.
            SEQAN_ASSERT_GT(it->virtualPosition + it->length, posEnd);
            TSize tmpDelta = delta;
            delta += posEnd - it->virtualPosition;
            it->physicalPosition += posEnd - it->virtualPosition;
            if (it->segmentSource == SOURCE_ORIGINAL)
                it->physicalOriginPosition = it->physicalPosition;
            it->length -= posEnd - it->virtualPosition;
            // We update this entry manually.
            it->virtualPosition -= tmpDelta;
            beginShiftPos = rmBeginPos + 1;
        }
        erase(tree._journalNodes, rmBeginPos, rmEndPos);
    }

    // Perform left-shift of the virtual positions.
    for (TIter it = begin(tree._journalNodes, Standard()) + beginShiftPos; it != end(tree._journalNodes, Standard()); ++it) {
        SEQAN_ASSERT_GEQ(it->virtualPosition, delta);
        it->virtualPosition -= delta;
    }
    // Perform update of physical host positions.
    if (beginShiftPos <= 3)
        beginShiftPos = 1;
    else
        beginShiftPos -= 2;
    for (TIter it = begin(tree._journalNodes, Standard()) + beginShiftPos; it != end(tree._journalNodes, Standard()); ++it) {
        if (it != begin(tree._journalNodes, Standard()) && it->segmentSource == SOURCE_PATCH)
            it->physicalOriginPosition = (it - 1)->physicalOriginPosition;
    }
//    std::cerr << __FILE__ << ":" << __LINE__ << " -- " << tree << std::endl;

    SEQAN_ASSERT(_checkSortedArrayTree(tree));
}


template <typename TNode, typename TJournalSpec, typename TPos>
inline
typename Position<typename Cargo<TNode>::Type >::Type
virtualToHostPosition(JournalEntries<TNode, TJournalSpec> const & journalEntries,
                      TPos pos)
{
    typedef JournalEntries<TNode, TJournalSpec> TJournalEntries;
    typedef typename Iterator<TJournalEntries const>::Type TIterator;

    TIterator it = findInJournalEntries(journalEntries, pos);
    // The easiest case is to find a segment from the original sequence that
    // contains the position.  This case also works if we hit the entry after
    // the last one in a segment from the original sequence.
    if (value(it).segmentSource == SOURCE_ORIGINAL) {
        TPos offset = pos - value(it).virtualPosition;
        return value(it).physicalPosition + offset;
    }
    // The harder case is to find a segment from the patch sequence.  We first
    // try to find an original segment right of it, if this fails a segment
    // left of it.  If this fails, 0 is returned.
    SEQAN_ASSERT_EQ(value(it).segmentSource, SOURCE_PATCH);
    TIterator it2 = it;
    for (++it; it != end(journalEntries, Standard()); ++it) {
        if (value(it).segmentSource == SOURCE_ORIGINAL)
            return value(it).physicalPosition;
    }
    it = it2;
    // Search left;
    while (true) {
        if (value(it).segmentSource == SOURCE_ORIGINAL)
            return value(it).physicalPosition + value(it).length;
        if (it == begin(journalEntries, Standard()))
            return 0;
        --it;
    }

    SEQAN_ASSERT_FAIL("Should never reach here!");
    return 0;
}

//---------------------------------------------------
// function hostToVirtualPosition()
//---------------------------------------------------

template <typename TCargo, typename TPos>
inline
typename Position<TCargo>::Type
hostToVirtualPosition(JournalEntries<TCargo, SortedArray> const & journalEntries, TPos const & hostPos)
{
    // std::cerr << journalEntries << "\n";
    typedef JournalEntries<TCargo, SortedArray> const TJournalEntries;
    typedef typename Iterator<TJournalEntries>::Type TIterator;
    typedef typename Position<TCargo>::Type TCargoPos;
    typedef typename Size<TCargo>::Type TCargoSize;
    typedef JournalEntryLtByPhysicalOriginPos<TCargoPos, TCargoSize> TComp;

    if (empty(journalEntries._journalNodes))
        return 0;

    TCargoPos result = journalEntries._originalStringLength;

    // Find upper bound based on the physical position in host string.
    //
    // After the call to upper_bound(), it can only point to a patch node if there are no more original nodes left.
    TCargo refCargo;
    refCargo.physicalOriginPosition = hostPos;
    TIterator it = std::upper_bound(begin(journalEntries._journalNodes, Standard()),
                                    end(journalEntries._journalNodes, Standard()),
                                    refCargo,
                                    TComp());

    // As stated above, if we end up at a patch node then there are no more original nodes left.
    if (it != end(journalEntries._journalNodes, Standard()) && it->segmentSource == SOURCE_PATCH)
        return result;
    // std::cerr << *it << "\n";

    // Memoize found virtual position if any as result if we end up right of the segment we go to below.
    if (it != end(journalEntries._journalNodes, Standard()))  // INVARIANT: Is original node.
        result = it->virtualPosition;

    if (it == begin(journalEntries._journalNodes, Standard()))
        return result;

    // Go to the next original node left of the current position to get the projection position.
    --it;
    // std::cerr << *it << "\n";
    while (it != begin(journalEntries._journalNodes, Standard()) && it->segmentSource == SOURCE_PATCH)
    {
        --it;
        // std::cerr << *it << "\n";
    }
    // If we end up at the beginning, finding no original node then return result.
    if (it == begin(journalEntries._journalNodes, Standard()) && it->segmentSource == SOURCE_PATCH)
        return result;
    SEQAN_ASSERT_NEQ(it->segmentSource, SOURCE_PATCH);
    // If hostPos is not in the found node then return result.
    if ((TCargoPos)hostPos >= (TCargoPos)(it->physicalPosition + it->length))
        return result;

    // Otherwise, compute virtual position.
    return it->virtualPosition + (hostPos - it->physicalPosition);
}

template <typename TNode, typename TJournalSpec, typename TPos>
inline
bool
isGapInHost(JournalEntries<TNode, TJournalSpec> const & journalEntries,
            TPos pos)
{
    typedef JournalEntries<TNode, TJournalSpec> const TJournalEntries;
    typedef typename Iterator<TJournalEntries>::Type TIterator;

    TIterator it = findInJournalEntries(journalEntries, pos);
    return value(it).segmentSource == SOURCE_PATCH;
}

// ----------------------------------------------------------------------------
// Function clear
// ----------------------------------------------------------------------------

template <typename TNode>
inline void
clear(JournalEntries<TNode, SortedArray> & journalEntries)
{
    clear(journalEntries._journalNodes);
    journalEntries._originalStringLength = 0u;
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

template <typename TNode>
inline bool
empty(JournalEntries<TNode, SortedArray> const & journalEntries)
{
    return empty(journalEntries._journalNodes);
}

}  // namespace seqan

#endif  // SEQAN_SEQUENCE_JOURNALED_JOURNAL_ENTRIES_SORTED_ARRAY_H_
