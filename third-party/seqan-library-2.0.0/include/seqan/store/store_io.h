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

#ifndef SEQAN_HEADER_STORE_IO_H
#define SEQAN_HEADER_STORE_IO_H

/* IOREV
 *
 * _doc_
 *
 *
 * if this file is about the amos file format why isn't it named accordingly?
 *
 * altogether it is unclear why sequence io is in file/ but store io is in
 * store/
 *
 */



namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// File tags
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

struct TagAmos_;
typedef Tag<TagAmos_> const Amos;


//////////////////////////////////////////////////////////////////////////////
// Auxillary functions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn FragmentStore#getClrRange
 * @brief Get the "clear" range of a read alignment.
 *
 * The clear range of a read alignment is the range of the part of the alignmetn that is not clipped.
 *
 * @signature void getClrRange(store, alignEl, begClr, endClr);
 *
 * @param[in]  store    The FragmentStore to work on.
 * @param[in]  alignEl  The @link AlignedReadStoreElement @endlink to work on.
 * @param[out] begClr   Begin of the clear range.
 * @param[out] endClr   End of the clear range.
 */


template <typename TSpec, typename TConfig, typename TPos, typename TGapAnchor, typename TSpecAlign, typename TBeginClr, typename TEndClr>
inline void
getClrRange(FragmentStore<TSpec, TConfig> const& fragStore,
            AlignedReadStoreElement<TPos, TGapAnchor, TSpecAlign> const& alignEl,
            TBeginClr& begClr,        // Out-parameter: left / begin position of the clear range
            TEndClr& endClr)        // Out-parameter: right / end position of the clear range
{
    typedef FragmentStore<TSpec, TConfig> TFragmentStore;
    typedef typename Size<TFragmentStore>::Type TSize;
    typedef typename Iterator<String<TGapAnchor> const, Standard>::Type TGapIter;

    TSize lenRead = length(fragStore.readSeqStore[alignEl.readId]);
    TGapIter itGap = begin(alignEl.gaps, Standard());
    TGapIter itGapEnd = end(alignEl.gaps, Standard());

    // Any gaps or clipped characters?
    if (itGap == itGapEnd) {
        begClr = 0;
        endClr = lenRead;
    } else {
        // Begin clear range
        begClr = (itGap->gapPos == 0) ? itGap->seqPos : 0;
        // End clear range
        --itGapEnd;
        if (static_cast<TSize>(itGapEnd->seqPos) != lenRead) endClr = lenRead;
        else {
            int diff = (itGap != itGapEnd) ? (*(itGapEnd - 1)).gapPos - (*(itGapEnd-1)).seqPos : 0;
            int newDiff = itGapEnd->gapPos - itGapEnd->seqPos;
            endClr = (newDiff < diff) ? lenRead - (diff - newDiff) : lenRead;
        }
    }

    // For reverse reads adapt clear ranges
    if (alignEl.beginPos > alignEl.endPos) {
        TBeginClr tmp = begClr;
        begClr = lenRead - endClr;
        endClr = lenRead - tmp;
    }
}



//////////////////////////////////////////////////////////////////////////////
// Read / Write of AMOS message files (*.afg)
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn FragmentStore#read
 * @brief Read the contents of a FragmentStore from a file.
 *
 * @signature int read(file, store, tag);
 *
 * @param[in,out] file  The @link StreamConcept @endlink to read from.
 * @param[in,out] store The FragmentStore to append to.
 * @param[in]     tag   The format to read from.  Can be Amos or Sam.
 *
 * @return int 0 in the case of success, non-0 value in case of errors.
 */

template <typename TKey, typename TValue, typename TIter>
void readAmosKeyValue(TKey &key, TValue &value, TIter &reader)
{
    clear(key);
    clear(value);

    // read "KEY:VALUE" pair
    readUntil(key, reader, EqualsChar<':'>());
    skipOne(reader, EqualsChar<':'>());
    readUntil(value, reader, IsWhitespace());
}

template<typename TFile, typename TSpec, typename TConfig>
inline int
read(FragmentStore<TSpec, TConfig>& fragStore,
     TFile & file,
     Amos)
{
    // Basic types
    typedef FragmentStore<TSpec, TConfig> TFragmentStore;
    typedef typename Id<TFragmentStore>::Type TId;
    typedef typename Size<TFragmentStore>::Type TSize;
    //typedef typename Value<TFile>::Type TValue;

    // All fragment store element types
    typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigElement;
    typedef typename Value<typename TFragmentStore::TLibraryStore>::Type TLibraryStoreElement;
    typedef typename Value<typename TFragmentStore::TMatePairStore>::Type TMatePairElement;
    typedef typename Value<typename TFragmentStore::TReadStore>::Type TReadStoreElement;
    typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;
    typedef typename TFragmentStore::TReadSeq TReadSeq;
    typedef typename TFragmentStore::TContigSeq TContigSeq;

    // All maps to mirror file ids to our ids
    typedef std::map<TId, TSize> TIdMap;
    // The following maps the library/fragment/read id from AMOS into the fragment store's ids.
    TIdMap libIdMap;
    TIdMap frgIdMap;
    TIdMap readIdMap;
    // For all paired reads (inferred from FRG), a mapping from the AMOS read id to the paired match id of the one
    // alignment read from the AMOS file.  Note that this has the assumption that there is only one alignment per read
    // in the AMOS file.
    TIdMap readToPairMatchId;
    // The id of the next pair match.
    unsigned nextPairMatchId = 0;

    typename DirectionIterator<TFile, Input>::Type reader(file);

    // Parse the file and convert the internal ids
    if (atEnd(reader))
        return 1;

    CharString  blockIdentifier, fieldIdentifier, eid, buffer;
    TReadSeq    readSeq;
    CharString  contigSeq;
    CharString  qual;
    TId _id;

    while (!atEnd(reader))
    {
        // New block?
        if (value(reader) == '{')
        {
            skipOne(reader, EqualsChar<'{'>());

            clear(blockIdentifier);
            readUntil(blockIdentifier, reader, NotFunctor<IsAlpha>());
            skipLine(reader);

            // reset id and eid
            _id = TReadStoreElement::INVALID_ID;
            clear(eid);

            // Library block
            if (blockIdentifier == "LIB")
            {
                // The LIB block contains the mean template length/insert size and standard deviation thereof.  Besides
                // the numeric id and (iid) a (external) name of the library (eid).
                TLibraryStoreElement libEl;
                while (!atEnd(reader) && value(reader) != '}')
                {
                    // read "KEY:VALUE" pair
                    readAmosKeyValue(fieldIdentifier, buffer, reader);

                    if (fieldIdentifier == "iid")
                        _id = lexicalCast<TId>(buffer);
                    else if (fieldIdentifier == "eid")
                        eid = buffer;
                    else if (fieldIdentifier == "mea")
                        libEl.mean = lexicalCast<double>(buffer);
                    else if (fieldIdentifier == "std")
                        libEl.std = lexicalCast<double>(buffer);

                    // Always skip over the remainder of the line.
                    skipUntil(reader, NotFunctor<IsWhitespace>());
                }
                skipOne(reader, EqualsChar<'}'>());

                // Insert library information into the fragmentstore.
                libIdMap.insert(std::make_pair(_id, length(fragStore.libraryStore)));
                appendValue(fragStore.libraryStore, libEl, Generous());
                appendValue(fragStore.libraryNameStore, eid, Generous());
            }
            // Fragment block
            else if (blockIdentifier == "FRG")
            {
                // The FRG block contains the parent library (lib), links to the paired sequencing reads (rds), internal
                // id (iid) and external id (eid).
                bool foundRds = false;
                TMatePairElement matePairEl;
                while (!atEnd(reader) && value(reader) != '}')
                {
                    // read "KEY:VALUE" pair
                    readAmosKeyValue(fieldIdentifier, buffer, reader);

                    if (fieldIdentifier == "iid")
                        _id = lexicalCast<TId>(buffer);
                    else if (fieldIdentifier == "eid")
                        eid = buffer;
                    else if (fieldIdentifier == "lib")
                        matePairEl.libId = lexicalCast<TId>(buffer);
                    else if (fieldIdentifier == "rds")
                    {
                        foundRds = true;
                        size_t comma = findFirst(buffer, EqualsChar<','>());
                        matePairEl.readId[0] = lexicalCast<TId>(prefix(buffer, comma));
                        matePairEl.readId[1] = lexicalCast<TId>(suffix(buffer, comma + 1));

                        // Store mapping to pair match id.
                        readToPairMatchId[matePairEl.readId[0]] = nextPairMatchId;
                        readToPairMatchId[matePairEl.readId[1]] = nextPairMatchId;
                        nextPairMatchId++;
                    }
                    // Always skip over the remainder of the line.
                    skipUntil(reader, NotFunctor<IsWhitespace>());
                }
                skipOne(reader, EqualsChar<'}'>());

                // Only insert valid mate pairs
                if (foundRds)
                {
                    frgIdMap.insert(std::make_pair(_id, length(fragStore.matePairStore)));
                    appendValue(fragStore.matePairStore, matePairEl, Generous());
                    appendValue(fragStore.matePairNameStore, eid, Generous());
                }
            }
            // Read block
            else if (blockIdentifier == "RED")
            {
                // If matePairId is not updated, this yields to a singleton read below.
                clear(readSeq);
                clear(qual);
                TId matePairId = TReadStoreElement::INVALID_ID;
                while (!atEnd(reader) && value(reader) != '}')
                {
                    // read "KEY:VALUE" pair
                    readAmosKeyValue(fieldIdentifier, buffer, reader);

                    if (fieldIdentifier == "iid")
                        _id = lexicalCast<TId>(buffer);
                    else if (fieldIdentifier == "eid")
                        eid = buffer;
                    else if (fieldIdentifier == "frg")
                        matePairId = lexicalCast<TId>(buffer);
                    else if (fieldIdentifier == "seq")
                    {
                        readUntil(readSeq, reader, EqualsChar<'.'>(), IsWhitespace());
                        skipOne(reader, EqualsChar<'.'>());
                    }
                    else if (fieldIdentifier == "qlt")
                    {
                        // TODO(holtgrew): Problematic if . part of qualities, need to count chars.
                        readUntil(qual, reader, EqualsChar<'.'>(), IsWhitespace());
                        skipOne(reader, EqualsChar<'.'>());
                    }
                    // Always skip over the remainder of the line.
                    skipUntil(reader, NotFunctor<IsWhitespace>());
                }
                // Set quality
                assignQualities(readSeq, qual);
                skipOne(reader, EqualsChar<'}'>());

                // Insert the read
                readIdMap.insert(std::make_pair(_id, length(fragStore.readStore)));
                appendRead(fragStore, readSeq, matePairId);
                appendValue(fragStore.readNameStore, eid, Generous());
            }
            // Contig block
            else if (blockIdentifier == "CTG")
            {
                clear(contigSeq);
                clear(qual);
                TContigElement contigEl;
                TSize fromAligned = length(fragStore.alignedReadStore);
                while (!atEnd(reader) && value(reader) != '}')
                {
                    // Are we entering a TLE block
                    if (value(reader) == '{')
                    {
                        // skip "{TLE\n"
                        skipOne(reader, EqualsChar<'{'>());
                        skipOne(reader, EqualsChar<'T'>());
                        skipOne(reader, EqualsChar<'L'>());
                        skipOne(reader, EqualsChar<'E'>());
                        skipUntil(reader, NotFunctor<IsWhitespace>());

                        TAlignedElement alignEl;
                        typedef typename TFragmentStore::TContigPos TContigPos;
                        TContigPos offsetPos = 0;
                        TContigPos clr1 = 0;
                        TContigPos clr2 = 0;
                        String<TContigPos> gaps;
                        while (!atEnd(reader) && value(reader) != '}')
                        {
                            // read "KEY:VALUE" pair
                            readAmosKeyValue(fieldIdentifier, buffer, reader);

                            if (fieldIdentifier == "src")
                                alignEl.readId = lexicalCast<TId>(buffer);
                            else if (fieldIdentifier == "off")
                            {
                                if (buffer != "-")
                                    offsetPos = lexicalCast<TContigPos>(buffer);
                                else
                                    offsetPos = 0;
                            }
                            else if (fieldIdentifier == "clr")
                            {
                                size_t comma = findFirst(buffer, EqualsChar<','>());
                                clr1 = lexicalCast<TContigPos>(prefix(buffer, comma));
                                clr2 = lexicalCast<TContigPos>(suffix(buffer, comma + 1));
                            }
                            else if (fieldIdentifier == "gap")
                            {
                                skipUntil(reader, NotFunctor<IsWhitespace>());
                                while (!atEnd(reader) && value(reader) != '.')
                                {
                                    clear(buffer);
                                    readUntil(buffer, reader, IsWhitespace());
                                    appendValue(gaps, lexicalCast<TContigPos>(buffer));
                                    skipUntil(reader, NotFunctor<IsWhitespace>());
                                }
                                skipOne(reader, EqualsChar<'.'>());
                            }
                            // Always skip over the remainder of the line.
                            skipUntil(reader, NotFunctor<IsWhitespace>());
                        }
                        skipOne(reader, EqualsChar<'}'>());
                        skipUntil(reader, NotFunctor<IsWhitespace>());

                        // Get the length of the read
                        TId readId = (readIdMap.find(alignEl.readId))->second;
                        TSize lenRead = length(value(fragStore.readSeqStore, readId));

                        // Create the gap anchors
                        typedef typename TFragmentStore::TContigGapAnchor TContigGapAnchor;
                        int offset = 0;
                        if ((clr1 < clr2) && (clr1>0)) offset = clr1;
                        else if ((clr1 > clr2) && (clr1 < static_cast<TContigPos>(lenRead))) offset = lenRead - clr1;
                        int diff = -1 * (int) (offset);
                        // Clipped begin
                        if (offset != 0) appendValue(alignEl.gaps, TContigGapAnchor(offset, 0), Generous());
                        // Internal gaps
                        typedef typename Iterator<String<TContigPos>, Standard>::Type TPosIter;
                        TPosIter posIt = begin(gaps, Standard());
                        TPosIter posItEnd = end(gaps, Standard());
                        TContigPos lastGap = 0;
                        TSize gapLen = 0;
                        TSize totalGapLen = 0;
                        for(;posIt!=posItEnd; goNext(posIt)) {
                            if (gapLen == 0) {
                                ++gapLen; ++totalGapLen;
                                ++diff;
                                lastGap = value(posIt);
                            }
                            else if (lastGap == value(posIt)) {
                                ++gapLen; ++totalGapLen;
                                ++diff;
                            }
                            else {
                                appendValue(alignEl.gaps, TContigGapAnchor(offset + lastGap, offset + lastGap + diff), Generous());
                                gapLen = 1; ++totalGapLen;
                                lastGap = value(posIt);
                                ++diff;
                            }
                        }
                        if (gapLen > 0) appendValue(alignEl.gaps, TContigGapAnchor(offset + lastGap, offset + lastGap + diff), Generous());
                        // Clipped end
                        if ((clr1 < clr2) && (clr2 < static_cast<TContigPos>(lenRead))) {
                            diff -= (lenRead - clr2);
                            appendValue(alignEl.gaps, TContigGapAnchor(lenRead, lenRead + diff), Generous());
                        } else if ((clr1 > clr2) && (clr2 > 0)) {
                            diff -= clr2;
                            appendValue(alignEl.gaps, TContigGapAnchor(lenRead, lenRead + diff), Generous());
                        }

                        // Set begin and end position
                        if (clr1 < clr2) {
                            alignEl.beginPos = offsetPos;
                            alignEl.endPos = offsetPos + totalGapLen + (clr2 - clr1);
                        } else {
                            alignEl.beginPos = offsetPos + totalGapLen + (clr1 - clr2);
                            alignEl.endPos = offsetPos;
                        }

                        // Append new align fragment, note: contigId must still be set
                        alignEl.id = length(fragStore.alignedReadStore);
                        appendValue(fragStore.alignedReadStore, alignEl, Generous());
                    }
                    else
                    {
                        // read "KEY:VALUE" pair
                        readAmosKeyValue(fieldIdentifier, buffer, reader);

                        if (fieldIdentifier == "iid")
                            _id = lexicalCast<TId>(buffer);
                        else if (fieldIdentifier == "eid")
                            eid = buffer;
                        else if (fieldIdentifier == "seq")
                        {
                            readUntil(contigSeq, reader, EqualsChar<'.'>(), IsWhitespace());
                            skipOne(reader, EqualsChar<'.'>());
                        }
                        else if (fieldIdentifier == "qlt")
                        {
                            // TODO(holtgrew): Problematic if . part of qualities, need to count chars.
                            readUntil(qual, reader, EqualsChar<'.'>(), IsWhitespace());
                            skipOne(reader, EqualsChar<'.'>());
                        }
                        // Always skip over the remainder of the line.
                        skipUntil(reader, NotFunctor<IsWhitespace>());
                    }
                }
                // Set quality
                skipOne(reader, EqualsChar<'}'>());

                // Create the gap anchors
                char gapChar = gapValue<char>();
                typedef typename Iterator<CharString, Standard>::Type TStringIter;
                TStringIter seqIt = begin(contigSeq, Standard());
                TStringIter seqItEnd = end(contigSeq, Standard());
                TStringIter qualIt = begin(qual, Standard());
                typedef typename TFragmentStore::TReadPos TPos;
                typedef typename TFragmentStore::TContigGapAnchor TContigGapAnchor;
                TPos ungappedPos = 0;
                TPos gappedPos = 0;
                bool gapOpen = false;
                for (; seqIt != seqItEnd; goNext(seqIt), goNext(qualIt), ++gappedPos)
                {
                    if (value(seqIt) == gapChar)
                    {
                        gapOpen = true;
                    }
                    else
                    {
                        if (gapOpen)
                        {
                            appendValue(contigEl.gaps, TContigGapAnchor(ungappedPos, gappedPos), Generous());
                            gapOpen = false;
                        }
                        typename Value<TContigSeq>::Type letter = getValue(seqIt);
                        assignQualityValue(letter, getValue(qualIt));
                        appendValue(contigEl.seq, letter, Generous());
                        ++ungappedPos;
                    }
                }
                if (gapOpen)
                    appendValue(contigEl.gaps, TContigGapAnchor(ungappedPos, gappedPos), Generous());

                // Set the contigId in all aligned reads
                TSize toAligned = length(fragStore.alignedReadStore);
                TId newContigId = length(fragStore.contigStore);
                for (; fromAligned < toAligned; ++fromAligned)
                    fragStore.alignedReadStore[fromAligned].contigId = newContigId;

                // Insert the contig
                appendValue(fragStore.contigStore, contigEl, Generous());
                appendValue(fragStore.contigNameStore, eid, Generous());
            }
            // Unknown block identifier
            else
            {
                // skip until closing bracket
                while (!atEnd(reader) && value(reader) != '}')
                    skipLine(reader);
            }
        }
        else
        {
            // not a block - skip line
            skipLine(reader);
        }
    }

    // Renumber all ids
    typedef typename TIdMap::const_iterator TIdMapIter;
    typedef typename Iterator<typename TFragmentStore::TMatePairStore>::Type TMateIter;
    TMateIter mateIt = begin(fragStore.matePairStore);
    TMateIter mateItEnd = end(fragStore.matePairStore);
    for(;mateIt != mateItEnd; goNext(mateIt)) {
        if (mateIt->libId != TMatePairElement::INVALID_ID) {
            TIdMapIter libIdPos = libIdMap.find(mateIt->libId);
            if (libIdPos != libIdMap.end())
                mateIt->libId = libIdPos->second;
            else
                mateIt->libId = TMatePairElement::INVALID_ID;
        }
        if (mateIt->readId[0] != TMatePairElement::INVALID_ID) {
            TIdMapIter readIdPos = readIdMap.find(mateIt->readId[0]);
            if (readIdPos != readIdMap.end())
                mateIt->readId[0] = readIdPos->second;
            else
                mateIt->readId[0] = TMatePairElement::INVALID_ID;
        }
        if (mateIt->readId[1]!= TMatePairElement::INVALID_ID) {
            TIdMapIter readIdPos = readIdMap.find(mateIt->readId[1]);
            if (readIdPos != readIdMap.end())
                mateIt->readId[1] = readIdPos->second;
            else
                mateIt->readId[0] = TMatePairElement::INVALID_ID;
        }
    }

    // Copy data from frgIdMap into the matePairId members of the readStore.
    typedef typename Iterator<typename TFragmentStore::TReadStore>::Type TReadIter;
    TReadIter readIt = begin(fragStore.readStore);
    TReadIter readItEnd = end(fragStore.readStore);
    for (;readIt != readItEnd; goNext(readIt))
    {
        if (readIt->matePairId != TReadStoreElement::INVALID_ID)
        {
            TIdMapIter mateIdPos = frgIdMap.find(readIt->matePairId);
            if (mateIdPos != frgIdMap.end())
                readIt->matePairId = mateIdPos->second;
            else
                readIt->matePairId = TReadStoreElement::INVALID_ID;
        }
    }

    // Copy data from readIdMap into the pairMatchId entries of the alignedReadStore.
    typedef typename Iterator<typename TFragmentStore::TAlignedReadStore>::Type TAlignIter;
    TAlignIter alignIt = begin(fragStore.alignedReadStore);
    TAlignIter alignItEnd = end(fragStore.alignedReadStore);
    for (;alignIt != alignItEnd; goNext(alignIt))
    {
        if (alignIt->readId != TAlignedElement::INVALID_ID)
        {
            TIdMapIter readIdPos = readIdMap.find(alignIt->readId);
            if (readIdPos != readIdMap.end())
            {
                //SEQAN_ASSERT(readToPairMatchId.find(alignIt->readId) != readToPairMatchId.end());
                if (readToPairMatchId.find(alignIt->readId) != readToPairMatchId.end())
                    alignIt->pairMatchId = readToPairMatchId[alignIt->readId];
                alignIt->readId = readIdPos->second;
            }
            else
            {
                alignIt->readId = TAlignedElement::INVALID_ID;
            }
        }
    }

    return 0;
}


//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn FragmentStore#write
 * @brief Write the contents of a FragmentStore to a file.
 *
 * @signature int write(file, store, tag);
 *
 * @param[in,out] file  The @link StreamConcept @endlink to write to.
 * @param[in]     store The FragmentStore to write to the file.
 * @param[in]     tag   The format to write out.  Types: Sam or Amos.
 *
 * @return int 0 in case of success, 1 in case of errors.
 */

template <typename TTarget, typename TSpec, typename TConfig>
inline void
write(TTarget & target,
      FragmentStore<TSpec, TConfig> & fragStore,
      Amos)
{
    // Basic types
    typedef FragmentStore<TSpec, TConfig> TFragmentStore;
    //typedef typename Id<TFragmentStore>::Type TId;
    typedef typename Size<TFragmentStore>::Type TSize;
    //typedef typename Value<TFile>::Type TValue;

    // All fragment store element types
    //typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigElement;
    //typedef typename Value<typename TFragmentStore::TLibraryStore>::Type TLibraryStoreElement;
    typedef typename Value<typename TFragmentStore::TMatePairStore>::Type TMatePairElement;
    typedef typename Value<typename TFragmentStore::TReadStore>::Type TReadStoreElement;
    //typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;

    typename DirectionIterator<TTarget, Output>::Type iter = directionIterator(target, Output());

    // Write Header
    write(iter, "{UNV\niid:1\neid:seqan\ncom:\nafg file created with SeqAn\n.\n}\n");

    // Write Libraries
    typedef typename Iterator<typename TFragmentStore::TLibraryStore, Standard>::Type TLibIter;
    TLibIter libIt = begin(fragStore.libraryStore, Standard() );
    TLibIter libItEnd = end(fragStore.libraryStore, Standard() );
    for(TSize idCount = 0;libIt != libItEnd; goNext(libIt), ++idCount)
    {
        write(iter, "{LIB\niid:");
        appendNumber(iter, idCount + 1);
        writeValue(iter, '\n');
        if (!empty(fragStore.libraryNameStore) && !empty(fragStore.libraryNameStore[idCount]))
        {
            write(iter, "eid:");
            write(iter, fragStore.libraryNameStore[idCount]);
            writeValue(iter, '\n');
        }
        write(iter, "{DST\nmea:");
        appendNumber(iter, libIt->mean);
        writeValue(iter, '\n');
        write(iter, "std:");
        appendNumber(iter, libIt->std);
        writeValue(iter, '\n');
        write(iter, "}\n}\n");
    }

    // Write Fragments / mate pairs
    typedef typename Iterator<typename TFragmentStore::TMatePairStore, Standard>::Type TMateIter;
    TMateIter mateIt = begin(fragStore.matePairStore, Standard() );
    TMateIter mateItEnd = end(fragStore.matePairStore, Standard() );
    for(TSize idCount = 0;mateIt != mateItEnd; goNext(mateIt), ++idCount)
    {
        write(iter, "{FRG\niid:");
        appendNumber(iter, idCount + 1);
        writeValue(iter, '\n');
        if (!empty(fragStore.matePairNameStore) && !empty(fragStore.matePairNameStore[idCount]))
        {
            write(iter, "eid:");
            write(iter, fragStore.matePairNameStore[idCount]);
            writeValue(iter, '\n');
        }
        write(iter, "lib:");
        appendNumber(iter, mateIt->libId + 1);
        writeValue(iter, '\n');
        if (mateIt->readId[0] != TMatePairElement::INVALID_ID && mateIt->readId[1] != TMatePairElement::INVALID_ID)
        {
            write(iter, "rds:");
            appendNumber(iter, mateIt->readId[0] + 1);
            writeValue(iter, ',');
            appendNumber(iter, mateIt->readId[1] + 1);
            writeValue(iter, '\n');
        }
        write(iter, "}\n");
    }

    // Get clear ranges
    typedef Pair<typename TFragmentStore::TReadPos, typename TFragmentStore::TReadPos> TClrRange;
    String<TClrRange> clrRange;
    resize(clrRange, length(fragStore.readStore), TClrRange(0,0));
    typedef typename Iterator<typename TFragmentStore::TAlignedReadStore, Standard>::Type TAlignIter;
    TAlignIter alignIt = begin(fragStore.alignedReadStore, Standard() );
    TAlignIter alignItEnd = end(fragStore.alignedReadStore, Standard() );
    for(;alignIt != alignItEnd; goNext(alignIt))
    {
        typename TFragmentStore::TReadPos begClr = 0;
        typename TFragmentStore::TReadPos endClr = 0;
        getClrRange(fragStore, value(alignIt), begClr, endClr);
        clrRange[alignIt->readId] = TClrRange(begClr, endClr);
    }

    // Write reads
    typedef typename Iterator<typename TFragmentStore::TReadStore, Standard>::Type TReadIter;
    TReadIter readIt = begin(fragStore.readStore, Standard() );
    TReadIter readItEnd = end(fragStore.readStore, Standard() );
    for(TSize idCount = 0;readIt != readItEnd; ++readIt, ++idCount)
    {
        // Skip reads without a name.
        if (length(value(fragStore.readNameStore, idCount)) == 0u)
            continue;
        write(iter, "{RED\niid:");
        appendNumber(iter, idCount + 1);
        writeValue(iter, '\n');
        if (!empty(fragStore.readNameStore) && !empty(fragStore.readNameStore[idCount]))
        {
            write(iter, "eid:");
            write(iter, fragStore.readNameStore[idCount]);
            writeValue(iter, '\n');
        }
        write(iter, "seq:\n");
        writeWrappedString(iter, fragStore.readSeqStore[idCount], 60);
        write(iter, ".\nqlt:\n");

        typedef typename Value<typename TFragmentStore::TReadSeqStore>::Type TReadSeq;
        typedef QualityExtractor<typename Value<TReadSeq>::Type> TQualityExtractor;
        ModifiedString<TReadSeq const, ModView<TQualityExtractor> > quals(fragStore.readSeqStore[idCount]);
        writeWrappedString(iter, quals, 60);
        write(iter,  ".\n");
        if (readIt->matePairId != TReadStoreElement::INVALID_ID)
        {
            write(iter, "frg:");
            appendNumber(iter, readIt->matePairId + 1);
            writeValue(iter, '\n');
        }
        if (clrRange[idCount].i1 != clrRange[idCount].i2)
        {
            write(iter, "clr:");
            appendNumber(iter, clrRange[idCount].i1);
            writeValue(iter, ',');
            appendNumber(iter, clrRange[idCount].i2);
            writeValue(iter, '\n');
        }
        write(iter, "}\n");
    }

    // Sort aligned reads according to contigId
    sortAlignedReads(fragStore.alignedReadStore, SortContigId());

    // Write Contigs
    typedef typename Iterator<typename TFragmentStore::TContigStore, Standard>::Type TContigIter;
    TContigIter contigIt = begin(fragStore.contigStore, Standard() );
    TContigIter contigItEnd = end(fragStore.contigStore, Standard() );
    alignIt = begin(fragStore.alignedReadStore);
    alignItEnd = end(fragStore.alignedReadStore);
    for(TSize idCount = 0;contigIt != contigItEnd; goNext(contigIt), ++idCount)
    {
        write(iter, "{CTG\niid:");
        appendNumber(iter, idCount + 1);
        writeValue(iter, '\n');
        if (!empty(fragStore.contigNameStore) && !empty(fragStore.contigNameStore[idCount]))
        {
            write(iter, "eid:");
            write(iter, fragStore.contigNameStore[idCount]);
            writeValue(iter, '\n');
        }
        String<char> qlt;
        write(iter, "seq:\n");
        typedef typename Iterator<typename TFragmentStore::TContigSeq>::Type TContigIter;
        TContigIter seqContigIt = begin(contigIt->seq);
        TContigIter seqContigItEnd = end(contigIt->seq);
        typedef typename Iterator<String<typename TFragmentStore::TContigGapAnchor> >::Type TGapsIter;
        TGapsIter itGaps = begin(contigIt->gaps);
        TGapsIter itGapsEnd = end(contigIt->gaps);
        int diff = 0;
        char gapChar = gapValue<char>();
        typename TFragmentStore::TContigPos mySeqPos = 0;
        TSize k = 0;
        for(;itGaps != itGapsEnd; goNext(itGaps))
        {
            while (mySeqPos < itGaps->seqPos)
            {

                if ((k % 60 == 0) && (k != 0))
                    writeValue(iter, '\n');
                ++k;
                writeValue(iter, value(seqContigIt));
                char c = ' ';
                convertQuality(c, getQualityValue(value(seqContigIt)));
                appendValue(qlt, c, Generous());
                goNext(seqContigIt);++mySeqPos;
            }
            for(int i = 0; i < ((int) itGaps->gapPos - (int) itGaps->seqPos) - diff; ++i)
            {
                if ((k % 60 == 0) && (k != 0))
                    writeValue(iter, '\n');
                ++k;
                writeValue(iter, gapChar);
                appendValue(qlt, '0', Generous());
            }
            diff = (itGaps->gapPos - itGaps->seqPos);
        }
        for(;seqContigIt != seqContigItEnd; goNext(seqContigIt))
        {
            if ((k % 60 == 0) && (k != 0))
                writeValue(iter, '\n');
            ++k;
            writeValue(iter, value(seqContigIt));
            char c = ' ';
            convertQuality(c, getQualityValue(value(seqContigIt)));
            appendValue(qlt, c, Generous());
        }
        write(iter, "\n.\nqlt:\n");
        writeWrappedString(iter, qlt, 60);
        write(iter, ".\n");

        while (alignIt != alignItEnd && idCount < alignIt->contigId)
            goNext(alignIt);

        for (; alignIt != alignItEnd && idCount == alignIt->contigId; goNext(alignIt))
        {
            write(iter, "{TLE\nsrc:");
            appendNumber(iter, alignIt->readId + 1);
            writeValue(iter, '\n');
            typedef typename Iterator<String<typename TFragmentStore::TReadGapAnchor> >::Type TReadGapsIter;
            TReadGapsIter itGaps = begin(alignIt->gaps);
            TReadGapsIter itGapsEnd = end(alignIt->gaps);

            // Create the gaps string and the clear ranges
            typename TFragmentStore::TReadPos lenRead = length(value(fragStore.readSeqStore, alignIt->readId));
            TSize clr1 = 0;
            TSize clr2 = lenRead;
            // Create first clear range
            if (itGaps != itGapsEnd && itGaps->gapPos == 0)
                clr1 = itGaps->seqPos;
            int diff = clr1;
            String<unsigned int> gaps;
            for (; itGaps != itGapsEnd; goNext(itGaps))
            {
                for (int i = 0; i + itGaps->seqPos < diff + itGaps->gapPos; ++i)
                    appendValue(gaps, itGaps->seqPos - clr1, Generous());
                // Clipped sequence
                if (diff + itGaps->gapPos < itGaps->seqPos)
                    clr2 = lenRead + diff + itGaps->gapPos - itGaps->seqPos;
                diff = (int)itGaps->seqPos - (int)itGaps->gapPos;
            }
            if (alignIt->beginPos > alignIt->endPos)
            {
                clr1 = lenRead - clr1;
                clr2 = lenRead - clr2;
            }
            write(iter, "off:");
            appendNumber(iter, std::min(alignIt->beginPos, alignIt->endPos));
            write(iter, "\nclr:");
            appendNumber(iter, clr1);
            writeValue(iter, ',');
            appendNumber(iter, clr2);
            writeValue(iter, '\n');
            if (!empty(gaps))
            {
                write(iter, "gap:\n");
                for (TSize z = 0; z < length(gaps); ++z)
                {
                    appendNumber(iter, gaps[z]);
                    writeValue(iter, '\n');
                }
                write(iter, ".\n");
            }
            write(iter, "}\n");
        }
        write(iter, "}\n");
    }
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn FragmentStore#writeContigs
 * @brief Write contigs from FragmentStore into a @link StreamConcept @endlink.
 *
 * @signature bool writeContigs(file, store, tag);
 *
 * @param[in,out] file  The @link StreamConcept @endlink to write to.
 * @param[in]     store The FragmentStore to write contigs of.
 * @param[in]     tag   A tag for the sequence format.
 *
 * @return bool true on success, false on errors.
 */

template <typename TSpec, typename TFSSpec, typename TFSConfig>
bool writeContigs(FormattedFile<Fastq, Output, TSpec> & file, FragmentStore<TFSSpec, TFSConfig> & store)
{
//IOREV _doc_
    for (unsigned i = 0; i < length(store.contigNameStore); ++i)
        writeRecord(file, store.contigNameStore[i], store.contigStore[i].seq);
    return true;
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn FragmentStore#loadContigs
 * @brief Load contigs into a FragmentStore.
 *
 * @signature bool loadContigs(store, fileName[, loadSeqs]);
 * @signature bool loadContigs(store, fileNameList[, loadSeqs]);
 *
 * @param[in,out] store        The FragmentStore to append the contigs to.
 * @param[in]     fileName     A @link CharString @endlink with the name of the file to load.
 * @param[in]     fileNameList A @link StringSet @endlink of @link CharString @endlink with a list of file names to
 *                             load.
 * @param[in]     loadSeqs     A <tt>bool</tt> indicating whether to load lazily.  If <tt>true</tt> then sequences are
 *                             loaded immediately.  If <tt>false</tt>, an emptycontig with a reference to the file is
 *                             created.  Its sequence can be loaded on demand by @link FragmentStore#lockContig
 *                             @endlink and @link FragmentStore#loadContig @endlink.
 *
 * @return bool true in case of success and false in case of error.
 */

template <typename TFSSpec, typename TFSConfig>
bool loadContigs(FragmentStore<TFSSpec, TFSConfig> &store, StringSet<CharString> const &fileNameList, bool loadSeqs)
{
//IOREV _nodoc_ although there is dddoc, there is no entry in html-doc
    typedef FragmentStore<TFSSpec, TFSConfig>            TFragmentStore;
    typedef typename TFragmentStore::TContigStore        TContigStore;
    typedef typename TFragmentStore::TContigFileStore    TContigFileStore;
    typedef typename Value<TContigStore>::Type            TContig;
    typedef typename Value<TContigFileStore>::Type        TContigFile;

    SeqFileIn seqFile;
    CharString meta, seq;

    for (unsigned f = 0; f < length(fileNameList); ++f)
    {
        if (!open(seqFile, toCString(fileNameList[f])))
            return false;

        TContigFile contigFile;
        contigFile.fileName = fileNameList[f];
        contigFile.firstContigId = length(store.contigStore);
        appendValue(store.contigFileStore, contigFile, Generous());

        while (!atEnd(seqFile))
        {
            resize(store.contigStore, length(store.contigStore) + 1, Generous());

            TContig & contig = back(store.contigStore);
            contig.usage = 0;
            contig.fileId = length(store.contigFileStore) - 1;
            contig.fileBeginPos = position(seqFile);

            if (loadSeqs)
                readRecord(meta, contig.seq, seqFile);
            else
                readRecord(meta, seq, seqFile);

            contig.fileEndPos = position(seqFile);

            cropAfterFirst(meta, IsWhitespace());
            appendValue(store.contigNameStore, meta, Generous());
        }
        close(seqFile);
    }
    return true;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFSSpec, typename TFSConfig>
bool loadContigs(FragmentStore<TFSSpec, TFSConfig> &store, CharString const &fileName, bool loadSeqs)
{
    StringSet<CharString> fileNames;
    appendValue(fileNames, fileName);
    return loadContigs(store, fileNames, loadSeqs);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFSSpec, typename TFSConfig, typename TFileNames>
bool loadContigs(FragmentStore<TFSSpec, TFSConfig> &store, TFileNames const &fileNames)
{
    return loadContigs(store, fileNames, true);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn FragmentStore#loadContig
 * @brief Manually load a contig sequence.
 *
 * @signature bool loadContig(store, contigId);
 *
 * @param[in,out] store    The FragmentStore to load the contig for.
 * @param[in]     contigId The id of the contig that was created earlier by @link FragmentStore#loadContigs @endlink.
 *
 * @return bool true on success, false on failure.
 */

template <typename TSpec, typename TConfig, typename TId>
bool loadContig(FragmentStore<TSpec, TConfig> &store, TId _id)
{
    typedef FragmentStore<TSpec, TConfig>                TFragmentStore;
    typedef typename TFragmentStore::TContigStore        TContigStore;
    typedef typename TFragmentStore::TContigFileStore    TContigFileStore;
    typedef typename Value<TContigStore>::Type            TContig;
    typedef typename Value<TContigFileStore>::Type        TContigFile;

    if ((TId)length(store.contigStore) <= _id) return false;
    TContig &contig = store.contigStore[_id];

    if (contig.fileId >= length(store.contigFileStore)) return false;

    TContigFile &contigFile = store.contigFileStore[contig.fileId];
    CharString meta;                                // dummy (seq name is already in the name store)

    SeqFileIn seqFile(toCString(contigFile.fileName));
    setPosition(seqFile, contig.fileBeginPos);
    readRecord(meta, contig.seq, seqFile);            // read contig sequence

    return true;
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn FragmentStore#lockContig
 * @brief Locks a contig sequence from being removed.
 *
 * This function increases the contig usage counter by 1 and ensures that the contig sequence is loaded.
 *
 * @signature bool lockContig(store, contigId);
 *
 * @param[in,out] store    The FragmentStore to lock the contig for.
 * @param[in]     contigId The id of the contig that was created earlier by @link FragmentStore#loadContigs @endlink.
 *
 * @return bool true on success, false on failure.
 */

template <typename TSpec, typename TConfig, typename TId>
bool lockContig(FragmentStore<TSpec, TConfig> &store, TId _id)
{
    typedef FragmentStore<TSpec, TConfig>                TFragmentStore;
    typedef typename TFragmentStore::TContigStore        TContigStore;
    typedef typename Value<TContigStore>::Type            TContig;

    if ((TId)length(store.contigStore) <= _id) return false;
    TContig &contig = store.contigStore[_id];

    if (contig.usage++ > 0 || !empty(contig.seq)) return true;
    return loadContig(store, _id);
}

/*!
 * @fn FragmentStore#unlockContig
 * @brief Removes a previous contig lock.
 *
 * This function decreases the contig usage counter by 1.
 *
 * @signature bool unlockContig(store, contigId);
 *
 * @param[in,out] store    The FragmentStore to unlock the contig for.
 * @param[in]     contigId The id of the contig that was created earlier by @link FragmentStore#loadContigs @endlink.
 *
 * @return bool true on success, false on failure.
 */

template <typename TSpec, typename TConfig, typename TId>
bool unlockContig(FragmentStore<TSpec, TConfig> &store, TId _id)
{
    if ((TId)length(store.contigStore) <= _id) return false;
    --store.contigStore[_id].usage;
    return true;
}

/*!
 * @fn FragmentStore#unlockAndFreeContig
 * @brief Removes a previous contig lock and clears the sequence if no further lock exists.
 *
 * This function decreases the contig usage counter by 1 and frees the sequences' memory if the counter equals 0.
 *
 * @signature bool unlockContig(store, contigId);
 *
 * @param[in,out] store    The FragmentStore to unlock the contig for.
 * @param[in]     contigId The id of the contig that was created earlier by @link FragmentStore#loadContigs @endlink.
 *
 * @return bool true on success, false on failure.
 */

template <typename TSpec, typename TConfig, typename TId>
bool unlockAndFreeContig(FragmentStore<TSpec, TConfig> &store, TId _id)
{
    typedef FragmentStore<TSpec, TConfig>                TFragmentStore;
    typedef typename TFragmentStore::TContigStore        TContigStore;
    typedef typename Value<TContigStore>::Type            TContig;

    if ((TId)length(store.contigStore) <= _id) return false;
    TContig &contig = store.contigStore[_id];

    if (--contig.usage == 0 && contig.fileId < length(store.contigFileStore))
    {
        typename TContig::TContigSeq emptySeq;
        swap(contig.seq, emptySeq);
        return true;
    }
    return false;
}

/*!
 * @fn FragmentStore#lockContigs
 * @brief Locks all contig sequences from being remove.
 *
 * @signature bool lockContigs(store);
 *
 * @param[in,out] store    The FragmentStore to lock the contigs for.
 *
 * @return bool true in case of success, false in case of errors.
 */

template <typename TSpec, typename TConfig>
bool lockContigs(FragmentStore<TSpec, TConfig> &store)
{
    bool result = true;
    for (unsigned _id = 0; _id < length(store.contigStore); ++_id)
        result &= lockContig(store, _id);
    return result;
}

/*!
 * @fn FragmentStore#unlockContigs
 * @brief Unlocks all contig sequences.
 *
 * @signature bool unlockContigs(store);
 *
 * @param[in,out] store    The FragmentStore to unlock the contigs for.
 *
 * @return bool true in case of success, false in case of errors.
 */

template <typename TSpec, typename TConfig>
bool unlockContigs(FragmentStore<TSpec, TConfig> &store)
{
    bool result = true;
    for (unsigned _id = 0; _id < length(store.contigStore); ++_id)
        result &= unlockContig(store, _id);
    return result;
}

/*!
 * @fn FragmentStore#unlockAndFreeContigs
 * @brief Unlocks all contig sequences and clears sequences without lock.
 *
 * @signature bool unlockAndFreeContigs(store);
 *
 * @param[in,out] store    The FragmentStore to unlock the contigs for.
 *
 * @return bool true in case of success, false in case of errors.
 */

template <typename TSpec, typename TConfig>
bool unlockAndFreeContigs(FragmentStore<TSpec, TConfig> &store)
{
    bool result = true;
    for (unsigned _id = 0; _id < length(store.contigStore); ++_id)
        result &= unlockAndFreeContig(store, _id);
    return result;
}


//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn FragmentStore#loadReads
 * @brief Loads reads into FragmentStore
 *
 * When two file names are given thent he files are expected to containt he same number of reads and reads with the same
 * index are assumed to be mate pairs.  Mate pairs are stored internally in an "interleaved mode": a read is read from
 * each file before reading the next one.
 *
 * @signature bool loadReads(store, fileName);
 * @signature bool loadReads(store, fileNameL, fileNameR);
 *
 * @param[in,out] store     The FragmentStore to append the reads to.
 * @param[in]     fileName  Path to single-end read file.
 * @param[in]     fileNameL Path to left read file in case of paired reads.
 * @param[in]     fileNameR Path to right read file in case of paired reads.
 *
 * @return bool true in case of success, false in case of errors.
 */

template <typename TSpec, typename TConfig, typename TFileName>
bool loadReads(FragmentStore<TSpec, TConfig> &store, TFileName &fileName)
{
    SEQAN_ASSERT_EQ(length(store.readStore), length(store.readSeqStore));
    SEQAN_ASSERT_EQ(length(store.readStore), length(store.readNameStore));

    typedef typename FragmentStore<TSpec, TConfig>::TReadStore TReadStore;
    typedef typename Value<TReadStore>::Type TReadStoreElement;

    SeqFileIn seqFile;
    if (!open(seqFile, toCString(fileName)))
        return false;

    readRecords(store.readNameStore, store.readSeqStore, seqFile);
    resize(store.readStore, length(store.readSeqStore), TReadStoreElement());
    return true;
}


template <typename TSpec, typename TConfig, typename TFileName>
bool loadReads(FragmentStore<TSpec, TConfig> & store, TFileName & fileNameL, TFileName & fileNameR)
{
    SEQAN_ASSERT_EQ(length(store.readStore), length(store.readSeqStore));
    SEQAN_ASSERT_EQ(length(store.readStore), length(store.readNameStore));

    typedef typename FragmentStore<TSpec, TConfig>::TReadSeq TReadSeq;

    SeqFileIn seqFileL, seqFileR;
    if (!open(seqFileL, toCString(fileNameL)) || !open(seqFileR, toCString(fileNameR)))
        return false;

    // Read in sequences
    TReadSeq seq[2];
    CharString meta[2];

    while (!atEnd(seqFileL) && !atEnd(seqFileR))
    {
        readRecord(meta[0], seq[0], seqFileL);
        readRecord(meta[1], seq[1], seqFileR);
        appendMatePair(store, seq[0], seq[1], meta[0], meta[1]);
    }
    return true;
}

}  // namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
