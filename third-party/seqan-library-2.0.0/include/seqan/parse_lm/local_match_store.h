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

// TODO(holtgrew): Allow the storage of additional fields.
/* This could look as follows:

   The values of the additional fields are stored as strings.  Either we store name/value pairs or, more efficiently,
   store a stringset for each match, with the same number of entries.  An additional String/Map maps from field name to
   index in this string sets.
 */

#ifndef INCLUDE_SEQAN_PARSE_LM_LOCAL_MATCH_STORE_H_
#define INCLUDE_SEQAN_PARSE_LM_LOCAL_MATCH_STORE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @class LocalMatch
 * @headerfile <seqan/parse_lm.h>
 * @brief Stores information about local matches.
 *
 * @signature template <typename TId, typename TPosition>
 *            class LocalMatch;
 *
 * @tparam TId       Type to use for subject/query id.
 * @tparam TPosition Type to use for storing positions.
 *
 * Matches on the reverse-complement are encoded by the begin position being greater than the end position.
 *
 * Sequence names are not stored in LocalMatch objects but in the @link LocalMatchStore @endlink they belong to.
 *
 * @see LocalMatchStore
 *
 *
 * @var TId LocalMatch::queryId
 * @brief The id of the query.
 *
 * @var TPosition LocalMatch::queryBeginPos
 * @brief Begin position of local match in the query.
 *
 * @var TPosition LocalMatch::queryEndPos
 * @brief End position of local match in the query.
 *
 * @var TPosition LocalMatch::subjectBeginPos
 * @brief Begin position of local match in the subject.
 *
 * @var TPosition LocalMatch::subjectEndPos
 * @brief End position of local match in the subject.
 *
 * @var TId LocalMatch::subjectId
 * @brief The id of the subject.
 */

template <typename TId, typename TPosition>
class LocalMatch
{
public:
    TId id;
    TId subjectId;
    TPosition subjectBeginPos;
    TPosition subjectEndPos;
    TId queryId;
    TPosition queryBeginPos;
    TPosition queryEndPos;

    LocalMatch() :
            id(MaxValue<TId>::VALUE),
            subjectId(MaxValue<TId>::VALUE),
            subjectBeginPos(MaxValue<TPosition>::VALUE),
            subjectEndPos(MaxValue<TPosition>::VALUE),
            queryId(MaxValue<TId>::VALUE),
            queryBeginPos(MaxValue<TPosition>::VALUE),
            queryEndPos(MaxValue<TPosition>::VALUE)
    {}

    LocalMatch(TId id_,
               TId subjectId_,
               TPosition subjectBeginPos_,
               TPosition subjectEndPos_,
               TId queryId_,
               TPosition queryBeginPos_,
               TPosition queryEndPos_) :
            id(id_),
            subjectId(subjectId_),
            subjectBeginPos(subjectBeginPos_),
            subjectEndPos(subjectEndPos_),
            queryId(queryId_),
            queryBeginPos(queryBeginPos_),
            queryEndPos(queryEndPos_)
    {}

    bool operator==(LocalMatch const & other) const
    {
        return id == other.id && subjectId == other.subjectId && subjectBeginPos == other.subjectBeginPos &&
                subjectEndPos == other.subjectEndPos && queryId == other.queryId &&
                queryBeginPos == other.queryBeginPos && queryEndPos == other.queryEndPos;
    }
};

/*!
 * @class LocalMatchStoreConfig
 * @headerfile <seqan/parse_lm.h>
 * @brief Stores information about local matches.
 *
 * @signature template <typename TSpec>
 *            struct LocalMatchStoreConfig;
 *
 * @tparam TSpec Specializing type.
 *
 *
 * @typedef LocalMatchStoreConfig::TId;
 * @brief Type to use for ids.
 * @signature typedef unsigned LocalMatchStoreConfig::TId;
 *
 * @typedef LocalMatchStoreConfig::TPosition;
 * @brief The type to use for positions.
 * @signature typedef (...) LocalMatchStoreConfig::TPosition;
 */

template <typename TSpec>
struct LocalMatchStoreConfig
{
    typedef unsigned TId;
    typedef typename Position<Dna5String>::Type TPosition;
};

/*!
 * @class LocalMatchStore
 * @headerfile <seqan/parse_lm.h>
 * @brief Stores information about local matches.
 *
 * @signature template <typename TSpec = void, TConfig = LocalMatchStoreConfig<TSpec> >
 *            struct LocalMatchStore;
 *
 * @tparam TSpec   Specialization tag.
 * @tparam TConfig Configuration class.
 *
 * The LocalMatchStore provides information about matches.  Similar to the @link FragmentStore @endlink, the
 * information is split into multiple sub stores.  Each sub store stores a part of the information.
 *
 * The LocalMatchStore#matchStore stores the information about a match.  The LocalMatchStore#sequenceNameStore stores
 * the sequence names.  These both sub stores are "core stores", they are filled with consistent information, i.e. for
 * each sequence id in the matchStore, there has to be a valid entry in the sequenceNameStore.
 *
 * The LocalMatchStore#cigarStore optionally stores CIGAR strings for the matches.  Its entries are referenced by
 * <tt>matchStore[i].id</tt>.
 *
 * @section Examples
 *
 * Read Lastz matches from a link RecordReader and then print them to stdout.
 *
 * @code{.cpp}
 * // Read local alignments from record reader.  Note that in real-world code,
 * // you should have error handling.
 * LocalMatchStore<> lmStore;
 * while (!atEnd(recordReader))
 *     readRecord(lmStore, recordReader, StellarGff());
 *
 * // Print local alignment information to stdout.
 * std::cout << "# Reverse strandness is indicated by end < begin\n"
 *           << "#db\tdb_beg\tdb_end\tquery\tq_beg\tq_end\n";
 * for (unsigned i = 0; i < length(lmStore.matchStore); ++i)
 *     std::cout << lmStore.sequenceNameStore[lmStore.matchStore[i].queryId] << "\t"
 *               << lmStore.matchStore[i].queryBeginPos << "\t"
 *               << lmStore.matchStore[i].queryEndPos << "\t"
 *               << lmStore.sequenceNameStore[lmStore.matchStore[i].subjectId] << "\t"
 *               << lmStore.matchStore[i].subjectBeginPos << "\t"
 *               << lmStore.matchStore[i].subjectEndPos << "\n";
 * @endcode
 * @see LocalMatch
 *
 * @var TMatchStore LocalMatchStore::matchStore
 * @brief String storing the LocalMatch local matches.
 *
 * @var TStringSet LocalMatchStore::sequenceNameStore
 * @brief StringSet storing the sequence names.
 *
 * @var TCigarString LocalMatchStore::cigarStore
 * @brief String storing the CIGAR strings.
 */

/*!
 * @fn LocalMatchStore#readRecord
 *
 * @headerfile <seqan/parse_lm.h>
 *
 * @brief Read Lastz "general" format record.
 *
 * @signature int readRecord(store, reader, tag);
 *
 * @param[in,out] store  LocalMatchStore object to read into.
 * @param[in,out] stream SinglePassRecordReader to read from.
 * @param[in]     tag    The tag for selecting the format, one of BlastnTabular, LastzGeneral, and StellarGff.
 *
 * @return int 0 on success, non-0 on errors and EOF
 */

template <typename TSpec=void, typename TConfig=LocalMatchStoreConfig<TSpec> >
class LocalMatchStore
{
public:
    // ----------------------------------------------------------------------------
    // Typedefs
    // ----------------------------------------------------------------------------

    typedef typename TConfig::TId TId;
    typedef typename TConfig::TPosition TPosition;

    typedef LocalMatch<TId, TPosition> TLocalMatch;
    typedef String<TLocalMatch> TMatchStore;

    typedef String<CigarElement<> > TCigarString;
    typedef String<TCigarString> TCigarStore;

    typedef StringSet<CharString> TNameStore;
    typedef NameStoreCache<TNameStore> TNameStoreCache_;

    // ----------------------------------------------------------------------------
    // Member Variables
    // ----------------------------------------------------------------------------

    TNameStore       sequenceNameStore;
    TNameStoreCache_ _sequenceNameStoreCache;

    TMatchStore       matchStore;
    TCigarStore       cigarStore;

    LocalMatchStore() :
            _sequenceNameStoreCache(sequenceNameStore)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function registerSequenceName
// ----------------------------------------------------------------------------

template <typename TLocalMatchStore>
inline void
registerSequenceName(TLocalMatchStore & store,
                     CharString const & sequenceName)
{
    unsigned id = 0;
    if (!getIdByName(store.sequenceNameStore, sequenceName, id, store._sequenceNameStoreCache))
    {
        id = length(store.sequenceNameStore);
        appendName(store.sequenceNameStore, sequenceName, store._sequenceNameStoreCache);
    }
}

// ----------------------------------------------------------------------------
// Function appendLocalMatch
// ----------------------------------------------------------------------------

/*!
 * @fn LocalMatchStore#appendLocalMatch
 * @brief Append a new local match to a @link LocalMatchStore @endlink
 *
 * @signature void appendLocalMatchStore(store, subjectId, subjectBeginPos, subjectEndPos, queryId, queryBeginPos, queryEndPos);
 * @signature void appendLocalMatchStore(store, subjectName, subjectBeginPos, subjectEndPos, queryName, queryBeginPos, queryEndPos, cigarStringBuffer);
 *
 * @param[in,out] store             The LocalMatchStore to add the local match to.
 * @param[in]     subjectId         Numeric subject sequence identifier, @link IntegerConcept @endlink.
 * @param[in]     subjectName       The textual name of the query sequence, @link CharString @endlink.
 * @param[in]     subjectBegin      The begin position of the match in the subject, @link IntegerConcept @endlink.
 * @param[in]     subjectEnd        The end position of the match in the subject, @link IntegerConcept @endlink.
 * @param[in]     queryId           Numeric query sequence identifier, @link IntegerConcept @endlink.
 * @param[in]     queryName         The textual name of the query sequence, @link CharString @endlink.
 * @param[in]     queryBegin        The begin position of the match in the query, @link IntegerConcept @endlink.
 * @param[in]     queryEnd          The end position of the match in the query, @link IntegerConcept @endlink.
 * @param[in]     cigarStringBuffer Buffer with the cigar string of the local alignment,  @link CharString @endlink.
 *
 * Matches on the reverse-complement are encoded by the begin position being greater than the begin position.
 */

template <typename TLocalMatchStore, typename TPosition>
inline void
appendLocalMatch(TLocalMatchStore & store,
                 unsigned const & subjectId,
                 TPosition const subjectBeginPos,
                 TPosition const subjectEndPos,
                 unsigned const & queryId,
                 TPosition const queryBeginPos,
                 TPosition const queryEndPos)
{
    typedef typename TLocalMatchStore::TLocalMatch TLocalMatch;

    SEQAN_ASSERT_LT(subjectId, length(store.sequenceNameStore));
    SEQAN_ASSERT_LT(queryId, length(store.sequenceNameStore));

    TLocalMatch localMatch(length(store.matchStore),
                           subjectId,
                           subjectBeginPos,
                           subjectEndPos,
                           queryId,
                           queryBeginPos,
                           queryEndPos);
    appendValue(store.matchStore, localMatch);
}

template <typename TLocalMatchStore, typename TPosition>
inline void
appendLocalMatch(TLocalMatchStore & store,
                 CharString const & subjectName,
                 TPosition const subjectBeginPos,
                 TPosition const subjectEndPos,
                 CharString const & queryName,
                 TPosition const queryBeginPos,
                 TPosition const queryEndPos)
{
    typedef typename TLocalMatchStore::TId         TId;

    // Get id for subject and query sequences;  Insert sequences into name stores/caches if not already there.
    TId subjectId = 0;
    if (!getIdByName(store.sequenceNameStore, subjectName, subjectId, store._sequenceNameStoreCache))
    {
        subjectId = length(store.sequenceNameStore);
        appendName(store.sequenceNameStore, subjectName, store._sequenceNameStoreCache);
    }
    TId queryId = 0;
    if (!getIdByName(store.sequenceNameStore, queryName, queryId, store._sequenceNameStoreCache))
    {
        queryId = length(store.sequenceNameStore);
        appendName(store.sequenceNameStore, queryName, store._sequenceNameStoreCache);
    }

    appendLocalMatch(store, subjectId, subjectBeginPos, subjectEndPos, queryId, queryBeginPos, queryEndPos);
}

template <typename TLocalMatchStore, typename TPosition>
inline void
appendLocalMatch(TLocalMatchStore & store,
                 CharString const & subjectName,
                 TPosition const subjectBeginPos,
                 TPosition const subjectEndPos,
                 CharString const & queryName,
                 TPosition const queryBeginPos,
                 TPosition const queryEndPos,
                 CharString const & cigarStringBuffer)
{
    // Append local match.
    appendLocalMatch(store, subjectName, subjectBeginPos, subjectEndPos, queryName, queryBeginPos, queryEndPos);
    // Make space for CIGAR string.
    resize(store.cigarStore, back(store.matchStore).id + 1);
    // TODO(holtgrew): Something can go wrong when parsing CIGAR string, need return value?
    // Parse out cigar string.
    /*
    typedef Stream<CharArray<char const *> > TCharArrayStream;
    TCharArrayStream cigarStream(&cigarStringBuffer[0], &cigarStringBuffer[0] + length(cigarStringBuffer));
    RecordReader<TCharArrayStream, SinglePass<> > recordReader(cigarStream);
    */
    DirectionIterator<CharString const, Input>::Type iter = begin(cigarStringBuffer);
    CharString numBuf;
    while (!atEnd(iter))
    {
        // Get number string into buffer.
        clear(numBuf);
        readUntil(numBuf, iter, NotFunctor<IsDigit>());
        unsigned num = 0;
        lexicalCast<unsigned>(num, numBuf);

        // Read operation char and advance record reader.
        char op;
        readOne(op, iter);
        // Append CIGAR element to CIGAR string.
        appendValue(back(store.cigarStore), CigarElement<>(op, num));
    }
}

}  // namespace seqan

#endif  // INCLUDE_SEQAN_PARSE_LM_LOCAL_MATCH_STORE_H_
