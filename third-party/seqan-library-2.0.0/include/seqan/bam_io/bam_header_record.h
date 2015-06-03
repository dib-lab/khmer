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
// BamHeaderRecord class, supporting types and functions for accessing tags
// in headers.
// ==========================================================================

#ifndef INCLUDE_SEQAN_BAM_IO_BAM_HEADER_RECORD_H_
#define INCLUDE_SEQAN_BAM_IO_BAM_HEADER_RECORD_H_

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Enum BamHeaderRecordType
// ----------------------------------------------------------------------------

/*!
 * @enum BamHeaderRecordType
 * @headerfile <seqan/bam_io.h>
 * @brief Enumeration for the header record type.
 *
 * @signature enum BamHeaderRecordType;
 *
 * @val BamHeaderRecordType BAM_HEADER_FIRST = 0;
 * @brief Is the first header (HD).
 *
 * @val BamHeaderRecordType BAM_HEADER_REFERENCE = 1;
 * @brief Is a reference (SQ) header.
 *
 * @val BamHeaderRecordType BAM_HEADER_READ_GROUP = 2;
 * @brief Is a read group (RG) header.
 *
 * @val BamHeaderRecordType BAM_HEADER_PROGRAM = 3;
 * @brief Is a program (PG) header.
 *
 * @val BamHeaderRecordType BAM_HEADER_COMMENT = 4;
 * @brief Is a comment (CO) header.
 */

enum BamHeaderRecordType
{
    BAM_HEADER_FIRST       = 0,
    BAM_HEADER_REFERENCE   = 1,
    BAM_HEADER_READ_GROUP  = 2,
    BAM_HEADER_PROGRAM     = 3,
    BAM_HEADER_COMMENT     = 4
};

// ----------------------------------------------------------------------------
// Enum BamSortOrder
// ----------------------------------------------------------------------------

/*!
 * @enum BamSortOrder
 * @headerfile <seqan/bam_io.h>
 * @brief Enumeration for the header record type.
 *
 * @signature enum BamSortOrder;
 *
 * @val BamSortOrder BAM_SORT_UNKNOWN = 0;
 * @brief BAM file sort order is unknown.
 *
 * @val BamSortOrder BAM_SORT_UNSORTED = 1;
 * @brief BAM file is unsorted.
 *
 * @val BamSortOrder BAM_SORT_QUERYNAME = 2;
 * @brief BAM file is sorted by query name;
 *
 * @val BamSortOrder BAM_SORT_COORDINATE = 3;
 * @brief BAM file is sorted by coordinate.
 */

enum BamSortOrder
{
    BAM_SORT_UNKNOWN    = 0,
    BAM_SORT_UNSORTED   = 1,
    BAM_SORT_QUERYNAME  = 2,
    BAM_SORT_COORDINATE = 3
};

// ----------------------------------------------------------------------------
// Class BamHeaderRecord
// ----------------------------------------------------------------------------

/*!
 * @class BamHeaderRecord
 * @headerfile <seqan/bam_io.h>
 * @brief Represents a header entry in a SAM file or the header section of the BAM header.
 *
 * @signature class BamHeaderRecord;
 *
 * @section Remarks
 *
 * Comment records are stored with one tag where the key is empty and the value is the comment.
 */

/*!
 * @fn BamHeaderRecord::BamHeaderRecord
 * @brief Constructor.
 * @signature BamHeaderRecord::BamRecord();
 *
 * Only the default constructor is provided.
 */

/*!
 * @typedef BamHeaderRecord::TTagName
 * @brief Type of the tag keys (@link CharString @endlink).
 * @signature BamHeaderRecord::TTagName;
 *
 * @typedef BamHeaderRecord::TTagValue
 * @brief Type of the tag values (@link CharString @endlink).
 * @signature BamHeaderRecord::TTagValue;
 *
 * @typedef BamHeaderRecord::TTag
 * @brief Type of the tag keys (@link Pair @endlink of <tt>TTagName</tt> and <tt>TTagValue</tt>).
 * @signature BamHeaderRecord::TTag;
 *
 * @typedef BamHeaderRecord::TTags
 * @brief Type of the tags string (@link AllocString @endlink of <tt>TTag</tt>).
 * @signature BamHeaderRecord::TTags;
 *
 * @var BamHeaderRecordType BamHeaderRecord::type;
 * @brief Type of the record.
 *
 * @var TRecordString BamHeaderRecord::tags;
 * @brief The header record's tags, of type @link BamHeaderRecord::TTags @endlink.
 */

class BamHeaderRecord
{
public:
    typedef CharString                  TTagName;
    typedef CharString                  TTagValue;
    typedef Pair<TTagName, TTagValue>   TTag;
    typedef String<TTag>                TTags;

    BamHeaderRecordType type;
    TTags tags;

    BamHeaderRecord() : type(BAM_HEADER_FIRST) {}
};

// ----------------------------------------------------------------------------
// Class BamHeader
// ----------------------------------------------------------------------------

/*!
 * @class BamHeader
 * @headerfile <seqan/bam_io.h>
 * @implements FormattedFileHeaderConcept
 * @signature typedef String<BamHeaderRecord> BamHeader;
 * @brief Represent the information of the BAM header.
 *
 * @see BamFileIn
 * @see BamFileOut
 */

typedef String<BamHeaderRecord> BamHeader;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function std::swap()
// ----------------------------------------------------------------------------

inline void
swap(BamHeaderRecord &a, BamHeaderRecord &b)
{
    std::swap(a.type, b.type);
    swap(a.tags, b.tags);
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

/*!
 * @fn BamHeaderRecord::clear
 * @brief Clear BamHeaderRecord object.
 *
 * @signature void clear(record);
 *
 * @param[in,out] record The record to clear.
 */

inline void
clear(BamHeaderRecord & record)
{
    clear(record.tags);
}

// ----------------------------------------------------------------------------
// Function findTagKey()
// ----------------------------------------------------------------------------

template <typename TKeyName>
inline bool
findTagKey(unsigned & idx, TKeyName const & key, BamHeaderRecord const & record)
{
    typedef BamHeaderRecord::TTags const TTags;
    typedef Iterator<TTags, Standard>::Type TIterator;

    idx = 0;
    TIterator itEnd = end(record.tags, Standard());
    for (TIterator it = begin(record.tags, Standard()); it != itEnd; goNext(it), ++idx)
        if (it->i1 == key)
            return true;

    return false;
}

// ----------------------------------------------------------------------------
// Function getTagValue()
// ----------------------------------------------------------------------------

/*!
 * @fn BamHeaderRecord::getTagValue
 * @brief Return tag value from a @link BamHeaderRecord @endlink.
 *
 * @signature bool getTagValue(tagValue, idx, record);
 * @signature bool getTagValue(tagValue, key, record);
 *
 * @param[out] tagValue The @link CharString @endlink to write the tag value to.
 * @param[in]  idx      An integer with the index of the tag in the header record.
 * @param[in]  key      A two-letter sequence with the key of the tag in the header record.
 *
 * @return bool true in case the value could be retrieved, false otherwise.
 *
 * @section Example
 *
 * @code{.cpp}
 * CharString tagValue;
 * bool keyFound = getTagValue(tagValue, "SN", record);
 * @endcode
 */

template <typename TId>
SEQAN_FUNC_ENABLE_IF(
    IsInteger<TId>,
    bool)
inline getTagValue(CharString & value, TId idx, BamHeaderRecord const & record)
{
    if ((unsigned)idx >= length(record.tags))
        return false;
    value = record.tags[idx].i2;
    return true;
}

template <typename TKeyName>
SEQAN_FUNC_ENABLE_IF(
    IsSequence<TKeyName>,
    bool)
inline getTagValue(CharString & value, TKeyName const & key, BamHeaderRecord const & record)
{
    unsigned idx = 0;
    if (!findTagKey(idx, key, record))
        return false;
    return getTagValue(value, idx, record);
}

/*!
 * @fn BamHeaderRecord::setTagValue
 * @brief Sets the value of a tag in a BamHeaderRecord.
 *
 * @signature void setTagValue(idx, value, record);
 * @signature void setTagValue(key, value, record);
 *
 * @param[in,out] idx    The index of the tag in the header record to set the value for.
 * @param[in]     key    The name of the tag (two-letter sequence) to set.
 * @param[in]     record The header record to set the value for.
 */

// TODO(holtgrew): Parameter order!

template <typename TId>
SEQAN_FUNC_ENABLE_IF(
    IsInteger<TId>,
    void)
inline setTagValue(TId idx, CharString const & value, BamHeaderRecord & record)
{
    if (idx >= length(record.tags))
        return;
    record.tags[idx].i2 = value;
}

template <typename TKeyName>
SEQAN_FUNC_ENABLE_IF(
    IsSequence<TKeyName>,
    void)
inline setTagValue(TKeyName const & key, CharString const & value, BamHeaderRecord & record)
{
    unsigned idx = 0;
    if (!findTagKey(idx, key, record))
    {
        idx = length(record.tags);
        resize(record.tags, idx + 1);
        record.tags[idx].i1 = key;
    }
    setTagValue(idx, value, record);
}

// ----------------------------------------------------------------------------
// Function getSortOrder()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Document me! Test me!


inline bool
searchRecord(unsigned & recordIdx,
             BamHeader const & header,
             BamHeaderRecordType recordType,
             unsigned startIdx)
{
    for (recordIdx = startIdx; recordIdx < length(header); ++recordIdx)
        if (header[recordIdx].type == recordType)
            return true;
    return false;
}

inline bool
searchRecord(unsigned & recordIdx,
             BamHeader const & header,
             BamHeaderRecordType recordType)
{
    return searchRecord(recordIdx, header, recordType, 0);
}


struct BamHeaderRecordTypeLess
{
    bool operator() (BamHeaderRecord const &a, BamHeaderRecord const &b) const
    {
        return a.type < b.type;
    }
};

struct BamHeaderRecordEqual
{
    bool operator() (BamHeaderRecord const &a, BamHeaderRecord const &b) const
    {
        return a.type == b.type && a.tags == b.tags;
    }
};


inline void
removeDuplicates(BamHeader & header)
{
    BamHeaderRecordTypeLess less;
    BamHeaderRecordEqual pred;

    std::stable_sort(begin(header, Standard()), end(header, Standard()), less);

    for (size_t uniqueBegin = 0, uniqueEnd = 1; uniqueEnd < length(header);)
    {
        if (less(header[uniqueBegin], header[uniqueEnd]))
            uniqueBegin = uniqueEnd;

        size_t j;
        for (j = uniqueBegin; j < uniqueEnd; ++j)
            if (pred(header[j], header[uniqueEnd]))
            {
                erase(header, uniqueEnd);
                break;
            }

        if (j == uniqueEnd)
            ++uniqueEnd;
    }
}

inline BamSortOrder
getSortOrder(BamHeader const & header)
{
    CharString soString;
    for (unsigned recIdx = 0; searchRecord(recIdx, header, BAM_HEADER_FIRST, recIdx); ++recIdx)
    {
        if (getTagValue(soString, "SO", header[recIdx]))
        {
            if (soString == "unsorted")
                return BAM_SORT_UNSORTED;
            else if (soString == "")
                return BAM_SORT_QUERYNAME;
            else if (soString == "")
                return BAM_SORT_COORDINATE;
            else
                return BAM_SORT_UNKNOWN;
        }
    }
    return BAM_SORT_UNKNOWN;
}

inline void
setSortOrder(BamHeader & header, BamSortOrder sortOrder)
{
    char const * soString;
    switch (sortOrder)
    {
        case BAM_SORT_UNSORTED:
            soString = "unsorted";
            break;

        case BAM_SORT_QUERYNAME:
            soString = "queryname";
            break;

        case BAM_SORT_COORDINATE:
            soString = "coordinate";
            break;

        default:
            soString = "unknown";
    }

    bool notFound = true;
    for (unsigned recIdx = 0; searchRecord(recIdx, header, BAM_HEADER_FIRST, recIdx); ++recIdx)
    {
        unsigned idx = 0;
        if (findTagKey(idx, "SO", header[recIdx]))
        {
            notFound = false;
            setTagValue(idx, soString, header[recIdx]);
        }
    }

    if (notFound)
    {
        unsigned recIdx = 0;
        if (searchRecord(recIdx, header, BAM_HEADER_FIRST))
        {
            setTagValue("SO", soString, header[recIdx]);
        }
        else
        {
            BamHeaderRecord rec;
            rec.type = BAM_HEADER_FIRST;
            setTagValue("VN", "1.4", rec);
            setTagValue("SO", soString, rec);
            insert(header, 0, rec);
        }
    }
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_BAM_IO_BAM_HEADER_RECORD_H_
