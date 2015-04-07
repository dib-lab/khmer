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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// BamHeaderRecord class, supporting types and functions for accessing tags
// in headers.
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_BAM_IO_BAM_HEADER_RECORD_H_
#define CORE_INCLUDE_SEQAN_BAM_IO_BAM_HEADER_RECORD_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Enum.BamHeaderRecordType
..cat:BAM I/O
..summary:Enumeration for the header record type.
..signature:BamHeaderRecordType
..value.BAM_HEADER_FIRST:@Class.BamHeaderRecord@ is of type $@HD$
..value.BAM_HEADER_REFERENCE:@Class.BamHeaderRecord@ is of type $@SQ$
..value.BAM_HEADER_READ_GROUP:@Class.BamHeaderRecord@ is of type $@RG$
..value.BAM_HEADER_PROGRAM:@Class.BamHeaderRecord@ is of type $@PG$
..value.BAM_HEADER_COMMENT:@Class.BamHeaderRecord@ is of type $@CO$
..include:seqan/bam_io.h
*/

enum BamHeaderRecordType
{
    BAM_HEADER_FIRST       = 0,
    BAM_HEADER_REFERENCE   = 1,
    BAM_HEADER_READ_GROUP  = 2,
    BAM_HEADER_PROGRAM     = 3,
    BAM_HEADER_COMMENT     = 4
};

/**
.Enum.BamSortOrder
..cat:BAM I/O
..summary:Enumeration for the header record type.
..signature:BamSortOrder
..value.BAM_SORT_UNKNOWN:BAM file sort order is unknown.
..value.BAM_SORT_UNSORTED:BAM file is unsorted.
..value.BAM_SORT_QUERYNAME:BAM file is sorted by query name.
..value.BAM_SORT_COORDINATE:BAM file is sorted by coordinates.
..include:seqan/bam_io.h
*/

enum BamSortOrder
{
    BAM_SORT_UNKNOWN    = 0,
    BAM_SORT_UNSORTED   = 1,
    BAM_SORT_QUERYNAME  = 2,
    BAM_SORT_COORDINATE = 3
};

/**
.Class.BamHeaderRecord
..cat:BAM I/O
..summary:Represents a header entry in a SAM file or the header section of the BAM header.
..signature:BamHeaderRecord
..remarks:Comment records are stored with one tag where the key is empty and the value is the comment.
..include:seqan/bam_io.h

.Memfunc.BamHeaderRecord#BamHeaderRecord
..class:Class.BamHeaderRecord
..signature:BamHeaderRecord()
..summary:Constructor.
..remarks:Only the default constructor is provided.

.Typedef.BamHeaderRecord#TTagName
..class:Class.BamHeaderRecord
..summary:Type of the tag keys.

.Typedef.BamHeaderRecord#TTagValue
..class:Class.BamHeaderRecord
..summary:Type of the tag values.

.Typedef.BamHeaderRecord#TTag
..class:Class.BamHeaderRecord
..summary:@Class.Pair@ to use for storing tags.

.Typedef.BamHeaderRecord#TTags
..class:Class.BamHeaderRecord
..summary:Type of the string of tag @Class.Pair|Pairs@.

.Memvar.BamHeaderRecord#type
..summary:Type of the record.
..class:Class.BamHeaderRecord
..type:Enum.BamHeaderRecordType

.Memvar.BamHeaderRecord#tags
..summary:The header record's tags.
..class:Class.BamHeaderRecord
..type:Typedef.BamHeaderRecord#TTags
*/

class BamHeaderRecord
{
public:
    typedef CharString TTagName;
    typedef CharString TTagValue;
    typedef Pair<TTagName, TTagValue> TTag;
    typedef String<TTag> TTags;

    BamHeaderRecordType type;
    String<Pair<TTagName, TTagValue> > tags;

    BamHeaderRecord() {}
};

/**
.Class.BamHeader
..cat:BAM I/O
..summary:Stores the information of the BAM header.
..signature:BamHeader
..see:Class.BamHeaderRecord
..include:seqan/bam_io.h

.Memfunc.BamHeader#BamHeader
..class:Class.BamHeader
..signature:BamHeader()
..summary:Constructor.
..remarks:Only the default constructor is provided.

.Memvar.BamHeader#sequenceInfos
..class:Class.BamHeader
..summary:String of $(seqid, length)$ with reference name / length information.
..type:nolink:$String<Pair<CharString, __int32> >$

.Memvar.BamHeader#records
..class:Class.BamHeader
..summary:String of @Class.BamHeaderRecord|BamHeaderRecords@.
..type:nolink:$String<BamHeaderRecord>$
*/

class BamHeader
{
public:
    typedef Pair<CharString, __int32> TSequenceInfo;
    
    String<Pair<CharString, __int32> > sequenceInfos;
    String<BamHeaderRecord> records;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

///.Function.clear.param.object.type:Class.BamHeaderRecord

inline void
clear(BamHeaderRecord & record)
{
    clear(record.tags);
}

// ----------------------------------------------------------------------------
// Function findTagKey()
// ----------------------------------------------------------------------------

/**
.Function.BamHeaderRecord#findTagKey
..cat:BAM I/O
..summary:Find a tag's key of a @Class.BamHeaderRecord@.
..signature:findTagKey(idx, key, record)
..param.idx:The index of the found key is stored here.
...type:nolink:$unsigned$
..param.key:The name of the tag key whose position is to be stored in $idx$.
...type:Shortcut.CharString
..param.record:The record to query.
...type:Class.BamHeaderRecord
..returns:$bool$, indicating whether the key could be found.
..include:seqan/bam_io.h
..example.code:
unsigned myIdx = 0;
bool keyFound = findTagKey(myIdx, "SN", record);
*/

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

/**
.Function.BamHeaderRecord#getTagValue
..cat:BAM I/O
..summary:Return tag value from a @Class.BamHeaderRecord@ or @Class.BamTagsDict@.
..signature:getTagValue(tagValue, idx, record)
..signature:getTagValue(tagValue, key, record)
..param.tagValue:The tag's value is stored here.
...type:Shortcut.CharString
..param.idx:The index of the tag whose value is to be retrieved.
...type:nolink:$unsigned$
..param.key:The name of tag whose value is to be retrieved.
...type:Shortcut.CharString
..param.record:The record to query.
...type:Class.BamHeaderRecord
..returns:$bool$, indicating whether the value could be retrieved, always $true$ if $idx$ is given.
..include:seqan/bam_io.h
..example.code:
CharString tagValue;
bool keyFound = getTagValue(tagValue, "SN", record);
..see:Function.BamHeaderRecord#findTagKey
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

// TODO(holtgrew): Parameter order!

/**
.Function.BamHeaderRecord#setTagValue
..cat:BAM I/O
..summary:Set tag value of a @Class.BamHeaderRecord@.
..signature:setTagValue(idx, tagValue, record)
..signature:setTagValue(key, tagValue, record)
..param.idx:The index of the tag whose value should be set.
...type:nolink:$unsigned$
..param.key:The name of tag whose value should be set.
...type:Shortcut.CharString
..param.tagValue:The new tag value.
...type:Shortcut.CharString
..param.record:The record to query.
...type:Class.BamHeaderRecord
..include:seqan/bam_io.h
..example.code:
setTagValue("SN", "chr1", record);
..see:Function.BamHeaderRecord#findTagKey
*/

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
    for (recordIdx = startIdx; recordIdx < length(header.records); ++recordIdx)
    {
        if (header.records[recordIdx].type == recordType)
            return true;
    }
    return false;
}

inline bool
searchRecord(unsigned & recordIdx,
             BamHeader const & header,
             BamHeaderRecordType recordType)
{
    return searchRecord(recordIdx, header, recordType, 0);
}

inline BamSortOrder
getSortOrder(BamHeader const & header)
{
    CharString soString;    
    for (unsigned recIdx = 0; searchRecord(recIdx, header, BAM_HEADER_FIRST, recIdx); ++recIdx)
    {
        if (getTagValue(soString, "SO", header.records[recIdx]))
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
    for (unsigned recIdx = 0; searchRecord(recIdx, header, BAM_HEADER_FIRST, recIdx); ++recIdx)
    {
        unsigned idx = 0;
        if (findTagKey(idx, "SO", header.records[recIdx]))
        {
            CharString soString;
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
            setTagValue(idx, soString, header.records[recIdx]);
        }
    }
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

///.Function.clear.param.object.type:Class.BamHeader

inline void
clear(BamHeader & header)
{
    clear(header.sequenceInfos);
    clear(header.records);
}


}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BAM_IO_BAM_HEADER_RECORD_H_
