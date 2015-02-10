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
// The class BamAlignmentRecord, flag checking methods, flag constants.
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_BAM_IO_BAM_RECORD_H_
#define CORE_INCLUDE_SEQAN_BAM_IO_BAM_RECORD_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

class BamAlignmentRecord;
inline void clear(BamAlignmentRecord & record);

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Enum.BamFlags
..cat:BAM I/O
..signature:BamFlags
..summary:Shortcuts to the bitmask flags for BAM/SAM files.
..value.BAM_FLAG_MULTIPLE:$0x0001$ Template has multiple fragments in sequencing.
..value.BAM_FLAG_ALL_PROPER:$0x0002$ All fragments have been aligned properly.
..value.BAM_FLAG_UNMAPPED:$0x0004$ This fragment could not be mapped.
..value.BAM_FLAG_NEXT_UNMAPPED:$0x0008$ Next fragment is unmapped.
..value.BAM_FLAG_RC:$0x0010$ This fragment is reverse-complemented.
..value.BAM_FLAG_NEXT_RC:$0x0020$ Next fragment is reverse-complemented.
..value.BAM_FLAG_FIRST:$0x0040$ This fragment is the first one in its template.
..value.BAM_FLAG_LAST:$0x0080$ This fragment is the last one in its template.
..value.BAM_FLAG_SECONDARY:$0x0100$ This is a secondary alignment.
..value.BAM_FLAG_QC_NO_PASS:$0x0200$ Does not pass quality controls.
..value.BAM_FLAG_DUPLICATE:$0x0400$ PCR or optical duplicate.
..remarks:Also see the SAM standard on these flags for more explanation.
..include:seqan/bam_io.h
*/

enum BamFlags
{
    BAM_FLAG_MULTIPLE      = 0x0001,
    BAM_FLAG_ALL_PROPER    = 0x0002,
    BAM_FLAG_UNMAPPED      = 0x0004,
    BAM_FLAG_NEXT_UNMAPPED = 0x0008,
    BAM_FLAG_RC            = 0x0010,
    BAM_FLAG_NEXT_RC       = 0x0020,
    BAM_FLAG_FIRST         = 0x0040,
    BAM_FLAG_LAST          = 0x0080,
    BAM_FLAG_SECONDARY     = 0x0100,
    BAM_FLAG_QC_NO_PASS    = 0x0200,
    BAM_FLAG_DUPLICATE     = 0x0400
};

/**
.Class.BamAlignmentRecord
..cat:BAM I/O
..summary:Represent a record from a BAM/SAM file.
..remarks:While also used to represent SAM records, called $BamAlignmentRecord$ since the data directly reflects a BAM record (0-based positions, identify references by ids, not names, tags stored in BAM format.)
..include:seqan/bam_io.h
..see:Enum.BamFlags

.Memfunc.BamAlignmentRecord#BamAlignmentRecord
..class:Class.BamAlignmentRecord
..summary:Constructor.
..signature:BamAlignmentRecord()
..remarks:Only the default constructor is provided.

.Memvar.BamAlignmentRecord#INVALID_POS
..class:Class.BamAlignmentRecord
..summary:Static member with invalid/sentinel position value.
..type:nolink:$__uint32$

.Memvar.BamAlignmentRecord#INVALID_REFID
..class:Class.BamAlignmentRecord
..summary:Static member with invalid/sentinel reference id (-1 as in BAM/SAM).
..type:nolink:$__int32$

.Memvar.BamAlignmentRecord#INVALID_LEN
..class:Class.BamAlignmentRecord
..summary:Static member with invalid/sentinel position value.
..type:nolink:$__int32$

.Memvar.BamAlignmentRecord#qName
..class:Class.BamAlignmentRecord
..summary:The read/query name.
..type:Shortcut.CharString

.Memvar.BamAlignmentRecord#flag
..class:Class.BamAlignmentRecord
..summary:The flag of this mapping, see @Enum.BamFlags@ for flag constants and the $hasFlag*$ functions.
..type:nolink:$__uint16$

.Memvar.BamAlignmentRecord#rID
..class:Class.BamAlignmentRecord
..summary:ID of reference for this fragment mapping (0-based, $INVALID_REFID$ for '*').
..type:nolink:$__int32$

.Memvar.BamAlignmentRecord#beginPos
..class:Class.BamAlignmentRecord
..summary:The position of this fragment mapping (0-based, $INVALID_POS$ for '*').
..type:nolink:$__int32$

.Memvar.BamAlignmentRecord#mapQ
..class:Class.BamAlignmentRecord
..summary:The mapping quality (255 for '*').
..type:nolink:$__uint8$

.Memvar.BamAlignmentRecord#bin
..class:Class.BamAlignmentRecord
..summary:The bin of the alignment, automatically computed when writing BAM.
..type:nolink:$__uint16$

.Memvar.BamAlignmentRecord#cigar
..class:Class.BamAlignmentRecord
..summary:The CIGAR string as string of @Class.CigarElement@ objects (empty for '*').
..type:nolink:$String<CigarElement<> >$

.Memvar.BamAlignmentRecord#rNextId
..class:Class.BamAlignmentRecord
..summary:ID of reference for next fragment mapping (0-based, $INVALID_REFID$ for '*')
..type:nolink:$__int32$

.Memvar.BamAlignmentRecord#pNext
..class:Class.BamAlignmentRecord
..summary:Position of next fragment mapping (0-based, $INVALID_POS$ for '*')
..type:nolink:$__uint32$

.Memvar.BamAlignmentRecord#tLen
..class:Class.BamAlignmentRecord
..summary:The inferred template size ($INVALID_LEN$ for '*')
..type:nolink:$__int32$

.Memvar.BamAlignmentRecord#seq
..class:Class.BamAlignmentRecord
..summary:The sequence string (empty for '*').
..type:Shortcut.CharString

.Memvar.BamAlignmentRecord#qual
..class:Class.BamAlignmentRecord
..summary:String with Phred scores (as in SAM file, empty for '*').
..type:Shortcut.CharString

.Memvar.BamAlignmentRecord#tags
..class:Class.BamAlignmentRecord
..summary:Raw BAM tag string, use @Class.BamTagsDict@ for comfortable access.
..type:Shortcut.CharString
*/

class BamAlignmentRecord
{
public:
    static __int32 const INVALID_POS = 2147483647;  // TODO(holtgrew): Should be MaxValue<__int32>::VALUE, but that is not a constant expression :(
    static __int32 const INVALID_REFID = -1;  // TODO(holtgrew): Rename to ...REF_ID.
    static __int32 const INVALID_LEN = 2147483647;
    static __uint32 const INVALID_QID = 4294967295u;  // TODO(holtgrew): Undocumented as of yet.

    __uint32 _qID;  // TODO(holtgrew): Undocumented as of yet.
    CharString qName;
    __uint16 flag;
    __int32 rID;
    __int32 beginPos;
    __uint8 mapQ;
    __uint16 bin;
    String<CigarElement<> > cigar;
    __int32 rNextId;
    __int32 pNext;
    __int32 tLen;
    CharString seq;
    CharString qual;
    CharString tags;  // raw tags in BAM format

    BamAlignmentRecord() : _qID(MaxValue<unsigned>::VALUE) { clear(*this); }
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

///.Function.clear.param.object.type:Class.BamAlignmentRecord
///.Function.clear.class:Class.BamAlignmentRecord

inline void
clear(BamAlignmentRecord & record)
{
    clear(record.qName);
    record._qID = MaxValue<__uint32>::VALUE;
    record.rID = BamAlignmentRecord::INVALID_REFID;
    record.beginPos = BamAlignmentRecord::INVALID_POS;
    record.mapQ = 255;
    record.bin = 0;
    clear(record.cigar);
    record.rNextId = BamAlignmentRecord::INVALID_REFID;
    record.pNext = BamAlignmentRecord::INVALID_POS;
    record.tLen = BamAlignmentRecord::INVALID_LEN;
    clear(record.seq);
    clear(record.qual);
    clear(record.tags);
}

// ----------------------------------------------------------------------------
// Function hasFlagMultiple()
// ----------------------------------------------------------------------------

/**
.Function.hasFlagMultiple
..class:Class.BamAlignmentRecord
..cat:BAM I/O
..summary:Return true if a @Class.BamAlignmentRecord@ has the "multiple" flag set.
..signature:hasFlagMultiple(record)
..param.record:The record to query.
...type:Class.BamAlignmentRecord
..returns:$bool$, indicating the flag's status.
..include:seqan/bam_io.h
..see:Function.hasFlagAllProper
..see:Function.hasFlagUnmapped
..see:Function.hasFlagNextUnmapped
..see:Function.hasFlagRC
..see:Function.hasFlagNextRC
..see:Function.hasFlagFirst
..see:Function.hasFlagLast
..see:Function.hasFlagSecondary
..see:Function.hasFlagQCNoPass
..see:Function.hasFlagDuplicate
*/

inline bool
hasFlagMultiple(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_MULTIPLE) == BAM_FLAG_MULTIPLE;
}

// ----------------------------------------------------------------------------
// Function hasFlagAllProper()
// ----------------------------------------------------------------------------

/**
.Function.hasFlagAllProper
..class:Class.BamAlignmentRecord
..cat:BAM I/O
..summary:Return true if a @Class.BamAlignmentRecord@ has the "all properly aligned" flag set.
..signature:hasFlagAllProper(record)
..param.record:The record to query.
...type:Class.BamAlignmentRecord
..returns:$bool$, indicating the flag's status.
..include:seqan/bam_io.h
..see:Function.hasFlagMultiple
..see:Function.hasFlagUnmapped
..see:Function.hasFlagNextUnmapped
..see:Function.hasFlagRC
..see:Function.hasFlagNextRC
..see:Function.hasFlagFirst
..see:Function.hasFlagLast
..see:Function.hasFlagSecondary
..see:Function.hasFlagQCNoPass
..see:Function.hasFlagDuplicate
*/

inline bool
hasFlagAllProper(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_ALL_PROPER) == BAM_FLAG_ALL_PROPER;
}

// ----------------------------------------------------------------------------
// Function hasFlagUnmapped()
// ----------------------------------------------------------------------------

/**
.Function.hasFlagUnmapped
..class:Class.BamAlignmentRecord
..cat:BAM I/O
..summary:Return true if a @Class.BamAlignmentRecord@ has the "fragment unmapped" flag set.
..signature:hasFlagUnmapped(record)
..param.record:The record to query.
...type:Class.BamAlignmentRecord
..returns:$bool$, indicating the flag's status.
..include:seqan/bam_io.h
..see:Function.hasFlagMultiple
..see:Function.hasFlagAllProper
..see:Function.hasFlagNextUnmapped
..see:Function.hasFlagRC
..see:Function.hasFlagNextRC
..see:Function.hasFlagFirst
..see:Function.hasFlagLast
..see:Function.hasFlagSecondary
..see:Function.hasFlagQCNoPass
..see:Function.hasFlagDuplicate
*/

inline bool
hasFlagUnmapped(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_UNMAPPED) == BAM_FLAG_UNMAPPED;
}

// ----------------------------------------------------------------------------
// Function hasFlagNextUnmapped()
// ----------------------------------------------------------------------------

/**
.Function.hasFlagNextUnmapped
..class:Class.BamAlignmentRecord
..cat:BAM I/O
..summary:Return true if a @Class.BamAlignmentRecord@ has the "next fragment unmapped" flag set.
..signature:hasFlagNextUnmapped(record)
..param.record:The record to query.
...type:Class.BamAlignmentRecord
..returns:$bool$, indicating the flag's status.
..include:seqan/bam_io.h
..see:Function.hasFlagMultiple
..see:Function.hasFlagAllProper
..see:Function.hasFlagUnmapped
..see:Function.hasFlagRC
..see:Function.hasFlagNextRC
..see:Function.hasFlagFirst
..see:Function.hasFlagLast
..see:Function.hasFlagSecondary
..see:Function.hasFlagQCNoPass
..see:Function.hasFlagDuplicate
*/

inline bool
hasFlagNextUnmapped(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_NEXT_UNMAPPED) == BAM_FLAG_NEXT_UNMAPPED;
}

// ----------------------------------------------------------------------------
// Function hasFlagRC()
// ----------------------------------------------------------------------------

/**
.Function.hasFlagRC
..class:Class.BamAlignmentRecord
..cat:BAM I/O
..cat:BAM I/O
..summary:Return true if a @Class.BamAlignmentRecord@ has the "reverse-complemented" flag set.
..signature:hasFlagRC(record)
..param.record:The record to query.
...type:Class.BamAlignmentRecord
..returns:$bool$, indicating the flag's status.
..include:seqan/bam_io.h
..see:Function.hasFlagMultiple
..see:Function.hasFlagAllProper
..see:Function.hasFlagUnmapped
..see:Function.hasFlagNextUnmapped
..see:Function.hasFlagNextRC
..see:Function.hasFlagFirst
..see:Function.hasFlagLast
..see:Function.hasFlagSecondary
..see:Function.hasFlagQCNoPass
..see:Function.hasFlagDuplicate
*/
inline bool
hasFlagRC(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_RC) == BAM_FLAG_RC;
}

// ----------------------------------------------------------------------------
// Function hasFlagNextRC()
// ----------------------------------------------------------------------------

/**
.Function.hasFlagNextRC
..class:Class.BamAlignmentRecord
..cat:BAM I/O
..summary:Return true if a @Class.BamAlignmentRecord@ has the "next fragment reverse-complemented" flag set.
..signature:hasFlagNextRC(record)
..param.record:The record to query.
...type:Class.BamAlignmentRecord
..returns:$bool$, indicating the flag's status.
..include:seqan/bam_io.h
..see:Function.hasFlagMultiple
..see:Function.hasFlagAllProper
..see:Function.hasFlagUnmapped
..see:Function.hasFlagNextUnmapped
..see:Function.hasFlagRC
..see:Function.hasFlagFirst
..see:Function.hasFlagLast
..see:Function.hasFlagSecondary
..see:Function.hasFlagQCNoPass
..see:Function.hasFlagDuplicate
*/

inline bool
hasFlagNextRC(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_NEXT_RC) == BAM_FLAG_NEXT_RC;
}

// ----------------------------------------------------------------------------
// Function hasFlagFirst()
// ----------------------------------------------------------------------------

/**
.Function.hasFlagFirst
..class:Class.BamAlignmentRecord
..cat:BAM I/O
..summary:Return true if a @Class.BamAlignmentRecord@ has the "first fragment of template" flag set.
..signature:hasFlagFirst(record)
..param.record:The record to query.
...type:Class.BamAlignmentRecord
..returns:$bool$, indicating the flag's status.
..include:seqan/bam_io.h
..see:Function.hasFlagMultiple
..see:Function.hasFlagAllProper
..see:Function.hasFlagUnmapped
..see:Function.hasFlagNextUnmapped
..see:Function.hasFlagRC
..see:Function.hasFlagNextRC
..see:Function.hasFlagLast
..see:Function.hasFlagSecondary
..see:Function.hasFlagQCNoPass
..see:Function.hasFlagDuplicate
*/

inline bool
hasFlagFirst(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_FIRST) == BAM_FLAG_FIRST;
}

// ----------------------------------------------------------------------------
// Function hasFlagLast()
// ----------------------------------------------------------------------------

/**
.Function.hasFlagLast
..class:Class.BamAlignmentRecord
..cat:BAM I/O
..summary:Return true if a @Class.BamAlignmentRecord@ has the "last fragment of template" flag set.
..signature:hasFlagLast(record)
..param.record:The record to query.
...type:Class.BamAlignmentRecord
..returns:$bool$, indicating the flag's status.
..include:seqan/bam_io.h
..see:Function.hasFlagMultiple
..see:Function.hasFlagAllProper
..see:Function.hasFlagUnmapped
..see:Function.hasFlagNextUnmapped
..see:Function.hasFlagRC
..see:Function.hasFlagNextRC
..see:Function.hasFlagFirst
..see:Function.hasFlagSecondary
..see:Function.hasFlagQCNoPass
..see:Function.hasFlagDuplicate
*/

inline bool
hasFlagLast(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_LAST) == BAM_FLAG_LAST;
}

// ----------------------------------------------------------------------------
// Function hasFlagSecondary()
// ----------------------------------------------------------------------------

/**
.Function.hasFlagSecondary
..class:Class.BamAlignmentRecord
..cat:BAM I/O
..summary:Return true if a @Class.BamAlignmentRecord@ has the "secondary alignment" flag set.
..signature:hasFlagSecondary(record)
..param.record:The record to query.
...type:Class.BamAlignmentRecord
..returns:$bool$, indicating the flag's status.
..include:seqan/bam_io.h
..see:Function.hasFlagMultiple
..see:Function.hasFlagAllProper
..see:Function.hasFlagUnmapped
..see:Function.hasFlagNextUnmapped
..see:Function.hasFlagRC
..see:Function.hasFlagNextRC
..see:Function.hasFlagFirst
..see:Function.hasFlagLast
..see:Function.hasFlagQCNoPass
..see:Function.hasFlagDuplicate
*/

inline bool
hasFlagSecondary(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_SECONDARY) == BAM_FLAG_SECONDARY;
}

// ----------------------------------------------------------------------------
// Function hasFlagQCNoPass()
// ----------------------------------------------------------------------------

/**
.Function.hasFlagQCNoPass
..class:Class.BamAlignmentRecord
..cat:BAM I/O
..summary:Return true if a @Class.BamAlignmentRecord@ has the "does not pass quality controls" flag set.
..signature:hasFlagQCNoPass(record)
..param.record:The record to query.
...type:Class.BamAlignmentRecord
..returns:$bool$, indicating the flag's status.
..include:seqan/bam_io.h
..see:Function.hasFlagMultiple
..see:Function.hasFlagAllProper
..see:Function.hasFlagUnmapped
..see:Function.hasFlagNextUnmapped
..see:Function.hasFlagRC
..see:Function.hasFlagNextRC
..see:Function.hasFlagFirst
..see:Function.hasFlagLast
..see:Function.hasFlagSecondary
..see:Function.hasFlagDuplicate
*/

inline bool
hasFlagQCNoPass(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_QC_NO_PASS) == BAM_FLAG_QC_NO_PASS;
}

// ----------------------------------------------------------------------------
// Function hasFlagDuplicate()
// ----------------------------------------------------------------------------

/**
.Function.hasFlagDuplicate
..class:Class.BamAlignmentRecord
..cat:BAM I/O
..summary:Return true if a @Class.BamAlignmentRecord@ has the "PCR or optical duplicate" flag set.
..signature:hasFlagDuplicate(record)
..param.record:The record to query.
...type:Class.BamAlignmentRecord
..returns:$bool$, indicating the flag's status.
..include:seqan/bam_io.h
..see:Function.hasFlagMultiple
..see:Function.hasFlagAllProper
..see:Function.hasFlagUnmapped
..see:Function.hasFlagNextUnmapped
..see:Function.hasFlagRC
..see:Function.hasFlagNextRC
..see:Function.hasFlagFirst
..see:Function.hasFlagLast
..see:Function.hasFlagSecondary
..see:Function.hasFlagQCNoPass
*/

inline bool
hasFlagDuplicate(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_DUPLICATE) == BAM_FLAG_DUPLICATE;
}

// ----------------------------------------------------------------------------
// Function getAlignmentLengthInRef()
// ----------------------------------------------------------------------------

/**
.Function.getAlignmentLengthInRef
..class:Class.BamAlignmentRecord
..cat:BAM I/O
..summary:Returns length of @Class.BamAlignmentRecord@'s projection in reference.
..signature:getAlignmentLengthInRef(record)
..param.record:The record to query.
...type:Class.BamAlignmentRecord
..returns:$unsigned$, the alignment length in the reference.
..include:seqan/bam_io.h
*/

inline unsigned
getAlignmentLengthInRef(BamAlignmentRecord const & record)
{
    unsigned l = 0;
    _getLengthInRef(record.cigar, l);
    return l;
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BAM_IO_BAM_RECORD_H_
