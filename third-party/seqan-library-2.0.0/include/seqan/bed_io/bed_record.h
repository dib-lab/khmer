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

#ifndef INCLUDE_SEQAN_BED_IO_BED_RECORD_H_
#define INCLUDE_SEQAN_BED_IO_BED_RECORD_H_

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct Bed3_;
typedef Tag<Bed3_> Bed3;

struct Bed4_;
typedef Tag<Bed4_> Bed4;

struct Bed5_;
typedef Tag<Bed5_> Bed5;

struct Bed6_;
typedef Tag<Bed6_> Bed6;

struct Bed12_;
typedef Tag<Bed12_> Bed12;

// ----------------------------------------------------------------------------
// Class BedRecord
// ----------------------------------------------------------------------------

/*!
 * @class BedRecord
 * @implements FormattedFileRecordConcept
 * @headerfile <seqan/bed_io.h>
 * @brief Data structure for storing BED records.
 *
 * @signature template <typename TSpec>
 *            class BedRecord;
 *
 * @tparam TSpec The specialization to use.  Default: <tt>Bed12</tt>.
 *
 * BED files allow the easy representation of intervals on the genome.  Originally, they were designed for tracks in the
 * UCSC genome browser. The original format has 12 columns but often variants using fewer columns with interpreted data
 * are used and the rest is kept as application dependent data.
 *
 * The BedRecord class allows for storing BED records. The various subclasses provide access to 3, 4, 5, 6, or 12 fields
 * of the BED format. For example, a <tt>BedRecord&lt;Bed&gt;></tt> has members variables for the first 5 columns of a BED
 * file. The remaining data is stored as the @link CharString @endlink member variable <tt>data</tt>.
 *
 * @section Remarks
 *
 * The <tt>ref</tt> field is the name of the reference as loaded from the BED file. The <tt>rID</tt> field can be used
 * to store a numeric reference id.
 *
 * Note that while the BED file format is 1-based, the coordinates in the BedRecord are 0-based.
 *
 * @var CharString BedRecord::ref;
 * @brief Name of the interval's reference name.
 *
 * @var __int32 BedRecord::beginPosition;
 * @brief Begin position on the reference.
 *
 * @var __int32 BedRecord::rID;
 * @brief Numeric id of the interval's reference (<tt>__int32</tt>, defaults to <tt>INVALID_REFID</tt>).
 *
 * @var __int32 BedRecord::INVALID_REFID;
 * @brief Constant for invalid references.
 * @signature static const __int32 BedRecord::INVALID_REFID = -1;
 *
 * @var __int32 BedRecord::endPosition;
 * @brief End position on the reference.
 *
 * @var __int32 BedRecord::INVALID_POS;
 * @brief Constant for invalid positions.
 * @signature static const __int32 BedRecord::INVALID_POS = -1;
 *
 * @var CharString BedRecord::data;
 * @brief Any data after the last position.
 */

/*!
 * @fn BedRecord::BedRecord
 * @brief Constructor.
 *
 * @signature BedRecord::BedRecord();
 */

template <typename TSpec = Bed12>
class BedRecord;

// ----------------------------------------------------------------------------
// Class Bed3 BedRecord
// ----------------------------------------------------------------------------

/*!
 * @class Bed3Record
 * @extends BedRecord
 * @headerfile <seqan/bed_io.h>
 * @brief BedRecord with 3 fields.
 *
 * @signature template <>
 *            class BedRecord<Bed3>;
 *
 * This BedRecord specialization stores the first three fields (ref, beginPos, endPos) of a BED file.
 */

template <>
class BedRecord<Bed3>
{
public:
    static const int INVALID_REFID = -1;
    static const int INVALID_POS = -1;

    // The chromosome name.
    CharString ref;
    // The id of the chromosome, -1 if not translated.
    __int32 rID;
    // The start position.
    __int32 beginPos;
    // The end position;
    __int32 endPos;
    // The remaining data from the file, unparsed.
    CharString data;

    BedRecord() : rID(INVALID_REFID), beginPos(INVALID_POS), endPos(INVALID_POS)
    {}

    void _clear()
    {
        rID = INVALID_REFID;
        beginPos = INVALID_POS;
        endPos = INVALID_POS;
        clear(ref);
        clear(data);
    }
};

// ----------------------------------------------------------------------------
// Class Bed4 BedRecord
// ----------------------------------------------------------------------------

/*!
 * @class Bed4Record
 * @extends Bed3Record
 * @headerfile <seqan/bed_io.h>
 * @brief BedRecord with 4 fields.
 *
 * @signature template <>
 *            class BedRecord<Bed3>;
 *
 * This BedRecord specialization stores the first four fields (ref, beginPos, endPos, name) of a BED file.
 *
 * @var CharString Bed4Record::name;
 * @brief The name of the interval (@link CharString @endlink).
 */

template <>
class BedRecord<Bed4> : public BedRecord<Bed3>
{
public:
    // The name of the feature.
    CharString name;

    BedRecord() : BedRecord<Bed3>()
    {}

    void _clear()
    {
        BedRecord<Bed3>::_clear();
        clear(name);
    }
};

// ----------------------------------------------------------------------------
// Class Bed5 BedRecord
// ----------------------------------------------------------------------------

/*!
 * @class Bed5Record
 * @extends Bed4Record
 * @headerfile <seqan/bed_io.h>
 *
 * @brief BedRecord with 5 fields.
 *
 * @signature template <>
 *            class BedRecord<Bed5>;
 *
 * This BedRecord specialization stores the first five fields (ref, beginPos, endPos, name, score) of a BED file.
 *
 * @var CharString Bed5Record::score;
 * @brief The score of the interval (stored as @link CharString @endlink to allow more flexible annotation).
 *
 * @section Remarks
 *
 * Storing the score as a @link CharString @endlink is provided for compatibility with bedtools.
 */

template <>
class BedRecord<Bed5> : public BedRecord<Bed4>
{
public:
    // The score of the feature, stored as CharString for compatibility with bedtools.
    CharString score;

    BedRecord() : BedRecord<Bed4>()
    {}

    void _clear()
    {
        BedRecord<Bed4>::_clear();
        clear(score);
    }
};

// ----------------------------------------------------------------------------
// Class Bed6 BedRecord
// ----------------------------------------------------------------------------

/*!
 * @class Bed6Record
 * @extends Bed5Record
 * @headerfile <seqan/bed_io.h>
 *
 * @brief BedRecord with 6 fields.
 *
 * @signature template <>
 *            class BedRecord<Bed6>;
 *
 * This BedRecord specialization stores the first six fields (ref, beginPos, endPos, name, score, strand) of a BED file.
 *
 * @var char Bed6Record::strand;
 * @brief The strand of the interval (stored as <tt>char</tt>, one of <tt>.</tt>, '-', and <tt>+</tt>).
 *
 * Defaults to '.'.
 */

template <>
class BedRecord<Bed6> : public BedRecord<Bed5>
{
public:
    // The strand of the feature.  One of '.', '-', and '+'.
    char strand;

    BedRecord() : BedRecord<Bed5>(), strand('.')
    {}

    void _clear()
    {
        BedRecord<Bed5>::_clear();
        strand = '.';
    }
};

// ----------------------------------------------------------------------------
// Class BedRgb
// ----------------------------------------------------------------------------

/*!
 * @class BedRgb
 * @headerfile <seqan/bed_io.h>
 * @brief RGB color for @link Bed12Record @endlink.
 *
 * @signature class BedRgb;
 *
 * @var __int32 BedRgb::red;
 * @brief Red value of RGB color (default is <tt>0</tt>).
 *
 * @var __int32 BedRgb::green;
 * @brief Green value of RGB color (default is <tt>0</tt>).
 *
 * @var __int32 BedRgb::blue;
 * @brief Blue value of RGB color (default is <tt>0</tt>).
 */

/*!
 * @fn BedRgb::BedRgb
 * @brief Default constructor and initialization of integer RGB values.
 *
 * @signature BedRgb::BedRgb();
 * @signature BedRgb::BedRgb(red, green, blue);
 *
 * @param[in] blue  __int32 blue value <tt>0-255</tt> (defaults to <tt>0</tt>).
 * @param[in] green __int32 green value <tt>0-255</tt> (defaults to <tt>0</tt>).
 * @param[in] red   __int32 red value <tt>0-255</tt> (defaults to <tt>0</tt>).
 */

class BedRgb
{
public:
    __int32 red, green, blue;

    BedRgb() : red(0), green(0), blue(0)
    {}

    BedRgb(int red, int green, int blue) : red(red), green(green), blue(blue)
    {}

    bool operator==(BedRgb const & other) const
    {
        return red == other.red && green == other.green && blue == other.blue;
    }

    bool operator!=(BedRgb const & other) const
    {
        return !(*this == other);
    }
};

// ----------------------------------------------------------------------------
// Class Bed12 BedRecord
// ----------------------------------------------------------------------------

/*!
 * @class Bed12Record
 * @extends BedRecord
 * @headerfile <seqan/bed_io.h>
 * @brief BedRecord with 12 fields.
 *
 * @signature template <>
 *            class BedRecord<Bed12>;
 *
 * This @link BedRecord @endlink specialization stores all fields of a BED file.
 *
 * @var __int32 Bed12Record::itemRgb;
 * @brief RGB color of item (@link BedRgb @endlink).
 *
 * @var __int32 Bed12Record::blockCount;
 * @brief The number of blocks.
 *
 * @var TIntString Bed12Record::blockBegins;
 * @brief The begin positions of the blocks (@link AllocString @endlink of <tt>__int32</tt>).
 *
 * @var TIntString Bed12Record::blockSizes;
 * @brief The sizes of the blocks (@link AllocString @endlink of <tt>__int32</tt>).
 *
 * @var __int32 Bed12Record::thickBegin;
 * @brief The begin position of thick drawing.
 *
 * @var __int32 Bed12Record::thickEnd;
 * @brief The end position of thick drawing.
 */

template <>
class BedRecord<Bed12> : public BedRecord<Bed6>
{
public:
    // The starting position of thick line for feature.
    __int32 thickBegin;
    // The end position of thick line for feature.
    __int32 thickEnd;
    // The color of the item.
    BedRgb itemRgb;
    // The number of blocks/exons for the feature.
    __int32 blockCount;
    // The list of block size.
    String<__int32> blockSizes;
    // List of block starts.
    String<__int32> blockBegins;

    BedRecord() : BedRecord<Bed6>(), thickBegin(INVALID_POS), thickEnd(INVALID_POS), blockCount(0)
    {}

    void _clear()
    {
        BedRecord<Bed6>::_clear();
        thickBegin = INVALID_POS;
        thickEnd = INVALID_POS;
        blockCount = 0;
        itemRgb = BedRgb();
        clear(blockSizes);
        clear(blockBegins);
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

/*!
 * @fn BedRecord#clear
 * @brief Reset BED record to state after default initialization.
 *
 * @signature void clear(record);
 *
 * @param[in,out] record BedRecord to reset.
 */

template <typename TSpec>
void clear(BedRecord<TSpec> & record)
{
    record._clear();
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_BED_IO_BED_RECORD_H_
