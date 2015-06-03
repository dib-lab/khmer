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
// Class for reading/writing files in BED format.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_BED_IO_BED_FILE_H_
#define SEQAN_INCLUDE_SEQAN_BED_IO_BED_FILE_H_

namespace seqan {

// ============================================================================
// Typedefs
// ============================================================================

// ----------------------------------------------------------------------------
// Type BedFileIn
// ----------------------------------------------------------------------------

/*!
 * @class BedFileIn
 * @signature typedef FormattedFile<Bed, Input> BedFileIn;
 * @extends FormattedFileIn
 * @headerfile <seqan/bed_io.h>
 * @brief Class for reading BED files.
 *
 * @see BedRecord
 */

typedef FormattedFile<Bed, Input>   BedFileIn;

// ----------------------------------------------------------------------------
// Type BedFileOut
// ----------------------------------------------------------------------------

/*!
 * @class BedFileOut
 * @signature typedef FormattedFile<Bed, Output> BedFileOut;
 * @extends FormattedFileOut
 * @headerfile <seqan/bed_io.h>
 * @brief Class for writing BED files.
 *
 * @see BedRecord
 */

typedef FormattedFile<Bed, Output>  BedFileOut;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction FormattedFileContext
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TStorageSpec>
struct FormattedFileContext<FormattedFile<Bed, TDirection, TSpec>, TStorageSpec>
{
    typedef CharString Type;
};

// ----------------------------------------------------------------------------
// Metafunction FileFormats
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec>
struct FileFormat<FormattedFile<Bed, TDirection, TSpec> >
{
    typedef Bed Type;
};

// ----------------------------------------------------------------------------
// Function readRecord(); BedRecord
// ----------------------------------------------------------------------------

// convient BedFile variant
template <typename TRecordSpec, typename TSpec>
inline void
readRecord(BedRecord<TRecordSpec> & record, FormattedFile<Bed, Input, TSpec> & file)
{
    readRecord(record, context(file), file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function write(); BedRecord
// ----------------------------------------------------------------------------

template <typename TSpec, typename TRecordSpec>
inline void
writeRecord(FormattedFile<Bed, Output, TSpec> & file, BedRecord<TRecordSpec> & record)
{
    writeRecord(file.iter, record, file.format);
}

}  // namespace seqan

#endif // SEQAN_INCLUDE_SEQAN_BED_IO_BED_FILE_H_
