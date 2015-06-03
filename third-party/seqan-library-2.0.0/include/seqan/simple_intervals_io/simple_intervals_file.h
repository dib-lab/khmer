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
// FormattedFile for .intervals format
// ==========================================================================

#ifndef INCLUDE_SEQAN_SIMPLE_INTERVALS_IO_SIMPLE_INTERVALS_FILE_H_
#define INCLUDE_SEQAN_SIMPLE_INTERVALS_IO_SIMPLE_INTERVALS_FILE_H_

#include <seqan/simple_intervals_io/simple_intervals_io.h>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Typedef SimpleIntervalsFileIn
// ----------------------------------------------------------------------------

/*!
 * @class SimpleIntervalsFileIn
 * @extends FormattedFileIn
 * @headerfile <seqan/ucsc_io.h>
 * @brief @link FormattedFile @endlink for reading .intervals files.
 *
 * @signature typedef FormattedFile<SimpleIntervals, Input> SimpleIntervalsFileIn;
 */
typedef FormattedFile<SimpleIntervals, Input>   SimpleIntervalsFileIn;

// ----------------------------------------------------------------------------
// Typedef SimpleIntervalsFileOut
// ----------------------------------------------------------------------------

/*!
 * @class SimpleIntervalsFileInOut
 * @extends FormattedFileOut
 * @headerfile <seqan/ucsc_io.h>
 * @brief @link FormattedFile @endlink for reading .intervals files.
 *
 * @signature typedef FormattedFile<SimpleIntervals, Output> SimpleIntervalsFileOut;
 */
typedef FormattedFile<SimpleIntervals, Output> SimpleIntervalsFileOut;


// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction MagicHeader
// ----------------------------------------------------------------------------

template <typename T>
struct MagicHeader<SimpleIntervals, T> : public MagicHeader<Nothing, T> {};

// ----------------------------------------------------------------------------
// Metafunction FileExtensions
// ----------------------------------------------------------------------------

template <typename T>
struct FileExtensions<SimpleIntervals, T>
{
    static char const * VALUE[1];  // default is one extension
};

template <typename T>
char const * FileExtensions<SimpleIntervals, T>::VALUE[1] =
{
    ".intervals"      // default output extension
};

// ----------------------------------------------------------------------------
// Metafunction FormattedFileContext
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TStorageSpec>
struct FormattedFileContext<FormattedFile<SimpleIntervals, TDirection, TSpec>, TStorageSpec>
{
    typedef SimpleIntervalsIOContext Type;
};

// ----------------------------------------------------------------------------
// Metafunction FileFormats
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec>
struct FileFormat<FormattedFile<SimpleIntervals, TDirection, TSpec> >
{
    typedef SimpleIntervals Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecord()
// ----------------------------------------------------------------------------

/*!
 * @fn SimpleIntervalsFileIn#readRecord
 * @brief Reading records from a SimpleIntervalsFileIn.
 *
 * @signature void readRecord(record, file);
 *
 * @param[out] record The resulting @link GenomicRegion @endlink.
 * @param[out] file   The SimpleIntervalsFileIn to read from.
 *
 * @throw IOError in case of problems.
 */
template <typename TSpec>
void readRecord(GenomicRegion & record, FormattedFile<SimpleIntervals, Input, TSpec> & file)
{
    readRecord(record, context(file), file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function write()
// ----------------------------------------------------------------------------

/*!
 * @fn SimpleIntervalsFileIn#writeRecord
 * @brief Writing records to a SimpleIntervalsFileIn.
 *
 * @signature void readRecord(record, file);
 *
 * @param[out] file   The SimpleIntervalsFileIn to write to.
 * @param[out] record The @link GenomicRegion @endlink to write out.
 *
 * @throw IOError in case of problems.
 */
template <typename TSpec>
void writeRecord(FormattedFile<SimpleIntervals, Output, TSpec> & file, GenomicRegion & record)
{
    writeRecord(file.iter, context(file), record, file.format);
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_SIMPLE_INTERVALS_IO_SIMPLE_INTERVALS_FILE_H_
