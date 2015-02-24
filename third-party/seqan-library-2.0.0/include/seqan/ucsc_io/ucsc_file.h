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
// Class for reading/writing files in UCSC format.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_UCSC_IO_UCSC_FILE_H_
#define SEQAN_INCLUDE_SEQAN_UCSC_IO_UCSC_FILE_H_

namespace seqan {

// ============================================================================
// Typedefs
// ============================================================================

typedef Tag<Ucsc_<> > Ucsc;


// ----------------------------------------------------------------------------
// Typedef UcscFileIn
// ----------------------------------------------------------------------------

/*!
 * @class UcscFileIn
 * @signature typedef FormattedFile<Ucsc, Input> UcscFileIn;
 * @extends FormattedFileIn
 * @headerfile <seqan/ucsc_io.h>
 * @brief Class for reading UCSC <tt>knownGenes.txt</tt> and <tt>knownIsoforms.txt</tt> files.
 *
 * @see UcscRecord
 */

typedef FormattedFile<Ucsc, Input>   UcscFileIn;

// ----------------------------------------------------------------------------
// Typedef UcscFileOut
// ----------------------------------------------------------------------------

/*!
 * @class UcscFileInOut
 * @signature typedef FormattedFile<Ucsc, Output> UcscFileOut;
 * @extends FormattedFileOut
 * @headerfile <seqan/ucsc_io.h>
 * @brief Class for writing UCSC <tt>knownGenes.txt</tt> and <tt>knownIsoforms.txt</tt> files.
 *
 * @see UcscRecord
 */

typedef FormattedFile<Ucsc, Output> UcscFileOut;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction FormattedFileContext
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TStorageSpec>
struct FormattedFileContext<FormattedFile<Ucsc, TDirection, TSpec>, TStorageSpec>
{
    typedef UcscIOContext Type;
};

// ----------------------------------------------------------------------------
// Metafunction FileFormats
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec>
struct FileFormat<FormattedFile<Ucsc, TDirection, TSpec> >
{
    typedef TagSelector<
                TagList<UcscKnownGene,
                TagList<UcscKnownIsoforms
                > >
            > Type;
};

// ----------------------------------------------------------------------------
// Function readRecord(); UcscRecord
// ----------------------------------------------------------------------------

// support for dynamically chosen file formats
template <typename TForwardIter>
inline void
readRecord(UcscRecord & /* record */,
           UcscIOContext & /* context */,
           TForwardIter & /* iter */,
           TagSelector<> const & /* format */)
{
    SEQAN_FAIL("UcscFileIn: File format not specified.");
}

template <typename TForwardIter, typename TTagList>
inline void
readRecord(UcscRecord & record,
           UcscIOContext & context,
           TForwardIter & iter,
           TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormat;

    if (isEqual(format, TFormat()))
        readRecord(record, context, iter, TFormat());
    else
        readRecord(record, context, iter, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

template <typename TSpec>
void readRecord(UcscRecord & record, FormattedFile<Ucsc, Input, TSpec> & file)
{
    readRecord(record, context(file), file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function writeRecord(); UcscRecord
// ----------------------------------------------------------------------------

// support for dynamically chosen file formats
template <typename TTarget>
inline void
writeRecord(TTarget & /* target */,
            UcscRecord const & /* record */,
            TagSelector<> const & /* format */)
{
    SEQAN_FAIL("UcscFileOut: File format not specified.");
}

template <typename TTarget, typename TTagList>
inline void
writeRecord(TTarget & target,
            UcscRecord const & record,
            TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormat;

    if (isEqual(format, TFormat()))
        writeRecord(target, record, TFormat());
    else
        writeRecord(target, record, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

template <typename TSpec>
inline void
writeRecord(FormattedFile<Ucsc, Output, TSpec> & file, UcscRecord const & record)
{
    writeRecord(file.iter, record, file.format);
}

}  // namespace seqan

#endif // SEQAN_INCLUDE_SEQAN_UCSC_IO_UCSC_FILE_H_
