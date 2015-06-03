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
// Class for reading/writing files in Vcf format.
// ==========================================================================
// TODO(weese:) add Bcf I/O and integrate it

#ifndef SEQAN_VCF_IO_VCF_FILE_H_
#define SEQAN_VCF_IO_VCF_FILE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

struct Vcf_;
typedef Tag<Vcf_> Vcf;

struct Bcf_;
typedef Tag<Bcf_> Bcf;

// ============================================================================
// Typedefs
// ============================================================================

// ----------------------------------------------------------------------------
// Typedef VcfFileIn
// ----------------------------------------------------------------------------

/*!
 * @class VcfFileIn
 * @signature typedef FormattedFile<Vcf, Input> VcfFileIn;
 * @extends FormattedFileIn
 * @headerfile <seqan/vcf_io.h>
 * @brief Class for reading VCF files.
 *
 * @see VcfHeader
 * @see VcfRecord
 */

typedef FormattedFile<Vcf, Input>   VcfFileIn;

// ----------------------------------------------------------------------------
// Typedef VcfFileOut
// ----------------------------------------------------------------------------

/*!
 * @class VcfFileOut
 * @signature typedef FormattedFile<Vcf, Output> VcfFileOut;
 * @extends FormattedFileOut
 * @headerfile <seqan/vcf_io.h>
 * @brief Class for writing VCF files.
 *
 * @see VcfHeader
 * @see VcfRecord
 */

typedef FormattedFile<Vcf, Output>  VcfFileOut;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Class MagicHeader
// ----------------------------------------------------------------------------

template <typename T>
struct MagicHeader<Vcf, T>
{
    static unsigned char const VALUE[16];
};

template <typename T>
unsigned char const MagicHeader<Vcf, T>::VALUE[16] =
{
    '#', '#', 'f', 'i', 'l', 'e', 'f', 'o', 'r', 'm', 'a', 't', '=', 'V', 'C', 'F'  // VCF's magic header
};

template <typename T>
struct MagicHeader<Bcf, T>
{
    static unsigned char const VALUE[5];
};

template <typename T>
unsigned char const MagicHeader<Bcf, T>::VALUE[5] = { 'B', 'C', 'F', '\2', '\1' };  // BCF2's magic header

// ----------------------------------------------------------------------------
// Class FileExtensions
// ----------------------------------------------------------------------------

template <typename T>
struct FileExtensions<Vcf, T>
{
    static char const * VALUE[1];    // default is one extension
};

template <typename T>
char const * FileExtensions<Vcf, T>::VALUE[1] =
{
    ".vcf"     // default output extension
};

template <typename T>
struct FileExtensions<Bcf, T>
{
    static char const * VALUE[1];    // default is one extension
};

template <typename T>
char const * FileExtensions<Bcf, T>::VALUE[1] =
{
    ".bcf"     // default output extension
};

// ----------------------------------------------------------------------------
// Metafunction FormattedFileContext
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TStorageSpec>
struct FormattedFileContext<FormattedFile<Vcf, TDirection, TSpec>, TStorageSpec>
{
    typedef StringSet<CharString>                                   TNameStore;
    typedef NameStoreCache<TNameStore>                              TNameStoreCache;
    typedef VcfIOContext<TNameStore, TNameStoreCache, TStorageSpec> Type;
};

// ----------------------------------------------------------------------------
// Metafunction FileFormats
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec>
struct FileFormat<FormattedFile<Vcf, TDirection, TSpec> >
{
// TODO(weese:) Enable this, as soon as someone implements BCF

//#if SEQAN_HAS_ZLIB
//    typedef TagSelector<
//                TagList<Bcf,
//                TagList<Vcf
//                > >
//            > Type;
//#else
    typedef Vcf Type;
//#endif
};

// --------------------------------------------------------------------------
// Function _mapBamFormatToCompressionFormat()
// --------------------------------------------------------------------------

inline BgzfFile
_mapFileFormatToCompressionFormat(Bcf)
{
    return BgzfFile();
}

// ----------------------------------------------------------------------------
// Function readHeader(); VcfHeader
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
readHeader(VcfHeader & header, FormattedFile<Vcf, Input, TSpec> & file)
{
    readHeader(header, context(file), file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function readRecord(); VcfRecord
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
readRecord(VcfRecord & record, FormattedFile<Vcf, Input, TSpec> & file)
{
    readRecord(record, context(file), file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function writeHeader(); VcfHeader
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
writeHeader(FormattedFile<Vcf, Output, TSpec> & file, VcfHeader & header)
{
    writeHeader(file.iter, header, context(file), file.format);
}

// ----------------------------------------------------------------------------
// Function writeRecord(); VcfRecord
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
writeRecord(FormattedFile<Vcf, Output, TSpec> & file, VcfRecord & record)
{
    writeRecord(file.iter, record, context(file), file.format);
}

}  // namespace seqan

#endif // SEQAN_VCF_IO_VCF_FILE_H_
