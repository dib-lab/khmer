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
// Basic definitions for the stream module
// ==========================================================================

#ifndef SEQAN_STREAM_STREAM_BASE_
#define SEQAN_STREAM_STREAM_BASE_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Enums
// ============================================================================

// --------------------------------------------------------------------------
// Compression Type Tags
// --------------------------------------------------------------------------

/*!
 * @defgroup FileCompressionTags File Compression Tags
 * @brief Tags for describing file compression formats.
 */

/*!
 * @tag FileCompressionTags#GZFile
 * @headerfile <seqan/stream.h>
 * @brief File compression using the popular <a href="http://gzip.org">gzip</a> format.
 * @signature typedef Tag<GZFile_> GZFile;
 */

struct GZFile_;
typedef Tag<GZFile_> GZFile;

/*!
 * @tag FileCompressionTags#BgzFile
 * @headerfile <seqan/stream.h>
 * @signature typedef Tag<BgzfFile_> BgzfFile;
 * @brief File compression using the BGZF (Block GZip Format).
 *
 * The file format is described in the <a href="http://samtools.github.io/hts-specs/SAMv1.pdf">SAM file format
 * description</a>.
 */

struct BgzfFile_;
typedef Tag<BgzfFile_> BgzfFile;

/*!
 * @tag FileCompressionTags#BZ2File
 * @headerfile <seqan/stream.h>
 *
 * @brief File compression using the popular <a href="http://bzip.org">bzip2</a> format.
 *
 * @signature typedef Tag<BZ2File_> BZ2File;
 */

struct BZ2File_;
typedef Tag<BZ2File_> BZ2File;

// --------------------------------------------------------------------------
// MagicHeader
// --------------------------------------------------------------------------

template <typename T>
struct MagicHeader<GZFile, T>
{
    static char const VALUE[3];
};

template <typename T>
char const MagicHeader<GZFile, T>::VALUE[3] = { 0x1f, '\x8b', 0x08 };  // gzip's magic number


template <typename T>
struct MagicHeader<BZ2File, T>
{
    static char const VALUE[3];
};

template <typename T>
char const MagicHeader<BZ2File, T>::VALUE[3] = { 0x42, 0x5a, 0x68 };  // bzip2's magic number

// --------------------------------------------------------------------------
// FileExtensions
// --------------------------------------------------------------------------

template <typename T>
struct FileExtensions<GZFile, T>
{
    static char const * VALUE[2];
};

template <typename T>
char const * FileExtensions<GZFile, T>::VALUE[2] =
{
    ".gz",      // default output extension
    ".z"
};


template <typename T>
struct FileExtensions<BgzfFile, T>
{
    static char const * VALUE[2];
};

template <typename T>
char const * FileExtensions<BgzfFile, T>::VALUE[2] =
{
    ".bgzf",      // default output extension
    ".bam"        // BAM files are bgzf compressed

    // if you add extensions here, extend getBasename() below
};


template <typename T>
struct FileExtensions<BZ2File, T>
{
    static char const * VALUE[2];
};

template <typename T>
char const * FileExtensions<BZ2File, T>::VALUE[2] =
{
    ".bz2",      // default output extension
    ".bz"
};

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function _isPipe()
// --------------------------------------------------------------------------

inline bool
_isPipe(const char * fileName)
{
#ifdef PLATFORM_WINDOWS
    struct _stat buf;
    if (_stat(fileName, &buf) == 0)
        if ((buf.st_mode & _S_IFMT) == _S_IFCHR)
            return true;
#else
    struct stat buf;
    if (stat(fileName, &buf) == 0)
        if ((buf.st_mode & S_IFMT) == S_IFIFO ||
            (buf.st_mode & S_IFMT) == S_IFCHR)
            return true;
#endif
    return false;
}

// --------------------------------------------------------------------------
// Function guessFormat()
// --------------------------------------------------------------------------

// read first bytes of a file/stream and compare with file format's magic header
template <typename TStream, typename TFormat_>
inline bool
guessFormatFromStream(TStream &istream, Tag<TFormat_>)
{
    typedef Tag<TFormat_> TFormat;

    SEQAN_ASSERT(istream.good());

    if ((char *)MagicHeader<TFormat>::VALUE == NULL)
        return true;

    bool match = true;

    // check magic header
    unsigned i;
    for (i = 0; i != sizeof(MagicHeader<TFormat>::VALUE) / sizeof(char); ++i)
    {
        int c = (int)istream.get();
        if (c != (unsigned char)MagicHeader<TFormat>::VALUE[i])
        {
            match = false;
            if (c != EOF)
                ++i;
            break;
        }
    }

    // unget all read characters
    for (; i > 0; --i)
        istream.unget();

    SEQAN_ASSERT(istream.good());

    return match;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_STREAM_BASE_
