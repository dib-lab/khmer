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
// TODO(esiragusa): tests

#ifndef SEQAN_STREAM_GUESS_FORMAT_
#define SEQAN_STREAM_GUESS_FORMAT_

namespace seqan {

// --------------------------------------------------------------------------
// Class FileExtensions
// --------------------------------------------------------------------------

template <typename TFormat, typename T>
struct FileExtensions;

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function guessFormatFromStream()
// --------------------------------------------------------------------------
// Base case: we get here if the file format could not be determined.

template <typename TFileSeq>
inline bool guessFormatFromStream(TFileSeq &, TagSelector<> &)
{
    return false;
}

// --------------------------------------------------------------------------
// Function guessFormatFromStream(TagSelector)
// --------------------------------------------------------------------------

template <typename TFileSeq, typename TTagList>
inline bool guessFormatFromStream(TFileSeq &seq, TagSelector<TTagList> &format)
{
    typedef typename TTagList::Type TFormat;

    if (value(format) == -1 || isEqual(format, TFormat()))
    {
        // if tagId is set to -1 (auto-detect) or the current format (TFormatTag) then test for TFormatTag format
        if (guessFormatFromStream(seq, TFormat()))
        {
            assign(format, TFormat());
            return true;
        }
    }
    return guessFormatFromStream(seq, static_cast<typename TagSelector<TTagList>::Base &>(format));
}

// --------------------------------------------------------------------------
// Function guessFormatFromFilename()
// --------------------------------------------------------------------------
// Base case: we get here if the file format could not be determined.

template <typename TFilename>
inline bool guessFormatFromFilename(TFilename const &, TagSelector<>)
{
    return false;
}

template <typename TFilename>
inline typename Prefix<TFilename const>::Type
getBasename(TFilename const & fname, TagSelector<> const &)
{
    return prefix(fname, length(fname));
}

// --------------------------------------------------------------------------
// Function guessFormatFromFilename(TagSelector)
// --------------------------------------------------------------------------

template <typename TFilename, typename TTagList>
inline bool guessFormatFromFilename(TFilename const &fname, TagSelector<TTagList> &format)
{
    typedef typename TTagList::Type TFormat;

    if (value(format) == -1 || isEqual(format, TFormat()))
    {
        // if tagId is set to -1 (auto-detect) or the current format (TFormatTag) then test for TFormatTag format
        if (guessFormatFromFilename(fname, TFormat()))
        {
            assign(format, TFormat());
            return true;
        }
    }
    return guessFormatFromFilename(fname, static_cast<typename TagSelector<TTagList>::Base &>(format));
}

template <typename TFilename, typename TTagList>
inline typename Prefix<TFilename const>::Type
getBasename(TFilename const & fname, TagSelector<TTagList> const &format)
{
    typedef typename TTagList::Type TFormat;

    if (isEqual(format, TFormat()))
        return getBasename(fname, TFormat());
    else
        return getBasename(fname, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

// --------------------------------------------------------------------------
// Function guessFormatFromFilename(Tag)
// --------------------------------------------------------------------------

/*!
 * @fn guessFormatFromFilename
 * @headerfile <seqan/file.h>
 * @brief Guesses a file format from a sequence file name.
 *
 * @signature bool guessFormatFromFilename(fileName, formatTag);
 *
 * @param[in] fileName  A filename of a sequence file.
 * @param[in] formatTag A file format tag.
 *
 * @return bool <tt>true</tt> if the format represented by <tt>formatTag</tt> was recognized in <tt>fileName</tt>.
 */

template <typename TFilename, typename TFormat_>
inline bool guessFormatFromFilename(TFilename const & fileName, Tag<TFormat_> /*formatTag*/)
{
    typedef typename Value<TFilename>::Type                                     TValue;
    typedef ModifiedString<TFilename const, ModView<FunctorLowcase<TValue> > >    TLowcase;
    typedef Tag<TFormat_>                                                       TFormat;

    TLowcase lowcaseFileName(fileName);
    for (unsigned i = 0; i < sizeof(FileExtensions<TFormat>::VALUE) / sizeof(char*); ++i)
        if (endsWith(lowcaseFileName, lowerString(FileExtensions<TFormat>::VALUE[i])))
            return true;

    return false;
}

template <typename TFilename, typename TFormat_>
inline typename Prefix<TFilename const>::Type
getBasename(TFilename const & fileName, Tag<TFormat_> const & /*formatTag*/)
{
    typedef typename Value<TFilename>::Type                                     TValue;
    typedef ModifiedString<TFilename const, ModView<FunctorLowcase<TValue> > >    TLowcase;
    typedef Tag<TFormat_>                                                       TFormat;

    TLowcase lowcaseFileName(fileName);
    for (unsigned i = 0; i < sizeof(FileExtensions<TFormat>::VALUE) / sizeof(char*); ++i)
        if (endsWith(lowcaseFileName, lowerString(FileExtensions<TFormat>::VALUE[i])))
            return prefix(fileName, length(fileName) - length(FileExtensions<TFormat>::VALUE[i]));

    return prefix(fileName, length(fileName));
}


// --------------------------------------------------------------------------
// Function _getFileExtensions()
// --------------------------------------------------------------------------

template <typename TStringSet, typename TFormat_>
inline void _getFileExtensions(TStringSet &stringSet, Tag<TFormat_> /*formatTag*/, bool primaryExtensionOnly = false)
{
    typedef Tag<TFormat_> TFormat;
    if (primaryExtensionOnly)
    {
        appendValue(stringSet, FileExtensions<TFormat>::VALUE[0]);    // first is primary
    }
    else
    {
        for (unsigned i = 0; i < sizeof(FileExtensions<TFormat>::VALUE) / sizeof(char*); ++i)
            appendValue(stringSet, FileExtensions<TFormat>::VALUE[i]);
    }
}

template <typename TStringSet, typename TTag>
inline void _getFileExtensions(TStringSet &stringSet, TagList<TTag, void> const /*formatTag*/, bool primaryExtensionOnly = false)
{
    _getFileExtensions(stringSet, TTag(), primaryExtensionOnly);
}

template <typename TStringSet, typename TTag, typename TSubList>
inline void _getFileExtensions(TStringSet &stringSet, TagList<TTag, TSubList> const /*formatTag*/, bool primaryExtensionOnly = false)
{
    _getFileExtensions(stringSet, TTag(), primaryExtensionOnly);
    _getFileExtensions(stringSet, TSubList(), primaryExtensionOnly);
}

template <typename TStringSet, typename TTagList>
inline void _getFileExtensions(TStringSet &stringSet, TagSelector<TTagList> const /*formatTag*/, bool primaryExtensionOnly = false)
{
    _getFileExtensions(stringSet, TTagList(), primaryExtensionOnly);
}

} // namespace seqan

#endif  // #ifndef SEQAN_STREAM_GUESS_FORMAT_
