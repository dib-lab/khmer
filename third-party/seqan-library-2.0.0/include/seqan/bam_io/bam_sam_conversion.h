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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Code to convert between SAM and BAM format tags (textual <-> binary).
// ==========================================================================

#ifndef INCLUDE_SEQAN_BAM_IO_BAM_SAM_CONVERSION_H_
#define INCLUDE_SEQAN_BAM_IO_BAM_SAM_CONVERSION_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function assignTagsSamToBam()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TBuffer>
struct AppendTagsSamToBamOneTagHelper_
{
    TTarget     &target;
    TBuffer     buffer;
    char        typeC;

    AppendTagsSamToBamOneTagHelper_(TTarget &target, TBuffer buffer, char typeC):
        target(target),
        buffer(buffer),
        typeC(typeC)
    {}

    template <typename Type>
    bool operator() (Type)
    {
        if (BamTypeChar<Type>::VALUE != typeC)
            return false;

        appendRawPod(target, lexicalCast<Type>(buffer));
        return true;
    }
};

template <typename TTarget, typename TForwardIter>
void _appendTagsSamToBamOneTag(TTarget & target, TForwardIter & iter, CharString & buffer)
{
    write(target, iter, 2);

    clear(buffer);
    write(buffer, iter, 3);
    SEQAN_ASSERT_EQ(buffer[0], ':');
    SEQAN_ASSERT_EQ(buffer[2], ':');

    char typeC = buffer[1];
    appendValue(target, typeC);

    switch (typeC)
    {
        case 'Z':
        case 'H':
            // BAM string
            // TODO(holtgrew): Could test on even length in case of 'H'.
            readUntil(target, iter, OrFunctor<IsTab, IsNewline>());
            appendValue(target, '\0');
            break;

        case 'B':
        {
            // BAM array

            // Read type.
            readOne(typeC, iter);
            appendValue(target, typeC);

            // Read array contents.
            clear(buffer);
            readUntil(buffer, iter, OrFunctor<IsTab, IsNewline>());

            size_t len = length(buffer);

            // Count number of entries (== number of commas after type character).
            __uint32 nEntries = 0;
            for (size_t i = 0; i != len; ++i)
                if (buffer[i] == ',')
                    ++nEntries;

            // Write out array length.
            appendRawPod(target, (__uint32)nEntries);

            // Write out array values.
            size_t startPos = 1;
            for (unsigned i = 0; i < nEntries; ++i)
            {
                SEQAN_ASSERT_LT(startPos, len);

                // search end of current entry
                size_t endPos = startPos;
                for (; endPos < len; ++endPos)
                    if (buffer[endPos] == ',')
                    {
                        buffer[endPos] = '\0';
                        break;
                    }

                AppendTagsSamToBamOneTagHelper_<TTarget, char*> func(target,
                                                                     toCString(buffer) + startPos,
                                                                     typeC);
                if (!tagApply(func, BamTagTypes()))
                    SEQAN_ASSERT_FAIL("Invalid tag type: %c!", typeC);

                startPos = endPos + 1;
            }
            break;
        }

        default:
        {
            // BAM simple value
            clear(buffer);
            readUntil(buffer, iter, OrFunctor<IsTab, IsNewline>());

            AppendTagsSamToBamOneTagHelper_<TTarget, CharString&> func(target, buffer, typeC);
            if (!tagApply(func, BamTagTypes()))
                SEQAN_ASSERT_FAIL("Invalid tag type: %c!", typeC);
        }
    }
}

/*!
 * @fn assignTagsSamToBam
 * @headerfile <seqan/bam_io.h>
 * @brief Assign tags in SAM format to tags in BAM format.
 *
 * @signature void assignTagsBamToSam(bamTags, samTags);
 *
 * @param[out] bamTags A sequence of <tt>char</tt> (e.g. @link CharString @endlink) for the target BAM tags.
 * @param[in]  samTags A sequence of <tt>char</tt> (e.g. @link CharString @endlink) for the source SAM tags.
 *
 * @see assignTagsBamToSam
 */

template <typename TTarget, typename TSource>
void appendTagsSamToBam(TTarget & target, TSource const & source)
{
    // Handle case of empty source sequence.
    if (empty(source))
        return;

    typedef typename Iterator<TSource const, Rooted>::Type TSourceIter;
    TSourceIter it = begin(source, Rooted());

    CharString buffer;

    while (!atEnd(it))
    {
        if (value(it) == '\t')
            skipOne(it);

        _appendTagsSamToBamOneTag(target, it, buffer);
    }
}

template <typename TTarget, typename TSource>
void assignTagsSamToBam(TTarget & target, TSource const & source)
{
    clear(target);
    appendTagsSamToBam(target, source);
}

// ----------------------------------------------------------------------------
// Function assignTagsBamToSam()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TSourceIter>
struct AssignTagsBamToSamOneTagHelper_
{
    TTarget     &target;
    TSourceIter &it;
    char        typeC;

    AssignTagsBamToSamOneTagHelper_(TTarget &target, TSourceIter &it, char typeC):
        target(target),
        it(it),
        typeC(typeC)
    {}

    template <typename Type>
    bool operator() (Type)
    {
        if (BamTypeChar<Type>::VALUE != typeC)
            return false;

        appendNumber(target, reinterpret_cast<Type const &>(*it));
        it += sizeof(Type);
        return true;
    }

    bool operator() (char)
    {
        if (BamTypeChar<char>::VALUE != typeC)
            return false;

        writeValue(target, getValue(it));
        ++it;
        return true;
    }
};

template <typename TTarget, typename TSourceIter>
void _appendTagsBamToSamOneTag(TTarget & target, TSourceIter & it)
{
    // Copy tag name.
    SEQAN_ASSERT_NOT(atEnd(it));
    writeValue(target, *it++);
    SEQAN_ASSERT_NOT(atEnd(it));
    writeValue(target, *it++);

    // Add ':'.
    writeValue(target, ':');

    char typeC = *it++;
    char c = FunctorLowcase<char>()(typeC);

    // The only integer type supported is a 32bit signed int (SAM Format Spec, 28 Feb 2014, Section 1.5)
    // This sucks as this projection is not identically reversible
    if (c == 'c' || c == 's' || c == 'i')
        writeValue(target, 'i');
    else
        writeValue(target, typeC);

    // Add ':'.
    writeValue(target, ':');

    switch (typeC)
    {
        case 'Z':
        case 'H':
            // BAM string
            SEQAN_ASSERT_NOT(atEnd(it));
            while (*it != '\0')
            {
                writeValue(target, *it);
                ++it;
                SEQAN_ASSERT_NOT(atEnd(it));
            }
            ++it;
            break;

        case 'B':
        {
            // BAM array
            typeC = *it++;
            writeValue(target, typeC);
            AssignTagsBamToSamOneTagHelper_<TTarget, TSourceIter> func(target, it, typeC);

            // Read array length.
            union {
                char raw[4];
                unsigned len;
            } tmp;
            for (unsigned i = 0; i < 4; ++i)
            {
                SEQAN_ASSERT_NOT(atEnd(it));
                tmp.raw[i] = *it++;
            }
            for (unsigned i = 0; i < tmp.len; ++i)
            {
                writeValue(target, ',');
                if (!tagApply(func, BamTagTypes()))
                    SEQAN_ASSERT_FAIL("Invalid tag type: %c!", typeC);
            }
            break;
        }

        default:
        {
            // BAM simple value
            AssignTagsBamToSamOneTagHelper_<TTarget, TSourceIter> func(target, it, typeC);
            if (!tagApply(func, BamTagTypes()))
                SEQAN_ASSERT_FAIL("Invalid tag type: %c!", typeC);
        }
    }
}

/*!
 * @fn assignTagsBamToSam
 * @headerfile <seqan/bam_io.h>
 * @brief Assign tags in BAM format to tags in SAM format.
 *
 * @signature void assignTagsBamToSam(samTags, bamTags);
 *
 * @param[out] samTags A sequence of <tt>char</tt> (e.g. @link CharString @endlink) for the target SAM tags.
 * @param[in]  bamTags A sequence of <tt>char</tt> (e.g. @link CharString @endlink) for the source BAM tags.
 *
 * @see assignTagsSamToBam
 */

template <typename TTarget, typename TSource>
inline void
appendTagsBamToSam(TTarget & target, TSource const & source)
{
    if (empty(source))
        return;

    typedef typename Iterator<TSource const, Rooted>::Type TSourceIter;
    TSourceIter it = begin(source, Rooted());

    while (true)
    {
        _appendTagsBamToSamOneTag(target, it);
        if (atEnd(it))
            return;
        writeValue(target, '\t');
    }
}

template <typename TTarget, typename TSource>
void assignTagsBamToSam(TTarget & target, TSource const & source)
{
    clear(target);
    appendTagsBamToSam(target, source);
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_BAM_IO_BAM_SAM_CONVERSION_H_
