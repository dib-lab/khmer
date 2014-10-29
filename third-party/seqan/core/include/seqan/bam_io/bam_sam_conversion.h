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
// Code to convert between SAM and BAM format tags (textual <-> binary).
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_BAM_IO_BAM_SAM_CONVERSION_H_
#define CORE_INCLUDE_SEQAN_BAM_IO_SAM_TAGS_TO_BAM_TAGS_H_

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

template <typename TTarget, typename TRecordReader>
void _assignTagsSamToBamOneTag(TTarget & target, TRecordReader & reader, CharString & buffer)
{
    SEQAN_ASSERT_NOT(atEnd(reader));
    int res = readNChars(target, reader, 2);  // Read tag name.
    (void)res;  // If run without assertions.
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_NOT(atEnd(reader));
    
    clear(buffer);
    res = readNChars(buffer, reader, 3);  // Read ':<type>:'.
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(buffer[0], ':');
    SEQAN_ASSERT_EQ(buffer[2], ':');
    char t = buffer[1];
    appendValue(target, t);

    SEQAN_ASSERT_NOT(atEnd(reader));
    
    switch (t)
    {
    case 'A':
        clear(buffer);
        res = readNChars(target, reader, 1);
        SEQAN_ASSERT_EQ(res, 0);
        break;
    case 'i':
        {
            clear(buffer);
            res = readUntilTabOrLineBreak(buffer, reader);
            SEQAN_ASSERT(res == 0 || res == EOF_BEFORE_SUCCESS);
            __int32 x = 0;
            bool b = lexicalCast2<__int32>(x, buffer);
            (void)b;
            SEQAN_ASSERT(b);
            char const * ptr = reinterpret_cast<char const *>(&x);
            for (int i = 0; i < 4; ++i, ++ptr)
                appendValue(target, *ptr);
        }
        break;
    case 'f':
        {
            clear(buffer);
            res = readUntilTabOrLineBreak(buffer, reader);
            SEQAN_ASSERT(res == 0 || res == EOF_BEFORE_SUCCESS);
            float x = 0;
            bool b = lexicalCast2<float>(x, buffer);
            (void)b;
            SEQAN_ASSERT(b);
            char const * ptr = reinterpret_cast<char const *>(&x);
            for (int i = 0; i < 4; ++i, ++ptr)
                appendValue(target, *ptr);
        }
        break;
    case 'H':
    case 'Z':
        {
            // TODO(holtgrew): Could test on even length in case of 'H'.
            res = readUntilTabOrLineBreak(target, reader);
            SEQAN_ASSERT(res == 0 || res == EOF_BEFORE_SUCCESS);
            appendValue(target, '\0');
        }
        break;
    case 'B':
        {
            CharString buffer2; // TODO(holtgrew): Also give from outside.
            
            // Read type.
            clear(buffer);
            res = readNChars(buffer, reader, 1);
            SEQAN_ASSERT_EQ(res, 0);
            char t2 = back(buffer);
            appendValue(target, t2);
            
            // Read array contents.
            clear(buffer);
            res = readUntilTabOrLineBreak(buffer, reader);
            SEQAN_ASSERT(res == 0 || res == EOF_BEFORE_SUCCESS);
            typename Iterator<CharString, Rooted>::Type it, it2;
            // Search first non-comma position.
            it = begin(buffer, Rooted());
            for (;!atEnd(it) && *it == ','; ++it)
                continue;
            // Count number of entries.
            __int32 nEntries = !atEnd(it);  // At least one if array not empty.
            for (it2 = it; !atEnd(it2); ++it2)
                nEntries += (*it2 == ',');
            // Write out array length to result.
            char const * ptr = reinterpret_cast<char const *>(&nEntries);
            for (int i = 0; i < 4; ++i, ++ptr)
                appendValue(target, *ptr);
                
            // Now, write out the arrays, depending on the entry type.
            // TODO(holtgrew): Whee, this could be a bit more compact...
            switch (t2)
            {
            case 'c':
                for (int i = 0; i < nEntries; ++i)
                {
                    SEQAN_ASSERT_NOT(atEnd(it));
                    clear(buffer2);
                    for (; !atEnd(it) && *it != ',' && *it != '\t'; goNext(it))
                        appendValue(buffer2, *it);
                    __int16 x = 0;  // short to avoid textual interpretation in lexicalCast<> below.
                    bool b = lexicalCast2<__int16>(x, buffer2);
                    (void)b;
                    SEQAN_ASSERT(b);
                    appendValue(target, static_cast<__int8>(x));
                    if (!atEnd(it) && *it == ',')
                        goNext(it);  // Skip ','.
                    else
                        break;  // End of field or end of string.
                }
                break;
            case 'C':
                for (int i = 0; i < nEntries; ++i)
                {
                    SEQAN_ASSERT_NOT(atEnd(it));
                    clear(buffer2);
                    for (; !atEnd(it) && *it != ',' && *it != '\t'; goNext(it))
                        appendValue(buffer2, *it);
                    __uint16 x = 0;  // short to avoid textual interpretation in lexicalCast<> below.
                    bool b = lexicalCast2<__uint16>(x, buffer2);
                    (void)b;
                    SEQAN_ASSERT(b);
                    appendValue(target, static_cast<__int8>(x));
                    if (!atEnd(it) && *it == ',')
                        goNext(it);  // Skip ','.
                    else
                        break;  // End of field or end of string.
                }
                break;
            case 's':
                for (int i = 0; i < nEntries; ++i)
                {
                    SEQAN_ASSERT_NOT(atEnd(it));
                    clear(buffer2);
                    for (; !atEnd(it) && *it != ',' && *it != '\t'; goNext(it))
                        appendValue(buffer2, *it);
                    __int16 x = 0;
                    bool b = lexicalCast2<__int16>(x, buffer2);
                    (void)b;
                    SEQAN_ASSERT(b);
                    char const * ptr = reinterpret_cast<char const *>(&x);  // write out byte-wise
                    for (int i = 0; i < 2; ++i, ++ptr)
                        appendValue(target, *ptr);
                    if (!atEnd(it) && *it == ',')
                        goNext(it);  // Skip ','.
                    else
                        break;  // End of field or end of string.
                }
                break;
            case 'S':
                for (int i = 0; i < nEntries; ++i)
                {
                    SEQAN_ASSERT_NOT(atEnd(it));
                    clear(buffer2);
                    for (; !atEnd(it) && *it != ',' && *it != '\t'; goNext(it))
                        appendValue(buffer2, *it);
                    __uint16 x = 0;
                    bool b = lexicalCast2<__uint16>(x, buffer2);
                    (void)b;
                    SEQAN_ASSERT(b);
                    char const * ptr = reinterpret_cast<char const *>(&x);  // write out byte-wise
                    for (int i = 0; i < 2; ++i, ++ptr)
                        appendValue(target, *ptr);
                    if (!atEnd(it) && *it == ',')
                        goNext(it);  // Skip ','.
                    else
                        break;  // End of field or end of string.
                }
                break;
            case 'i':
                for (int i = 0; i < nEntries; ++i)
                {
                    SEQAN_ASSERT_NOT(atEnd(it));
                    clear(buffer2);
                    for (; !atEnd(it) && *it != ',' && *it != '\t'; goNext(it))
                        appendValue(buffer2, *it);
                    __int32 x = 0;
                    bool b = lexicalCast2<__int32>(x, buffer2);
                    (void)b;
                    SEQAN_ASSERT(b);
                    char const * ptr = reinterpret_cast<char const *>(&x);  // write out byte-wise
                    for (int i = 0; i < 4; ++i, ++ptr)
                        appendValue(target, *ptr);
                    if (!atEnd(it) && *it == ',')
                        goNext(it);  // Skip ','.
                    else
                        break;  // End of field or end of string.
                }
                break;
            case 'I':
                for (int i = 0; i < nEntries; ++i)
                {
                    SEQAN_ASSERT_NOT(atEnd(it));
                    clear(buffer2);
                    for (; !atEnd(it) && *it != ',' && *it != '\t'; goNext(it))
                        appendValue(buffer2, *it);
                    __uint32 x = 0;
                    bool b = lexicalCast2<__uint32>(x, buffer2);
                    (void)b;
                    SEQAN_ASSERT(b);
                    char const * ptr = reinterpret_cast<char const *>(&x);  // write out byte-wise
                    for (int i = 0; i < 4; ++i, ++ptr)
                        appendValue(target, *ptr);
                    if (!atEnd(it) && *it == ',')
                        goNext(it);  // Skip ','.
                    else
                        break;  // End of field or end of string.
                }
                break;
            case 'f':
                for (int i = 0; i < nEntries; ++i)
                {
                    SEQAN_ASSERT_NOT(atEnd(it));
                    clear(buffer2);
                    for (; !atEnd(it) && *it != ',' && *it != '\t'; goNext(it))
                        appendValue(buffer2, *it);
                    float x = 0;
                    bool b = lexicalCast2<float>(x, buffer2);
                    (void)b;
                    SEQAN_ASSERT(b);
                    char const * ptr = reinterpret_cast<char const *>(&x);  // write out byte-wise
                    for (int i = 0; i < 4; ++i, ++ptr)
                        appendValue(target, *ptr);
                    if (!atEnd(it) && *it == ',')
                        goNext(it);  // Skip ','.
                    else
                        break;  // End of field or end of string.
                }
                break;
            default:
                SEQAN_FAIL("Invalid array type: %c!", t2);
            }
        }
        break;
    default:
        SEQAN_ASSERT_FAIL("Invalid tag type: %c!", t);
    }
}

/**
.Function.assignTagsSamToBam
..cat:BAM I/O
..summary:Assign tags in SAM format to tags in BAM format.
..signature:assignTagsSamToBam(bamTags, samTags)
..param.bamTags:Destination BAM tags.
...type:Shortcut.CharString
..param.samTags:Source SAM tags.
...type:Shortcut.CharString
..returns:$void$
..include:seqan/bam_io.h
*/

template <typename TTarget, typename TSource>
void assignTagsSamToBam(TTarget & target, TSource & source)
{
    // Handle case of empty source sequence.
    if (empty(source))
        clear(target);

    typedef typename Iterator<TSource, Standard>::Type TSourceIter;
    TSourceIter it = begin(source, Standard());
    TSourceIter itEnd = end(source, Standard());
    
    typedef Stream<CharArray<char const *> > TStream;
    typedef RecordReader<TStream, SinglePass<> > TRecordReader;
    
    TStream stream(it, itEnd);
    TRecordReader reader(stream);
    
    CharString buffer;

    while (!atEnd(reader))
    {
        if (value(reader) == '\t')
            goNext(reader);
        SEQAN_ASSERT_NOT(atEnd(reader));
        
        _assignTagsSamToBamOneTag(target, reader, buffer);
    }
}

// ----------------------------------------------------------------------------
// Function assignTagsBamToSam()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TSourceIter>
void _assignTagsBamToSamOneTag(TTarget & target, TSourceIter & it, std::stringstream & ss)
{
    // Copy tag name.
    SEQAN_ASSERT_NOT(atEnd(it));
    appendValue(target, *it++);
    SEQAN_ASSERT_NOT(atEnd(it));
    appendValue(target, *it++);
    unsigned char t = *it;
    
    // Add ':'.
    appendValue(target, ':');
    
    // Add type.
    SEQAN_ASSERT_NOT(atEnd(it));
    if (*it == 'c' || *it == 'C' || *it == 's' || *it == 'S' || *it == 'i' || *it == 'I')
        appendValue(target, 'i');
    else
        appendValue(target, *it);
    ++it;

    // Add ':'.
    appendValue(target, ':');

    // Convert the payload, depending on the field's type.
    
    ss.str("");
    
    switch (t)
    {
    case 'A':
        appendValue(target, *it++);
        break;
    case 'c':
        {
            SEQAN_ASSERT_NOT(atEnd(it));
            __int8 x = *it++;
            ss << static_cast<int>(x);  // Cast to prevent printing as textual char.
            append(target, ss.str());
        }
        break;
    case 'C':
        {
            SEQAN_ASSERT_NOT(atEnd(it));
            __uint8 x = *it++;
            ss << static_cast<unsigned>(x);  // Cast to prevent printing as textual char.
            append(target, ss.str());
        }
        break;
    case 's':
        {
            __int16 x = 0;
            char * ptr = reinterpret_cast<char *>(&x);
            for (int i = 0; i < 2; ++i)
            {
                SEQAN_ASSERT_NOT(atEnd(it));
                *ptr++ = *it++;
            }
            ss << x;
            append(target, ss.str());
        }
        break;
    case 'S':
        {
            __uint16 x = 0;
            char * ptr = reinterpret_cast<char *>(&x);
            for (int i = 0; i < 2; ++i)
            {
                SEQAN_ASSERT_NOT(atEnd(it));
                *ptr++ = *it++;
            }
            ss << x;
            append(target, ss.str());
        }
        break;
    case 'i':
        {
            __int32 x = 0;
            char * ptr = reinterpret_cast<char *>(&x);
            for (int i = 0; i < 4; ++i)
            {
                SEQAN_ASSERT_NOT(atEnd(it));
                *ptr++ = *it++;
            }
            ss << x;
            append(target, ss.str());
        }
        break;
    case 'I':
        {
            __uint32 x = 0;
            char * ptr = reinterpret_cast<char *>(&x);
            for (int i = 0; i < 4; ++i)
            {
                SEQAN_ASSERT_NOT(atEnd(it));
                *ptr++ = *it++;
            }
            ss << x;
            append(target, ss.str());
        }
        break;
    case 'f':
        {
            float x = 0;
            char * ptr = reinterpret_cast<char *>(&x);
            for (int i = 0; i < 4; ++i)
            {
                SEQAN_ASSERT_NOT(atEnd(it));
                *ptr++ = *it++;
            }
            ss << x;
            append(target, ss.str());
        }
        break;
    case 'Z':
        {
            while (*it != '\0')
            {
                SEQAN_ASSERT_NOT(atEnd(it));
                appendValue(target, *it++);
            }
            SEQAN_ASSERT_NOT(atEnd(it));
            it++;
        }
        break;
    case 'H':
        {
            while (*it != '\0')
            {
                SEQAN_ASSERT_NOT(atEnd(it));
                appendValue(target, *it++);
            }
            SEQAN_ASSERT_NOT(atEnd(it));
            it++;
        }
        break;
    case 'B':
        {
            // Read type.
            char t2 = *it++;
            appendValue(target, t2);
            // Read array length.
            __int32 x = 0;
            char * ptr = reinterpret_cast<char *>(&x);
            for (int i = 0; i < 4; ++i)
            {
                SEQAN_ASSERT_NOT(atEnd(it));
                *ptr++ = *it++;
            }
            // Depending on t2, read array.
            // TODO(holtgrew): Whee, this could be a bit more compact...
            switch (t2)
            {
            case 'c':
                for (__int32 i = 0; i < x; ++i)
                {
                    appendValue(target, ',');
                    ss.str("");
                    __int8 y = *it++;
                    ss << static_cast<int>(y);  // Cast to prevent textual interpretation.
                    append(target, ss.str());
                }
                break;
            case 'C':
                for (__int32 i = 0; i < x; ++i)
                {
                    appendValue(target, ',');
                    ss.str("");
                    __uint8 y = *it++;
                    ss << static_cast<int>(y);  // Cast to prevent textual interpretation.
                    append(target, ss.str());
                }
                break;
            case 's':
                for (__int32 i = 0; i < x; ++i)
                {
                    appendValue(target, ',');
                    ss.str("");
                    __int16 y = 0;
                    char * ptr = reinterpret_cast<char *>(&y);
                    for (int i = 0; i < 2; ++i)
                    {
                        SEQAN_ASSERT_NOT(atEnd(it));
                        *ptr++ = *it++;
                    }
                    ss << y;
                    append(target, ss.str());
                }
                break;
            case 'S':
                for (__int32 i = 0; i < x; ++i)
                {
                    appendValue(target, ',');
                    ss.str("");
                    __uint16 y = 0;
                    char * ptr = reinterpret_cast<char *>(&y);
                    for (int i = 0; i < 2; ++i)
                    {
                        SEQAN_ASSERT_NOT(atEnd(it));
                        *ptr++ = *it++;
                    }
                    ss << y;
                    append(target, ss.str());
                }
                break;
            case 'i':
                for (__int32 i = 0; i < x; ++i)
                {
                    appendValue(target, ',');
                    ss.str("");
                    __int32 y = 0;
                    char * ptr = reinterpret_cast<char *>(&y);
                    for (int i = 0; i < 4; ++i)
                    {
                        SEQAN_ASSERT_NOT(atEnd(it));
                        *ptr++ = *it++;
                    }
                    ss << y;
                    append(target, ss.str());
                }
                break;
            case 'I':
                for (__int32 i = 0; i < x; ++i)
                {
                    appendValue(target, ',');
                    ss.str("");
                    __uint32 y = 0;
                    char * ptr = reinterpret_cast<char *>(&y);
                    for (int i = 0; i < 4; ++i)
                    {
                        SEQAN_ASSERT_NOT(atEnd(it));
                        *ptr++ = *it++;
                    }
                    ss << y;
                    append(target, ss.str());
                }
                break;
            case 'f':
                for (__int32 i = 0; i < x; ++i)
                {
                    appendValue(target, ',');
                    ss.str("");
                    float y = 0;
                    char * ptr = reinterpret_cast<char *>(&y);
                    for (int i = 0; i < 4; ++i)
                    {
                        SEQAN_ASSERT_NOT(atEnd(it));
                        *ptr++ = *it++;
                    }
                    ss << y;
                    append(target, ss.str());
                }
                break;
            default:
                SEQAN_FAIL("Invalid array type: %c!", t2);
            }
        }
        break;
    default:
        SEQAN_ASSERT_FAIL("Invalid tag type: %c!", t);
    }
}

/**
.Function.assignTagsBamToSam
..cat:BAM I/O
..summary:Assign tags in BAM format to tags in SAM format.
..signature:assignTagsSamToBam(bamTags, samTags)
..param.samTags:Destination SAM tags.
...type:Shortcut.CharString
..param.bamTags:Source BAM tags.
...type:Shortcut.CharString
..returns:$void$
..include:seqan/bam_io.h
..see:Function.assignTagsSamToBam
*/

template <typename TTarget, typename TSource>
void assignTagsBamToSam(TTarget & target, TSource const & source)
{
    // Handle case of empty source sequence.
    if (empty(source))
        clear(target);

    std::stringstream ss;
    clear(target);
    
    typedef typename Iterator<TSource const, Rooted>::Type TSourceIter;
    TSourceIter it = begin(source, Rooted());
    
    bool first = true;
    while (!atEnd(it))
    {
        if (!first)
            appendValue(target, '\t');
        first = false;
        _assignTagsBamToSamOneTag(target, it, ss);
    }
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BAM_IO_SAM_TAGS_TO_BAM_TAGS_H_
