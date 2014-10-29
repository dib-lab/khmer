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
// Code for read/write access to BAM tag dicts.
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_BAM_IO_BAM_TAGS_DICT_H_
#define CORE_INCLUDE_SEQAN_BAM_IO_BAM_TAGS_DICT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Class.BamTagsDict
..cat:BAM I/O
..cat:Fragment Store
..signature:BamTagsDict
..summary:Indexes start positions of BAM tags in a @Shortcut.CharString@ and provides a dict-like API.
..example.code:
CharString samStr = "AA:Z:value1\tAB:Z:value2\tAC:i:30";
CharString bamStr;
assignSamToBam(bamStr, samStr);
BamTagsDict tags(bamStr);
std::cerr << length(tags) << std::endl;  // #=> "3"
for (unsigned i = 0; i < length(tags); ++i)
{
    std::cerr << getTagKey(tags, i) << " -> " << getTagValue(tags, i) << std::endl;
    if (getTagValue(tags, i)[0] == 'i')  // is 32 bit integer
    {
        __int32 x = 0;
        bool res = extractTagValue(x, tags, i);
        SEQAN_ASSERT_MSG(res, "Not a valid integer at pos %u!", i);
        std::cerr << "     " << x << std::endl;
    }
}
// #=> "AA -> Zvalue1"
// #=> "AB -> Zvalue2"
// #-> "AC -> i<binary representation of 30>"
#  #-> "      30"
..include:seqan/bam_io.h

.Memfunc.BamTagsDict#BamTagsDict
..class:Class.BamTagsDict
..signature:BamTagsDict()
..summary:Constructor.
..remarks:Only the default constructor is provided.
*/

class BamTagsDict
{
public:
    Holder<CharString> _host;
    String<unsigned> _positions;

    BamTagsDict() {}

    BamTagsDict(CharString & tags) : _host(tags) {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Host
// ----------------------------------------------------------------------------

template <>
struct Host<BamTagsDict>
{
    typedef CharString Type;
};

template <>
struct Host<BamTagsDict const>
{
    typedef CharString const Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function host()
// ----------------------------------------------------------------------------

inline Host<BamTagsDict>::Type &
host(BamTagsDict & bamTags)
{
    return value(bamTags._host);
}

inline Host<BamTagsDict const>::Type &
host(BamTagsDict const & bamTags)
{
    return value(bamTags._host);
}

// ----------------------------------------------------------------------------
// Function hasIndex()
// ----------------------------------------------------------------------------

/**
.Function.BamTagsDict#hasIndex
..class:Class.BamTagsDict
..cat:Fragment Store
..summary:Return $true$ if @Class.BamTagsDict@ has an index.
..signature:hasIndex(bamTags)
..param.bamTags:SAM Tags to query
...type:Class.BamTagsDict
..returns:$bool$
..include:<seqan/store_ex.h>
*/

inline bool
hasIndex(BamTagsDict const & bamTags)
{
    return length(bamTags._positions) != 0u;
}

inline bool
hasIndex(BamTagsDict & bamTags)
{
    return hasIndex(const_cast<BamTagsDict const &>(bamTags));
}

// ----------------------------------------------------------------------------
// Function getBamTypeSize()
// ----------------------------------------------------------------------------

// Return sizeof() of the type identified with the given char.  Returns -2 if not
// valid, -1 if of variable length.

/**
.Function.getBamTypeSize
..class:Class.BamTagsDict
..cat:BAM I/O
..signature:getBamTypeSize(c)
..summary:Return size of the type identified by $c$.
..param.c:The BAM type identifier
..returns:$int$ with the $sizeof()$ of the type, -1 for variable sized types, -2 for invalid parameters.
..include:seqan/bam_io.h
*/

inline int
getBamTypeSize(char c)
{
    switch (c)
    {
        case 'A':
            return 1;
        case 'f':
            return 4;
        case 'Z':
        case 'H':
        case 'B':
            return -1;
        case 'c':
        case 'C':
            return 1;
        case 's':
        case 'S':
            return 2;
        case 'i':
        case 'I':
            return 4;
    }
    return -2;
}

// ----------------------------------------------------------------------------
// Function buildIndex()
// ----------------------------------------------------------------------------

/**
.Function.BamTagsDict#buildIndex
..class:Class.BamTagsDict
..cat:Fragment Store
..summary:Build index for a @Class.BamTagsDict@ object.
..signature:buildIndex(bamTags)
..param.bamTags:SAM Tags to build index for.
...type:Class.BamTagsDict
..returns:$void$
..include:<seqan/bam_io.h>
*/

inline void
buildIndex(BamTagsDict & bamTags)
{
    typedef Host<BamTagsDict>::Type TCharString;
    typedef Iterator<TCharString, Rooted>::Type TCharStringIter;

    clear(bamTags._positions);
    if (empty(value(bamTags._host)))
        return;  // Done.

    appendValue(bamTags._positions, 0);
    for (TCharStringIter it = begin(host(bamTags)); !atEnd(it);)
    {
        it += 2;
        char c = *it;
        if (c == 'H' || c == 'Z')
        {
            while (!atEnd(it) && *it != '\0')
                ++it;
            ++it;
        }
        else if (c == 'B')
        {
            ++it;
            c = *it;
            ++it;
            __uint32 len = 0;
            memcpy(&len, &*it, 4);
            it += 4;
            it += len * getBamTypeSize(c);
        }
        else
        {
            ++it;
            it += getBamTypeSize(c);
        }

        appendValue(bamTags._positions, position(it));
    }
    // if (!empty(value(bamTags._host)))
    //     appendValue(bamTags._positions, length(host(bamTags)) + 1);  // +1 since there is not tab at the end
}

// ----------------------------------------------------------------------------
// Function setHost()
// ----------------------------------------------------------------------------

inline Holder<CharString> &
_dataHost(BamTagsDict & bamTags)
{
    return bamTags._host;
}

inline void
setHost(BamTagsDict & me, CharString & host_)
{
    SEQAN_CHECKPOINT;
	setValue(_dataHost(me), host_);
    clear(me._positions);
}

inline void
setHost(BamTagsDict & me, CharString const & host_)
{
    SEQAN_CHECKPOINT;
	setValue(_dataHost(me), host_);
    clear(me._positions);
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

///.Function.length.param.object.type:Class.BamTagsDict

inline unsigned
length(BamTagsDict const & tags)
{
    if (empty(value(tags._host)))
        return 0;
    if (!hasIndex(tags))
        buildIndex(const_cast<BamTagsDict &>(tags));
    return length(tags._positions) - 1;
}

// ----------------------------------------------------------------------------
// Function getTagType()
// ----------------------------------------------------------------------------

/**
.Function.BamTagsDict#getTagType
..class:Class.BamTagsDict
..cat:BAM I/O
..signature:getTagType(tagsDict, idx)
..summary:Get key of a tag by index.
..param.tagsDict:The @Class.BamTagsDict@ to retrieve data from.
..param.idx:Index of the tag whose key to retrieve.
..returns:$char$, the SAM/BAM identifier of the type.
..include:seqan/bam_io.h
*/

template <typename TPos>
inline char
getTagType(BamTagsDict & tags, TPos idx)
{
    if (!hasIndex(tags))
        buildIndex(tags);
    return host(tags)[tags._positions[idx] + 2];
}

// ----------------------------------------------------------------------------
// Function getTagKey()
// ----------------------------------------------------------------------------

/**
.Function.BamTagsDict#getTagKey
..class:Class.BamTagsDict
..cat:BAM I/O
..signature:getTagKey(tagsDict, idx)
..summary:Return key of a tag by index.
..param.tagsDict:The @Class.BamTagsDict@ to retrieve data from.
...type:Class.BamTagsDict
..param.idx:Index of the tag whose key to retrieve.
..returns:Infix of the underlying string.
..remarks:See @Class.BamTagsDict@ for an example.
..include:seqan/bam_io.h
*/

template <typename TPos>
inline Infix<Host<BamTagsDict>::Type>::Type
getTagKey(BamTagsDict & tags, TPos idx)
{
    if (!hasIndex(tags))
        buildIndex(tags);
    return infix(host(tags), tags._positions[idx], tags._positions[idx] + 2);
}

template <typename TPos>
inline Infix<Host<BamTagsDict const>::Type>::Type
getTagKey(BamTagsDict const & tags, TPos idx)
{
    return getTagKey(const_cast<BamTagsDict &>(tags), idx);
}

// ----------------------------------------------------------------------------
// Function findTagKey()
// ----------------------------------------------------------------------------

/**
.Function.BamTagsDict#findTagKey
..summary:Find a tag by its key for a @Class.BamTagsDict@ object.
..class:Class.BamTagsDict
..signature:findTagKey(idx, tagsDict, name)
..param.idx:Index of the tag with the given key.
...type:nolink:$unsigned$
..param.tagsDict:The @Class.BamTagsDict@ to retrieve data from.
..param.name:Name of the key to find.
...type:Shortcut.CharString
..returns:$bool$, indicating whether such a key could be found.
..include:seqan/bam_io.h
*/

inline bool
findTagKey(unsigned & idx, BamTagsDict & tags, CharString const & name)
{
    for (idx = 0; idx < length(tags); ++idx)
        if (getTagKey(tags, idx) == name)
            return true;
    return false;
}

inline bool
findTagKey(unsigned & idx, BamTagsDict const & tags, CharString const & name)
{
    return findTagKey(idx, const_cast<BamTagsDict &>(tags), name);
}

// ----------------------------------------------------------------------------
// Function getTagValue()
// ----------------------------------------------------------------------------

/**
.Function.BamTagsDict#getTagValue
..class:Class.BamTagsDict
..cat:BAM I/O
..summary:Return the value of a tag by its index in the @Class.BamTagsDict@.
..signature:getTagValue(tagsDict, idx)
..param.tagsDict:The @Class.BamTagsDict@ to retrieve data from.
...type:Class.BamTagsDict
..param.idx:Index of the tag whose value to retrieve.
..returns:@Shortcut.CharString@ with the raw tags data.
..remarks:Note that you will get $<type char> + payload$ in case of @Class.BamTagsDict@.
..remarks:See @Class.BamTagsDict@ for an example.
..include:seqan/bam_io.h
*/

template <typename TIdx>
inline CharString
getTagValue(BamTagsDict & tags, TIdx idx)
{
    if (!hasIndex(tags))
        buildIndex(tags);

    // TODO(holtgrew): Can't we use positions to speed this up?
    
    typedef typename Position<CharString>::Type TPos;
    TPos beginPos = tags._positions[idx] + 2;
    TPos endPos = beginPos + 1;
    
    char theType = getTagType(tags, idx);
    if (theType == 'Z' || theType == 'H')
    {
        typedef typename Iterator<CharString, Rooted>::Type TIterator;
        TIterator it = begin(host(tags), Rooted()) + beginPos + 1;
        for (; !atEnd(it) && *it != '\0'; goNext(it))
            endPos += 1;
        endPos += 1;
    }
    else if (theType == 'B')
    {
        __uint32 len = 0;
        memcpy(&len, &host(tags)[tags._positions[idx]] + 4, 4);
        char c = host(tags)[tags._positions[idx] + 3];
        int typeSize = getBamTypeSize(c);
        SEQAN_ASSERT_GT(typeSize, 0);
        endPos += 5 + len * typeSize;
    }
    else
    {
        endPos += getBamTypeSize(theType);
    }
    
    return infix(host(tags), beginPos, endPos);
}

template <typename TPos>
inline CharString //Infix<Host<BamTagsDict const>::Type>::Type
getTagValue(BamTagsDict const & tags, TPos idx)
{
    return getValue(const_cast<BamTagsDict &>(tags), idx);
}

// ----------------------------------------------------------------------------
// Function extractTagValue()
// ----------------------------------------------------------------------------

/**
.Function.BamTagsDict#extractTagValue
..class:Class.BamTagsDict
..cat:BAM I/O
..signature:extractTagValue(dest, tags, idx)
..summary:Extract and cast "atomic" value from tags string with index $idx$.
..param.dest:The variable to write the value to.
...remarks:The value is first copied in a variable of the type indicated in the BAM file. Then it is cast into the type of $dest$.
..param.tags:@Class.BamTagsDict@ object.
...type:Class.BamTagsDict
..params.idx:Index of the tag in the tag list.
..returns:$bool$, indicating the success.
..remarks:The function only works for atomic types such as $int$, not for $char*$ or arrays.
..remarks:See @Class.BamTagsDict@ for an example.
..see:Function.BamTagsDict#getTagValue
..include:seqan/bam_io.h
*/

template <typename TDest, typename TIdx>
inline bool
extractTagValue(TDest & dest, BamTagsDict & tags, TIdx idx)
{
    if (!hasIndex(tags))
        buildIndex(tags);

    char typeC = host(tags)[tags._positions[idx] + 2];
    if (typeC == 'A')
    {
        char x = 0;
        char * ptr = reinterpret_cast<char *>(&x);
        memcpy(ptr, &host(tags)[tags._positions[idx] + 3], 1);
        dest = static_cast<TDest>(x);
    }
    else if (typeC == 'c')
    {
        __int8 x = 0;
        char * ptr = reinterpret_cast<char *>(&x);
        memcpy(ptr, &host(tags)[tags._positions[idx] + 3], 1);
        dest = static_cast<TDest>(x);
    }
    else if (typeC == 'C')
    {
        __uint8 x = 0;
        char * ptr = reinterpret_cast<char *>(&x);
        memcpy(ptr, &host(tags)[tags._positions[idx] + 3], 1);
        dest = static_cast<TDest>(x);
    }
    else if (typeC == 's')
    {
        __int16 x = 0;
        char * ptr = reinterpret_cast<char *>(&x);
        memcpy(ptr, &host(tags)[tags._positions[idx] + 3], 2);
        dest = static_cast<TDest>(x);
    }
    else if (typeC == 'S')
    {
        __uint16 x = 0;
        char * ptr = reinterpret_cast<char *>(&x);
        memcpy(ptr, &host(tags)[tags._positions[idx] + 3], 2);
        dest = static_cast<TDest>(x);
    }
    else if (typeC == 'i')
    {
        __int32 x = 0;
        char * ptr = reinterpret_cast<char *>(&x);
        memcpy(ptr, &host(tags)[tags._positions[idx] + 3], 4);
        dest = static_cast<TDest>(x);
    }
    else if (typeC == 'I')
    {
        __uint32 x = 0;
        char * ptr = reinterpret_cast<char *>(&x);
        memcpy(ptr, &host(tags)[tags._positions[idx] + 3], 4);
        dest = static_cast<TDest>(x);
    }
    else if (typeC == 'f')
    {
        float x = 0;
        char * ptr = reinterpret_cast<char *>(&x);
        memcpy(ptr, &host(tags)[tags._positions[idx] + 3], 4);
        dest = static_cast<TDest>(x);
    }
    else // variable sized type or invald
    {
        return false;
    }
    return true;
}


// ----------------------------------------------------------------------------
// Function getBamTypeChar()
// ----------------------------------------------------------------------------

/**
.Function.getBamTypeChar
..class:Class.BamTagsDict
..cat:BAM I/O
..summary:Return char identifying the type of the atomic argument.
..signature:getBamTypeChar<T>()
..param.T:The type to get the BAM char for.
..returns:$char$ describing the BAM type. One of $ACcSsIifZ$.
..remarks:Note that this function is defined for the $__int16$, $__uint16$ etc. but not for the types $short$, $int$ etc. An exception are 8-bit characters/char, where it is defined for $__int8$, $__uint8$, and $char$ unless $char$ is equal to one of the other two types. This is important when used in @Function.BamTagsDict#setTagValue@ etc. since BAM gives type chars for printable characters, signed 8-bit numbers and unsigned 8-bit numbers.
..remarks:If $__int8$ and $__uint8$ are not identical to $char$, we can make this decision from the type, otherwise we cannot and we will give the integer types a higher precedence.
..remarks:In your programs, this should not make any difference, only the written SAM/BAM will differ.
..include:seqan/bam_io.h
*/

template <typename T>
inline char getBamTypeChar()
{
	if (IsSameType<T, __int8>::Type::VALUE)
		return 'C';
	if (IsSameType<T, __uint8>::Type::VALUE)
		return 'c';
	if (IsSameType<T, char>::Type::VALUE)
		return 'A';
	if (IsSameType<T, __int16>::Type::VALUE)
		return 's';
	if (IsSameType<T, __uint16>::Type::VALUE)
		return 'S';
	if (IsSameType<T, __int32>::Type::VALUE)
		return 'i';
	if (IsSameType<T, __int32>::Type::VALUE)
		return 'I';
	if (IsSameType<T, float>::Type::VALUE)
		return 'f';
    if (IsSameType<T, char *>::Type::VALUE || IsSameType<T, char const *>::Type::VALUE)
        return 'Z';
	else
		return '\0';
}

// ----------------------------------------------------------------------------
// Function setTagValue()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Test me!

/**
.Function.BamTagsDict#setTagValue
..class:Class.BamTagsDict
..cat:BAM I/O
..summary:Set the value of a tag through a @Class.BamTagsDict@.
..signature:setTagValue(tags, key, val[, typeC])
..param.tags:The dict to modify.
...type:Class.BamTagsDict
..param.key:The key of the tag.
...type:Shortcut.CharString
...remarks:Must be a string of length 2.
..param.val:The value to set the the tag to.
..param.typeC:BAM type char to use.
...type:nolink:By default, the type is inflected using @Function.getBamTypeChar@.
...remarks:For portability (so the generated files are the same on all platforms), use a signed/unsigned qualified type for $val$ or give $typeC$. Also see the remarks for @Function.getBamTypeChar@.
..returns:$bool$ indicating the success. This function can fail if the key is not a valid tag id (e.g. does not have length 2) or if the type of $val$ is not an atomic value or a string (anything but $char *$, $char const *$, a character, integer or float type is invalid).
..see:Function.getBamTypeChar
..remarks:Note that $setTagValue$ does not cast the type, so $typeC$ only influences the type character written out but $val$ is written out in binary without modification.
..include:seqan/bam_io.h
..example.text:An example setting some atomic tag values.
..example.code:
CharString rawTagsText;
BamTagsDict tags(rawTagsText);
setTagValue(tags, "XA", 9);    // int
setTagValue(tags, "XB", 9u);   // unsigned int
setTagValue(tags, "XC", 'X');  // char
..example.text:If $char$ is equal to $__int8$ or $__uint8$ then the last line produces an entry with type 'c' or 'C'. To make sure that the type char 'A' (for "printable character") is written to the file, give it explicitely:
..example.code:
setTagValue(tags, "XC", 'X', 'A');  // Overrwrite XC, enforce type 'printable character'.
..example.text:Note that on most systems $int$s have a width of 32 bytes, but the C++ standard leaves this open. For all types but characters, you should not give an explicit type char but use one of the types with explicit width and signed/unsigned qualifier such as $__int32$, $__uint32$ etc.
..example.code:
// The following is not recommended since the type of $x$ is not "unsigned 32 bit int."
__int32 x = -1;
setTagValue(tags, "XB", x, 'I');
// Instead, explicitely use an unsigned type if you need one.  Note that your compiler
// might warn you about assigning -1 to an unsigned variable so you know that you are
// probably doing something unintended.
__uint32 y = -1;
setTagValue(tags, "XB", y);

// Do not do this!
setTagValue(tags, "XA", 9, 'f');    // BOGUS since 9 is not a floating point number.
*/

// Convert "atomic" value to BAM tag.  Return whether val was atomic.
template <typename T>
bool _toBamTagValue(CharString & result, T const & val, char const typeC)
{
    appendValue(result, typeC);

    if (typeC == 'A' || typeC == 'c' || typeC == 'C')
    {
        resize(result, length(result) + 1);
        char * dst = reinterpret_cast<char *>(&result[0]) + length(result) - 1;
        char const * src = reinterpret_cast<char const *>(&val);
        memcpy(dst, src, 1);
    }
    else if (typeC == 's' || typeC == 'S')
    {
        resize(result, length(result) + 2);
        char * dst = reinterpret_cast<char *>(&result[0]) + length(result) - 2;
        char const * src = reinterpret_cast<char const *>(&val);
        memcpy(dst, src, 2);
    }
    else if (typeC == 'i' || typeC == 'I' || typeC == 'f')
    {
        resize(result, length(result) + 4);
        char * dst = reinterpret_cast<char *>(&result[0]) + length(result) - 4;
        char const * src = reinterpret_cast<char const *>(&val);
        memcpy(dst, src, 4);
    }
    else if (typeC == 'Z')
    {
        unsigned oldSize = length(result);
        unsigned valLen = length(val) + 1;
        resize(result, length(result) + valLen);
        char * dst = reinterpret_cast<char *>(&result[0] + oldSize);
        char const * src = reinterpret_cast<char const *>(val);
        memcpy(dst, src, valLen);
        *(dst + valLen - 1) = '\0';
    }
    else // non-string and variable sized type or invald
    {
        return false;
    }
    return true;
}

// Sets an atomic value in a BamTagsDict.
// Returns true successful, can fail if val not atomic or key is not a valid tag id (2 chars).

template <typename T>
inline bool
setTagValue(BamTagsDict & tags, CharString const & key, T const & val, char const typeC)
{
    if (!hasIndex(tags))
        buildIndex(tags);

    // Build value to insert/append.
    if (length(key) != 2u)
        return false;
    CharString bamTagVal;
    // append(bamTagVal, key);
    if (!_toBamTagValue(bamTagVal, val, typeC))
        return false;
    
    unsigned idx = 0;
    if (findTagKey(idx, tags, key))
    {
        // TODO(holtgrew): Speed this up with positions?
        CharString tmp;
        tmp = getTagValue(tags, idx);
        replace(host(tags), tags._positions[idx] + 2, tags._positions[idx] + 2 + length(tmp), bamTagVal);
    }
    else
    {
        append(host(tags), key);
        append(host(tags), bamTagVal);
    }

    // Remove index and return success.
    clear(tags._positions);  // Also necessary when appending?
    return true;
}

template <typename T>
inline bool
setTagValue(BamTagsDict & tags, CharString const & key, T const & val)
{
    return setTagValue(tags, key, val, getBamTypeChar<T>());
}

// ----------------------------------------------------------------------------
// Function eraseTag()
// ----------------------------------------------------------------------------

/**
.Function.BamTagsDict#eraseTag
..class:Class.BamTagsDict
..summary:Erase tag from @Class.BamTagsDict@.
..cat:BAM I/O
..signature:eraseTag(tagsDict, key)
..param.tags:The dict to erase from.
...type:Class.BamTagsDict
..param.key:The key of the entry to remove.
...type:Shortcut.CharString
..returns:$bool$, indicating whether the key was present.
..include:seqan/bam_io.h
 */

inline bool
eraseTag(BamTagsDict & tags, CharString const & key)
{
    if (!hasIndex(tags))
        buildIndex(tags);
    
    unsigned idx = 0;
    if (!findTagKey(idx, tags, key))
        return false;

    // TODO(holtgrew): Speed this up with positions?
    CharString tmp;
    tmp = getTagValue(tags, idx);
    erase(host(tags), tags._positions[idx], tags._positions[idx + 1]);
    
    return true;
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BAM_IO_BAM_TAGS_DICT_H_
