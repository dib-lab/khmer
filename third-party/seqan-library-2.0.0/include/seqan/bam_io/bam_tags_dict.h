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
// Code for read/write access to BAM tag dicts.
// ==========================================================================

#ifndef INCLUDE_SEQAN_BAM_IO_BAM_TAGS_DICT_H_
#define INCLUDE_SEQAN_BAM_IO_BAM_TAGS_DICT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

class BamTagsDict;

inline bool hasIndex(BamTagsDict const & bamTags);
inline void buildIndex(BamTagsDict const & bamTags);

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

template <>
struct Size<BamTagsDict>
{
    typedef unsigned Type;
};

template <>
struct Position<BamTagsDict>
{
    typedef unsigned Type;
};

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @class BamTagsDict
 * @headerfile <seqan/bam_io.h>
 * @brief Indexes start positions of BAM tags in a @link CharString @endlink and provides a dict-like API.
 *
 * @signature class BamTagsDict;
 *
 * @section Example
 *
 * @include demos/bam_io/bam_tags_dict.cpp
 *
 * Output is:
 *
 * @include demos/bam_io/bam_tags_dict.cpp.stdout
 *
 * @see getBamTypeSize
 * @see getBamTypeChar
 */

/*!
 * @fn BamTagsDict::BamTagsDict
 * @brief Constructor
 *
 * @signature BamTagsDict::BamTagsDict();
 */

class BamTagsDict
{
public:
    typedef Host<BamTagsDict>::Type TBamTagsSequence;
    typedef Position<TBamTagsSequence>::Type TPos;

    Holder<TBamTagsSequence> _host;
    mutable String<TPos> _positions;

    BamTagsDict() {}

    explicit
    BamTagsDict(TBamTagsSequence & tags) : _host(tags) {}

    template <typename TPos>
    inline Infix<Host<BamTagsDict const>::Type>::Type
    operator[] (TPos pos) const
    {
        if (!hasIndex(*this))
            buildIndex(*this);
        return infix(host(*this), _positions[pos], _positions[pos + 1]);
    }
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

/*!
 * @fn BamTagsDict#hasIndex
 * @brief Returns whether the BamTagsDict has an index.
 *
 * @signature bool hasIndex(dict);
 *
 * @param[in] dict The @link BamTagsDict @endlink to query.
 *
 * @return bool <tt>true</tt> if <tt>dict</tt> has an index and <tt>false</tt> otherwise.
 */

inline bool
hasIndex(BamTagsDict const & bamTags)
{
    return !empty(bamTags._positions) || empty(host(bamTags));
}

// ----------------------------------------------------------------------------
// Function getBamTypeSize()
// ----------------------------------------------------------------------------

/*!
 * @fn getBamTypeSize
 * @headerfile <seqan/bam_io.h>
 * @brief Return size of the type identified by a type char.
 *
 * @signature int getBamTypeSize(c);
 *
 * @param[in] c A <tt>char</tt> that identifies a type.
 *
 * @return int The size of the type in bytes, -1 for variable-length types, -2 for invalid paramters.
 *
 * @see BamTagsDict
 * @see getBamTypeChar
 */

// Return sizeof() of the type identified with the given char.  Returns -2 if not
// valid, -1 if of variable length.

struct GetBamTypeSizeHelper_
{
    int     &resultSize;
    char    typeC;

    GetBamTypeSizeHelper_(int &resultSize, char typeC) :
        resultSize(resultSize),
        typeC(typeC)
    {}

    template <typename Type>
    bool operator() (Type) const
    {
        if (BamTypeChar<Type>::VALUE != typeC)
            return false;

        resultSize = sizeof(Type);
        return true;
    }
};


inline int
getBamTypeSize(char c)
{
    switch (toUpperValue(c))
    {
        case 'Z':
        case 'H':
        case 'B':
            return -1;

        default:
            int result = -2;
            GetBamTypeSizeHelper_ func(result, c);
            tagApply(func, BamTagTypes());
            return result;
    }
}

// ----------------------------------------------------------------------------
// Function buildIndex()
// ----------------------------------------------------------------------------

/*!
 * @fn BamTagsDict#buildIndex
 * @brief Build index for a BamTagsDict  object.
 *
 * @signature void buildIndex(bamTags);
 *
 * @param[in,out] bamTags The BamTagsDict object to build the index for.
 */

inline void
buildIndex(BamTagsDict const & bamTags)
{
    typedef Host<BamTagsDict>::Type TTagString;
    typedef Iterator<TTagString const, Standard>::Type TIter;

    clear(bamTags._positions);
    if (empty(host(bamTags)))
        return;  // Done.

    appendValue(bamTags._positions, 0);
    TIter itBegin = begin(host(bamTags), Standard());
    TIter itEnd = end(host(bamTags), Standard());
    for (TIter it = itBegin; it != itEnd; )
    {
        SEQAN_ASSERT(it < itEnd);
        // skip tag name (e.g. "NM")
        it += 2;

        // get tag type (e.g. 'I')
        char c = *(it++);
        if (c == 'H' || c == 'Z')
        {
            // skip string and its end-of-string marker
            while (*it != '\0')
            {
                ++it;
                SEQAN_ASSERT(it != itEnd);
            }
            ++it;
        }
        else if (c == 'B')
        {
            // skip array of PODs
            c = *(it++);
            union {
                char raw[4];
                __uint32 len;
            } tmp;
            arrayCopyForward(it, it + 4, tmp.raw);
            it += 4 + tmp.len * getBamTypeSize(c);
        }
        else
        {
            // skip POD type (e.g. byte, int)
            it += getBamTypeSize(c);
        }
        appendValue(bamTags._positions, it - itBegin);
    }
    // if (!empty(value(bamTags._host)))
    //     appendValue(bamTags._positions, length(host(bamTags)) + 1);  // +1 since there is not tab at the end
}

// ----------------------------------------------------------------------------
// Function _dataHost()
// ----------------------------------------------------------------------------

inline Holder<Host<BamTagsDict>::Type> &
_dataHost(BamTagsDict & bamTags)
{
    return bamTags._host;
}

inline Holder<Host<BamTagsDict>::Type> const &
_dataHost(BamTagsDict const & bamTags)
{
    return bamTags._host;
}

// ----------------------------------------------------------------------------
// Function setHost()
// ----------------------------------------------------------------------------

template <typename THost>
inline void
setHost(BamTagsDict & me, THost & host_)
{
    SEQAN_CHECKPOINT;
    setValue(_dataHost(me), host_);
    clear(me._positions);
}

template <typename THost>
inline void
setHost(BamTagsDict & me, THost const & host_)
{
    SEQAN_CHECKPOINT;
    setValue(_dataHost(me), host_);
    clear(me._positions);
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

/*!
 * @fn BamTagsDict#length
 * @brief Returns the number of entries in a BamTagsDict.
 *
 * @signature unsigned length(tagsDict);
 *
 * @param[in] tagsDict The BamTagsDict object to query for its length.
 *
 * @return TSize The number of entries in the BamTagsDict. <tt>TSize</tt> is the result of
 *               <tt>Size&lt;BamTagsDict&gt;::Type</tt>.
 */

inline Size<BamTagsDict>::Type
length(BamTagsDict const & tags)
{
    if (empty(host(tags)))
        return 0;
    if (!hasIndex(tags))
        buildIndex(tags);
    return length(tags._positions) - 1;
}

// ----------------------------------------------------------------------------
// Function getTagType()
// ----------------------------------------------------------------------------

/*!
 * @fn BamTagsDict#getTagType
 * @brief Returns the tag type char for an entry of a BamTagsDict.
 *
 * @signature char getTagType(tags, id);
 *
 * @param[in] tags The BamTagsDict to query.
 * @param[in] id   The id of the tag for which to determine the type. See @link BamTagsDict#findTagKey @endlink.
 *
 * @return char A <tt>char</tt> that identifies the tag type.
 */

template <typename TId>
inline char
getTagType(BamTagsDict const & tags, TId id)
{
    return tags[id][2];
}

// ----------------------------------------------------------------------------
// Function getTagKey()
// ----------------------------------------------------------------------------

/*!
 * @fn BamTagsDict#getTagKey
 * @brief Return key of a tag by index.
 *
 * @signature TKey getTagKey(tagsDict, id);
 *
 * @param[in] tagsDict The BamTagsDict to query.
 * @param[in] id       The index of the dict entry.
 *
 * @return TKey An infix of a @link CharString @endlink.  Will be a two-character char sequence.
 */

template <typename TId>
inline Infix<Host<BamTagsDict const>::Type>::Type
getTagKey(BamTagsDict const & tags, TId id)
{
    return prefix(tags[id], 2);
}

// ----------------------------------------------------------------------------
// Function findTagKey()
// ----------------------------------------------------------------------------

/*!
 * @fn BamTagsDict#findTagKey
 * @brief Find a tag by its key for a @link BamTagsDict @endlink object.
 *
 * @signature bool findTagKey(id, tagsDict, key);
 *
 * @param[out] id       The id of the found tag.
 * @param[in]  tagsDict The BamTagsDict to query.
 * @param[in]  key      The key to query for: @link CharString @endlink.
 *
 * @return bool <tt>true</tt> if the key could be found and <tt>false</tt> otherwise.
 */

template <typename TId, typename TKey>
inline bool
findTagKey(TId & id, BamTagsDict const & tags, TKey const & key)
{
    for (id = 0; id < (TId)length(tags); ++id)
        if (getTagKey(tags, id) == key)
            return true;
    return false;
}

// ----------------------------------------------------------------------------
// Function extractTagValue()
// ----------------------------------------------------------------------------

/*!
 * @fn BamTagsDict#extractTagValue
 * @brief Extract and cast "atomic" value from tags string with index <tt>id</tt>.
 *
 * @signature bool extractTagValue(dest, tags, id)
 *
 * @param[out] dest The variable to write the value to.The value is first copied in a variable of the type indicated in
 *                  the BAM file. Then it is cast into the type of <tt>dest</tt>.
 *
 * @param[in] tags The BamTagsDict object to query.
 * @param[in] id   The id of the tag to extract the value from. See @link BamTagsDict#findTagKey @endlink.
 *
 * @return bool <tt>true</tt> if the value could be extracted.
 *
 * @section Remarks
 *
 * The function only works for atomic types such as <tt>int</tt>, not for <tt>char*</tt> or arrays.
 *
 * See @link BamTagsDict @endlink for an example.
 */

template <typename TResultType, typename TIter>
struct ExtractTagValueHelper_
{
    TResultType &result;
    TIter rawIter;
    char typeC;

    ExtractTagValueHelper_(TResultType &result, char typeC, TIter rawIter) :
        result(result),
        rawIter(rawIter),
        typeC(typeC)
    {}

    template <typename Type>
    bool operator() (Type) const
    {
        if (BamTypeChar<Type>::VALUE != typeC)
            return false;

        union {
            char raw[sizeof(Type)];
            Type i;
        } tmp;

        arrayCopyForward(rawIter, rawIter + sizeof(Type), tmp.raw);
        result = static_cast<TResultType>(tmp.i);
        return true;
    }
};

template <typename TResultValue, typename TId>
SEQAN_FUNC_ENABLE_IF(Is<IntegerConcept<TResultValue> >, bool)
extractTagValue(TResultValue & val, BamTagsDict const & tags, TId id)
{
    typedef Infix<Host<BamTagsDict const>::Type>::Type TInfix;
    typedef Iterator<TInfix, Standard>::Type TIter;

    TInfix inf = tags[id];
    if (length(inf) < 4 || inf[2] == 'Z')
        return false;

    TIter it = begin(inf, Standard()) + 2;
    char typeC = getValue(it++);
    ExtractTagValueHelper_<TResultValue, TIter> func(val, typeC, it);

    return tagApply(func, BamTagTypes());
}

template <typename TResultValue, typename TId>
SEQAN_FUNC_ENABLE_IF(IsSequence<TResultValue>, bool)
extractTagValue(TResultValue & val, BamTagsDict const & tags, TId id)
{
    typedef Infix<Host<BamTagsDict const>::Type>::Type TInfix;

    TInfix inf = tags[id];
    if (length(inf) < 4 || inf[2] != 'Z')
        return false;

    val = infix(inf, 3, length(inf) - 1);
    return true;
}

// ----------------------------------------------------------------------------
// Function getBamTypeChar()
// ----------------------------------------------------------------------------

/*!
 * @fn getBamTypeChar
 * @headerfile <seqan/bam_io.h>
 * @brief Return char identifying the type of the argument type.
 *
 * @signature char getBamTypeChar<T>();
 * @signature BamTypeChar<T>::VALUE
 *
 * @tparam T The type to query for its type char.
 *
 * @section Remarks
 *
 * Note that this function is defined for the <tt>__int16</tt>, <tt>__uint16</tt> etc. but not for the types
 * <tt>short</tt>, <tt>int</tt> etc. An exception are 8-bit characters/char, where it is defined for <tt>__int8</tt>,
 * <tt>__uint8</tt>, and <tt>char</tt> unless <tt>char</tt> is equal to one of the other two types. This is important
 * when used in @link BamTagsDict#setTagValue @endlink etc. since BAM gives type chars for printable characters, signed
 * 8-bit numbers and unsigned 8-bit numbers.
 *
 * If <tt>__int8</tt> and <tt>__uint8</tt> are not identical to <tt>char</tt>, we can make this decision from the type,
 * otherwise we cannot and we will give the integer types a higher precedence.
 *
 * In your programs, this should not make any difference, only the written SAM/BAM will differ.
 *
 * @see BamTagsDict
 * @see getBamTypeSize
 */

template <typename TValue>
struct BamTypeChar<TValue const> :
    BamTypeChar<TValue> {};

template <typename T>
inline char getBamTypeChar(T const &)
{
    return BamTypeChar<T>::VALUE;
}

// ----------------------------------------------------------------------------
// Function setTagValue()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Test me!

/*!
 * @fn BamTagsDict#setTagValue
 *
 * @headerfile <seqan/bam_io.h>
 *
 * @brief Set the value of a tag through a @link BamTagsDict @endlink.
 *
 * @signature bool setTagValue(tags, key, val[, typeC]);
 *
 * @param[in,out] tags  The BamTagsDict to modify.
 * @param[in]     key   The key of the tag. Must be a sequence of length 2.
 * @param[in]     val   The value to set the tag to.
 * @param[in]     typeC BAM type char to use. For portability (so the generated files are the same on all platforms), use
 *                      a signed/unsigned qualified type for <tt>val</tt> or give <tt>typeC</tt>. Also see the remarks
 *                      for @link getBamTypeChar @endlink. Types: getBamTypeChar@.
 *
 * @return bool <tt>true</tt> on success, <tt>false</tt> on failure.  This function can fail if the key is not a valid tag id (e.g. does
 *              not have length 2) or if the type of <tt>val</tt> is not an atomic value or a string (anything but
 *              <tt>char *</tt>, <tt>char const *</tt>, a character, integer or float type is invalid).
 *
 * @section Remarks
 *
 * Note that <tt>setTagValue</tt> does not cast the type, so <tt>typeC</tt> only influences the type character written
 * out but <tt>val</tt> is written out in binary without modification.
 *
 * @section Examples
 *
 * An example setting some atomic tag values.
 *
 * @code{.cpp}
 * CharString rawTagsText;
 * BamTagsDict tags(rawTagsText);
 * setTagValue(tags, "XA", 9);    // int
 * setTagValue(tags, "XB", 9u);   // unsigned int
 * setTagValue(tags, "XC", 'X');  // char
 * @endcode
 *
 * If <tt>char</tt> is equal to <tt>__int8</tt> or <tt>__uint8</tt> then the last line produces an entry with type 'c'
 * or 'C'. To make sure that the type char 'A' (for "printable character") is written to the file, give it explicitely:
 *
 * @code{.cpp}
 * setTagValue(tags, "XC", 'X', 'A');  // Overrwrite XC, enforce type 'printable character'.
 * @endcode
 *
 * Note that on most systems <tt>int</tt>s have a width of 32 bytes, but the C++ standard leaves this open. For all
 * types but characters, you should not give an explicit type char but use one of the types with explicit width and
 * signed/unsigned qualifier such as <tt>__int32</tt>, <tt>__uint32</tt> etc.
 *
 * @code{.cpp}
 * // The following is not recommended since the type of <tt>x</tt> is not "unsigned 32 bit int."
 * __int32 x = -1;
 * setTagValue(tags, "XB", x, 'I');
 * // Instead, explicitely use an unsigned type if you need one.  Note that your compiler
 * // might warn you about assigning -1 to an unsigned variable so you know that you are
 * // probably doing something unintended.
 * __uint32 y = -1;
 * setTagValue(tags, "XB", y);
 *
 * // Do not do this!
 * setTagValue(tags, "XA", 9, 'f');    // BOGUS since 9 is not a floating point number.
 * @endcode
 *
 * @see getBamTypeChar
 */

template <typename TBamValueSequence, typename TValue>
struct ToBamTagValueHelper_
{
    TBamValueSequence &result;
    TValue val;
    char typeC;

    ToBamTagValueHelper_(TBamValueSequence &result, char typeC, TValue val) :
        result(result),
        val(val),
        typeC(typeC)
    {}

    template <typename Type>
    bool operator() (Type) const
    {
        if (BamTypeChar<Type>::VALUE != typeC)
            return false;

        union {
            char raw[sizeof(Type)];
            Type i;
        } tmp;

        tmp.i = static_cast<Type>(val);
        append(result, toRange(&tmp.raw[0], &tmp.raw[sizeof(Type)]));
        return true;
    }
};

// Convert "atomic" value to BAM tag.  Return whether val was atomic.
template <typename TBamValueSequence, typename TValue>
SEQAN_FUNC_ENABLE_IF(Is<IntegerConcept<TValue> >, bool)
_toBamTagValue(TBamValueSequence & result, TValue const & val, char typeC)
{
    if (typeC == 'Z')
        return false;

    appendValue(result, typeC);
    ToBamTagValueHelper_<TBamValueSequence, TValue> func(result, typeC, val);
    if (tagApply(func, BamTagTypes()))
        return true;

    resize(result, length(result) - 1);
    return false;
}

template <typename TBamValueSequence, typename TValue>
SEQAN_FUNC_ENABLE_IF(IsSequence<TValue>, bool)
_toBamTagValue(TBamValueSequence & result, TValue const & val, char typeC)
{
    if (typeC != 'Z')
        return false;

    appendValue(result, typeC);
    append(result, val);
    appendValue(result, '\0');
    return true;
}


// Sets an atomic value in a BamTagsDict.
// Returns <tt>true</tt> successful, can fail if val not atomic or key is not a valid tag id (2 chars).

template <typename TKey, typename TValue>
inline bool
setTagValue(BamTagsDict & tags, TKey const & key, TValue const & val, char typeC)
{
    if (!hasIndex(tags))
        buildIndex(tags);

    // Build value to insert/append.
    if (length(key) != 2u)
        return false;

    Position<BamTagsDict>::Type id = 0;
    if (findTagKey(id, tags, key))
    {
        CharString bamTagVal;
        if (!_toBamTagValue(bamTagVal, val, typeC))
            return false;

        replace(host(tags), tags._positions[id] + 2, tags._positions[id + 1], bamTagVal);
        clear(tags._positions);
    }
    else
    {
        append(host(tags), key);
        if (!_toBamTagValue(host(tags), val, typeC))
        {
            resize(host(tags), length(host(tags)) - length(key));
            return false;
        }
        appendValue(tags._positions, length(host(tags)));
    }

    return true;
}

template <typename TKey, typename TValue>
inline bool
setTagValue(BamTagsDict & tags, TKey const & key, TValue const & val)
{
    return setTagValue(tags, key, val, BamTypeChar<TValue>::VALUE);
}

/*!
 * @fn BamTagsDict#appendTagValue
 *
 * @headerfile <seqan/bam_io.h>
 *
 * @brief Append a tag/value pair to a @link BamTagsDict @endlink.
 *
 * @signature bool appendTagValue(tags, key, val[, typeC]);
 *
 * @param[in,out] tags  The BamTagsDict to modify.
 * @param[in]     key   The key of the tag. Must be a sequence of length 2.
 * @param[in]     val   The value to set the tag to.
 * @param[in]     typeC BAM type char to use. For portability (so the generated files are the same on all platforms), use
 *                      a signed/unsigned qualified type for <tt>val</tt> or give <tt>typeC</tt>. Also see the remarks
 *                      for @link getBamTypeChar @endlink. Types: getBamTypeChar@.
 *
 * @return bool <tt>true</tt> on success, <tt>false</tt> on failure.  This function can fail if the key is not a valid tag id (e.g. does
 *              not have length 2) or if the type of <tt>val</tt> is not an atomic value or a string (anything but
 *              <tt>char *</tt>, <tt>char const *</tt>, a character, integer or float type is invalid).
 *
 * @section Remarks
 *
 * @link BamTagsDict#setTagValue @endlink behaves like <tt>appendTagValue</tt> if <tt>key</tt> was not part of <tt>tags</tt>
 * before. However, in this case <tt>appendTagValue</tt> is faster.
 */

template <typename TSequence, typename TKey, typename TValue>
inline bool
appendTagValue(TSequence & tags, TKey const & key, TValue const & val, char typeC)
{
    // Build value to insert/append.
    if (length(key) != 2u)
        return false;

    append(tags, key);
    return _toBamTagValue(tags, val, typeC);
}

template <typename TKey, typename TValue>
inline bool
appendTagValue(BamTagsDict & tags, TKey const & key, TValue const & val, char typeC)
{
    if (appendTagValue(host(tags), key, val, typeC))
    {
        appendValue(tags._positions, length(host(tags)));
        return true;
    }
    return false;
}


template <typename TDictOrString, typename TKey, typename TValue>
inline bool
appendTagValue(TDictOrString & tags, TKey const & key, TValue const & val)
{
    return appendTagValue(tags, key, val, BamTypeChar<TValue>::VALUE);
}


// ----------------------------------------------------------------------------
// Function eraseTag()
// ----------------------------------------------------------------------------

/*!
 * @fn BamTagsDict#eraseTag
 * @brief Erase a tag from BamTagsDict.
 *
 * @signature bool eraseTag(tagsDict, key);
 *
 * @param[in,out] tagsDict The BamTagsDict to erase the tag from.
 * @param[in]     key      The key of the tag to erase.
 *
 * @return bool <tt>true</tt> if the tag could be erased, <tt>false</tt> if the key wasn't present.
 */

template <typename TKey>
inline SEQAN_FUNC_DISABLE_IF(Is<IntegerConcept<TKey> >, bool)
eraseTag(BamTagsDict & tags, TKey const & key)
{
    if (!hasIndex(tags))
        buildIndex(tags);

    Position<BamTagsDict>::Type id = 0;
    if (!findTagKey(id, tags, key))
        return false;

    erase(host(tags), tags._positions[id], tags._positions[id + 1]);
    clear(tags._positions);
    return true;
}

template <typename TId>
inline SEQAN_FUNC_ENABLE_IF(Is<IntegerConcept<TId> >, bool)
eraseTag(BamTagsDict & tags, TId const & id)
{
    typedef typename Iterator<String<typename BamTagsDict::TPos>, Standard>::Type TIter;
    if (!hasIndex(tags))
        buildIndex(tags);

    typename BamTagsDict::TPos delta = tags._positions[id + 1] - tags._positions[id];
    erase(host(tags), tags._positions[id], tags._positions[id + 1]);
    erase(tags._positions, id);
    TIter it = begin(tags._positions, Standard()) + id;
    TIter itEnd = end(tags._positions, Standard());
    for (; it != itEnd; ++it)
        *it -= delta;
    return true;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_BAM_IO_BAM_TAGS_DICT_H_
