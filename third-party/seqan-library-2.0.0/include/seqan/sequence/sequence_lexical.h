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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Implementation of efficient lexical sequence comparison.
// ==========================================================================

#ifndef SEQAN_HEADER_LEXICAL_H
#define SEQAN_HEADER_LEXICAL_H


namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Switches for prefix ordering mode
//////////////////////////////////////////////////////////////////////////////

/*!
 * @defgroup PrefixOrderTags Prefix Order Tags
 * @brief Specify whether less-than or greather-than comparison is meant.
 *
 *
 * @tag PrefixOrderTags#TagPrefixLess
 * @headerfile <seqan/sequence.h>
 * @brief A prefix is smaller.
 *
 * @signature typedef Tag<TagPrefixLess_> const TagPrefixLess;
 *
 *
 * @tag PrefixOrderTags#TagPrefixGreater
 * @headerfile <seqan/sequence.h>
 * @brief A prefix is larger.
 *
 * @signature typedef Tag<TagPrefixGreater_> const TagPrefixGreater;
 */

struct TagPrefixLess_ {};
typedef Tag<TagPrefixLess_> const TagPrefixLess;

struct TagPrefixGreater_ {};
typedef Tag<TagPrefixGreater_> const TagPrefixGreater;


/*!
 * @mfn DefaultPrefixOrder
 * @headerfile <seqan/sequence.h>
 * @brief The default prefix order.
 *
 * @signature DefaultPrefixOrder<T>::Type;
 *
 * @tparam T The type to query for the prefix order.
 *
 * @return Type The prefix order tag type of T, see @link PrefixOrderTags @endlink.
 */

template <typename T>
struct DefaultPrefixOrder
{
    typedef TagPrefixLess Type;
};

//////////////////////////////////////////////////////////////////////////////
// Lexical
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class Lexical
 * @headerfile <seqan/sequence.h>
 * @brief Comparator for lexical comparison.
 *
 * @signature template <[typename TSpec]>
 *            class Lexical;
 *
 * @tparam TSpec The specializing size type, defaults to <tt>size_t</tt>.
 *
 * This class implement comparator objects that perform (lexical) comparisons between two sequences.  The result of the
 * comparison is stored in the data members of the instance and can be accessed by some functions, for example @link
 * Lexical#isLess @endlink or @link Lexical#isEqual @endlink.
 *
 * In most cases, there is no need for an explicite use of comparators, but sometimes this concept provide the
 * opportunity to speed up the code.
 *
 * @section Examples
 *
 * This program compares the strings <tt>str1</tt> and <tt>str2</tt>:
 *
 * @code{.cpp}
 * if (isLess(str1, str2)) //first comparison
 * {
 *     //str1 < str2
 * }
 * else if (isGreater(str1, str2)) //second comparison
 * {
 *     //str1 > str2
 * }
 * else
 * {
 *     //str == str2
 * }
 * @endcode
 *
 * Using a comparator, the same program only needs one comparison instead of two:
 *
 * @code{.cpp}
 * Lexical <> comparator(str1, str2); //comparison is executed here
 * if (isLess(comparator))
 * {
 *     //str1 < str2
 * }
 * else if (lexGreater(comparator))
 * {
 *     //str1 > str2
 * }
 * else
 * {
 *     //str == str2
 * }
 * @endcode
 *
 * The state of a default constructed <tt>Lexical</tt> instance is undefined until it is set by a call of @link
 * Lexical#compare @endlink.
 *
 * @see Comparator
 */

template <typename TSpec = size_t>
struct Lexical
{
public:
    typename Size<Lexical>::Type data_lcp;
    char data_compare;

public:
    Lexical()
    {
SEQAN_CHECKPOINT
    }

    template <typename TLeft, typename TRight>
    Lexical(TLeft const & left, TRight const & right)
    {
SEQAN_CHECKPOINT
        compare(*this, left, right);
    }

    Lexical(Lexical const & other):
        data_lcp(other.data_lcp),
        data_compare(other.data_compare)
    {
SEQAN_CHECKPOINT
    }

    Lexical & operator=(Lexical const & other)
    {
SEQAN_CHECKPOINT
        data_compare = other.data_compare;
        data_lcp = other.data_lcp;
        return *this;
    }

    ~Lexical() {}
//____________________________________________________________________________

    enum
    {
        EQUAL = 1,
        LESS = 2,
        GREATER = 4,
        LEFT_IS_PREFIX = 8,
        RIGHT_IS_PREFIX = 16
    };
};


//////////////////////////////////////////////////////////////////////////////
// Metafunctions
//////////////////////////////////////////////////////////////////////////////
// Comparator: returns object that can compare objects of type T

// TODO(holtgrew): Does this belong to the StringConcept?

/*!
 * @mfn Comparator
 * @headerfile <seqan/sequence.h>
 * @brief Type of comparator object
 *
 * @signature Comparator<T>::Type;
 *
 * @tparam T Type for which the comparator type is to be determined.
 *
 * @return Type the comparator type.
 *
 * Comparators are objects that can be used to compare other objects and store the result of comparisons.
 */

template <typename T>
struct Comparator
{
    typedef Lexical<typename Size<T>::Type> Type;
};

//////////////////////////////////////////////////////////////////////////////
// Size

template <typename TSpec>
struct Size<Lexical<TSpec> >
{
    typedef TSpec Type;
};

template <typename TSpec>
struct Size<Lexical<TSpec> const>
{
    typedef TSpec Type;
};

//////////////////////////////////////////////////////////////////////////////
// Spec

template <typename TSpec>
struct Spec<Lexical<TSpec> >
{
    typedef TSpec Type;
};

template <typename TSpec>
struct Spec<Lexical<TSpec> const>
{
    typedef TSpec Type;
};


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// compare
//////////////////////////////////////////////////////////////////////////////

// TODO(holtgrew): Does this belong to Lexical?

/*!
 * @fn Lexical#compare
 * @headerfile <seqan/sequence.h>
 * @brief Compares two objects.
 *
 * @signature void compare(comparator, left, right);
 *
 * @param[out] comparator Object that stores the results. Types: Lexical
 * @param[in]  left       The first objects.
 * @param[in]  right      The second objects that is compared to <tt>left</tt>.
 *
 * @see Comparator
 */

template <typename TSpec, typename TLeft, typename TRight>
inline void
compare_(Lexical<TSpec> & lexical,
         TLeft & left,
         TRight & right)
{
    typename Iterator<TLeft, Standard>::Type left_it = begin(left, Standard());
    typename Size<TLeft>::Type left_length = length(left);
    typename Iterator<TRight, Standard>::Type right_it = begin(right, Standard());
    typename Size<TRight>::Type right_length = length(right);

    if (left_length == right_length) lexical.data_compare = Lexical<TSpec>::EQUAL;
    else if (left_length < right_length) lexical.data_compare = Lexical<TSpec>::LEFT_IS_PREFIX;
    else
    {
        lexical.data_compare = Lexical<TSpec>::RIGHT_IS_PREFIX;
        left_length = right_length;
    }

    lexical.data_lcp = 0;
    for (lexical.data_lcp = 0; lexical.data_lcp < left_length; ++lexical.data_lcp)
    {
        if (*left_it < *right_it)
        {
            lexical.data_compare = Lexical<TSpec>::LESS;
            break;
        }
        if (*left_it > *right_it)
        {
            lexical.data_compare = Lexical<TSpec>::GREATER;
            break;
        }
        ++left_it;
        ++right_it;
    }
}
//////////////////////////////////////////////////////////////////////////////

template <typename TSpec, typename TLeft, typename TRight>
inline void
compare(Lexical<TSpec> & lexical,
        TLeft const & left,
        TRight const & right)
{
    compare_(lexical, left, right);
}

// TODO(holtgrew): Are these bugs present in currently supported VC++ versions or is this only a legacy issue?
//workaround for VC++ "const arrays" bug
template <typename TSpec, typename TLeftValue, typename TRight>
inline void
compare(Lexical<TSpec> & lexical,
        TLeftValue const * left,
        TRight const & right)
{
    compare_(lexical, left, right);
}
template <typename TSpec, typename TLeftValue, typename TRightValue>
inline void
compare(Lexical<TSpec> & lexical,
        TLeftValue const * left,
        TRightValue const * right)
{
    compare_(lexical, left, right);
}
template <typename TSpec, typename TLeft, typename TRightValue>
inline void
compare(Lexical<TSpec> & lexical,
        TLeft const & left,
        TRightValue const * right)
{
    compare_(lexical, left, right);
}

//////////////////////////////////////////////////////////////////////////////
// isEqual
//////////////////////////////////////////////////////////////////////////////

// TODO(holtgrew): Do we need the left/right specialization of this function here?

/*!
 * @fn Lexical#isEqual
 * @headerfile <seqan/sequence.h>
 * @brief Operator "==".
 *
 * @signature bool isEqual(left, right);
 * @signature bool isEqual(comparator);
 *
 * @param[in] left The first parameter.
 * @param[in] right The second parameter that is compared to <tt>left</tt>.
 * @param[in] comparator A comparator. Types: Lexical
 *
 * @return bool <tt>true</tt> if <tt>left</tt> equals <tt>right</tt>, <tt>false</tt> otherwise.
 */

template <typename TLeft, typename TRight >
inline bool
isEqual(TLeft const & left,
        TRight const & right)
{
SEQAN_CHECKPOINT
    return left == right;
}

template <typename TSpec>
inline bool
isEqual(Lexical<TSpec> const & _lex)
{
SEQAN_CHECKPOINT
    return (_lex.data_compare & Lexical<TSpec>::EQUAL);
}

//////////////////////////////////////////////////////////////////////////////
// isNotEqual
//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Lexical#isNotEqual
 * @headerfile <seqan/sequence.h>
 * @brief Operator "!=".
 *
 * @signature bool isNotEqual(left, right);
 * @signature bool isNotEqual(comparator);
 *
 * @param[in] left The first parameter.
 * @param[in] right The second parameter that is compared to <tt>left</tt>.
 * @param[in] comparator A comparator. Types: Lexical
 *
 * @return bool <tt>true</tt> if <tt>left</tt> does not equal <tt>right</tt>, <tt>false</tt> otherwise.
 */

template <typename TLeft, typename TRight >
inline bool
isNotEqual(TLeft const & left,
         TRight const & right)
{
SEQAN_CHECKPOINT
    return left != right;
}

template <typename TSpec>
inline bool
isNotEqual(Lexical<TSpec> const & _lex)
{
SEQAN_CHECKPOINT
    return !(_lex.data_compare & Lexical<TSpec>::EQUAL);
}

//////////////////////////////////////////////////////////////////////////////
// isLess
//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Lexical#isLess
 * @headerfile <seqan/sequence.h>
 * @brief Operator "&lt;".
 *
 * @signature bool isLess(left, right);
 * @signature bool isLess(comparator);
 *
 * @param[in] left The first parameter.
 * @param[in] right The second parameter that is compared to <tt>left</tt>.
 * @param[in] comparator A comparator. Types: Lexical
 *
 * @return bool <tt>true</tt> if <tt>left</tt> is less than <tt>right</tt>, <tt>false</tt> otherwise.
 */

template <typename TLeft, typename TRight, typename TPrefixOrder >
inline bool
isLess(TLeft const & left,
       TRight const & right,
       Tag<TPrefixOrder> const tag)
{
SEQAN_CHECKPOINT
    typename Comparator<TLeft>::Type _lex(left, right);
    return isLess(_lex, tag);
}
template <typename TLeft, typename TRight>
inline bool
isLess(TLeft const & left,
       TRight const & right)
{
SEQAN_CHECKPOINT
    return left < right;
}

template <typename TSpec>
inline bool
isLess(Lexical<TSpec> const & _lex,
       TagPrefixLess)
{
SEQAN_CHECKPOINT
   return (_lex.data_compare & (Lexical<TSpec>::LESS | Lexical<TSpec>::LEFT_IS_PREFIX)) != 0;
}
template <typename TSpec>
inline bool
isLess(Lexical<TSpec> const & _lex,
       TagPrefixGreater)
{
SEQAN_CHECKPOINT
   return (_lex.data_compare & (Lexical<TSpec>::LESS | Lexical<TSpec>::RIGHT_IS_PREFIX)) != 0;
}
template <typename TSpec>
inline bool
isLess(Lexical<TSpec> const & _lex)
{
SEQAN_CHECKPOINT
    return isLess(_lex, typename DefaultPrefixOrder< Lexical<TSpec> >::Type());
}

//////////////////////////////////////////////////////////////////////////////
// isLessOrEqual
//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Lexical#isLessOrEqual
 * @headerfile <seqan/sequence.h>
 * @brief Operator "&lt;=".
 *
 * @signature bool isLessOrEqual(left, right);
 * @signature bool isLessOrEqual(comparator);
 *
 * @param[in] left The first parameter.
 * @param[in] right The second parameter that is compared to <tt>left</tt>.
 * @param[in] comparator A comparator. Types: Lexical
 *
 * @return bool <tt>true</tt> if <tt>left</tt> is less than or equal to <tt>right</tt>, <tt>false</tt> otherwise.
 */

template <typename TLeft, typename TRight, typename TPrefixOrder >
inline bool
isLessOrEqual(TLeft const & left,
        TRight const & right,
        Tag<TPrefixOrder> const tag)
{
SEQAN_CHECKPOINT
    typename Comparator<TLeft>::Type _lex(left, right);
    return isLessOrEqual(_lex, tag);
}
template <typename TLeft, typename TRight>
inline bool
isLessOrEqual(TLeft const & left,
        TRight const & right)
{
SEQAN_CHECKPOINT
    return left <= right;
}

template <typename TSpec>
inline bool
isLessOrEqual(Lexical<TSpec> const & _lex,
        TagPrefixLess)
{
SEQAN_CHECKPOINT
   return (_lex.data_compare & (Lexical<TSpec>::LESS | Lexical<TSpec>::EQUAL | Lexical<TSpec>::LEFT_IS_PREFIX)) != 0;
}
template <typename TSpec>
inline bool
isLessOrEqual(Lexical<TSpec> const & _lex,
        TagPrefixGreater)
{
SEQAN_CHECKPOINT
   return (_lex.data_compare & (Lexical<TSpec>::LESS | Lexical<TSpec>::EQUAL | Lexical<TSpec>::RIGHT_IS_PREFIX)) != 0;
}
template <typename TSpec>
inline bool
isLessOrEqual(Lexical<TSpec> const & _lex)
{
SEQAN_CHECKPOINT
    return isLessOrEqual(_lex, typename DefaultPrefixOrder< Lexical<TSpec> >::Type());
}

//////////////////////////////////////////////////////////////////////////////
// isGreater
//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Lexical#isGreater
 * @headerfile <seqan/sequence.h>
 * @brief Operator "&gt;".
 *
 * @signature bool isGreater(left, right);
 * @signature bool isGreater(comparator);
 *
 * @param[in] left The first parameter.
 * @param[in] right The second parameter that is compared to <tt>left</tt>.
 * @param[in] comparator A comparator. Types: Lexical
 *
 * @return bool <tt>true</tt> if <tt>left</tt> is greater than <tt>right</tt>, <tt>false</tt> otherwise.
 */

template <typename TLeft, typename TRight, typename TPrefixOrder >
inline bool
isGreater(TLeft const & left,
        TRight const & right,
        Tag<TPrefixOrder> const tag)
{
SEQAN_CHECKPOINT
    typename Comparator<TLeft>::Type _lex(left, right);
    return isGreater(_lex, tag);
}
template <typename TLeft, typename TRight>
inline bool
isGreater(TLeft const & left,
        TRight const & right)
{
SEQAN_CHECKPOINT
    return left > right;
}

template <typename TSpec>
inline bool
isGreater(Lexical<TSpec> const & _lex,
        TagPrefixLess)
{
SEQAN_CHECKPOINT
   return (_lex.data_compare & (Lexical<TSpec>::GREATER | Lexical<TSpec>::RIGHT_IS_PREFIX)) != 0;
}
template <typename TSpec>
inline bool
isGreater(Lexical<TSpec> const & _lex,
        TagPrefixGreater)
{
SEQAN_CHECKPOINT
   return (_lex.data_compare & (Lexical<TSpec>::GREATER | Lexical<TSpec>::LEFT_IS_PREFIX)) != 0;
}
template <typename TSpec>
inline bool
isGreater(Lexical<TSpec> const & _lex)
{
SEQAN_CHECKPOINT
    return isGreater(_lex, typename DefaultPrefixOrder< Lexical<TSpec> >::Type());
}

//////////////////////////////////////////////////////////////////////////////
// isGreaterOrEqual
//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Lexical#isGreaterOrEqual
 * @headerfile <seqan/sequence.h>
 * @brief Operator "&gt;=".
 *
 * @signature bool isGreaterOrEqual(left, right);
 * @signature bool isGreaterOrEqual(comparator);
 *
 * @param[in] left The first parameter.
 * @param[in] right The second parameter that is compared to <tt>left</tt>.
 * @param[in] comparator A comparator. Types: Lexical
 *
 * @return bool <tt>true</tt> if <tt>left</tt> is greater than or equal to <tt>right</tt>, <tt>false</tt> otherwise.
 */

template <typename TLeft, typename TRight, typename TPrefixOrder >
inline bool
isGreaterOrEqual(TLeft const & left,
        TRight const & right,
        Tag<TPrefixOrder> const tag)
{
SEQAN_CHECKPOINT
    typename Comparator<TLeft>::Type _lex(left, right);
    return isGreaterOrEqual(_lex, tag);
}
template <typename TLeft, typename TRight>
inline bool
isGreaterOrEqual(TLeft const & left,
        TRight const & right)
{
SEQAN_CHECKPOINT
    return left >= right;
}

template <typename TSpec>
inline bool
isGreaterOrEqual(Lexical<TSpec> const & _lex,
        TagPrefixLess)
{
SEQAN_CHECKPOINT
   return (_lex.data_compare & (Lexical<TSpec>::GREATER | Lexical<TSpec>::EQUAL | Lexical<TSpec>::RIGHT_IS_PREFIX)) != 0;
}
template <typename TSpec>
inline bool
isGreaterOrEqual(Lexical<TSpec> const & _lex,
        TagPrefixGreater)
{
SEQAN_CHECKPOINT
   return (_lex.data_compare & (Lexical<TSpec>::GREATER | Lexical<TSpec>::EQUAL | Lexical<TSpec>::LEFT_IS_PREFIX)) != 0;
}
template <typename TSpec>
inline bool
isGreaterOrEqual(Lexical<TSpec> const & _lex)
{
SEQAN_CHECKPOINT
    return isGreaterOrEqual(_lex, typename DefaultPrefixOrder< Lexical<TSpec> >::Type());
}

//////////////////////////////////////////////////////////////////////////////
// isPrefix
//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Lexical#isPrefix
 * @headerfile <seqan/sequence.h>
 * @brief Test whether a sequence is the prefix of another sequence.
 *
 * @signature bool isPrefix(left, right);
 * @signature bool isPrefix(comparator);
 *
 * @param[in] left       The putative prefix.
 * @param[in] right      The second sequence.
 * @param[in] comparator A comparator. Types: Lexical
 *
 * @return bool <tt>true</tt> if <tt>left</tt> is a prefix of<tt>right</tt>, <tt>false</tt> otherwise.
 *
 * By definition, a sequence is a prefix of itself: <tt>isPrefix("abc", "abc")</tt> is <tt>true</tt>.
 */

template <typename TLeft, typename TRight >
inline bool
isPrefix(TLeft const & left,
        TRight const & right)
{
SEQAN_CHECKPOINT
    typename Comparator<TLeft>::Type _lex(left, right);
    return isPrefix(_lex);
}
template <typename TSpec>
inline bool
isPrefix(Lexical<TSpec> const & _lex)
{
SEQAN_CHECKPOINT
    return (_lex.data_compare & (Lexical<TSpec>::LEFT_IS_PREFIX | Lexical<TSpec>::EQUAL)) != 0;
}


//////////////////////////////////////////////////////////////////////////////
// hasPrefix
//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Lexical#hasPrefix
 * @headerfile <seqan/sequence.h>
 * @brief Test whether a sequence is the prefix of another sequence.
 *
 * @signature bool isPrefix(left, right);
 * @signature bool isPrefix(comparator);
 *
 * @param[in] left       The first sequence.
 * @param[in] right      The putative prefix.
 * @param[in] comparator A comparator. Types: Lexical
 *
 * @return bool <tt>true</tt> if <tt>left</tt> is a prefix of<tt>right</tt>, <tt>false</tt> otherwise.
 *
 * By definition, a sequence is a prefix of itself: <tt>hasPrefix("abc", "abc")</tt> is <tt>true</tt>.
 */

template <typename TLeft, typename TRight >
inline bool
hasPrefix(TLeft const & left,
        TRight const & right)
{
SEQAN_CHECKPOINT
    typename Comparator<TLeft>::Type _lex(left, right);
    return hasPrefix(_lex);
}
template <typename TSpec>
inline bool
hasPrefix(Lexical<TSpec> const & _lex)
{
SEQAN_CHECKPOINT
    return (_lex.data_compare & (Lexical<TSpec>::RIGHT_IS_PREFIX | Lexical<TSpec>::EQUAL)) != 0;
}

//////////////////////////////////////////////////////////////////////////////
// lcpLength
//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Lexical#lcpLength
 * @headerfile <seqan/sequence.h>
 * @brief Length of the longest common prefix.
 *
 * @signature TSize lcpLength(left, right);
 * @signature TSize lcpLength(comparator);
 *
 * @param[in] left       The first sequence.
 * @param[in] right      The second sequence.
 * @param[in] comparator A comparator. Types: Lexical
 *
 * @return TSize The length of the longest common prefix of <tt>left</tt> and <tt>right</tt>.  TSize is the Size type of
 *               the left size type.
 *
 * By definition, a sequence is a prefix of itself: <tt>hasPrefix("abc", "abc")</tt> is <tt>true</tt>.
 */

template <typename TLeft, typename TRight >
inline typename Size<TLeft>::Type
lcpLength(TLeft const & left, TRight const & right)
{
SEQAN_CHECKPOINT
    typename Comparator<TLeft>::Type _lex(left, right);
    return lcpLength(_lex);
}

template <typename TSpec>
inline typename Size< Lexical<TSpec> >::Type
lcpLength(Lexical<TSpec> const & _lex)
{
SEQAN_CHECKPOINT
    return _lex.data_lcp;
}

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
