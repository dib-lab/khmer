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

/**
.Tag.Prefix Order:
..cat:Sequences
..summary:Specify whether less-than or greather-than comparison is meant.
..tag.TagPrefixLess:A prefix is smaller.
...text:For example: $"abc" < "abcde"$.
..tag.TagPrefixGreater:A prefix is greater.
...text:For example: $"abc" > "abcde"$.
..remarks:The default for all comparison functions is $TagPrefixLess$.
..include:seqan/sequence.h
*/
struct TagPrefixLess_ {};
typedef Tag<TagPrefixLess_> const TagPrefixLess;

struct TagPrefixGreater_ {};
typedef Tag<TagPrefixGreater_> const TagPrefixGreater;


/**
.Metafunction.DefaultPrefixOrder:
..hidefromindex
..summary:The default prefix order.
..signature:DefaultPrefixOrder<T>::Type
..param.T:Type for which the prefix order is determined.
..returns.param.Type:Prefix order tag for type of $T$.
..see:Tag.Prefix Order
..include:seqan/sequence.h
*/
template <typename T>
struct DefaultPrefixOrder
{
	typedef TagPrefixLess Type;
};

//////////////////////////////////////////////////////////////////////////////
// Lexical
//////////////////////////////////////////////////////////////////////////////

/**
.Class.Lexical:
..cat:Basic
..summary:Comparator for lexical comparison.
..signature:Lexical<TSpec>
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...text:This type can be used for specializations of $Lexical$.
...remarks:$TSpec$ is by default interpreted as size-type.
...default:$size_t$
..remarks:
...text:This class implement comparator objects that perform (lexical) comparisons between two sequences.
The result of the comparison is stored in the data members of the instance and can be
accessed by some functions, for example @Function.isLess@ or @Function.isEqual@.
...text:In most cases, there is no need for an explicite use of comparators,
but sometimes this concept provide the opportunity to speed up the code.
..example:
...text:This program compares the strings $str1$ and $str2$:
...code:if (isLess(str1, str2)) //first comparison
{
	//str1 < str2
}
else if (isGreater(str1, str2)) //second comparison
{
	//str1 > str2
}
else
{
	//str == str2
}
...text:Using a comparator, the same program only needs one comparison instead of two:
...code:Lexical <> comparator(str1, str2); //comparison is executed here
if (isLess(comparator))
{
	//str1 < str2
}
else if (lexGreater(comparator))
{
	//str1 > str2
}
else
{
	//str == str2
}
...text:The state of a default constructed $Lexical$ instance is undefined until
it is set by a call of @Function.compare@.
..see:Metafunction.Comparator
..include:seqan/sequence.h
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

/**
.Metafunction.Comparator:
..cat:Basic
..summary:Type of comparator object
..signature:Comparator<T>::Type
..param.T:Type for which the comparator type is to be determined.
..returns.param.Type:Comparator type
..remarks:Comparators are objects that can be used to compare other objects and store the
result of comparisons.
..include:seqan/sequence.h
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

/**
.Function.compare:
..cat:Comparisons
..summary:Compares two objects.
..signature:compare(comparator, left, right)
..param.left:The first objects.
..param.right:The second objects that is compared to $left$.
..param.comparator:Object that stores the results.
...type:Class.Lexical
..see:Metafunction.Comparator
..include:seqan/sequence.h
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

/**
.Function.isEqual:
..cat:Comparisons
..class:Class.Lexical
..summary:Operator "==".
..signature:isEqual(left, right)
..signature:isEqual(comparator)
..param.left:The first parameter.
..param.right:The second parameter that is compared to $left$.
..param.comparator:A comparator.
...type:Class.Lexical
..returns:$true$ if $left$ equals $right$, $false$ otherwise.
..see:Metafunction.Comparator
..include:seqan/sequence.h
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

/**
.Function.isNotEqual:
..cat:Comparisons
..class:Class.Lexical
..summary:Operator "!=".
..signature:isNotEqual(left, right)
..signature:isNotEqual(comparator)
..param.left:The first parameter.
..param.right:The second parameter that is compared to $left$.
..param.comparator:A comparator.
...type:Class.Lexical
..returns:$true$ if $left$ is not equal to $right$, $false$ otherwise.
..see:Metafunction.Comparator
..include:seqan/sequence.h
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

/**
.Function.isLess:
..cat:Comparisons
..class:Class.Lexical
..summary:Operator "<".
..signature:isLess(left, right [, prefix_order_tag])
..signature:isLess(comparator)
..param.left:The first parameter.
..param.right:The second parameter that is compared to $left$.
..param.prefix_order_tag:Tag that specify whether prefixes are less or greater. (optional)
...text:If omitted, the default tag is determined by @Metafunction.DefaultPrefixOrder@ for the type of $left$.
...see:Tag.Prefix Order
..param.comparator:A comparator.
...type:Class.Lexical
..returns:$true$ if $left$ is less than $right$, $false$ otherwise.
..see:Metafunction.Comparator
..remarks:
...text:Sequences are compared in lexicographical order.
..see:Tag.Prefix Order
..see:Metafunction.DefaultPrefixOrder
..include:seqan/sequence.h
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

/**
.Function.isLessOrEqual:
..cat:Comparisons
..class:Class.Lexical
..summary:Operator "<=".
..signature:isLessOrEqual(left, right [, prefix_order_tag])
..signature:isLessOrEqual(comparator)
..param.left:The first parameter.
..param.right:The second parameter that is compared to $left$.
..param.prefix_order_tag:Tag that specify whether prefixes are less or greater. (optional)
...text:If omitted, the default tag is determined by @Metafunction.DefaultPrefixOrder@ for the type of $left$.
...see:Tag.Prefix Order
..param.comparator:A comparator.
...type:Class.Lexical
..returns:$true$ if $left$ is less than or equal to $right$, $false$ otherwise.
..see:Metafunction.Comparator
..remarks:
...text:Sequences are compared in lexicographical order.
..see:Tag.Prefix Order
..see:Metafunction.DefaultPrefixOrder
..include:seqan/sequence.h
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

/**
.Function.isGreater:
..cat:Comparisons
..class:Class.Lexical
..summary:Operator ">".
..signature:isGreater(left, right [, prefix_order_tag])
..signature:isGreater(comparator)
..param.left:The first parameter.
..param.right:The second parameter that is compared to $left$.
..param.prefix_order_tag:Tag that specify whether prefixes are less or greater. (optional)
...text:If omitted, the default tag is determined by @Metafunction.DefaultPrefixOrder@ for the type of $left$.
...see:Tag.Prefix Order
..param.comparator:A comparator.
...type:Class.Lexical
..returns:$true$ if $left$ is greater than $right$, $false$ otherwise.
..see:Metafunction.Comparator
..remarks:
...text:Sequences are compared in lexicographical order.
..see:Tag.Prefix Order
..see:Metafunction.DefaultPrefixOrder
..include:seqan/sequence.h
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

/**
.Function.isGreaterOrEqual:
..cat:Comparisons
..class:Class.Lexical
..summary:Operator ">=".
..signature:isGreaterOrEqual(left, right [, prefix_order_tag])
..signature:isGreaterOrEqual(comparator)
..param.left:The first parameter.
..param.right:The second parameter that is compared to $left$.
..param.prefix_order_tag:Tag that specify whether prefixes are less or greater. (optional)
...text:If omitted, the default tag is determined by @Metafunction.DefaultPrefixOrder@ for the type of $left$.
...see:Tag.Prefix Order
..param.comparator:A comparator.
...type:Class.Lexical
..returns:$true$ if $left$ is greater than or equal to $right$, $false$ otherwise.
..see:Metafunction.Comparator
..remarks:
...text:Sequences are compared in lexicographical order.
..see:Tag.Prefix Order
..see:Metafunction.DefaultPrefixOrder
..include:seqan/sequence.h
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

/**
.Function.isPrefix:
..cat:Comparisons
..class:Class.Lexical
..summary:Test whether a sequence is prefix of another sequence.
..signature:isPrefix(left, right)
..signature:isPrefix(comparator)
..param.left:The first sequence, the putative prefix.
..param.right:The second sequence.
..param.comparator:A comparator.
...type:Class.Lexical
..returns:$true$ if $left$ is a prefix of $right$, $false$ otherwise.
..see:Metafunction.Comparator
..remarks:By definition, the whole sequence is a prefix of itself too: $isPrefix("abc", "abc") == true$.
..include:seqan/sequence.h
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

/**
.Function.hasPrefix:
..cat:Comparisons
..class:Class.Lexical
..summary:Test whether a sequence is prefix of another sequence.
..signature:hasPrefix(left, right)
..signature:hasPrefix(comparator)
..param.left:The first sequence.
..param.right:The second sequence, the putative prefix.
..param.comparator:A comparator.
...type:Class.Lexical
..returns:$true$ if $right$ is a prefix of $left$, $false$ otherwise.
..see:Metafunction.Comparator
..see:Function.isPrefix
..remarks:By definition, the whole sequence is a prefix of itself too: $hasPrefix("abc", "abc") == true$.
..include:seqan/sequence.h
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

/**
.Function.lcpLength:
..summary:Length of longest common prefix.
..cat:Comparisons
..class:Class.Lexical
..signature:lcpLength(left, right)
..signature:lcpLength(comparator)
..param.left:The first sequence.
..param.right:The second sequence that is compared to $left$.
..param.comparator:A comparator.
...type:Class.Lexical
..returns:The length of the longest common prefix of $left$ and $right$.
..see:Metafunction.Comparator
..include:seqan/sequence.h
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
