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

#ifndef SEQAN_HEADER_MODIFIER_ALPHABET_H
#define SEQAN_HEADER_MODIFIER_ALPHABET_H

namespace seqan
{

//////////////////////////////////////////////////////////////////////////////

/**
.Class.ModifiedAlphabet:
..summary:Modifies value types.
..cat:Modifier
..signature:ModifiedAlphabet<TAlphabet, TSpec>
..param.TAlphabet:Original value type.
..param.TSpec:The modifier type.
...metafunction:Metafunction.Spec
...remarks:There is no default specialization.
..include:seqan/modifier.h
*/

template <typename THost, typename TSpec>
class ModifiedAlphabet;


//////////////////////////////////////////////////////////////////////////////
// sizes
//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec>
struct BitsPerValue<ModifiedAlphabet<THost, TSpec> >
    : BitsPerValue<THost>
{};

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec>
struct ValueSize<ModifiedAlphabet<THost, TSpec> >
    : ValueSize<THost>
{};

//////////////////////////////////////////////////////////////////////////////
// conversions
//////////////////////////////////////////////////////////////////////////////

template <typename TTarget, typename T, typename THost, typename TSpec>
inline typename Convert<TTarget, THost>::Type
convertImpl(Convert<TTarget, T> const convert_,
			ModifiedAlphabet<THost, TSpec> const & source_)
{
    SEQAN_CHECKPOINT;
	return convertImpl(convert_, static_cast<THost const &>(source_));
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec>
inline unsigned ordValue(ModifiedAlphabet<THost,TSpec> const &c) 
{
    SEQAN_CHECKPOINT;
	return ordValue(static_cast<THost const &>(c));
}

template <typename TStream, typename THost, typename TSpec>
TStream & operator<<(TStream & stream, ModifiedAlphabet<THost, TSpec> const & c)
{
    stream << convert<char>(c);
    return stream;
}

//////////////////////////////////////////////////////////////////////////////
// comparisons
//////////////////////////////////////////////////////////////////////////////


template <typename THost, typename TSpec, typename TRight>
struct CompareType<ModifiedAlphabet<THost, TSpec>, TRight>
{
	typedef typename CompareType<THost, TRight>::Type Type;
};


//////////////////////////////////////////////////////////////////////////////
// operator ==

template <typename THost, typename TSpec, typename TRight>
inline bool
operator==(ModifiedAlphabet<THost, TSpec> const & left_, 
           TRight const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}

template <typename TLeft, typename THost, typename TSpec>
inline bool
operator==(TLeft const & left_, 
           ModifiedAlphabet<THost, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}

template <typename TLeftHost, typename TLeftSpec, typename TRightHost, typename TRightSpec>
inline bool
operator==(ModifiedAlphabet<TLeftHost, TLeftSpec> const & left_, 
           ModifiedAlphabet<TRightHost, TRightSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<TLeftHost, TLeftSpec> TLeft;
	typedef ModifiedAlphabet<TRightHost, TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}

template <typename THost, typename TSpec>
inline bool
operator==(ModifiedAlphabet<THost, TSpec> const & left_, 
           ModifiedAlphabet<THost, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	return ordValue(left_) == ordValue(right_);
}

//SimpleType

template <typename TValue, typename TSpec, typename THost, typename TSpec2>
inline bool
operator==(SimpleType<TValue, TSpec> const & left_, 
           ModifiedAlphabet<THost, TSpec2> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}
template <typename TValue, typename TSpec, typename THost, typename TSpec2>
inline bool
operator==(ModifiedAlphabet<THost, TSpec2> const & left_,
           SimpleType<TValue, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}

//Proxy

template <typename TSpec, typename THost, typename TSpec2>
inline bool
operator==(Proxy<TSpec> const & left_, 
           ModifiedAlphabet<THost, TSpec2> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef Proxy<TSpec> TLeft;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}
template <typename TSpec, typename THost, typename TSpec2>
inline bool
operator==(ModifiedAlphabet<THost, TSpec2> const & left_,
			 Proxy<TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}

//////////////////////////////////////////////////////////////////////////////
// operator !=

template <typename THost, typename TSpec, typename TRight>
inline bool
operator!=(ModifiedAlphabet<THost, TSpec> const & left_, 
           TRight const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}

template <typename TLeft, typename THost, typename TSpec>
inline bool
operator!=(TLeft const & left_, 
           ModifiedAlphabet<THost, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}

template <typename TLeftHost, typename TLeftSpec, typename TRightHost, typename TRightSpec>
inline bool
operator!=(ModifiedAlphabet<TLeftHost, TLeftSpec> const & left_, 
			 ModifiedAlphabet<TRightHost, TRightSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<TLeftHost, TLeftSpec> TLeft;
	typedef ModifiedAlphabet<TRightHost, TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}

template <typename THost, typename TSpec>
inline bool
operator!=(ModifiedAlphabet<THost, TSpec> const & left_, 
           ModifiedAlphabet<THost, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	return ordValue(left_) != ordValue(right_);
}

//SimpleType

template <typename TValue, typename TSpec, typename THost, typename TSpec2>
inline bool
operator!=(SimpleType<TValue, TSpec> const & left_, 
           ModifiedAlphabet<THost, TSpec2> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}
template <typename TValue, typename TSpec, typename THost, typename TSpec2>
inline bool
operator!=(ModifiedAlphabet<THost, TSpec2> const & left_,
           SimpleType<TValue, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}

//Proxy

template <typename TSpec, typename THost, typename TSpec2>
inline bool
operator!=(Proxy<TSpec> const & left_, 
           ModifiedAlphabet<THost, TSpec2> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef Proxy<TSpec> TLeft;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}
template <typename TSpec, typename THost, typename TSpec2>
inline bool
operator!=(ModifiedAlphabet<THost, TSpec2> const & left_,
           Proxy<TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}

//////////////////////////////////////////////////////////////////////////////
// operator <=

template <typename THost, typename TSpec, typename TRight>
inline bool
operator<=(ModifiedAlphabet<THost, TSpec> const & left_, 
           TRight const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}

template <typename TLeft, typename THost, typename TSpec>
inline bool
operator<=(TLeft const & left_, 
           ModifiedAlphabet<THost, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}

template <typename TLeftHost, typename TLeftSpec, typename TRightHost, typename TRightSpec>
inline bool
operator<=(ModifiedAlphabet<TLeftHost, TLeftSpec> const & left_, 
			 ModifiedAlphabet<TRightHost, TRightSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<TLeftHost, TLeftSpec> TLeft;
	typedef ModifiedAlphabet<TRightHost, TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}

template <typename THost, typename TSpec>
inline bool
operator<=(ModifiedAlphabet<THost, TSpec> const & left_, 
           ModifiedAlphabet<THost, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	return ordValue(left_) <= ordValue(right_);
}

//SimpleType

template <typename TValue, typename TSpec, typename THost, typename TSpec2>
inline bool
operator<=(SimpleType<TValue, TSpec> const & left_, 
           ModifiedAlphabet<THost, TSpec2> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}
template <typename TValue, typename TSpec, typename THost, typename TSpec2>
inline bool
operator<=(ModifiedAlphabet<THost, TSpec2> const & left_,
           SimpleType<TValue, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}

//Proxy

template <typename TSpec, typename THost, typename TSpec2>
inline bool
operator<=(Proxy<TSpec> const & left_, 
           ModifiedAlphabet<THost, TSpec2> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef Proxy<TSpec> TLeft;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}
template <typename TSpec, typename THost, typename TSpec2>
inline bool
operator<=(ModifiedAlphabet<THost, TSpec2> const & left_,
           Proxy<TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}

//////////////////////////////////////////////////////////////////////////////
// operator <

template <typename THost, typename TSpec, typename TRight>
inline bool
operator<(ModifiedAlphabet<THost, TSpec> const & left_, 
          TRight const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}

template <typename TLeft, typename THost, typename TSpec>
inline bool
operator<(TLeft const & left_, 
          ModifiedAlphabet<THost, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}

template <typename TLeftHost, typename TLeftSpec, typename TRightHost, typename TRightSpec>
inline bool
operator<(ModifiedAlphabet<TLeftHost, TLeftSpec> const & left_, 
          ModifiedAlphabet<TRightHost, TRightSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<TLeftHost, TLeftSpec> TLeft;
	typedef ModifiedAlphabet<TRightHost, TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}

template <typename THost, typename TSpec>
inline bool
operator<(ModifiedAlphabet<THost, TSpec> const & left_, 
          ModifiedAlphabet<THost, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	return ordValue(left_) < ordValue(right_);
}

//SimpleType

template <typename TValue, typename TSpec, typename THost, typename TSpec2>
inline bool
operator<(SimpleType<TValue, TSpec> const & left_, 
          ModifiedAlphabet<THost, TSpec2> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}
template <typename TValue, typename TSpec, typename THost, typename TSpec2>
inline bool
operator<(ModifiedAlphabet<THost, TSpec2> const & left_,
          SimpleType<TValue, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}

//Proxy

template <typename TSpec, typename THost, typename TSpec2>
inline bool
operator<(Proxy<TSpec> const & left_, 
          ModifiedAlphabet<THost, TSpec2> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef Proxy<TSpec> TLeft;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}
template <typename TSpec, typename THost, typename TSpec2>
inline bool
operator<(ModifiedAlphabet<THost, TSpec2> const & left_,
          Proxy<TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}

//////////////////////////////////////////////////////////////////////////////
// operator >=

template <typename THost, typename TSpec, typename TRight>
inline bool
operator>=(ModifiedAlphabet<THost, TSpec> const & left_, 
           TRight const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}

template <typename TLeft, typename THost, typename TSpec>
inline bool
operator>=(TLeft const & left_, 
           ModifiedAlphabet<THost, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}

template <typename TLeftHost, typename TLeftSpec, typename TRightHost, typename TRightSpec>
inline bool
operator>=(ModifiedAlphabet<TLeftHost, TLeftSpec> const & left_, 
           ModifiedAlphabet<TRightHost, TRightSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<TLeftHost, TLeftSpec> TLeft;
	typedef ModifiedAlphabet<TRightHost, TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}

template <typename THost, typename TSpec>
inline bool
operator>=(ModifiedAlphabet<THost, TSpec> const & left_, 
           ModifiedAlphabet<THost, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	return ordValue(left_) >= ordValue(right_);
}

//SimpleType

template <typename TValue, typename TSpec, typename THost, typename TSpec2>
inline bool
operator>=(SimpleType<TValue, TSpec> const & left_, 
           ModifiedAlphabet<THost, TSpec2> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}
template <typename TValue, typename TSpec, typename THost, typename TSpec2>
inline bool
operator>=(ModifiedAlphabet<THost, TSpec2> const & left_,
           SimpleType<TValue, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}

//Proxy

template <typename TSpec, typename THost, typename TSpec2>
inline bool
operator>=(Proxy<TSpec> const & left_, 
           ModifiedAlphabet<THost, TSpec2> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef Proxy<TSpec> TLeft;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}
template <typename TSpec, typename THost, typename TSpec2>
inline bool
operator>=(ModifiedAlphabet<THost, TSpec2> const & left_,
           Proxy<TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}

//////////////////////////////////////////////////////////////////////////////
// operator >

template <typename THost, typename TSpec, typename TRight>
inline bool
operator>(ModifiedAlphabet<THost, TSpec> const & left_, 
          TRight const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}

template <typename TLeft, typename THost, typename TSpec>
inline bool
operator>(TLeft const & left_, 
          ModifiedAlphabet<THost, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}

template <typename TLeftHost, typename TLeftSpec, typename TRightHost, typename TRightSpec>
inline bool
operator>(ModifiedAlphabet<TLeftHost, TLeftSpec> const & left_, 
          ModifiedAlphabet<TRightHost, TRightSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<TLeftHost, TLeftSpec> TLeft;
	typedef ModifiedAlphabet<TRightHost, TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}

template <typename THost, typename TSpec>
inline bool
operator>(ModifiedAlphabet<THost, TSpec> const & left_, 
          ModifiedAlphabet<THost, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	return ordValue(left_) > ordValue(right_);
}

//SimpleType

template <typename TValue, typename TSpec, typename THost, typename TSpec2>
inline bool
operator>(SimpleType<TValue, TSpec> const & left_, 
          ModifiedAlphabet<THost, TSpec2> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}
template <typename TValue, typename TSpec, typename THost, typename TSpec2>
inline bool
operator>(ModifiedAlphabet<THost, TSpec2> const & left_,
          SimpleType<TValue, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}

//Proxy

template <typename TSpec, typename THost, typename TSpec2>
inline bool
operator>(Proxy<TSpec> const & left_, 
          ModifiedAlphabet<THost, TSpec2> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef Proxy<TSpec> TLeft;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}
template <typename TSpec, typename THost, typename TSpec2>
inline bool
operator>(ModifiedAlphabet<THost, TSpec2> const & left_,
          Proxy<TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename TTarget, typename TSource>
inline void
_initializeAlphabetConversionTable(TTarget *,
								   TSource const &)
{
    SEQAN_CHECKPOINT;
	//default: do nothing (because this array is not used)
	//define this function for each conversion table
}


template <typename TTarget, typename TSource>
struct AlphabetConversionTable_
{
	enum { SIZE = InternalValueSize_<TSource>::VALUE };
public:
	static TTarget * table;
	static TTarget * initialize()
	{
        SEQAN_CHECKPOINT;
        static TTarget table_store[SIZE];
		static bool _is_initialized = false;
		if (! _is_initialized)
		{
			_initializeAlphabetConversionTable(table_store, TSource());
		}
		_is_initialized = true;
		return table_store;
	}
};

template <typename TTarget, typename TSource>
TTarget * AlphabetConversionTable_<TTarget, TSource>::table = AlphabetConversionTable_<TTarget, TSource>::initialize();


//////////////////////////////////////////////////////////////////////////////

template <typename TTarget, typename TSource>
inline void
_initializeAlphabetOrdTable(TTarget *,
							TSource const &)
{
    SEQAN_CHECKPOINT;
	//default: do nothing (because this array is not used)
	//define this function for each conversion table
}


template <typename TSource>
struct AlphabetOrdTable_
{
	enum { SIZE = InternalValueSize_<TSource>::VALUE };
public:
	static unsigned * table;
	static unsigned * initialize()
	{
        SEQAN_CHECKPOINT;
        static unsigned table_store[SIZE];
		static bool _is_initialized = false;
		if (! _is_initialized)
		{
			_initializeAlphabetOrdTable(table_store, TSource());
		}
		_is_initialized = true;
		return table_store;
	}
};


template <typename TSource>
unsigned * AlphabetOrdTable_<TSource>::table = AlphabetOrdTable_<TSource>::initialize();

}  // namespace seqan

#endif
