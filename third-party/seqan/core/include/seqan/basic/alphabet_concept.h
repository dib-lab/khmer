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
// Author: David Weese <david.weese@fu-berlin.de>
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Concept definitions for alphabets.
// ==========================================================================

// SEQAN_NO_GENERATED_FORWARDS

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_ALPHABET_CONCEPT_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_ALPHABET_CONCEPT_H_

namespace seqan {

// ============================================================================
// Concepts for generic alphabets
// ============================================================================

/**
.Concept.AlphabetConcept
..cat:Alphabets
..baseconcept:Concept.AssignableConcept
..baseconcept:Concept.DefaultConstructibleConcept
..baseconcept:Concept.CopyConstructibleConcept
..summary:Natural container value.
..include:seqan/basic.h
..example.text:Valid expressions ($v$ is of type $T$):
..example.code:
unsigned bpv = BitsPerValue<T>::VALUE;

.Metafunction.BitsPerValue
..concept:Concept.AlphabetConcept
..cat:Basic
..summary:Number of bits needed to store a value.
..signature:BitsPerValue<T>::VALUE
..param.T:A class.
...type:Concept.AlphabetConcept
..returns.param.VALUE:Number of bits needed to store $T$.
...default:$sizeof<T> * 8$
..see:Metafunction.ValueSize
..include:seqan/basic.h
 */

// Forwards for Metafunctions and Functions.
template <typename T> struct BitsPerValue;

// minimal requirements for the alphabet of a String class
SEQAN_CONCEPT_REFINE(AlphabetConcept, (TValue), (Assignable)(DefaultConstructible)(CopyConstructible))
{
    typedef typename BitsPerValue<TValue>::Type TBitsPerValue;
    
    TValue val, val2;

    SEQAN_CONCEPT_USAGE(AlphabetConcept)
    {
        SEQAN_STATIC_ASSERT_MSG(BitsPerValue<TValue>::VALUE != 0, "Alphabet types must implement the BitsPerValue metafunction with non-zero value.");

        // assign must be available as an equivalent to '='
        assign(val, val2);
//      swap(val, val2);

        TBitsPerValue b = BitsPerValue<TValue>::VALUE;
        
        ignoreUnusedVariableWarning(b);
    }
};

// ============================================================================
// Concepts For Alphabets From The Mathematics Domain.
// ============================================================================

/**
.Concept.OrderedAlphabetConcept
..concept:Concept.OrderedAlphabetConcept
..cat:Alphabets
..summary:Totally strict ordered alphabet.
..baseconcept:Concept.AlphabetConcept
..baseconcept:Concept.ComparableConcept
..include:seqan/basic.h

.Function.operator<.concept:Concept.OrderedAlphabetConcept

.Metafunction.MaxValue
..concept:Concept.OrderedAlphabetConcept
..cat:Miscellaneous
..summary:Supremum for a given type.
..signature:MaxValue<T>::VALUE
..param.T:An ordered type.
...type:Concept.OrderedAlphabetConcept
..returns.param.VALUE:A value $sup$ for which holds: $sup >= i$ for all values $i$ of type $T$.
..remarks:Note tat
..see:Function.maxValue
..include:seqan/basic.h

.Metafunction.MinValue
..concept:Concept.OrderedAlphabetConcept
..cat:Miscellaneous
..summary:Infimum for a given type.
..signature:MinValue<T>::VALUE
..param.T:An ordered type.
...type:Concept.OrderedAlphabetConcept
..returns.param.VALUE:A value $inf$ for which holds: $inf <= i$ for all values $i$ of type $T$.
..remarks:Note tat
..see:Function.minValue
..include:seqan/basic.h

.Function.supremumValueImpl
..concept:Concept.OrderedAlphabetConcept
..hidefromindex
..cat:Alphabets
..summary:Implements @Function.maxValue@.
..signature:supremumValueImpl(value_pointer_tag)
..param.value_pointer_tag:A pointer that is used as a tag to specify the value type.
...type:Concept.OrderedAlphabetConcept
...remarks:The pointer needs not to point to a valid object, so it is possible to use a null pointer here.
..returns:A value $inf$ that holds: $inf >= i$ for all values $i$.
..remarks.text:This function implements @Function.maxValue@. 
It is recommended to use @Function.maxValue@ rather than $supremumValueImpl$.
..status:deprecated, will be removed in favour of @Metafunction.MaxValue@
..include:seqan/basic.h

.Function.maxValue
..concept:Concept.OrderedAlphabetConcept
..cat:Alphabets
..summary:Supremum for a given type.
..signature:maxValue<T>()
..param.T:An ordered type.
...type:Concept.OrderedAlphabetConcept
..returns:A value $inf$ that holds: $inf >= i$ for all values $i$ of type $T$.
..remarks:
The function is implemented in @Function.supremumValueImpl@. 
Do not specialize $maxValue$, specialize @Function.supremumValueImpl@ instead!
..see:Function.supremumValueImpl
..status:deprecated, will be removed in favour of @Metafunction.MaxValue@
..include:seqan/basic.h

// TODO(holtgrew): Rename to minValueImpl?

.Function.infimumValueImpl
..concept:Concept.OrderedAlphabetConcept
..hidefromindex
..cat:Alphabets
..summary:Implements @Function.minValue@.
..signature:infimumValueImpl(value_pointer_tag)
..param.value_pointer_tag:A pointer that is used as a tag to specify the value type.
...type:Concept.OrderedAlphabetConcept
...remarks:The pointer needs not to point to a valid object, so it is possible to use a null pointer here.
..returns:A value $inf$ that holds: $inf <= i$ for all values $i$.
..remarks.text:This function implements @Function.minValue@. 
It is recommended to use @Function.minValue@ rather than $infimumValueImpl$.
..status:deprecated, will be removed in favour of @Metafunction.MinValue@
..include:seqan/basic.h

.Function.minValue
..concept:Concept.OrderedAlphabetConcept
..cat:Alphabets
..summary:Infimum for a given type.
..signature:minValue<T>()
..param.T:An ordered type.
...type:Concept.OrderedAlphabetConcept
..returns:A value $inf$ that holds: $inf <= i$ for all values $i$ of type $T$.
..remarks.text:The function is implemented in @Function.infimumValueImpl@. 
Do not specialize $minValue$, specialize @Function.infimumValueImpl@ instead!
..see:Function.infimumValueImpl
..see:Function.maxValue
..status:deprecated, will be removed in favour of @Metafunction.MinValue@
..include:seqan/basic.h
*/

// Forwards for Metafunctions and Functions.
template <typename T> struct MinValue;
template <typename T> struct MaxValue;
template <typename T> T const & minValue();
template <typename T> T const & minValue(T);
template <typename T> T const & maxValue();
template <typename T> T const & maxValue(T);

SEQAN_CONCEPT_REFINE(OrderedAlphabetConcept, (TValue), (AlphabetConcept)(Comparable))
{
    TValue val;

    SEQAN_CONCEPT_USAGE(OrderedAlphabetConcept)
    {
        // type consistency checks
        sameType(minValue(val), val);
        sameType(minValue<TValue>(), val);
        sameType(MinValue<TValue>::VALUE, val);
        sameType(maxValue(val), val);
        sameType(maxValue<TValue>(), val);
        sameType(MaxValue<TValue>::VALUE, val);

        // TODO(holtgrew): This does not work in C++98, we need C++11 with constexpr.
        // TODO(holtgrew): Do these tests for each alphabet in runtime tests.
        // sanity checks
        // SEQAN_STATIC_ASSERT_MSG(MinValue<TValue>::VALUE <= MaxValue<TValue>::VALUE, "Minimal alphabet value must be less or equal to the maximal value.");

        // TODO(holtgrew): This does not work in C++98, we need C++11 with constexpr, cannot cast non-integral and non-enumeration types at compile time in C++98.
        // 0 must be an element of the alphabet, as we want to be able
        // to initialize a TValue variable to omit uninitialized warnings.
        // SEQAN_STATIC_ASSERT_MSG(MinValue<TValue>::VALUE <= static_cast<TValue>(0), "0 must be convertible to a valid alphabet value.");
        // SEQAN_STATIC_ASSERT_MSG(static_cast<TValue>(0) <= MaxValue<TValue>::VALUE, "0 must be convertible to a valid alphabet value.");
    }
};

/**
.Concept.FiniteOrderedAlphabetConcept
..cat:Alphabets
..summary:An type that is of finite domain and totally ordered and thus has a minimum and maximum value.
..baseconcept:Concept.OrderedAlphabetConcept
..include:seqan/basic.h

.Function.ordValue
..concept:Concept.FiniteOrderedAlphabetConcept
..summary:Maps an alphabet 1-to-1 to the interval [0..ValueSize).
..cat:Alphabets
..signature:ordValue(value)
..param.value:Arbitrary character value.
...type:Class.SimpleType
..returns:An unsigned value (result of @Metafunction.Size@$<typeof(value)>$ between 0 and @Metafunction.ValueSize@ of the type of value.
..note:This function first converts value to its unsigned value type and after that to an $unsigned int$.
You can't use $(unsigned int)c$ for a character $c$ as on some systems $char$ is signed and a $-1$ would be mapped to $0xffffffff$ instead of $0x000000ff$.
..include:seqan/sequence.h

.Function.valueSize
..concept:Concept.FiniteOrderedAlphabetConcept
..cat:Alphabets
..summary:Returns size of an alphabet.
..signature:valueSize<T>()
..param.T:Type to query for value size.
..returns:Number of values in type $T$.
...type:Metafunction.ValueSize
..include:seqan/basic.h
..see:Metafunction.ValueSize

.Metafunction.ValueSize:
..concept:Concept.FiniteOrderedAlphabetConcept
..cat:Basic
..summary:Number of different values a value type object can have.
..signature:ValueSize<T>::VALUE
..param.T:A class.
...type:Concept.FiniteOrderedAlphabetConcept
..returns.param.VALUE:Value size of $T$.
..remarks:
This function is only defined for integral types like $unsigned$, $int$, or @Spec.Dna@.
For floating point numbers and the 64 bit types $__int64$ and $__uint64$, it returns 0 since there is no standard compliant way to return the number of values for these types.
..remarks:
Note that you cannot get pointers or references to $ValueSize<T>::VALUE$ in your program.
You can use @Function.valueSize@ in your programs without problems, though.
When you get problems in your tests, use the "unary plus" workaround from the examples section.
..example.text:The temporary assignment workaround.
..example.code:
SEQAN_ASSERT_EQ(ValueSize<bool>::VALUE, 2u);    // Linker error.
SEQAN_ASSERT_EQ(+ValueSize<bool>::VALUE, 2u);  // OK
SEQAN_ASSERT_EQ(valueSize<bool>(), 2u);         // OK
..see:Function.valueSize
..see:Metafunction.Value
..include:seqan/basic.h
 */

// Forwards for Metafunctions and Functions.
template <typename T> struct ValueSize;
template <typename T> typename ValueSize<T>::Type valueSize();
// Forwards for Metafunctions and Functions.
template <typename TValue> typename ValueSize<TValue>::Type ordValue(TValue const & c);

SEQAN_CONCEPT_REFINE(FiniteOrderedAlphabetConcept, (TValue), (OrderedAlphabetConcept))
{
    typedef typename ValueSize<TValue>::Type TSize;

    TValue  val;
    TSize   size;

    SEQAN_CONCEPT_ASSERT((UnsignedIntegerConcept<TSize>));

    SEQAN_CONCEPT_USAGE(FiniteOrderedAlphabetConcept)
    {
        // a finite alphabet must be countable
        sameType(ordValue(val), size);
        sameType(valueSize<TValue>(), size);
        sameType(ValueSize<TValue>::VALUE, size);

        // alphabet must be non-empty
        SEQAN_STATIC_ASSERT_MSG(static_cast<TSize>(0) < ValueSize<TValue>::VALUE, "Alphabet size be greater than zero.");
        
        // convert integer to alphabet value
        val = 0;
        val = size;
        TValue val2(0);
        TValue val3(size);

        ignoreUnusedVariableWarning(val2);
        ignoreUnusedVariableWarning(val3);
    }
};

// ============================================================================
// Concepts For Alphabets From The Bioinformatics Domain.
// ============================================================================

/**
.Concept.AlphabetWithGapsConcept
..cat:Alphabets
..baseconcept:Concept.AlphabetConcept
..summary:An alphabet that includes a specific gap character.
..include:seqan/basic.h

.Function.gapValue
..concept:Concept.AlphabetWithGapsConcept
..cat:Alphabets
..cat:Alignments
..summary:Return the "gap" value from an alphabet.
..signature:gapValue<T>()
..param.T:The alphabet type to query the "gap" value from.
...type:Concept.AlphabetWithGapsConcept
..returns:The gap character.
..remarks.text:The function is implemented in @Function.gapValueImpl@.
Do not specialize $gapValue$, specialize @Function.gapValueImpl@ instead!
..see:Function.gapValueImpl
..include:seqan/basic.h

.Function.gapValueImpl
..concept:Concept.AlphabetWithGapsConcept
..hidefromindex
..cat:Alphabets
..cat:Alignments
..summary:Implements @Function.gapValue@.
..signature:gapValueImpl(valuePointerTag)
..param.valuePointerTag:A pointer that is used as a tag to specify the value type.
....type:Concept.AlphabetWithGapsConcept
...remarks:The pointer needs not to point to a valid object, so it is possible to use a null pointer here.
..returns:A gap character.
..see:Function.gapValue
..remarks.text:This function implements @Function.gapValue@.
It is recommended to use @Function.gapValue@ rather than $gapValueImpl$.
..include:seqan/basic.h
*/

// Forwards for Metafunctions and Functions.
template <typename T> T gapValue();
template <typename T> T gapValueImpl(T *);

SEQAN_CONCEPT_REFINE(AlphabetWithGapsConcept, (TValue), (AlphabetConcept))
{
    TValue val;

    SEQAN_CONCEPT_USAGE(AlphabetWithGapsConcept)
    {
        // Test the availability and return type of gapValue() and gapValueImpl().
        sameType(gapValue<TValue>(), val);
        sameType(gapValueImpl<TValue>(static_cast<TValue *>(0)), val);
    }
};

/**
.Concept.AlphabetWithUnknownValueConcept
..cat:Alphabets
..baseconcept:Concept.AlphabetConcept
..summary:An alphabet which includes a specific "unknown" character.
..include:seqan/basic.h

.Function.unknownValueImpl
..concept:Concept.AlphabetWithUnknownValueConcept
..hidefromindex
..cat:Alphabets
..summary:Implements @Function.unknownValue@.
..signature:gapValueImpl(valuePointerTag)
..param.valuePointerTag:A pointer that is used as a tag to specify the value type.
...type:Concept.AlphabetWithUnknownValueConcept
...remarks:The pointer needs not to point to a valid object, so it is possible to use a null pointer here.
..returns:A "unknown" character.
..see:Function.unknownValue
..remarks.text:This function implements @Function.unknownValue@.
It is recommended to use @Function.gapValue@ rather than $gapValueImpl$.
..include:seqan/basic.h

.Function.unknownValue
..concept:Concept.AlphabetWithUnknownValueConcept
..cat:Alphabets
..summary:Return the "unknown" value from an alphabet.
..signature:unknownValue<T>()
..param.T:The alphabet type to query the "unknown" value from.
...type:Concept.AlphabetWithUnknownValueConcept
..returns:The "unknown" value.
 */

// Forwards for Metafunctions and Functions.
template <typename T> T unknownValue();
template <typename T> T unknownValueImpl(T *);

SEQAN_CONCEPT_REFINE(AlphabetWithUnknownValueConcept, (TValue), (AlphabetConcept))
{
    TValue val;

    SEQAN_CONCEPT_USAGE(AlphabetWithUnknownValueConcept)
    {
        // Test the availability and return type of unknownValue() and unknownValueImpl().
        sameType(unknownValue<TValue>(), val);
        sameType(unknownValueImpl<TValue>(static_cast<TValue *>(0)), val);
    }
};

/**
.Concept.AlphabetWithQualitiesConcept
..cat:Alphabets
..baseconcept:Concept.AlphabetConcept
..summary:An alphabet where qualities can be attached to the characters.
..include:seqan/basic.h

.Metafunction.QualityValueSize
..concept:Concept.AlphabetWithQualitiesConcept
..cat:Alphabets
..summary:Return the number of quality values in characters from alphabet with qualities.
..signature:QualityValueSize<TAlphabet>::VALUE
..param.TAlphabet:The alphabet to query for its value size.
...type:Concept.AlphabetWithQualitiesConcept
..returns:The cardinality of the set of qualities.
..include:seqan/basic.h

.Metafunction.HasQualities
..concept:Concept.AlphabetWithQualitiesConcept
..cat:Alphabets
..summary:Return whether the given type stores qualities besides the alphabet.
..signature:HasQualities<TAlphabet>::VALUE
..signature:HasQualities<TAlphabet>::Type
..param.TAlphabet:The alphabe to query.
..returns:$true$, $false$, $True$, or $False$.
..include:seqan/basic.h

.Function.getQualityValue
..concept:Concept.AlphabetWithQualitiesConcept
..cat:Alphabets
..signature:getQualityValue(c)
..summary:Returns the quality of a character from an alphabet with integrated quality, e.g. the quality associated with a specified element from a sequence.
..param.c:Character to retrieve the quality from.
...type:Concept.AlphabetWithQualitiesConcept
..returns:Quality value of $c$.
...type:nolink:int
...remarks:The quality value is an integral value between 0 and 62 (inclusive).
..see:Function.assignQualityValue
..see:Function.convertQuality
..include:seqan/basic.h
..example.code:
    String<Dna5Q> seq = "TATA";
    // Assign quality value to first 'T' in sequence seq
    assignQualityValue(seq[0], 35);
    // Print quality value of first 'T', and default quality value of first 'A'
    std::cout << getQualityValue(seq[0]) << std::endl; // Defined as 35
    std::cout << getQualityValue(seq[1]) << std::endl; // Default value 60


.Function.assignQualityValue
..concept:Concept.AlphabetWithQualitiesConcept
..cat:Alphabets
..signature:assignQualityValue(c, q)
..summary:Assigns quality to a character from an alphabet with integrated quality, e.g. to a specified element from a sequence.
..param.c:Target character to assign quality to.
...type:Concept.AlphabetWithQualitiesConcept
..param.q:Quality to assign to the character.
...type:nolink:int
...type:nolink:char
...remarks:The quality value is an integral value between 0 and 62 (inclusive).
..remarks:If $q$ is a $char$ then $'!'$ is subtracted from $q$. This is useful for ASCII encoded PHRED scores.
..see:Function.getQualityValue
..see:Function.convertQuality
..include:seqan/basic.h
 */

// TODO(holtgrew): What about different quality types? Guess scaling? Look at how other packages do this.

SEQAN_CONCEPT_REFINE(AlphabetWithQualitiesConcept, (TValue), (AlphabetConcept))
{
    TValue val;

    SEQAN_CONCEPT_USAGE(AlphabetWithQualitiesConcept)
    {
        // TODO(holtgrew): Write me!
    }
};


}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_ALPHABET_CONCEPT_H_
