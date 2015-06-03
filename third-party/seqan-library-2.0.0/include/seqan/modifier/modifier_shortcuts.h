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

#ifndef SEQAN_HEADER_MODIFIER_SHORTCUTS_H
#define SEQAN_HEADER_MODIFIER_SHORTCUTS_H

namespace seqan
{

// ==========================================================================
// Shortcuts for Modified Strings.
// ==========================================================================



/*!
 * @typedef DnaStringComplement
 * @headerfile <seqan/modifier.h>
 * @brief Modifier for the complement of DnaString.
 *
 * @signature typedef ModifiedString<DnaString, ModView<FunctorComplementDna> > > DnaStringComplement;
 */

/*!
 * @typedef Dna5StringComplement
 * @headerfile <seqan/modifier.h>
 * @brief Modifier for the complement of Dna5String.
 *
 * @signature typedef ModifiedString<Dna5String, ModView<FunctorComplementDna> > > DnaStringComplement;
 */

/*!
 * @typedef RnaStringComplement
 * @headerfile <seqan/modifier.h>
 * @brief Modifier for the complement of RnaString.
 *
 * @signature typedef ModifiedString<RnaString, ModView<FunctorComplementDna> > > DnaStringComplement;
 */

/*!
 * @typedef Rna5StringComplement
 * @headerfile <seqan/modifier.h>
 * @brief Modifier for the complement of Rna5String.
 *
 * @signature typedef ModifiedString<Rna5String, ModView<FunctorComplementDna> > > DnaStringComplement;
 */

/*!
 * @typedef DnaStringReverse
 * @headerfile <seqan/modifier.h>
 * @brief Modifier for the reverse of a DnaString.
 *
 * @signature typedef ModifiedString<DnaString, ModReverse> DnaStringReverse;
 */

/*!
 * @typedef Dna5StringReverse
 * @headerfile <seqan/modifier.h>
 * @brief Modifier for the reverse of a Dna5String.
 *
 * @signature typedef ModifiedString<Dna5String, ModReverse> DnaStringReverse;
 */

/*!
 * @typedef RnaStringReverse
 * @headerfile <seqan/modifier.h>
 * @brief Modifier for the reverse of a RnaString.
 *
 * @signature typedef ModifiedString<RnaString, ModReverse> DnaStringReverse;
 */

/*!
 * @typedef Rna5StringReverse
 * @headerfile <seqan/modifier.h>
 * @brief Modifier for the reverse of a Rna5String.
 *
 * @signature typedef ModifiedString<Rna5String, ModReverse> DnaStringReverse;
 */

/*!
 * @typedef DnaStringReverseComplement
 * @headerfile <seqan/modifier.h>
 * @brief Modifier for the reverse-complement of a DnaString.
 *
 * @signature typedef ModifiedString<ModifiedString<DnaString, ModView<FunctorComplement<Dna5> >, ModReverse> DnaStringReverseComplement;
 */

/*!
 * @typedef Dna5StringReverseComplement
 * @headerfile <seqan/modifier.h>
 * @brief Modifier for the reverse-complement of a Dna5String.
 *
 * @signature typedef ModifiedString<ModifiedString<Dna5String, ModView<FunctorComplement<Dna55> >, ModReverse> Dna5StringReverseComplement;
 */

/*!
 * @typedef RnaStringReverseComplement
 * @headerfile <seqan/modifier.h>
 * @brief Modifier for the reverse-complement of a RnaString.
 *
 * @signature typedef ModifiedString<ModifiedString<RnaString, ModView<FunctorComplement<Rna5> >, ModReverse> RnaStringReverseComplement;
 */

/*!
 * @typedef Rna5StringReverseComplement
 * @headerfile <seqan/modifier.h>
 * @brief Modifier for the reverse-complement of a Rna5String.
 *
 * @signature typedef ModifiedString<ModifiedString<Rna5String, ModView<FunctorComplement<Rna55> >, ModReverse> Rna5StringReverseComplement;
 */

typedef ModView<FunctorComplement<Dna> >    ModComplementDna;
typedef ModView<FunctorComplement<Dna5> >    ModComplementDna5;
typedef ModView<FunctorComplement<Rna> >    ModComplementRna;
typedef ModView<FunctorComplement<Rna5> >    ModComplementRna5;

template <typename THost>
struct ReverseString
{
    typedef ModifiedString<THost, ModReverse> Type;
};

template <typename THost>
struct ComplementString
{
    typedef ModifiedString<THost, ModView<FunctorComplement<typename Value<THost>::Type> > > Type;
};

template <typename THost>
struct ReverseComplementString
{
    typedef typename ReverseString<typename ComplementString<THost>::Type>::Type Type;
};


// outdated shortcuts (better use metafunctions above)
typedef ComplementString<DnaString>::Type           DnaStringComplement;
typedef ComplementString<Dna5String>::Type          Dna5StringComplement;
typedef ComplementString<RnaString>::Type           RnaStringComplement;
typedef ComplementString<Rna5String>::Type          Rna5StringComplement;

typedef ReverseString<DnaString>::Type              DnaStringReverse;
typedef ReverseString<Dna5String>::Type             Dna5StringReverse;
typedef ReverseString<RnaString>::Type              RnaStringReverse;
typedef ReverseString<Rna5String>::Type             Rna5StringReverse;

typedef ReverseComplementString<DnaString>::Type    DnaStringReverseComplement;
typedef ReverseComplementString<Dna5String>::Type    Dna5StringReverseComplement;
typedef ReverseComplementString<RnaString>::Type    RnaStringReverseComplement;
typedef ReverseComplementString<Rna5String>::Type    Rna5StringReverseComplement;

// --------------------------------------------------------------------------
// Function complement()
// --------------------------------------------------------------------------

/*!
 * @fn complement
 * @headerfile <seqan/modifier.h>
 * @brief Complement a sequence or a StringSet in-place.
 *
 * @signature void complement(sequence);
 * @signature void complement(stringSet);
 *
 * @param[in,out] stringSet The StringSet to complement.
 * @param[in,out] sequence  The ContainerConcept to complement.
 *
 * Complementing only works for alphabets where a <tt>FunctorComplement&lt;T&gt;</tt> is implemented.  SeqAn comes with
 * such definitions for Dna, Dna5, Rna, Rna5.
 */

template <typename TSequence >
inline void complement(TSequence & sequence)
{
    convert(sequence, FunctorComplement<typename Value<TSequence>::Type>());
}

template <typename TSequence >
inline void complement(TSequence const & sequence)
{
    convert(sequence, FunctorComplement<typename Value<TSequence>::Type>());
}

// --------------------------------------------------------------------------
// Function complement()
// --------------------------------------------------------------------------

template < typename TSequence, typename TSpec >
inline void complement(StringSet<TSequence, TSpec> & stringSet)
{
    unsigned seqCount = length(stringSet);
    for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo)
        complement(stringSet[seqNo]);
}

template < typename TSequence, typename TSpec >
inline void complement(StringSet<TSequence, TSpec> const & stringSet)
{
    unsigned seqCount = length(stringSet);
    for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo)
        complement(stringSet[seqNo]);
}

// --------------------------------------------------------------------------
// Function complementString()
// --------------------------------------------------------------------------

template <typename THost>
inline ModifiedString<THost, ModView<FunctorComplement<typename Value<THost>::Type> > >
complementString(THost & host)
{
    return ModifiedString<THost, ModView<FunctorComplement<typename Value<THost>::Type> > >(host);
}

template <typename THost>
inline ModifiedString<THost const, ModView<FunctorComplement<typename Value<THost>::Type> > >
complementString(THost const & host)
{
    return ModifiedString<THost const, ModView<FunctorComplement<typename Value<THost>::Type> > >(host);
}

// --------------------------------------------------------------------------
// Function reverseComplement()
// --------------------------------------------------------------------------

/*!
 * @fn reverseComplement
 * @headerfile <seqan/modifier.h>
 * @brief Reverse-complement a sequence or a StringSet in-place.
 *
 * @signature void reverseComplement(sequence);
 * @signature void reverseComplement(stringSet);
 *
 * @param[in,out] stringSet The StringSet to complement.
 * @param[in,out] sequence  The ContainerConcept to complement.
 *
 * StringSet objects are reverse-complemented element-wise, i.e. the order of the sequences stays the same but the
 * sequences in the StringSet themselves are reverse-complemented.
 *
 * Complementing only works for alphabets where a <tt>FunctorComplement&lt;T&gt;</tt> is implemented.  SeqAn comes with
 * such definitions for Dna, Dna5, Rna, Rna5.
 */


template < typename TSequence, typename TParallelTag >
inline void reverseComplement(TSequence & sequence, Tag<TParallelTag> parallelTag)
{
    complement(sequence);
    reverse(sequence, parallelTag);
}

// TODO(holtgrew): How is doing anything in-place on a const value possible?
// (weese:) it is possible for rvalue references like temporary Segments/ModifiedStrings
template < typename TSequence, typename TParallelTag >
inline void reverseComplement(TSequence const & sequence, Tag<TParallelTag> parallelTag)
{
    complement(sequence);
    reverse(sequence, parallelTag);
}

template < typename TSequence, typename TSpec, typename TParallelTag >
inline void reverseComplement(StringSet<TSequence, TSpec> & stringSet, Tag<TParallelTag>)
{
    int seqCount = length(stringSet);
    SEQAN_OMP_PRAGMA(parallel for if(IsSameType<Tag<TParallelTag>, Parallel>::VALUE))
    for(int seqNo = 0; seqNo < seqCount; ++seqNo)
        reverseComplement(stringSet[seqNo], Serial());
}

// TODO(holtgrew): How is doing anything in-place on a const value possible?
// (weese:) it is possible for rvalue references like temporary Segments/ModifiedStrings
template < typename TSequence, typename TSpec, typename TParallelTag >
inline void reverseComplement(StringSet<TSequence, TSpec> const & stringSet, Tag<TParallelTag>)
{
    int seqCount = length(stringSet);
    SEQAN_OMP_PRAGMA(parallel for if(IsSameType<Tag<TParallelTag>, Parallel>::VALUE))
    for(int seqNo = 0; seqNo < seqCount; ++seqNo)
        reverseComplement(stringSet[seqNo], Serial());
}

template <typename TText>
inline void reverseComplement(TText & text)
{
    reverseComplement(text, Serial());
}

template <typename TText>
inline void reverseComplement(TText const & text)
{
    reverseComplement(text, Serial());
}

// --------------------------------------------------------------------------
// Function reverseComplementString()
// --------------------------------------------------------------------------

template <typename THost>
inline typename ReverseComplementString<THost>::Type
reverseComplementString(THost & host)
{
    typedef typename ReverseComplementString<THost>::Type TRevComp;
    return TRevComp(host);
}

template <typename THost>
inline typename ReverseComplementString<THost const>::Type
reverseComplementString(THost const & host)
{
    typedef typename ReverseComplementString<THost const>::Type TRevComp;
    return TRevComp(host);
}

// --------------------------------------------------------------------------
// Function toLower()
// --------------------------------------------------------------------------

/*!
 * @fn toLower
 * @headerfile <seqan/modifier.h>
 * @brief Convert characters in sequence or StringSet to lower case in-place.
 *
 * @signature void toLower(sequence);
 * @signature void toLower(stringSet);
 *
 * @param[in,out] sequence  The ContainerConcept to convert to lower-case.
 * @param[in,out] stringSet The StringSet to convert to lower-case.
 */

template < typename TSequence >
inline void toLower(TSequence & sequence)
{
    convert(sequence, FunctorLowcase<typename Value<TSequence>::Type>());
}

// TODO(holtgrew): How is doing anything in-place on a const value possible?
// (weese:) it is possible for rvalue references like temporary Segments/ModifiedStrings
template < typename TSequence >
inline void toLower(TSequence const & sequence)
{
    convert(sequence, FunctorLowcase<typename Value<TSequence>::Type>());
}

template < typename TSequence, typename TSpec >
inline void toLower(StringSet<TSequence, TSpec> & stringSet)
{
    unsigned seqCount = length(stringSet);
    for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo)
        toLower(stringSet[seqNo]);
}

// TODO(holtgrew): How is doing anything in-place on a const value possible?
// (weese:) it is possible for rvalue references like temporary Segments/ModifiedStrings
template < typename TSequence, typename TSpec >
inline void toLower(StringSet<TSequence, TSpec> const & stringSet)
{
    unsigned seqCount = length(stringSet);
    for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo)
        toLower(stringSet[seqNo]);
}

// --------------------------------------------------------------------------
// Function toUpper()
// --------------------------------------------------------------------------

/*!
 * @fn toUpper
 * @headerfile <seqan/modifier.h>
 * @brief Convert characters in sequence or StringSet to upper case in-place.
 *
 * @signature void toUpper(sequence);
 * @signature void toUpper(stringSet);
 *
 * @param[in,out] sequence  The ContainerConcept to convert to upper-case.
 * @param[in,out] stringSet The StringSet to convert to upper-case.
 */

template < typename TSequence >
inline void toUpper(TSequence & sequence)
{
    convert(sequence, FunctorUpcase<typename Value<TSequence>::Type>());
}

// TODO(holtgrew): How is doing anything in-place on a const value possible?
// (weese:) it is possible for rvalue references like temporary Segments/ModifiedStrings
template < typename TSequence >
inline void toUpper(TSequence const & sequence)
{
    convert(sequence, FunctorUpcase<typename Value<TSequence>::Type>());
}

template < typename TSequence, typename TSpec >
inline void toUpper(StringSet<TSequence, TSpec> & stringSet)
{
    unsigned seqCount = length(stringSet);
    for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo)
        toUpper(stringSet[seqNo]);
}

// TODO(holtgrew): How is doing anything in-place on a const value possible?
// (weese:) it is possible for rvalue references like temporary Segments/ModifiedStrings
template < typename TSequence, typename TSpec >
inline void toUpper(StringSet<TSequence, TSpec> const & stringSet)
{

    unsigned seqCount = length(stringSet);
    for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo)
        toUpper(stringSet[seqNo]);
}

// --------------------------------------------------------------------------
// Function lowerString()
// --------------------------------------------------------------------------

template <typename THost>
inline ModifiedString<THost const, ModView<FunctorLowcase<typename Value<THost>::Type> > >
lowerString(THost const & host)
{
    return ModifiedString<THost const, ModView<FunctorLowcase<typename Value<THost>::Type> > >(host);
}

// --------------------------------------------------------------------------
// Function upperString()
// --------------------------------------------------------------------------

template <typename THost>
inline ModifiedString<THost const, ModView<FunctorUpcase<typename Value<THost>::Type> > >
upperString(THost const & host)
{
    return ModifiedString<THost const, ModView<FunctorUpcase<typename Value<THost>::Type> > >(host);
}

}

#endif
