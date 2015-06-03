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
// Shortcuts for certain often-used string types.
// ==========================================================================

#ifndef SEQAN_HEADER_SEQUENCE_SHORTCUTS_H
#define SEQAN_HEADER_SEQUENCE_SHORTCUTS_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

/*!
 * @typedef CharString
 * @headerfile <seqan/sequence.h>
 * @brief An AllocString with char alphabet.
 *
 * You can efficiently cast a CharString to <tt>char *</tt> using <tt>toCString</tt>.
 *
 * @signature typedef String<char> CharString;
 */

typedef String<char, Alloc<void> > CharString;

//____________________________________________________________________________

/*!
 * @typedef CharIterator
 * @headerfile <seqan/sequence.h>
 * @brief An iterator overa a @link CharString @endlink.
 *
 * @signature typedef Iterator<CharString, Rooted>::Type CharIterator;
 */

typedef Iterator<CharString, Rooted>::Type CharIterator;

//////////////////////////////////////////////////////////////////////////////

/*!
 * @typedef UnicodeString
 * @headerfile <seqan/sequence.h>
 * @brief An AllocString with wchar_t alphabet.
 *
 * @signature typedef String<wchar_t> CharString;
 */

typedef String<wchar_t, Alloc<void> > UnicodeString;

//____________________________________________________________________________

typedef Iterator<UnicodeString, Rooted>::Type UnicodeIterator;

//////////////////////////////////////////////////////////////////////////////

/*!
 * @typedef DnaString
 * @headerfile <seqan/sequence.h>
 * @brief An AllocString with Dna alphabet.
 *
 * @signature typedef String<Dna> CharString;
 */

typedef String<Dna, Alloc<void> > DnaString;

//____________________________________________________________________________

/*!
 * @typedef DnaIterator
 * @headerfile <seqan/sequence.h>
 * @brief A rooted iterator over a Dna.
 *
 * @signature typedef Iterator<DnaString, Rooted>::Type DnaIterator;
 */

typedef Iterator<DnaString, Rooted>::Type DnaIterator;

//////////////////////////////////////////////////////////////////////////////

/*!
 * @typedef Dna5String
 * @headerfile <seqan/sequence.h>
 * @brief An AllocString with Dna5 alphabet.
 *
 * @signature typedef String<Dna5> CharString;
 */

typedef String<Dna5, Alloc<void> > Dna5String;

//////////////////////////////////////////////////////////////////////////////

/*!
 * @typedef DnaQString
 * @headerfile <seqan/sequence.h>
 * @brief An AllocString with DnaQ alphabet.
 *
 * @signature typedef String<DnaQ> CharString;
 */

typedef String<DnaQ, Alloc<void> > DnaQString;

//////////////////////////////////////////////////////////////////////////////

/*!
 * @typedef Dna5QString
 * @headerfile <seqan/sequence.h>
 * @brief An AllocString with Dna5Q alphabet.
 *
 * @signature typedef String<Dna5Q> CharString;
 */

typedef String<Dna5Q, Alloc<void> > Dna5QString;

//____________________________________________________________________________

/*!
 * @typedef Dna5Iterator
 * @headerfile <seqan/sequence.h>
 * @brief A rooted iterator over a Dna5.
 *
 * @signature typedef Iterator<Dna5String, Rooted>::Type Dna5Iterator;
 */

typedef Iterator<Dna5String, Rooted>::Type Dna5Iterator;

//////////////////////////////////////////////////////////////////////////////

/*!
 * @typedef RnaString
 * @headerfile <seqan/sequence.h>
 * @brief An AllocString with Rna alphabet.
 *
 * @signature typedef String<Rna> CharString;
 */

typedef String<Rna, Alloc<void> > RnaString;

//____________________________________________________________________________

/*!
 * @typedef RnaIterator
 * @headerfile <seqan/sequence.h>
 * @brief A rooted iterator over a Rna.
 *
 * @signature typedef Iterator<RnaString, Rooted>::Type RnaIterator;
 */

typedef Iterator<RnaString, Rooted>::Type RnaIterator;

//////////////////////////////////////////////////////////////////////////////

/*!
 * @typedef Rna5String
 * @headerfile <seqan/sequence.h>
 * @brief An AllocString with Rna5 alphabet.
 *
 * @signature typedef String<Rna5> CharString;
 */

typedef String<Rna5, Alloc<void> > Rna5String;

//____________________________________________________________________________

/*!
 * @typedef Rna5Iterator
 * @headerfile <seqan/sequence.h>
 * @brief A rooted iterator over a Rna5.
 *
 * @signature typedef Iterator<Rna5String, Rooted>::Type Rna5Iterator;
 */

typedef Iterator<Rna5String, Rooted>::Type Rna5Iterator;

//////////////////////////////////////////////////////////////////////////////

/*!
 * @typedef IupacString
 * @headerfile <seqan/sequence.h>
 * @brief An AllocString with Iupac alphabet.
 *
 * @signature typedef String<Iupac> CharString;
 */

typedef String<Iupac, Alloc<void> > IupacString;

//____________________________________________________________________________

/*!
 * @typedef IupacIterator
 * @headerfile <seqan/sequence.h>
 * @brief A rooted iterator over a Iupac.
 *
 * @signature typedef Iterator<IupacString, Rooted>::Type IupacIterator;
 */

typedef Iterator<IupacString, Rooted>::Type IupacIterator;

//////////////////////////////////////////////////////////////////////////////

/*!
 * @typedef Peptide
 * @headerfile <seqan/sequence.h>
 * @brief An AllocString with AminoAcid alphabet.
 *
 * @signature typedef String<AminoAcid> Peptide;
 */

typedef String<AminoAcid, Alloc<void> > Peptide;

//____________________________________________________________________________

/*!
 * @typedef PeptideIterator
 * @headerfile <seqan/sequence.h>
 * @brief A rooted iterator over a Peptide.
 *
 * @signature typedef Iterator<Peptide, Rooted>::Type PeptideIterator;
 */

typedef Iterator<Peptide, Rooted>::Type PeptideIterator;

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
