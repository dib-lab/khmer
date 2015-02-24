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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// Code for Dna(5) to AminoAcid Translation
// ==========================================================================


#ifndef INCLUDE_SEQAN_TRANSLATION_TRANSLATION_H_
#define INCLUDE_SEQAN_TRANSLATION_TRANSLATION_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// -----------------------------------------------------------------------
// Enum TranslationFrames
// -----------------------------------------------------------------------

/*!
 * @enum TranslationFrames
 * @headerfile <seqan/translation.h>
 * @brief Class Enum with frames for @link translate @endlink()
 *
 * @signature enum class TranslationFrames : uint8_t { ... };
 *
 * Please not that this is part of the translation module which requires C++11.
 *
 * @val TranslationFrames SINGLE_FRAME = 0;
 * @brief Translate the sequence(s) "as is", n input sequences result in n output sequences.
 *
 * @val TranslationFrames WITH_REVERSE_COMPLEMENT = 1;
 * @brief Translate the sequence(s) as well as their reverse complements (n -> * 2n).
 *
 * @val TranslationFrames WITH_FRAME_SHIFTS = 2;
 * @brief Translate the sequence(s) as well as their shifted frames (n -> 3n).
 *
 * @val TranslationFrames SIX_FRAME = 3;
 * @brief Equals (WITH_REVERSE_COMPLEMENT | WITH_FRAME_SHIFTS); shifted frames of original and reverse complement are
 *        translated (n -> 6n).
 */

enum TranslationFrames
{
    SINGLE_FRAME             = 0,
    WITH_REVERSE_COMPLEMENT  = 1,
    WITH_FRAME_SHIFTS        = 2,
    SIX_FRAME                = 3
};

// -----------------------------------------------------------------------
// Tag Frames_ (internal)
// -----------------------------------------------------------------------

template <uint8_t num>
struct Frames_
{};

// ============================================================================
// Metafunctions
// ============================================================================

// -----------------------------------------------------------------------
// Metafunction ReverseComplement_
// -----------------------------------------------------------------------

// type of the reverse complement of a string, also works with infixes and
// ModifiedStrings

template <typename T>
struct ReverseComplement_
{
    typedef ModifiedString<
        ModifiedString<T, ModView<
                           FunctorComplement<
                               typename Value<T>::Type > > >, ModReverse> Type;
};

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function _ord()
// --------------------------------------------------------------------------

// returns ordValue of a DNA(5) or RNA(5) character
// for everything else (e.g. char) the character is converted to Dna5 first

template <typename T>
inline typename ValueSize<T>::Type
_ord(T const & c)
{
    return ordValue(Dna5(c));
}

inline ValueSize<Dna>::Type
_ord(Dna const & c)
{
    return ordValue(c);
}

inline ValueSize<Dna5>::Type
_ord(Dna5 const & c)
{
    return ordValue(c);
}

inline ValueSize<Rna>::Type
_ord(Rna const & c)
{
    return ordValue(c);
}

inline ValueSize<Rna5>::Type
_ord(Rna5 const & c)
{
    return ordValue(c);
}

// --------------------------------------------------------------------------
// Function _translateTriplet()
// --------------------------------------------------------------------------


template <typename TOrd, GeneticCodeSpec CODE_SPEC>
inline AminoAcid
_translateTriplet(TOrd const & c1,
                  TOrd const & c2,
                  TOrd const & c3,
                  GeneticCode<CODE_SPEC> const & /**/)
{
    return (( c1 > 3 ) || ( c2 > 3 ) || ( c3 > 3 ) )
            ? 'X'
            : TranslateTableDnaToAminoAcid_<
                GeneticCode<CODE_SPEC> >::VALUE[c1][c2][c3];
}

// --------------------------------------------------------------------------
// Function _translateString()
// --------------------------------------------------------------------------

template <typename TOutString, typename TInString, GeneticCodeSpec CODE_SPEC>
inline void
_translateString(TOutString & target,
                 TInString const & source,
                 GeneticCode<CODE_SPEC> const & /**/)
{
    SEQAN_ASSERT_EQ(length(source)/3, length(target));
    typedef typename Position<TInString>::Type TPos;

    for (TPos i = 0; i+2 < length(source); i+=3)
    {
        target[i/3] = _translateTriplet(_ord(value(source, i  )),
                                        _ord(value(source, i+1)),
                                        _ord(value(source, i+2)),
                                        GeneticCode<CODE_SPEC>());
    }
}

template <typename TOutString, typename TSpec, typename TInString, GeneticCodeSpec CODE_SPEC>
inline void
_translateString(Segment<TOutString, TSpec> SEQAN_FORWARD_RETURN target,
                 TInString const & source,
                 GeneticCode<CODE_SPEC> const & /**/)
{
    SEQAN_ASSERT_EQ(length(source)/3, length(target));
    typedef typename Position<TInString>::Type TPos;

    for (TPos i = 0; i+2 < length(source); i+=3)
    {
        target[i/3] = _translateTriplet(_ord(value(source, i  )),
                                        _ord(value(source, i+1)),
                                        _ord(value(source, i+2)),
                                        GeneticCode<CODE_SPEC>());
    }
}

// --------------------------------------------------------------------------
// Function _translateImplLoop()
// --------------------------------------------------------------------------

// single frame
template <typename TSpec1, typename TSpec2, typename TSpec3, typename TInString,
          GeneticCodeSpec CODE_SPEC>
inline void
_translateImplLoop(StringSet<String<AminoAcid, TSpec1>, TSpec2> & target,
                   unsigned const i,
                   StringSet<TInString, TSpec3> const & source,
                   GeneticCode<CODE_SPEC> const & /**/,
                   Frames_<1u> const & /**/)
{
    typedef GeneticCode<CODE_SPEC> TCode;
    _translateString(target[i], source[i], TCode());
}

// with reverse complement
template <typename TSpec1, typename TSpec2, typename TSpec3, typename TInString,
          GeneticCodeSpec CODE_SPEC>
inline void
_translateImplLoop(StringSet<String<AminoAcid, TSpec1>, TSpec2> & target,
                   unsigned const i,
                   StringSet<TInString, TSpec3> const & source,
                   GeneticCode<CODE_SPEC> const & /**/,
                   Frames_<2u> const & /**/)
{
    typedef typename Value<StringSet<TInString, TSpec3> const>::Type TVal;
    typedef typename ReverseComplement_<TVal>::Type TRevComp;
    typedef GeneticCode<CODE_SPEC> TCode;

    if (i % 2)
    {
        TRevComp revComp(value(source, i/2));
        _translateString(target[i], revComp, TCode());
    }
    else
    {
        _translateString(target[i], source[i/2], TCode());
    }
}

// three frame
template <typename TSpec1, typename TSpec2, typename TSpec3, typename TInString,
          GeneticCodeSpec CODE_SPEC>
inline void
_translateImplLoop(StringSet<String<AminoAcid, TSpec1>, TSpec2> & target,
                   unsigned const i,
                   StringSet<TInString, TSpec3> const & source,
                   GeneticCode<CODE_SPEC> const & /**/,
                   Frames_<3u> const & /**/)
{
    typedef GeneticCode<CODE_SPEC> TCode;
    _translateString(target[i], suffix(source[i/3], i % 3), TCode());
}

// six frame
template <typename TSpec1, typename TSpec2, typename TSpec3, typename TInString,
          GeneticCodeSpec CODE_SPEC>
inline void
_translateImplLoop(StringSet<String<AminoAcid, TSpec1>, TSpec2> & target,
                   unsigned const i,
                   StringSet<TInString, TSpec3> const & source,
                   GeneticCode<CODE_SPEC> const & /**/,
                   Frames_<6u> const & /**/)
{
    typedef typename Prefix<
        typename Value<
            StringSet<TInString, TSpec3> const>::Type>::Type TVal;
    typedef typename ReverseComplement_<TVal>::Type TRevComp;
    typedef GeneticCode<CODE_SPEC> TCode;

    if ((i % 6) > 2)
    {
        TRevComp revComp(prefix(value(source, i/6),
                                length(value(source,i/6)) - (i % 3)));
        _translateString(target[i], revComp, TCode());
    }
    else
    {
        _translateString(target[i], suffix(source[i/6], i % 3), TCode());
    }

}

// --------------------------------------------------------------------------
// Function _translateImplLoopOMPWrapper()
// --------------------------------------------------------------------------


template <typename TSource, typename TTarget, uint8_t frames,
          GeneticCodeSpec CODE_SPEC>
inline void
_translateImplLoopOMPWrapper(TTarget & target,
                             TSource const & source,
                             GeneticCode<CODE_SPEC> const & /**/,
                             Frames_<frames> const & /**/,
                             Parallel const & /**/)
{
    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic))
    for (int64_t i = 0; i < static_cast<int64_t>(length(target)); ++i)
        _translateImplLoop(target, i, source, GeneticCode<CODE_SPEC>(),
                           Frames_<frames>());
}

template <typename TSource, typename TTarget, uint8_t frames,
          GeneticCodeSpec CODE_SPEC>
inline void
_translateImplLoopOMPWrapper(TTarget & target,
                             TSource const & source,
                             GeneticCode<CODE_SPEC> const & /**/,
                             Frames_<frames> const & /**/,
                             Serial const & /**/)
{
    typedef typename Size<TTarget>::Type TPos;
    for (TPos i = 0; i < length(target); ++i)
        _translateImplLoop(target, i, source, GeneticCode<CODE_SPEC>(),
                           Frames_<frames>());
}

// --------------------------------------------------------------------------
// Function _translateImpl()
// --------------------------------------------------------------------------

// general case
template <typename TSpec1, typename TSpec2, typename TSpec3, typename TInString,
          typename TParallelism, GeneticCodeSpec CODE_SPEC, unsigned char n>
inline void
_translateImpl(StringSet<String<AminoAcid, TSpec1>, TSpec2> & target,
               StringSet<TInString, TSpec3> const & source,
               GeneticCode<CODE_SPEC> const & /**/,
               Frames_<n> const & /**/,
               TParallelism const & /**/)
{
    typedef typename Position<StringSet<TInString, TSpec3> >::Type TPos;

    resize(target, length(source) * n, Exact());
    for (TPos i = 0; i < length(target); ++i)
    {
        // current dnastring's length / 3 (3DNA -> 1 AA)
        TPos len = length(source[i/n]) / 3;
        // shorten for shifted frames
        if (( n > 2 ) && ( length(source[i/n]) % 3 ) < ( i%3 ))
            --len;
        resize(target[i], len, Exact());
    }

    _translateImplLoopOMPWrapper(target, source, GeneticCode<CODE_SPEC>(),
                                 Frames_<n>(),
                                 TParallelism());
}

// ConcatDirect specialization
template <typename TSpec1, typename TSpec3, typename TInString,
          typename TParallelism, GeneticCodeSpec CODE_SPEC, unsigned char n>
inline void
_translateImpl(StringSet<String<AminoAcid,
                                TSpec1>, Owner<ConcatDirect<> > > & target,
               StringSet<TInString, TSpec3> const & source,
               GeneticCode<CODE_SPEC> const & /**/,
               Frames_<n> const & /**/,
               TParallelism const & /**/)
{
    typedef typename Position<StringSet<TInString, TSpec3> >::Type TPos;

    resize(target.limits, length(source) * n + 1, Exact());
    target.limits[0] = 0;
    for (TPos i = 0; i+1 < length(target.limits); ++i)
    {
        // current dnastring's length / 3 (3DNA -> 1 AA)
        TPos len = length(source[i/n]) / 3;
        // shorten for shifted frames
        if (( n > 2 ) && ( length(source[i/n]) % 3 ) < ( i%3 ))
            --len;
        target.limits[i+1] = target.limits[i] + len;
    }

    resize(target.concat, back(target.limits), Exact());

    _translateImplLoopOMPWrapper(target, source, GeneticCode<CODE_SPEC>(),
                                 Frames_<n>(),
                                 TParallelism());
}

// --------------------------------------------------------------------------
// Function _translateInputWrap()
// --------------------------------------------------------------------------

// stringset to stringset
template <typename TSpec1, typename TSpec2, typename TSpec3, typename TInString,
          typename TParallelism, GeneticCodeSpec CODE_SPEC, unsigned char n>
inline void
_translateInputWrap(StringSet<String<AminoAcid, TSpec1>, TSpec2> & target,
                    StringSet<TInString, TSpec3> const & source,
                    GeneticCode<CODE_SPEC> const & /**/,
                    Frames_<n> const & /**/,
                    TParallelism const & /**/)
{
    _translateImpl(target, source, GeneticCode<CODE_SPEC>(), Frames_<n>(),
                   TParallelism());
}

// single string to stringset conversion
template <typename TSpec1, typename TSpec2, typename TInString,
          typename TParallelism, GeneticCodeSpec CODE_SPEC, unsigned char n>
inline void
_translateInputWrap(StringSet<String<AminoAcid, TSpec1>, TSpec2> & target,
                    TInString const & source,
                    GeneticCode<CODE_SPEC> const & /**/,
                    Frames_<n> const & /**/,
                    TParallelism const & /**/)
{
    StringSet<TInString, Dependent<> > set;
    appendValue(set, source);
    _translateImpl(target, set, GeneticCode<CODE_SPEC>(), Frames_<n>(),
                   TParallelism());
}

//bail out because multiple frames don't fit in one string
template <typename TSpec1, typename TInString, typename TParallelism,
          GeneticCodeSpec CODE_SPEC, unsigned char n>
inline void
_translateInputWrap(String<AminoAcid, TSpec1> & /**/,
                    TInString const & /**/,
                    GeneticCode<CODE_SPEC> const & /**/,
                    Frames_<n> const & /**/,
                    TParallelism const & /**/)
{
    SEQAN_FAIL("Implementation error, multiple frames selected, but only a "
               "singe target string.");
}

// single string to single string conversion
template <typename TSpec1, typename TInString, typename TParallelism,
          GeneticCodeSpec CODE_SPEC>
inline void
_translateInputWrap(String<AminoAcid, TSpec1> & target,
                    TInString const & source,
                    GeneticCode<CODE_SPEC> const & /**/,
                    Frames_<1> const & /**/,
                    TParallelism const & /**/)
{
    resize(target, length(source)/3, Exact());
    _translateString(target, source, GeneticCode<CODE_SPEC>());
}

// --------------------------------------------------------------------------
// Function translate()
// --------------------------------------------------------------------------

/*!
 * @fn translate
 * @headerfile <seqan/translation.h>
 * @brief translate sequences of Dna or Rna into amino acid alphabet, optionally with frames
 * @signature int translate(target, source[, frames][, geneticCode][, TParallelism])
 * @signature int translate(target, source[, frames][, geneticCodeSpec][, TParallelism])
 *
 * @param[out]      target      The amino acid sequence(s).  @link StringSet @endlink of @link AminoAcid @endlink
 *                              or @link String @endlink of @link AminoAcid @endlink if source is a single string
 *                              and frames is <tt>SINGLE_FRAME</tt>.
 * @param[in]       source      Source sequences @link String @endlink or @link StringSet @endlink.
 *                              If the value type is not Dna, Dna5, Rna, Rna5 then it is converted
 *                              to Dna5.
 * @param[in]       frame       The @link TranslationFrames @endlink, defaults to SINGLE_FRAME.
 * @param[in]       geneticCode The @link GeneticCode @endlink to use, defaults
 *                              to <tt>GeneticCode&lt;CANONICAL&gt;</tt>
 *                              (use to specify GeneticCode at compile-time)
 * @param[in]       geneticCodeSpec The @link GeneticCodeSpec @endlink to use
 *                              (use to specify GenetiCode at run-time)
 * @param[in]       TParallelism Whether to use SMP or not, see @link ParallelismTags @endlink .
 *
 * @return int 0 on success, and -1 on incompatible parameters (e.g. multiple frames but target type not StringSet).
 *
 * If OpenMP is supported by platform and TParallelism is not specified as
 * "Serial", translation will be parallelized. The only exception is when doing
 * single-frame translation of a single string -- which is never parallelized.
 *
 * The translation process is fastest when using ConcatDirect-StringSets for
 * both input and output StringSets and when not having to convert the alphabet
 * of the source, i.e. feeding AminoAcid-Strings, not CharStrings (although
 * the latter is possible).
 *
 * Please note that specifying the GeneticCode at compile time avoids having
 * unrequired conversion tables in memory. The implementation profits slightly
 * from having SEQAN_CXX11_STANDARD defined.
 * @section Example
 *
 * @code{.cpp}
 * StringSet<Dna5String> dnaSeqs;
 *
 * // do something that fills up dnaSeqs, e.g. read from file or assign
 *
 * StringSet<String<AminoAcid>, Owner<ConcatDirect<> > > aaSeqs;
 *
 * translate(aaSeqs, dnaSeqs, SIX_FRAME);
 *
 * // do something with the aaSeqs
 * @endcode
 *
 * @see TranslationFrames
 * @see GeneticCode
 */

template <typename TTarget, typename TSource, typename TParallelism,
          GeneticCodeSpec CODE_SPEC>
inline void
translate(TTarget & target,
          TSource const & source,
          TranslationFrames const frames,
          GeneticCode<CODE_SPEC> const & /**/,
          TParallelism const & /**/)
{
    typedef GeneticCode<CODE_SPEC> TCode;
    switch (frames)
    {
    case SINGLE_FRAME:
        return _translateInputWrap(target, source, TCode(), Frames_<1>(),
                                   TParallelism());
    case WITH_REVERSE_COMPLEMENT:
        return _translateInputWrap(target, source, TCode(), Frames_<2>(),
                                   TParallelism());
    case WITH_FRAME_SHIFTS:
        return _translateInputWrap(target, source, TCode(), Frames_<3>(),
                                   TParallelism());
    case SIX_FRAME:
        return _translateInputWrap(target, source, TCode(), Frames_<6>(),
                                   TParallelism());
    }
}

template <typename TTarget, typename TSource, GeneticCodeSpec CODE_SPEC>
inline void
translate(TTarget & target,
          TSource const & source,
          TranslationFrames const frames,
          GeneticCode<CODE_SPEC> const & /**/)
{
    return translate(target, source, frames, GeneticCode<CODE_SPEC>(),
                     Parallel());
}

template <typename TTarget, typename TSource, GeneticCodeSpec CODE_SPEC>
inline void
translate(TTarget & target,
          TSource const & source,
          GeneticCode<CODE_SPEC> const & /**/)
{
    return translate(target, source, SINGLE_FRAME,
                     GeneticCode<CODE_SPEC>(), Parallel());
}

template <typename TTarget, typename TSource>
inline void
translate(TTarget & target,
          TSource const & source,
          TranslationFrames const frames)
{
    return translate(target, source, frames, GeneticCode<>(), Parallel());
}

template <typename TTarget, typename TSource>
inline void
translate(TTarget & target,
          TSource const & source)
{
    return translate(target, source, SINGLE_FRAME,
                     GeneticCode<>(), Parallel());
}

template <typename TTarget, typename TSource, typename TParallelism>
inline void
translate(TTarget & target,
          TSource const & source,
          TranslationFrames const frames,
          TParallelism const & /**/)
{
    return translate(target, source, frames, GeneticCode<>(), TParallelism());
}

// -----------------------------------------------------------------------
// Function translate() with run-time spec selection
// -----------------------------------------------------------------------

template <typename TTarget, typename TSource, typename TParallelism,
          GeneticCodeSpec currentSpec, typename TRestList>
inline void
_translate(TTarget & /**/,
           TSource const & /**/,
           TranslationFrames const /**/,
           GeneticCodeSpec const & /**/,
           TagList<GeneticCode<currentSpec>, TRestList> const & /**/,
           TParallelism const & /**/,
           True const & /**/)
{
    SEQAN_FAIL("Invalid genetic code translation table selected.\n");
}

// forward declare because of double-recursion
template <typename TTarget, typename TSource, typename TParallelism,
          GeneticCodeSpec currentSpec, typename TRestList>
inline void
_translate(TTarget & target,
           TSource const & source,
           TranslationFrames const frames,
           GeneticCodeSpec const & geneticCodeSpec,
           TagList<GeneticCode<currentSpec>, TRestList> const & /**/,
           TParallelism const & /**/);

template <typename TTarget, typename TSource, typename TParallelism,
          GeneticCodeSpec currentSpec, typename TRestList>
inline void
_translate(TTarget & target,
           TSource const & source,
           TranslationFrames const frames,
           GeneticCodeSpec const & geneticCodeSpec,
           TagList<GeneticCode<currentSpec>, TRestList> const & /**/,
           TParallelism const & /**/,
           False const & /**/)
{
    return _translate(target, source, frames, geneticCodeSpec, TRestList(),
                      TParallelism());
}

template <typename TTarget, typename TSource, typename TParallelism,
          GeneticCodeSpec currentSpec, typename TRestList>
inline void
_translate(TTarget & target,
           TSource const & source,
           TranslationFrames const frames,
           GeneticCodeSpec const & geneticCodeSpec,
           TagList<GeneticCode<currentSpec>, TRestList> const & /**/,
           TParallelism const & /**/)
{
    typedef TagList<GeneticCode<currentSpec>, TRestList> TTagList;

    if (geneticCodeSpec == currentSpec)
        return translate(target, source, frames, GeneticCode<currentSpec>(),
                         TParallelism());

    return _translate(target, source, frames, geneticCodeSpec,
                             TTagList(), TParallelism(),
                             typename IsSameType<TRestList, void>::Type());
}

template <typename TTarget, typename TSource, typename TParallelism>
inline void
translate(TTarget & target,
          TSource const & source,
          TranslationFrames const frames,
          GeneticCodeSpec const & geneticCodeSpec,
          TParallelism const & /**/)
{
    return _translate(target, source, frames, geneticCodeSpec, GeneticCodes_(),
                      TParallelism());
}

// convenience
template <typename TTarget, typename TSource>
inline void
translate(TTarget & target,
          TSource const & source,
          TranslationFrames const frames,
          GeneticCodeSpec const & geneticCode)
{
    return translate(target, source, frames, geneticCode, Parallel());
}

}

#endif // INCLUDE_SEQAN_TRANSLATION_TRANSLATION_H_
