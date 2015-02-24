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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// Input/Output on FASTA and FASTQ files.
// ==========================================================================

#ifndef SEQAN_SEQ_IO_FASTA_FASTQ_H_
#define SEQAN_SEQ_IO_FASTA_FASTQ_H_

namespace seqan {

// ============================================================================
// Tags
// ============================================================================

// --------------------------------------------------------------------------
// Tag Fasta
// --------------------------------------------------------------------------

struct TagFasta_;
typedef Tag<TagFasta_> Fasta;

// --------------------------------------------------------------------------
// Tag Fastq
// --------------------------------------------------------------------------

struct TagFastq_;
typedef Tag<TagFastq_> Fastq;

// --------------------------------------------------------------------------
// Tag Raw
// --------------------------------------------------------------------------

struct Raw_;
typedef Tag<Raw_> Raw;


// ============================================================================
// Metafunctions
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction MagicHeader
// --------------------------------------------------------------------------

// TODO(weese:) The following defines makes the old guessFormat functions in file_format_mmap.h obsolete. Disable them!
template <typename T>
struct MagicHeader<Fasta, T>
{
    static char const VALUE[1];
};
template <typename T>
char const MagicHeader<Fasta, T>::VALUE[1] = { '>' };  // Fasta's first character


template <typename T>
struct MagicHeader<Fastq, T>
{
    static char const VALUE[1];
};
template <typename T>
char const MagicHeader<Fastq, T>::VALUE[1] = { '@' };  // Fastq's first character


template <typename T>
struct MagicHeader<Raw, T> :
    MagicHeader<Nothing, T> {};

// --------------------------------------------------------------------------
// Metafunction FileExtensions
// --------------------------------------------------------------------------

template <typename T>
struct FileExtensions<Fasta, T>
{
    static char const * VALUE[2];
};
template <typename T>
char const * FileExtensions<Fasta, T>::VALUE[2] =
{
    ".fa",      // default output extension
    ".fasta"
//    ".faa",     // FASTA Amino Acid file
//    ".ffn",     // FASTA nucleotide coding regions file
//    ".fna",     // FASTA Nucleic Acid file
//    ".frn"
};


template <typename T>
struct FileExtensions<Fastq, T>
{
    static char const * VALUE[2];
};
template <typename T>
char const * FileExtensions<Fastq, T>::VALUE[2] =
{
    ".fq",      // default output extension
    ".fastq"
};


template <typename T>
struct FileExtensions<Raw, T>
{
    static char const * VALUE[1];
};
template <typename T>
char const * FileExtensions<Raw, T>::VALUE[1] =
{
    ".txt"      // default output extension
//    ".seq"
};

// ----------------------------------------------------------------------------
// Class FastaIgnoreFunctor_
// ----------------------------------------------------------------------------

template <typename TAlphabet>
struct FastaIgnoreFunctor_
{
    typedef typename If< Or< IsSameType<TAlphabet, char>,
                         Or< IsSameType<TAlphabet, signed char>,
                             IsSameType<TAlphabet, unsigned char> > >,
                         IsNewline,                     // ignore only newline if the target alphabet is a char
                         IsWhitespace                   // ignore whitespace as well for all other alphabets
                     >::Type Type;
};

template <typename TAlphabet>
struct FastaIgnoreOrAssertFunctor_
{
    typedef typename FastaIgnoreFunctor_<TAlphabet>::Type               TIgnore;
    typedef AssertFunctor<IsInAlphabet<TAlphabet>, ParseError, Fasta>   TAsserter;

    typedef typename If< Or< IsSameType<TAlphabet, char>,
                         Or< IsSameType<TAlphabet, signed char>,
                             IsSameType<TAlphabet, unsigned char> > >,
                         TIgnore,                       // don't assert in case of char alphabets
                         OrFunctor<TIgnore, TAsserter>  // assert being part of the alphabet for other alphabets
                     >::Type Type;
};


// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class SequenceOutputOptions
// ----------------------------------------------------------------------------

/*!
 * @class SequenceOutputOptions
 * @headerfile <seqan/seq_io.h>
 * @signature struct SequenceOutputOptions
 * @brief Configuration for writing sequence (FASTA/FASTQ) files.
 *
 * This struct is used for the configuration of writing out FASTA and FASTQ files.
 *
 * @var int SequenceOutputOptions::lineLength;
 * @brief Length of the lines when writing out.
 *
 * Set to <tt>-1</tt> for default behaviour (no line break for FASTQ, line length of 70 for FASTA) and <tt>0</tt> for
 * disabling line breaks.
 *
 * @var bool SequenceOutputOptions::qualMeta;
 * @brief Whether or not to write the meta information into the <tt>"+"</tt> line before the qualities (interpreted for
 *        FASTQ only). Default is <tt>false</tt>.
 */

// TODO(holtgrew): Would it be worth having two/three shortcuts for "short reads" and "genomic sequence" and faster or can the compiler optimize the creation away?

struct SequenceOutputOptions
{
    int lineLength;
    bool qualMeta;

    explicit
    SequenceOutputOptions(int lineLength = -1, bool qualMeta = false) :
        lineLength(lineLength),
        qualMeta(qualMeta)
    {}
};

// ----------------------------------------------------------------------------
// Class QualityExtractor
// ----------------------------------------------------------------------------

template <typename TValue>
struct QualityExtractor : public std::unary_function<TValue, char>
{
    inline char operator()(TValue const & x) const
    {
        return '!' + static_cast<char>(getQualityValue(x));
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecord(TagSelector); Qualities inside seq
// ----------------------------------------------------------------------------

template <typename TIdString, typename TSeqString, typename TFwdIterator>
inline void
readRecord(TIdString & /* meta */, TSeqString & /* seq */, TFwdIterator & /* iter */,
           TagSelector<> const & /* format */)
{}

template <typename TIdString, typename TSeqString, typename TFwdIterator, typename TTagList>
inline void
readRecord(TIdString & meta, TSeqString & seq, TFwdIterator & iter, TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormat;

    if (isEqual(format, TFormat()))
        readRecord(meta, seq, iter, TFormat());
    else
        readRecord(meta, seq, iter, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

// ----------------------------------------------------------------------------
// Function readRecord(TagSelector); Qualities inside qual
// ----------------------------------------------------------------------------

template <typename TIdString, typename TSeqString, typename TQualString, typename TFwdIterator>
inline void
readRecord(TIdString & /* meta */, TSeqString & /* seq */, TQualString & /* qual */, TFwdIterator & /* iter */,
           TagSelector<> const & /* format */)
{}

template <typename TIdString, typename TSeqString, typename TQualString, typename TFwdIterator, typename TTagList>
inline void
readRecord(TIdString & meta, TSeqString & seq, TQualString & qual, TFwdIterator & iter, TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormat;

    if (isEqual(format, TFormat()))
        readRecord(meta, seq, qual, iter, TFormat());
    else
        readRecord(meta, seq, qual, iter, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

// ----------------------------------------------------------------------------
// Function readRecord(Raw);
// ----------------------------------------------------------------------------

template <typename TIdString, typename TSeqString, typename TFwdIterator>
inline void readRecord(TIdString & meta, TSeqString & seq, TFwdIterator & iter, Raw)
{
    typedef typename Value<TSeqString>::Type                        TAlphabet;
    typedef typename FastaIgnoreOrAssertFunctor_<TAlphabet>::Type   TIgnoreOrAssert;

    clear(meta);
    clear(seq);
    readUntil(seq, iter, IsNewline(), TIgnoreOrAssert());
    skipLine(iter);
}

// ----------------------------------------------------------------------------
// Function readRecord(Raw);
// ----------------------------------------------------------------------------

template <typename TIdString, typename TSeqString, typename TQualString, typename TFwdIterator>
inline void readRecord(TIdString & meta, TSeqString & seq, TQualString & qual, TFwdIterator & iter, Raw const & raw)
{
    clear(qual);
    readRecord(meta, seq, iter, raw);
}

// ----------------------------------------------------------------------------
// Function readRecord(Fasta); Qualities inside seq
// ----------------------------------------------------------------------------

template <typename TIdString, typename TSeqString, typename TFwdIterator>
inline void readRecord(TIdString & meta, TSeqString & seq, TFwdIterator & iter, Fasta)
{
    typedef typename Value<TSeqString>::Type                        TAlphabet;
    typedef typename FastaIgnoreOrAssertFunctor_<TAlphabet>::Type   TIgnoreOrAssert;
    typedef EqualsChar<'>'>                                         TFastaBegin;

    clear(meta);
    clear(seq);

    skipUntil(iter, TFastaBegin());     // forward to the next '>'
    skipOne(iter);                      // assert and skip '>'

    readLine(meta, iter);               // read Fasta id
    readUntil(seq, iter, TFastaBegin(), TIgnoreOrAssert()); // read Fasta sequence
}

// ----------------------------------------------------------------------------
// Function readRecord(Fasta); Qualities inside qual
// ----------------------------------------------------------------------------

template <typename TIdString, typename TSeqString, typename TQualString, typename TFwdIterator>
inline void readRecord(TIdString & meta, TSeqString & seq, TQualString & qual, TFwdIterator & iter, Fasta)
{
    clear(qual);
    readRecord(meta, seq, iter, Fasta());
}

// ----------------------------------------------------------------------------
// Function readRecord(Fastq); Qualities inside seq
// ----------------------------------------------------------------------------

template <typename TIdString, typename TSeqString, typename TFwdIterator>
inline void readRecord(TIdString & meta, TSeqString & seq, TFwdIterator & iter, Fastq)
{
    typedef typename Value<TSeqString>::Type                                TAlphabet;
    typedef typename FastaIgnoreFunctor_<TAlphabet>::Type                   TIgnore;
    typedef typename FastaIgnoreOrAssertFunctor_<TAlphabet>::Type           TIgnoreOrAssert;
    typedef EqualsChar<'@'>                                                 TFastqBegin;
    typedef EqualsChar<'+'>                                                 TQualsBegin;

    clear(meta);
    clear(seq);

    skipUntil(iter, TFastqBegin());     // forward to the next '@'
    skipOne(iter);                      // skip '@'

    readLine(meta, iter);               // read Fastq id

    readUntil(seq, iter, TQualsBegin(), TIgnoreOrAssert());     // read Fastq sequence
    skipOne(iter, TQualsBegin());       // assert and skip '+'
    skipLine(iter);                     // skip optional 2nd Fastq id

    // as '@' could be either part of the base-call qualities or mark the
    // beginning of the next record, we count the number of consumed quality
    // values instead of (1) reading only 1 line or (2) until the next '@'
    CountDownFunctor<NotFunctor<TIgnore> > qualCountDown(length(seq));

    if (HasQualities<TAlphabet>::VALUE)
    {
        // TODO: Replace this temporary by Modifier
        CharString qual;
        TIgnore qualIgnore;
        readUntil(qual, iter, qualCountDown, qualIgnore);  // read Fastq qualities
        assignQualities(seq, qual);
    }
    else
    {
        skipUntil(iter, qualCountDown);                     // skip Fastq qualities
    }
    skipUntil(iter, TFastqBegin());                         // forward to the next '@'
}

// ----------------------------------------------------------------------------
// Function readRecord(Fastq); Qualities inside qual
// ----------------------------------------------------------------------------

template <typename TIdString, typename TSeqString, typename TQualString, typename TFwdIterator>
inline void readRecord(TIdString & meta, TSeqString & seq, TQualString & qual, TFwdIterator & iter, Fastq)
{
    typedef typename Value<TSeqString>::Type                                TSeqAlphabet;
    typedef typename Value<TQualString>::Type                               TQualAlphabet;
    typedef typename FastaIgnoreOrAssertFunctor_<TSeqAlphabet>::Type        TSeqIgnoreOrAssert;
    typedef typename FastaIgnoreFunctor_<TQualAlphabet>::Type               TQualIgnore;
    typedef typename FastaIgnoreOrAssertFunctor_<TQualAlphabet>::Type       TQualIgnoreOrAssert;
    typedef EqualsChar<'@'>                                                 TFastqBegin;
    typedef EqualsChar<'+'>                                                 TQualsBegin;

    clear(meta);
    clear(seq);
    clear(qual);

    skipUntil(iter, TFastqBegin());     // forward to the next '@'
    skipOne(iter);                      // skip '@'

    readLine(meta, iter);               // read Fastq id

    readUntil(seq, iter, TQualsBegin(), TSeqIgnoreOrAssert());  // read Fastq sequence
    skipOne(iter, TQualsBegin());       // assert and skip '+'
    skipLine(iter);                     // skip optional 2nd Fastq id

    // as '@' could be either part of the base-call qualities or mark the
    // beginning of the next record, we count the number of consumed quality
    // values instead of (1) reading only 1 line or (2) until the next '@'
    CountDownFunctor<NotFunctor<TQualIgnore> > qualCountDown(length(seq));
    TQualIgnoreOrAssert qualIgnore;

    readUntil(qual, iter, qualCountDown, qualIgnore);  // read Fastq qualities
    skipUntil(iter, TFastqBegin());     // forward to the next '@'
}

// ----------------------------------------------------------------------------
// Function writeRecord(Raw); Qualities inside seq
// ----------------------------------------------------------------------------

template <typename TFwdIterator, typename TIdString, typename TSeqString>
inline void
writeRecord(TFwdIterator & iter, TIdString const & /* meta */, TSeqString const & seq, Raw const &,
            SequenceOutputOptions const & = SequenceOutputOptions())
{
    write(iter, seq);
    writeValue(iter, '\n');
}

// ----------------------------------------------------------------------------
// Function writeRecord(Raw); Qualities inside qual
// ----------------------------------------------------------------------------

template <typename TFwdIterator, typename TIdString, typename TSeqString, typename TQualString>
inline void
writeRecord(TFwdIterator & iter, TIdString const & /* meta */, TSeqString const & seq, TQualString const & /* qual */, Raw const &,
            SequenceOutputOptions const & = SequenceOutputOptions())
{
    write(iter, seq);
    writeValue(iter, '\n');
}


// ----------------------------------------------------------------------------
// Function writeRecord(Fasta);
// ----------------------------------------------------------------------------

template <typename TTarget, typename TIdString, typename TSeqString>
inline void
writeRecord(TTarget & target,
            TIdString const & meta,
            TSeqString const & seq,
            Fasta const & /*tag*/,
            SequenceOutputOptions const & options = SequenceOutputOptions())
{
    writeValue(target, '>');
    write(target, meta);
    writeValue(target, '\n');

    writeWrappedString(target, seq, (options.lineLength < 0)? 70 : options.lineLength); // 70bp wrapping, by default
}

template <typename TTarget, typename TIdString, typename TSeqString, typename TQualString>
inline void
writeRecord(TTarget & target,
            TIdString const & meta,
            TSeqString const & seq,
            TQualString const & /*tag*/,
            Fasta const & fasta,
            SequenceOutputOptions const & options = SequenceOutputOptions())
{
    writeRecord(target, meta, seq, fasta, options);
}


// ----------------------------------------------------------------------------
// Function writeRecord(Fastq); Qualities inside qual
// ----------------------------------------------------------------------------

template <typename TTarget, typename TIdString, typename TSeqString, typename TQualString>
inline void
writeRecord(TTarget & target,
            TIdString const & meta,
            TSeqString const & seq,
            TQualString const & qual,
            Fastq const & /*tag*/,
            SequenceOutputOptions const & options = SequenceOutputOptions())
{
    writeValue(target, '@');
    write(target, meta);
    writeValue(target, '\n');

    int lineLength = (options.lineLength < 0)? 0 : options.lineLength;  // no wrapping, by default
    writeWrappedString(target, seq, lineLength);

    writeValue(target, '+');
    if (options.qualMeta)
        write(target, meta);
    writeValue(target, '\n');

    writeWrappedString(target, qual, lineLength);
}

// ----------------------------------------------------------------------------
// Function writeRecord(Fastq); Qualities inside seq
// ----------------------------------------------------------------------------

template <typename TTarget, typename TIdString, typename TSeqString>
inline void
writeRecord(TTarget & target,
            TIdString const & meta,
            TSeqString const & seq,
            Fastq const & tag,
            SequenceOutputOptions const & options = SequenceOutputOptions())
{
    typedef QualityExtractor<typename Value<TSeqString>::Type> TQualityExtractor;
    ModifiedString<TSeqString const, ModView<TQualityExtractor> > quals(seq);
    writeRecord(target, meta, seq, quals, tag, options);
}


// ----------------------------------------------------------------------------
// Function writeRecord(TagSelector); Qualities inside seq
// ----------------------------------------------------------------------------

template <typename TFwdIterator, typename TIdString, typename TSeqString>
inline void
writeRecord(TFwdIterator & /* iter */, TIdString const & /* meta */, TSeqString const & /* seq */,
            TagSelector<> const & /* format */, SequenceOutputOptions const & /* options */)
{}

template <typename TFwdIterator, typename TIdString, typename TSeqString, typename TTagList>
inline void
writeRecord(TFwdIterator & iter, TIdString const & meta, TSeqString const & seq,
            TagSelector<TTagList> const & format, SequenceOutputOptions const & options = SequenceOutputOptions())
{
    typedef typename TTagList::Type TFormat;

    if (isEqual(format, TFormat()))
        writeRecord(iter, meta, seq, TFormat(), options);
    else
        writeRecord(iter, meta, seq, static_cast<typename TagSelector<TTagList>::Base const &>(format), options);
}

// ----------------------------------------------------------------------------
// Function writeRecord(TagSelector); Qualities inside qual
// ----------------------------------------------------------------------------

template <typename TFwdIterator, typename TIdString, typename TSeqString, typename TQualString>
inline void
writeRecord(TFwdIterator & /* iter */, TIdString const & /* meta */, TSeqString const & /* seq */, TQualString const & /* qual */,
            TagSelector<> const & /* format */, SequenceOutputOptions const & /* options */)
{}

template <typename TFwdIterator, typename TIdString, typename TSeqString, typename TQualString, typename TTagList>
inline void
writeRecord(TFwdIterator & iter, TIdString const & meta, TSeqString const & seq, TQualString const & qual,
            TagSelector<TTagList> const & format, SequenceOutputOptions const & options = SequenceOutputOptions())
{
    typedef typename TTagList::Type TFormat;

    if (isEqual(format, TFormat()))
        writeRecord(iter, meta, seq, qual, TFormat(), options);
    else
        writeRecord(iter, meta, seq, qual, static_cast<typename TagSelector<TTagList>::Base const &>(format), options);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQ_IO_FASTA_FASTQ_H_
