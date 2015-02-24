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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
//         David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Class for reading/writing files in Fasta or Fastq format.
// ==========================================================================

#ifndef SEQAN_SEQ_IO_SEQUENCE_FILE_H_
#define SEQAN_SEQ_IO_SEQUENCE_FILE_H_

namespace seqan {

// ============================================================================
// Classes
// ============================================================================

// ============================================================================
// Typedefs
// ============================================================================

/*!
 * @class SeqFileIn
 * @signature typedef FormattedFile<Fastq, Input> SeqFileIn;
 * @extends FormattedFileIn
 * @headerfile <seqan/seq_io.h>
 * @brief Class for reading RAW, FASTA, FASTQ, EMBL and GENBANK files containing unaligned sequences.
 */

typedef FormattedFile<Fastq, Input>     SeqFileIn;

/*!
 * @class SeqFileOut
 * @signature typedef FormattedFile<Fastq, Output> SeqFileOut;
 * @extends FormattedFileOut
 * @headerfile <seqan/seq_io.h>
 * @brief Class for writing RAW, FASTA, FASTQ, EMBL and GENBANK files containing unaligned sequences.
 */

typedef FormattedFile<Fastq, Output>    SeqFileOut;

// --------------------------------------------------------------------------
// Tag AutoSeqFormat
// --------------------------------------------------------------------------
// if TagSelector is set to -1, the file format is auto-detected

/*!
 * @class AutoSeqFormat
 * @extends TagSelector
 * @headerfile <seqan/file.h>
 * @brief Auto-detects and stores a file format.
 *
 * @signature typedef TagList<Fastq, TagList<Fasta, TagList<Raw> > > SeqFormats;
 * @signature typedef TagSelector<SeqFormat> AutoSeqFormat;
 */

typedef
    TagList<Fastq,
    TagList<Fasta,
    TagList<Embl,
    TagList<GenBank,
    TagList<Raw
    > > > > >
    SeqInFormats;

typedef
    TagList<Fastq,
    TagList<Fasta,
    TagList<Raw
    > > >
    SeqOutFormats;

typedef TagSelector<SeqInFormats>   SeqInFormat;
typedef TagSelector<SeqOutFormats>  SeqOutFormat;

// deprecated
typedef SeqInFormat AutoSeqFormat;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction SeqFileBuffer_
// ----------------------------------------------------------------------------

template <typename TSeqStringSet, typename TSpec>
struct SeqFileBuffer_
{
    typedef typename Value<TSeqStringSet>::Type     TSeqString;
    typedef typename Value<TSeqString>::Type        TSeqAlphabet;

    typedef String<TSeqAlphabet>                    Type;
};

// ----------------------------------------------------------------------------
// Metafunction SeqFileContext_
// ----------------------------------------------------------------------------

template <typename TDirection>
struct SeqFileContext_;

template <>
struct SeqFileContext_<Input>
{
    Tuple<CharString, 3>    buffer;
    Dna5QString             hybrid;
};

template <>
struct SeqFileContext_<Output>
{
    SequenceOutputOptions   options;
};

// ----------------------------------------------------------------------------
// Metafunction FormattedFileContext
// ----------------------------------------------------------------------------

template <typename TSpec, typename TDirection, typename TStorageSpec>
struct FormattedFileContext<FormattedFile<Fastq, TDirection, TSpec>, TStorageSpec>
{
    typedef SeqFileContext_<TDirection> Type;
};

// ----------------------------------------------------------------------------
// Metafunction FileFormats
// ----------------------------------------------------------------------------

template <typename TSpec>
struct FileFormat<FormattedFile<Fastq, Input, TSpec> >
{
    typedef TagSelector<SeqInFormats> Type;
};

template <typename TSpec>
struct FileFormat<FormattedFile<Fastq, Output, TSpec> >
{
    typedef TagSelector<SeqOutFormats> Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecord(); Qualities in seq string
// ----------------------------------------------------------------------------

/*!
 * @fn SeqFileIn#readRecord
 * @brief Read one @link FormattedFileRecordConcept @endlink from a @link SeqFileIn @endlink object.
 *
 * @signature void readRecord(meta, seq, qual, fileIn);
 *
 * @param[out] meta         The @link StringConcept @endlink object where to read the meta information into.
 * @param[out] seq          The @link StringConcept @endlink object where to read the sequence information into.
 * @param[out] qual         The @link StringConcept @endlink object where to read the quality information into.
 * @param[in,out] fileIn    The @link SeqFileIn @endlink object to read from.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 */

template <typename TSpec, typename TIdString, typename TSeqString>
inline SEQAN_FUNC_ENABLE_IF(And<Is<InputStreamConcept<typename FormattedFile<Fastq, Input, TSpec>::TStream> >,
                                Not<HasQualities<typename Value<TSeqString>::Type> > >, void)
readRecord(TIdString & meta, TSeqString & seq, FormattedFile<Fastq, Input, TSpec> & file)
{
    readRecord(meta, seq, file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function readRecord(); No qualities in seq string
// ----------------------------------------------------------------------------

template <typename TSpec, typename TIdString, typename TSeqString>
inline SEQAN_FUNC_ENABLE_IF(And<Is<InputStreamConcept<typename FormattedFile<Fastq, Input, TSpec>::TStream> >,
                                HasQualities<typename Value<TSeqString>::Type> >, void)
readRecord(TIdString & meta, TSeqString & seq, FormattedFile<Fastq, Input, TSpec> & file)
{
    readRecord(meta, seq, context(file).buffer[2], file.iter, file.format);
    assignQualities(seq, context(file).buffer[2]);
}

// ----------------------------------------------------------------------------
// Function readRecord(); With separate qualities
// ----------------------------------------------------------------------------

template <typename TSpec, typename TIdString, typename TSeqString, typename TQualString>
inline SEQAN_FUNC_ENABLE_IF(Is<InputStreamConcept<typename FormattedFile<Fastq, Input, TSpec>::TStream> >, void)
readRecord(TIdString & meta, TSeqString & seq, TQualString & qual, FormattedFile<Fastq, Input, TSpec> & file)
{
    readRecord(meta, seq, qual, file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function readRecords()
// ----------------------------------------------------------------------------

/*!
 * @fn SeqFileIn#readRecords
 * @brief Read many @link FormattedFileRecordConcept @endlink from a @link SeqFileIn @endlink object.
 * @signature void readRecords(metas, seqs, quals, fileIn, numRecord);
 * @see SeqFileIn#readRecord
 */

template <typename TPtrA, typename TPtrB>
inline void
swapPtr(TPtrA &a, TPtrB &b)
{
    union {
        TPtrA a;
        TPtrB b;
    } tmp1, tmp2;
    tmp1.a = a;
    tmp2.b = b;
    b = tmp1.b;
    a = tmp2.a;
}

template <typename TIdStringSet, typename TSeqStringSet, typename TSpec, typename TSize>
inline void readRecords(TIdStringSet & meta,
                        TSeqStringSet & seq,
                        FormattedFile<Fastq, Input, TSpec> & file,
                        TSize maxRecords)
{
    typedef typename SeqFileBuffer_<TSeqStringSet, TSpec>::Type TSeqBuffer;

    TSeqBuffer seqBuffer;

    // reuse the memory of context(file).buffer for seqBuffer (which has a different type but same sizeof(Alphabet))
    swapPtr(seqBuffer.data_begin, context(file).buffer[1].data_begin);
    swapPtr(seqBuffer.data_end, context(file).buffer[1].data_end);
    seqBuffer.data_capacity = context(file).buffer[1].data_capacity;

    for (; !atEnd(file) && maxRecords > 0; --maxRecords)
    {
        readRecord(context(file).buffer[0], seqBuffer, file);
        appendValue(meta, context(file).buffer[0]);
        appendValue(seq, seqBuffer);
    }

    swapPtr(seqBuffer.data_begin, context(file).buffer[1].data_begin);
    swapPtr(seqBuffer.data_end, context(file).buffer[1].data_end);
    context(file).buffer[1].data_capacity = seqBuffer.data_capacity;
    seqBuffer.data_capacity = 0;
}

// ----------------------------------------------------------------------------
// Function readRecords(); Without max records
// ----------------------------------------------------------------------------

template <typename TIdStringSet, typename TSeqStringSet, typename TSpec>
inline void readRecords(TIdStringSet & meta,
                        TSeqStringSet & seq,
                        FormattedFile<Fastq, Input, TSpec> & file)
{
    readRecords(meta, seq, file, MaxValue<__uint64>::VALUE);
}

// ----------------------------------------------------------------------------
// Function readRecords(); With separate qualities
// ----------------------------------------------------------------------------

template <typename TIdStringSet, typename TSeqStringSet, typename TQualStringSet, typename TSpec, typename TSize>
inline void readRecords(TIdStringSet & meta,
                        TSeqStringSet & seq,
                        TQualStringSet & qual,
                        FormattedFile<Fastq, Input, TSpec> & file,
                        TSize maxRecords)
{
    typedef typename SeqFileBuffer_<TSeqStringSet, TSpec>::Type TSeqBuffer;

    TSeqBuffer seqBuffer;

    // reuse the memory of context(file).buffer for seqBuffer (which has a different type but same sizeof(Alphabet))
    std::swap(reinterpret_cast<char* &>(seqBuffer.data_begin), context(file).buffer[1].data_begin);
    std::swap(reinterpret_cast<char* &>(seqBuffer.data_end), context(file).buffer[1].data_end);
    std::swap(seqBuffer.data_capacity, context(file).buffer[1].data_capacity);

    for (; !atEnd(file) && maxRecords > 0; --maxRecords)
    {
        readRecord(context(file).buffer[0], seqBuffer, context(file).buffer[2], file);
        appendValue(meta, context(file).buffer[0]);
        appendValue(seq, seqBuffer);
        appendValue(qual, context(file).buffer[2]);
    }

    std::swap(reinterpret_cast<char* &>(seqBuffer.data_begin), context(file).buffer[1].data_begin);
    std::swap(reinterpret_cast<char* &>(seqBuffer.data_end), context(file).buffer[1].data_end);
    std::swap(seqBuffer.data_capacity, context(file).buffer[1].data_capacity);
}

// ----------------------------------------------------------------------------
// Function readRecords(); With separate qualities; Without max records
// ----------------------------------------------------------------------------

template <typename TIdStringSet, typename TSeqStringSet, typename TQualStringSet, typename TSpec>
inline void readRecords(TIdStringSet & meta,
                        TSeqStringSet & seq,
                        TQualStringSet & qual,
                        FormattedFile<Fastq, Input, TSpec> & file)
{
    readRecords(meta, seq, qual, file, MaxValue<__uint64>::VALUE);
}

// ----------------------------------------------------------------------------
// Function writeRecord()
// ----------------------------------------------------------------------------

/*!
 * @fn SeqFileOut#writeRecord
 * @brief Write one @link FormattedFileRecordConcept @endlink into a @link SeqFileOut @endlink object.
 *
 * @signature void writeRecord(fileOut, meta, seq, qual);
 *
 * @param[in,out] fileOut   The @link SeqFileOut @endlink object to write into.
 * @param[in] meta          The @link StringConcept @endlink object where to read the meta information from.
 * @param[in] seq           The @link StringConcept @endlink object where to read the sequence information from.
 * @param[in] qual          The @link StringConcept @endlink object where to read the quality information from.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 */

template <typename TSpec, typename TIdString, typename TSeqString>
inline SEQAN_FUNC_ENABLE_IF(Is<OutputStreamConcept<typename FormattedFile<Fastq, Output, TSpec>::TStream> >, void)
writeRecord(FormattedFile<Fastq, Output, TSpec> & file,
            TIdString const & meta,
            TSeqString const & seq)
{
    writeRecord(file.iter, meta, seq, file.format, context(file).options);
}

// ----------------------------------------------------------------------------
// Function writeRecord(); With separate qualities
// ----------------------------------------------------------------------------

template <typename TSpec, typename TIdString, typename TSeqString, typename TQualString>
inline SEQAN_FUNC_ENABLE_IF(Is<OutputStreamConcept<typename FormattedFile<Fastq, Output, TSpec>::TStream> >, void)
writeRecord(FormattedFile<Fastq, Output, TSpec> & file,
            TIdString const & meta,
            TSeqString const & seq,
            TQualString const & qual)
{
    writeRecord(file.iter, meta, seq, qual, file.format, context(file).options);
}

// ----------------------------------------------------------------------------
// Function writeRecords()
// ----------------------------------------------------------------------------

/*!
 * @fn SeqFileOut#writeRecords
 * @brief Write many @link FormattedFileRecordConcept @endlink into a @link SeqFileOut @endlink object.
 * @signature void writeRecords(fileOut, metas, seqs, quals);
 * @see SeqFileOut#writeRecord
 */

template <typename TSpec, typename TIdStringSet, typename TSeqStringSet>
inline void
writeRecords(FormattedFile<Fastq, Output, TSpec> & file,
             TIdStringSet const & meta,
             TSeqStringSet const & seq)
{
    for (typename Size<TIdStringSet>::Type i = 0; i != length(seq); ++i)
        writeRecord(file, meta[i], seq[i]);
}

// ----------------------------------------------------------------------------
// Function writeRecords(); With separate qualities
// ----------------------------------------------------------------------------

template <typename TSpec, typename TIdStringSet, typename TSeqStringSet, typename TQualStringSet>
inline void
writeRecords(FormattedFile<Fastq, Output, TSpec> & file,
             TIdStringSet const & meta,
             TSeqStringSet const & seq,
             TQualStringSet const & qual)
{
    for (typename Size<TIdStringSet>::Type i = 0; i != length(seq); ++i)
        writeRecord(file, meta[i], seq[i], qual[i]);
}

}  // namespace seqan

#endif // SEQAN_SEQ_IO_SEQUENCE_FILE_H_
