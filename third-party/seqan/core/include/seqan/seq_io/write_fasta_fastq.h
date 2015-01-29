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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// Code for writing FASTA and FASTQ files
// ==========================================================================

#ifndef SEQAN_SEQ_IO_WRITE_H_
#define SEQAN_SEQ_IO_WRITE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class SequenceOutputOptions
// ----------------------------------------------------------------------------

/**
.Class.SequenceOutputOptions
..cat:Input/Output
..summary:Configuration for writing sequence (FASTA/FASTQ) files.
..description:
This $struct$ is used for the configuration of writing out FASTA and FASTQ files.
..include:seqan/seq_io.h

.Memvar.SequenceOutputOptions#lineLength
..class:Class.SequenceOutputOptions
..type:nolink:$int$
..summary:Length of the lines when writing out.
..description:Set to $-1$ for default behaviour (no line break for FASTQ, line length of 70 for FASTA) and $0$ for disabling line breaks.

.Memvar.SequenceOutputOptions#qualMeta
..class:Class.SequenceOutputOptions
..type:nolink:$bool$
..summary:Whether or not to write the meta information into the $"+"$ line before the qualities (interpreted for FASTQ only). Default is $false$.
*/

// TODO(holtgrew): Would it be worth having two/three shortcuts for "short reads" and "genomic sequence" and faster or can the compiler optimize the creation away?

struct SequenceOutputOptions
{
public:
    int lineLength;
    bool qualMeta;

    explicit
    SequenceOutputOptions(int lineLength = -1, bool qualMeta = false) : lineLength(lineLength), qualMeta(qualMeta)
    {}    
};

// ============================================================================
// Metafunctions
// ============================================================================


// ============================================================================
// Functions
// ============================================================================


// ----------------------------------------------------------------------------
// Function writeRecord()
// ----------------------------------------------------------------------------

/**
.Function.FASTA/FASTQ I/O#writeRecord
..summary:Write one FASTA or FASTQ record.
..signature:int writeRecord(stream, id, seq, tag[, options])
..signature:int writeRecord(stream, id, seq, quals, tag[, options])
..param.stream:The stream to write to.
...type:Concept.StreamConcept
..param.id:ID/Meta information line to write out.
...type:Concept.SequenceConcept
..param.seq:Sequence to write out.
...type:Concept.SequenceConcept
..param.quals:ASCII quality characters to write out.
...type:Concept.SequenceConcept
..param.tag:The format selector.
...type:nolink:$Fasta$, $Fastq$
..param.options:if not supplied defaults are chosen.
...type:Class.SequenceOutputOptions
..include:seqan/seq_io.h
*/

template <typename TStream, typename TIdString, typename TSeqString>
inline int
writeRecord(TStream & stream,
            TIdString const & meta,
            TSeqString const & seq,
            Fasta const & /*tag*/,
            SequenceOutputOptions const & options)
{
    int res = streamWriteChar(stream, '>');
    if (res)
        return res;

    if (streamWriteBlock(stream, &*begin(meta, Standard()), length(meta)) != length(meta))
        return 1;

    res = streamWriteChar(stream, '\n');
    if (res)
        return res;

    int lineLength = (options.lineLength >= 0) ? options.lineLength : 70;

    if (lineLength > 0)
    {
        // write stream character by character
        typename Iterator<TSeqString const, Standard>::Type it = begin(seq);
        typename Iterator<TSeqString const, Standard>::Type it_end = end(seq);
        for (int l = 0; it < it_end; ++it)
        {
            res = streamWriteChar(stream, (char)*it);
            if (res)
                return res;
            if (++l == lineLength)
            {
                res = streamWriteChar(stream, '\n');
                l = 0;
                if (res)
                    return res;
            }
        }
        if (res)
            return res;
    } else
    {
        for (typename Iterator<TSeqString const, Rooted>::Type it = begin(seq, Rooted()); !atEnd(it); ++it)
            res = streamWriteChar(stream, (char)*it);
        if (res)
            return res;
    }

    res = streamWriteChar(stream, '\n');
    return res;
}

// FASTA
template <typename TStream, typename TIdString, typename TSeqString>
inline int
writeRecord(TStream & stream,
            TIdString const & meta,
            TSeqString const & seq,
            Fasta const & /*tag*/)
{
    return writeRecord(stream, meta, seq, Fasta(), SequenceOutputOptions());
}


template <typename TStream, typename TSequence, typename TQualString>
inline int _writeRecordFastq(TStream & stream, TSequence const & seq, TQualString const &, SequenceOutputOptions const & options, True const & /*HasQualities<Value<TSequence>::Type>::VALUE*/)
{
    int res = 0;

    if (options.lineLength > 0)
    {
        // write stream character by character
        typename Iterator<TSequence const>::Type it = begin(seq);
        typename Iterator<TSequence const>::Type it_end = end(seq);
        for (int l = 0; it < it_end; ++it)
        {
            res = streamWriteChar(stream, static_cast<char>('!' + getQualityValue(*it)));
            if (res)
                return res;
            if (++l == options.lineLength)
            {
                res = streamWriteChar(stream, '\n');
                l = 0;
                if (res)
                    return res;
            }
        }
        if (res)
            return res;
    }
    else
    {
        typename Iterator<TSequence const>::Type it = begin(seq);
        typename Iterator<TSequence const>::Type it_end = end(seq);
        for (; it < it_end; ++it)
        {
            res = streamWriteChar(stream, static_cast<char>('!' + getQualityValue(*it)));
            if (res)
                return res;
        }
    }

    return 0;
}


template <typename TStream, typename TSequence, typename TQualString>
inline int _writeRecordFastq(TStream & stream, TSequence const & seq, TQualString const & qual, SequenceOutputOptions const & options, False const & /*HasQualities<Value<TSequence>::Type>::VALUE*/)
{
    int res = 0;

    if (empty(qual)) // we don't actually have qualities
    {
        if (options.lineLength > 0)
        {
            for (int i = 0, l = 0; i < (int)length(seq); ++i)
            {
                res = streamWriteChar(stream, char(126));
                if (res)
                    return res;
                if (++l == options.lineLength)
                {
                    res = streamWriteChar(stream, '\n');
                    l = 0;
                    if (res)
                        return res;
                }
            }
        } else
        {
            for (int i = 0; i < (int)length(seq); ++i)
            {
                res = streamWriteChar(stream, char(33 + 40));
                if (res)
                    return res;
            }
        }
    } else
    {
        if (options.lineLength > 0)
        {
            // write stream character by character
            typename Iterator<TQualString const>::Type it = begin(qual);
            typename Iterator<TQualString const>::Type it_end = end(qual);
            for (int l = 0; it < it_end; ++it)
            {
                res = streamWriteChar(stream, (char)*it);
                if (res)
                    return res;
                if (++l == options.lineLength)
                {
                    res = streamWriteChar(stream, '\n');
                    l = 0;
                    if (res)
                        return res;
                }
            }
            if (res)
                return res;
        } else
        {
            if (streamWriteBlock(stream, &*begin(qual, Standard()), length(qual)) != length(qual))
                return 1;
        }
    }
    return 0;
}

// FASTQ
template <typename TStream,
          typename TIdString,
          typename TQualString,
          typename TSeqString>
inline int
writeRecord(TStream & stream,
            TIdString const & meta,
            TSeqString const & seq,
            TQualString const & qual,
            Fastq const & /*tag*/,
            SequenceOutputOptions const & options)
{
    int res = streamWriteChar(stream, '@');
    if (res)
        return res;

    if (streamWriteBlock(stream, &*begin(meta, Standard()), length(meta)) != length(meta))
        return 1;

    res = streamWriteChar(stream, '\n');
    if (res)
        return res;

    if (options.lineLength > 0)
    {
        // write stream character by character
        typename Iterator<TSeqString const>::Type it = begin(seq);
        typename Iterator<TSeqString const>::Type it_end = end(seq);
        for (int l = 0; it < it_end; ++it)
        {
            res = streamWriteChar(stream, (char)*it);
            if (res)
                return res;
            if (++l == options.lineLength)
            {
                res = streamWriteChar(stream, '\n');
                l = 0;
                if (res)
                    return res;
            }
        }
        if (res)
            return res;
    } else
    {
        for (typename Iterator<TSeqString const, Rooted>::Type it = begin(seq, Rooted()); !atEnd(it); ++it)
            res = streamWriteChar(stream, (char)*it);
        if (res)
            return res;
    }

    if (streamWriteBlock(stream, "\n+", 2) != 2)
        return 1;

    if (options.qualMeta)
    {
        if (streamWriteBlock(stream, &*begin(meta, Standard()), length(meta)) != length(meta))
            return res;
    }

    res = streamWriteChar(stream, '\n');
    if (res)
        return res;

    res = _writeRecordFastq(stream, seq, qual, options, typename HasQualities<typename Value<TSeqString>::Type>::Type());
    if (res)
        return res;
    res = streamWriteChar(stream, '\n');
    return res;
}

// FASTQ and we have no qualities
template <typename TStream, typename TIdString, typename TQualString, typename TSeqString>
inline int
writeRecord(TStream & stream,
            TIdString const & meta,
            TSeqString const & seq,
            TQualString const & qual,
            Fastq const & /*tag*/)
{
    return writeRecord(stream, meta, seq, qual, Fastq(), SequenceOutputOptions());
}

// FASTQ and we don't have the qualities
template <typename TStream, typename TIdString, typename TSeqString>
inline int
writeRecord(TStream & stream,
            TIdString const & meta,
            TSeqString const & seq,
            Fastq const & /*tag*/,
            SequenceOutputOptions const & options)
{
    return writeRecord(stream, meta, seq, CharString(), Fastq(), options);
}
// FASTQ and we don't have the qualities
template <typename TStream,
          typename TIdString,
          typename TSeqString>
inline int
writeRecord(TStream & stream,
            TIdString const & meta,
            TSeqString const & seq,
            Fastq const & /*tag*/)
{
    return writeRecord(stream, meta, seq, CharString(), Fastq(), SequenceOutputOptions());
}



// ----------------------------------------------------------------------------
// Function write2()
// ----------------------------------------------------------------------------

/**
.Function.FASTA/FASTQ I/O#write2
..summary:Write FASTA or FASTQ records.
..signature:int write2(stream, ids, seqs, tag[, options])
..signature:int write2(stream, ids, seqs, quals, tag[, options])
..param.stream:The stream to write to.
...type:Concept.StreamConcept
..param.ids:IDs/Metainformation strings to write out.
...type:Class.StringSet
..param.seqs:Sequences to write out.
...type:Class.StringSet
..param.quals:ASCII quality characters to write out.
...type:Class.StringSet
..param.tag:The format selector.
...type:nolink:$Fasta$, $Fastq$
..param.options:if not supplied defaults are chosen.
...type:Class.SequenceOutputOptions
..include:seqan/seq_io.h
*/

// FASTA
template <typename TStream, typename TIdString, typename TIdSpec, typename TSeqString, typename TSeqSpec>
int write2(TStream & stream,
         StringSet<TIdString, TIdSpec> const & sequenceIds,
         StringSet<TSeqString, TSeqSpec> const & sequences,
         Fasta const & /*tag*/,
         SequenceOutputOptions const & options)
{
    if (length(sequenceIds) != length(sequences))
        return -1;

    typedef StringSet<TIdString, TIdSpec> const TIdSet;
    typedef StringSet<TSeqString, TSeqSpec> const TSeqSet;

    typename Iterator<TIdSet>::Type  itMeta     = begin(sequenceIds);
    typename Iterator<TIdSet>::Type  itMeta_end = end(sequenceIds);
    typename Iterator<TSeqSet>::Type itSeq      = begin(sequences);
//    typename Iterator<TSeqSet>::Type itSeq_end  = end(sequences);

    for (; itMeta != itMeta_end; ++itMeta, ++itSeq)
    {
        int res = writeRecord(stream, *itMeta, *itSeq, Fasta(), options);
        if (res)
            return res;
    }
    return 0;
}

template <typename TStream,
          typename TIdString, typename TIdSpec,
          typename TSeqString, typename TSeqSpec>
int write2(TStream & stream,
         StringSet<TIdString, TIdSpec> const & sequenceIds,
         StringSet<TSeqString, TSeqSpec> const & sequences,
         Fasta const & /*tag*/)
{
    return write2(stream, sequenceIds, sequences, Fasta(), SequenceOutputOptions());

}

// FASTQ
template <typename TStream, typename TIdString, typename TIdSpec, typename TSeqString, typename TSeqSpec,
          typename TQualString, typename TQualSpec>
int write2(TStream & stream,
         StringSet<TIdString, TIdSpec> const & sequenceIds,
         StringSet<TSeqString, TSeqSpec> const & sequences,
         StringSet<TQualString, TQualSpec> const & qualities,
         Fastq const & /*tag*/,
         SequenceOutputOptions const & options)
{
    if (length(sequenceIds) != length(sequences) ||
        length(qualities) != length(sequences))
        return -1;

    typedef StringSet<TIdString, TIdSpec> const TIdSet;
    typedef StringSet<TSeqString, TSeqSpec> const TSeqSet;
    typedef StringSet<TQualString, TSeqSpec> const TQualSet;

    typename Iterator<TIdSet>::Type   itMeta      = begin(sequenceIds);
    typename Iterator<TIdSet>::Type   itMeta_end  = end(sequenceIds);
    typename Iterator<TSeqSet>::Type  itSeq       = begin(sequences);
//    typename Iterator<TSeqSet>::Type  itSeq_end   = end(sequences);
    typename Iterator<TQualSet>::Type itQual      = begin(qualities);
    // typename Iterator<TQualSet>::Type itQual_end  = end(qualities);

    for (; itMeta != itMeta_end; ++itMeta, ++itSeq, ++itQual)
    {
        int res = writeRecord(stream,*itMeta, *itSeq, *itQual, Fastq(), options);
        if (res)
            return res;
    }
    return 0;
}

template <typename TStream, typename TIdString, typename TIdSpec, typename TSeqString, typename TSeqSpec,
          typename TQualString, typename TQualSpec>
int write2(TStream & stream,
         StringSet<TIdString, TIdSpec> const & sequenceIds,
         StringSet<TSeqString, TSeqSpec> const & sequences,
         StringSet<TQualString, TQualSpec> const & qualities,
         Fastq const & /*tag*/)
{
    return write2(stream, sequenceIds, sequences, qualities, Fastq(), SequenceOutputOptions());
}

template <typename TStream, typename TIdString, typename TIdSpec, typename TSeqString, typename TSeqSpec>
int write2(TStream & stream,
         StringSet<TIdString, TIdSpec> const & sequenceIds,
         StringSet<TSeqString, TSeqSpec> const & sequences,
         Fastq const & /*tag*/,
         SequenceOutputOptions const & options)
{
    typedef StringSet<TIdString, TIdSpec> const TIdSet;
    typedef StringSet<TSeqString, TSeqSpec> const TSeqSet;

    typename Iterator<TIdSet>::Type   itMeta      = begin(sequenceIds);
    typename Iterator<TIdSet>::Type   itMeta_end  = end(sequenceIds);
    typename Iterator<TSeqSet>::Type  itSeq       = begin(sequences);
//    typename Iterator<TSeqSet>::Type  itSeq_end   = end(sequences);

    for (; itMeta != itMeta_end; ++itMeta, ++itSeq)
    {
        int res = writeRecord(stream, *itMeta, *itSeq, Fastq(), options);
        if (res)
            return res;
    }
    return 0;
}

template <typename TStream, typename TIdString, typename TIdSpec, typename TSeqString, typename TSeqSpec>
int write2(TStream & stream,
         StringSet<TIdString, TIdSpec> const & sequenceIds,
         StringSet<TSeqString, TSeqSpec> const & sequences,
         Fastq const & /*tag*/)
{
    return write2(stream, sequenceIds, sequences, Fastq(), SequenceOutputOptions());
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQ_IO_WRITE_H_
