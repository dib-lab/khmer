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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

// TODO(holtgrew): Add EMBL support.
// TODO(holtgrew): Add function readFasta(SEQ, FILENAME).

#ifndef CORE_INCLUDE_SEQAN_SEQ_IO_SEQUENCE_SEQ_IO_H_
#define CORE_INCLUDE_SEQAN_SEQ_IO_SEQUENCE_SEQ_IO_H_

#include <memory>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ---------------------------------------------------------------------------
// Class SequenceStream
// ---------------------------------------------------------------------------

/**
.Class.SequenceStream
..cat:Input/Output
..summary:High-level reading and writing of sequences.
..description:
Building upon the more low-level sequence I/O functionality of SeqAn, this class provides easier to use I/O facilities.
Especially, the underlying @Class.Stream|stream layer@ and using @Class.RecordReader@s is hidden from the user.
This is achieved by using dynamic polymorphism which comes at some performance cost.
..remarks:Operation Mode
..remarks.text:
When reading, there are two operation modes:
Normal reading and reading of "persistent" records.
When reading in "persistent" mode, @Class.SequenceStream@ will scan over each record twice:
Once for determining its size and once for actually reading the sequences.
After the first pass, we can allocate a buffer of the exact size we need.
This can save memory up to a factor of two, at the cost of scanning each record twice.
Note that this is only possible for reading uncompressed files.
..remarks:File Format and File Type
..remarks.text:
The file type determines whether a file is stored as raw text or whether it is compressed.
Examples for file types are text files or $gzip$ compressed files ($FILE.gz$).
The file format considers the contents of the raw/decompressed file.
Examples for file formats are FASTA, FASTQ, or EMBL.
..remarks.text:
When reading, the file type and format are guessed from the file itself.
You do not have to specify any but you can force the @Class.SequenceStream@ to use the ones you provide.
When writing, you should specify a file type and format when constructing the @Class.SequenceStream@ object.
Otherwise, it will default to writing out raw-text FASTA files.
..example.text:
Read the sequence file (FASTA or FASTQ) from $argv[1]$ record by record.
The identifiers and sequences of the stream are printed to stdout.
See the documentation of @Function.SequenceStream#readRecord@, @Function.SequenceStream#readBatch@, and @Function.SequenceStream#readAll@ for more examples, including record-wise reading, reading in batches, and reading all records in a file.
..example.file:demos/seq_io/sequence_stream_read.cpp
..example.text:
Open a SequenceStream for writing and write two sequences to it.
..example.file:demos/seq_io/sequence_stream_write.cpp
..wiki:Tutorial/SimpleSeqIO|Simple Sequence I/O
..include:seqan/seq_io.h

.Memfunc.SequenceStream#SequenceStream
..summary:Constructor
..description:See documentation of @Class.SequenceStream@ for more information.
..class:Class.SequenceStream
..signature:SequenceStream()
..signature:SequenceStream(fileName[, operationMode[, format[, fileType]]])
..param.fileName:Path to the file to open.
...type:nolink:$char const *$
..param.operationMode:Mode to open the file in. Optional.
...default:@Enum.SequenceStream\colon\colonOperationMode.value.READ@
...type:Enum.SequenceStream\colon\colonOperationMode
..param.format:Mode to open the file in. Optional.
...type:Enum.SequenceStream\colon\colonFileFormat
...default:@Enum.SequenceStream\colon\colonFileFormat.value.AUTO_FORMAT@
..param.fileType:Mode to open the file in. Optional.
...type:Enum.SequenceStream\colon\colonFileType
...default:@Enum.SequenceStream\colon\colonFileType.value.AUTO_TYPE@

.Enum.SequenceStream\colon\colonOperationMode
..cat:Input/Output
..summary:Select the operation mode of a @Class.SequenceStream@.
..value.READ:Open stream for reading.
..value.READ_PERSISTENT:Open stream for reading, mark as "persisent reading". See @Class.SequenceStream@ for more information on the difference between normal and persistent reading.
..value.WRITE:Open stream for writing.
..include:seqan/seq_io.h

.Enum.SequenceStream\colon\colonFileFormat
..cat:Input/Output
..summary:Select the file format to read/write.
..description:The file format is the format of the possibly compressed content.
..value.AUTO_FORMAT:Auto-detect format from file content on reading and from the file name on writing. If Auto-detection fails, FASTA is used.
..value.FASTA:Force reading/writing of FASTA.
..value.FASTQ:Force reading/writing of FASTQ.
..include:seqan/seq_io.h

.Enum.SequenceStream\colon\colonFileType
..cat:Input/Output
..summary:Select the file type to read/write.
..description:The file type is the type of the file itself, i.e. plain text or compressed.
..value.AUTO_TYPE:Auto-detect format from file content on reading and from the file name on writing. If Auto-detection fails, $PLAIN_TEXT$ is used.
..value.PLAIN_TEXT:Force reading/writing of plain text.
..value.GZ:Force reading/writing with gzip compression.
..value.BZ2:Force reading/writing with bzip compression.
..include:seqan/seq_io.h
*/

class SequenceStream
{
public:
    // -----------------------------------------------------------------------
    // Operation Mode
    // -----------------------------------------------------------------------

    // This enum is used to select the operation mode.
    enum OperationMode
    {
        READ,
        READ_PERSISTENT,
        WRITE
    };

    // This enum is used to select the file format.
    enum FileFormat
    {
        FASTA,
        FASTQ,
        AUTO_FORMAT
    };

    // This enum is used to select the file type, compression if any.
    enum FileType
    {
        AUTO_TYPE,
        PLAIN_TEXT
#if SEQAN_HAS_ZLIB
        ,
        GZ
#endif  // #if SEQAN_HAS_ZLIB
#if SEQAN_HAS_BZIP2
        ,
        BZ2
#endif  // #if SEQAN_HAS_BZIP2
    };

    // -----------------------------------------------------------------------
    // Member Variables
    // -----------------------------------------------------------------------

    // The configuration for the output format options.
    SequenceOutputOptions outputOptions;

    CharString filename;
    OperationMode operationMode;
    bool _atEnd;
    bool _isGood;

#if __cplusplus <= 199711L
    // C++98
    std::auto_ptr<SequenceStreamImpl_> _impl;
#else  // #if __cplusplus <= 199711L
    // C++11
    std::unique_ptr<SequenceStreamImpl_> _impl;
#endif  // #if __cplusplus <= 199711L

    SeqIOFileType_::Type _fileType;
    SeqIOFileFormat_::Type _fileFormat;

    // -----------------------------------------------------------------------
    // Constructor
    // -----------------------------------------------------------------------

    SequenceStream() : _atEnd(false), _isGood(false), _fileType(SeqIOFileType_::FILE_TYPE_TEXT),
        _fileFormat(SeqIOFileFormat_::FILE_FORMAT_FASTA)
    {}

    SequenceStream(char const * filename,
                   OperationMode operationMode = READ,
                   FileFormat format = AUTO_FORMAT,
                   FileType fileType = AUTO_TYPE) :
        filename(filename), operationMode(operationMode), _atEnd(false), _isGood(true), _fileType(SeqIOFileType_::FILE_TYPE_TEXT),
        _fileFormat(SeqIOFileFormat_::FILE_FORMAT_FASTA)
    {
        _init(operationMode, format, fileType);
    }

    void _init(OperationMode operationMode, FileFormat format, FileType fileType)
    {
        // Translate from FileFormat to SeqIOFileFormat_::Type.
        switch (format)
        {
            case FASTQ:
                _fileFormat = SeqIOFileFormat_::FILE_FORMAT_FASTQ;
                break;
            case FASTA:
                _fileFormat = SeqIOFileFormat_::FILE_FORMAT_FASTA;
                break;
            case AUTO_FORMAT:
                _fileFormat = SeqIOFileFormat_::FILE_FORMAT_AUTO;
                break;
        }
        // Translate from FileFormat to SeqIOFileFormat_::Type.
        switch (fileType)
        {
            case AUTO_TYPE:
                _fileType = SeqIOFileType_::FILE_TYPE_AUTO;
                break;
            case PLAIN_TEXT:
                _fileType = SeqIOFileType_::FILE_TYPE_TEXT;
                break;
#if SEQAN_HAS_ZLIB
            case GZ:
                _fileType = SeqIOFileType_::FILE_TYPE_GZ;
                break;
#endif  // #if SEQAN_HAS_ZLIB
#if SEQAN_HAS_BZIP2
            case BZ2:
                _fileType = SeqIOFileType_::FILE_TYPE_BZ2;
                break;
#endif  // #if SEQAN_HAS_BZIP2
        }

        bool isRead = (operationMode != WRITE);
        bool hintDoublePass = (operationMode == READ_PERSISTENT);
        _impl.reset(new SequenceStreamImpl_(filename, _fileFormat, _fileType, isRead, hintDoublePass));
        // Copy out, possibly detected/adjusted file type and format.
        _fileType = _impl->_fileType;
        _fileFormat = _impl->_fileFormat;
        _isGood = _impl->_isGood && (_fileType != SeqIOFileType_::FILE_TYPE_ERROR) &&
                (_fileFormat != SeqIOFileFormat_::FILE_FORMAT_ERROR);
        _atEnd = _impl->_atEnd;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

/**
.Function.SequenceStream#open
..summary:Open or re-open a file using a SequenceStream.
..class:Class.SequenceStream
..signature:void open(seqStream, fileName[, operationMode[, format[, fileType]]])
..param.seqStream:The SequenceStream object to open.
...type:Class.SequenceStream
..param.fileName:Path to the file to open.
...type:nolink:$char const *$
..param.operationMode:Mode to open the file in. Optional.
...default:@Enum.SequenceStream\colon\colonOperationMode.value.READ@
...type:Enum.SequenceStream\colon\colonOperationMode
..param.format:Mode to open the file in. Optional.
...type:Enum.SequenceStream\colon\colonFileFormat
...default:@Enum.SequenceStream\colon\colonFileFormat.value.AUTO_FORMAT@
..param.fileType:Mode to open the file in. Optional.
...type:Enum.SequenceStream\colon\colonFileType
...default:@Enum.SequenceStream\colon\colonFileType.value.AUTO_TYPE@
*/

inline void open(SequenceStream & seqIO,
                 char const * filename,
                 SequenceStream::OperationMode operationMode = SequenceStream::READ,
                 SequenceStream::FileFormat format = SequenceStream::AUTO_FORMAT,
                 SequenceStream::FileType fileType = SequenceStream::AUTO_TYPE)
{
    seqIO.filename = filename;
    seqIO.operationMode = operationMode;
    seqIO._atEnd = false;
    seqIO._isGood = true;
    seqIO._fileType = SeqIOFileType_::FILE_TYPE_TEXT;
    seqIO._fileFormat = SeqIOFileFormat_::FILE_FORMAT_FASTA;

    seqIO._init(operationMode, format, fileType);
}

// ----------------------------------------------------------------------------
// Function close()
// ----------------------------------------------------------------------------

/**
.Function.SequenceStream#close
..class:Class.SequenceStream
..summary:Close the @Class.SequenceStream@.
..signature:void close(seqIO)
..param.seqIO:The @Class.SequenceStream@ object to close from.
...type:Class.SequenceStream
..include:seqan/seq_io.h
*/

inline void close(SequenceStream & seqIO)
{
    seqIO._impl->close();
}

// ----------------------------------------------------------------------------
// Function flush()
// ----------------------------------------------------------------------------

/**
.Function.SequenceStream#flush
..class:Class.SequenceStream
..summary:Write all data from SequenceStream to disk.
..signature:void close(seqIO)
..param.seqIO:The @Class.SequenceStream@ object to flush.
...type:Class.SequenceStream
..include:seqan/seq_io.h
*/

inline void flush(SequenceStream & seqIO)
{
    seqIO._impl->flush();
}

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

/**
.Function.SequenceStream#atEnd
..class:Class.SequenceStream
..summary:Check whether a @Class.SequenceStream@ is at the end of the file.
..signature:bool isGood(seqIO)
..param.seqIO:The @Class.SequenceStream@ object to read from.
...type:Class.SequenceStream
..returns:A $bool$, $true$ indicating that the $seqIO$ is at the end of the file, $false$ indicating otherwise.
..see:Function.SequenceStream#isGood
..include:seqan/seq_io.h
*/

inline bool atEnd(SequenceStream const & seqIO)
{
    return seqIO._atEnd;
}

// TODO(holtgrew): We'd rather only have the const variant.
inline bool atEnd(SequenceStream & seqIO)
{
    return seqIO._atEnd;
}

// ----------------------------------------------------------------------------
// Function isGood()
// ----------------------------------------------------------------------------

/**
.Function.SequenceStream#isGood
..class:Class.SequenceStream
..summary:Check whether a @Class.SequenceStream@ object is ready for reading.
..signature:bool isGood(seqIO)
..param.seqIO:The @Class.SequenceStream@ object to read from.
...type:Class.SequenceStream
..returns:A $bool$, $true$ indicating that the $seqIO$ is ready for reading, $false$ that there was an error or it is at the end of the file.
..see:Function.SequenceStream#atEnd
..include:seqan/seq_io.h
*/

inline bool isGood(SequenceStream const & seqIO)
{
    return seqIO._isGood;
}

// ----------------------------------------------------------------------------
// Function readRecord()
// ----------------------------------------------------------------------------

/**
.Function.SequenceStream#readRecord
..class:Class.SequenceStream
..summary:Read the next sequence record from @Class.SequenceStream@.
..signature:int readRecord(id, seq, seqIO)
..signature:int readRecord(id, seq, quals, seqIO)
..param.id:The identifier of the sequence is written here.
...type:Shortcut.CharString
..param.seq:The sequence of the record is written here.
...type:Class.String
..param.quals:The qualities of the sequence is written here. Optional.
...type:Shortcut.CharString
...remarks:If the sequence has no qualities, @Function.clear@ is called on $quals$ to indicate this.
..param.seqIO:The @Class.SequenceStream@ object to read from.
...type:Class.SequenceStream
..returns:An integer, $0$ on success, $1$ on errors.
...type:nolink:$int$
..example:Read the first sequence of a FASTA file.
..example.code:
int main()
{
    seqan::SequenceStream seqIO("in.fasta", seqan::SequenceStream::READ_SINGLE);
    seqan::CharString id;
    seqan::Dna5String seq;

    if (atEnd(seqIO))
    {
        std::cerr << "ERROR: File does not contain any sequences!\n";
        return 1;
    }
    int res = readRecord(id, seq, seqIO);
    if (res != 0)
    {
        std::cerr << "ERROR: Could not read first record!\n";
        return 1;
    }

    return 0;
}
..include:seqan/seq_io.h
*/

template <typename TId, typename TSequence, typename TQualities>
int readRecord(TId & id, TSequence & seq, TQualities & qual, SequenceStream & seqIO)
{
    int res = 0;

    if (seqIO._fileFormat == SeqIOFileFormat_::FILE_FORMAT_FASTA)
        res = seqIO._impl->readRecord(id, seq, qual, Fasta());
    else if (seqIO._fileFormat == SeqIOFileFormat_::FILE_FORMAT_FASTQ)
        res = seqIO._impl->readRecord(id, seq, qual, Fastq());
    else
        res = 1;

    seqIO._isGood = seqIO._impl->_isGood;
    seqIO._atEnd = seqIO._impl->_atEnd;
    return res;
}

template <typename TId, typename TSequence>
int readRecord(TId & id, TSequence & seq, SequenceStream & seqIO)
{
    int res = 0;

    if (seqIO._fileFormat == SeqIOFileFormat_::FILE_FORMAT_FASTA)
        res = seqIO._impl->readRecord(id, seq, Fasta());
    else if (seqIO._fileFormat == SeqIOFileFormat_::FILE_FORMAT_FASTQ)
        res = seqIO._impl->readRecord(id, seq, Fastq());
    else
        res = 1;

    seqIO._isGood = seqIO._impl->_isGood;
    seqIO._atEnd = seqIO._impl->_atEnd;
    return res;
}

// ----------------------------------------------------------------------------
// Function readBatch()
// ----------------------------------------------------------------------------

/**
.Function.SequenceStream#readBatch
..class:Class.SequenceStream
..summary:Read a given number of sequence records from @Class.SequenceStream@.
..signature:int readBatch(ids, seqs, seqIO, num)
..signature:int readBatch(ids, seqs, quals, seqIO, num)
..param.ids:The identifiers of the sequence are written here.
...type:nolink:@Class.StringSet@ of @Shortcut.CharString@.
..param.seq:The sequence of the record is written here.
...type:Class.StringSet
..param.quals:The qualities of the sequence is written here. Optional.
...type:nolink:@Class.StringSet@ of @Shortcut.CharString@.
...remarks:If the sequences have no qualities, as in FASTA files, the @Class.StringSet@ will contain empty strings.
..param.seqIO:The @Class.SequenceStream@ object to read from.
..returns:An integer, $0$ on success, $1$ on errors.
...type:nolink:$int$
..example:Read the first sequences of a FASTA file, up to ten.
..example.code:
int main()
{
    seqan::SequenceStream seqIO("in.fasta", seqan::SequenceStream::READ_BATCH);
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::Dna5String> seqs;

    int res = readBatch(ids, seqs, seqIO, 10);
    if (res != 0)
    {
        std::cerr << "ERROR: Could not read records!\n";
        return 1;
    }

    return 0;
}
..include:seqan/seq_io.h
*/

template <typename TId, typename TIdSpec, typename TSequence, typename TSeqSpec, typename TQualities,
          typename TQualSpec>
int readBatch(StringSet<TId, TIdSpec> & ids,
              StringSet<TSequence, TSeqSpec> & seqs,
              StringSet<TQualities, TQualSpec> & quals,
              SequenceStream & seqIO,
              unsigned num)
{
    int res = 0;

    if (seqIO._fileFormat == SeqIOFileFormat_::FILE_FORMAT_FASTA)
        res = seqIO._impl->readBatch(ids, seqs, quals, num, Fasta());
    else if (seqIO._fileFormat == SeqIOFileFormat_::FILE_FORMAT_FASTQ)
        res = seqIO._impl->readBatch(ids, seqs, quals, num, Fastq());
    else
        res = 1;

    seqIO._isGood = seqIO._impl->_isGood;
    seqIO._atEnd = seqIO._impl->_atEnd;
    return res;
}

template <typename TId, typename TIdSpec, typename TSequence, typename TSeqSpec>
int readBatch(StringSet<TId, TIdSpec> & ids,
              StringSet<TSequence, TSeqSpec> & seqs,
              SequenceStream & seqIO,
              unsigned num)
{
    int res = 0;

    if (seqIO._fileFormat == SeqIOFileFormat_::FILE_FORMAT_FASTA)
        res = seqIO._impl->readBatch(ids, seqs, num, Fasta());
    else if (seqIO._fileFormat == SeqIOFileFormat_::FILE_FORMAT_FASTQ)
        res = seqIO._impl->readBatch(ids, seqs, num, Fastq());
    else
        res = 1;

    seqIO._isGood = seqIO._impl->_isGood;
    seqIO._atEnd = seqIO._impl->_atEnd;
    return res;
}

// ----------------------------------------------------------------------------
// Function readAll()
// ----------------------------------------------------------------------------

/**
.Function.SequenceStream#readAll
..class:Class.SequenceStream
..summary:Read all sequence records from a @Class.SequenceStream@ object.
..signature:int readAll(ids, seqs, seqIO)
..signature:int readAll(ids, seqs, quals, seqIO)
..param.ids:The identifiers of the sequence are written here.
...type:nolink:@Class.StringSet@ of @Shortcut.CharString@.
..param.seq:The sequence of the record is written here.
...type:Class.StringSet
..param.quals:The qualities of the sequence is written here. Optional.
...type:nolink:@Class.StringSet@ of @Shortcut.CharString@.
...remarks:If the sequences have no qualities, as in FASTA files, the @Class.StringSet@ will contain empty strings.
..param.seqIO:The @Class.SequenceStream@ object to read from.
...type:Class.SequenceStream
..returns:An integer, $0$ on success, $1$ on errors.
...type:nolink:$int$
..example:Read the sequences of a FASTA file.
..example.code:
int main()
{
    seqan::SequenceStream seqIO("in.fasta", seqan::SequenceStream::READ_ALL);
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::Dna5String> seqs;

    int res = readAll(ids, seqs, seqIO);
    if (res != 0)
    {
        std::cerr << "ERROR: Could not read records!\n";
        return 1;
    }

    return 0;
}
..include:seqan/seq_io.h
*/

template <typename TId, typename TIdSpec, typename TSequence, typename TSeqSpec, typename TQualities,
          typename TQualSpec>
int readAll(StringSet<TId, TIdSpec> & ids,
            StringSet<TSequence, TSeqSpec> & seqs,
            StringSet<TQualities, TQualSpec> & quals,
            SequenceStream & seqIO)
{
    int res = 0;

    if (seqIO._fileFormat == SeqIOFileFormat_::FILE_FORMAT_FASTA)
        res = seqIO._impl->readAll(ids, seqs, quals, Fasta());
    else if (seqIO._fileFormat == SeqIOFileFormat_::FILE_FORMAT_FASTQ)
        res = seqIO._impl->readAll(ids, seqs, quals, Fastq());
    else
        res = 1;

    seqIO._isGood = seqIO._impl->_isGood;
    seqIO._atEnd = seqIO._impl->_atEnd;
    return res;
}

template <typename TId, typename TIdSpec, typename TSequence, typename TSeqSpec>
int readAll(StringSet<TId, TIdSpec> & ids,
            StringSet<TSequence, TSeqSpec> & seqs,
            SequenceStream & seqIO)
{
    int res = 0;

    if (seqIO._fileFormat == SeqIOFileFormat_::FILE_FORMAT_FASTA)
        res = seqIO._impl->readAll(ids, seqs, Fasta());
    else if (seqIO._fileFormat == SeqIOFileFormat_::FILE_FORMAT_FASTQ)
        res = seqIO._impl->readAll(ids, seqs, Fastq());
    else
        res = 1;

    seqIO._isGood = seqIO._impl->_isGood;
    seqIO._atEnd = seqIO._impl->_atEnd;
    return res;
}

// ----------------------------------------------------------------------------
// Function writeRecord()
// ----------------------------------------------------------------------------

/**
.Function.SequenceStream#writeRecord
..class:Class.SequenceStream
..summary:Write one sequence record from to a @Class.SequenceStream@ object.
..description:
The record is appended to the file if you have written out any previously.
When writing out @Spec.Dna5@, qualities are automatically taken from the sequence characters.
..signature:int writeRecord(seqIO, id, seq, [options])
..signature:int writeRecord(seqIO, id, seq, quals, [options])
..param.seqIO:The @Class.SequenceStream@ object to write to.
...type:Class.SequenceStream
..param.id:The identifier to write.
...type:Shortcut.CharString
..param.seq:The sequence to write.
...type:Class.String
..param.quals:The qualities to write out.
...type:Shortcut.CharString
...remarks:If the sequence has no qualities, @Function.clear@ is called on $quals$ to indicate this.
..param.options:The configuration for writing FASTA and FASTQ files.
...type:Class.SequenceOutputOptions
..returns:An integer, $0$ on success, $1$ on errors.
...type:nolink:$int$
..example:Write out two sequences to a FASTQ file.
..example.code:
int main()
{
    seqan::SequenceStream seqIO("in.fasta", seqan::SequenceStream::WRITE);
    seqan::StringSet<seqan::CharString> ids;
    appendValue(ids, "seq1");
    appendValue(ids, "seq2");
    seqan::StringSet<seqan::Dna5String> seqs;
    appendValue(seqs, "CGAT");
    appendValue(seqs, "TTTT");

    for (unsigned i = 0; i < length(ids); ++i)
    {
        int res = writeRecord(seqIO, ids[0], seqs[0]);
        if (res != 0)
        {
            std::cerr << "ERROR: Could not write records!\n";
            return 1;
        }
    }

    return 0;
}
..include:seqan/seq_io.h
*/

template <typename TId, typename TSequence, typename TQualities>
int writeRecord(SequenceStream & seqIO,
                TId const & id,
                TSequence const & seq,
                TQualities const & qual,
                SequenceOutputOptions const & options)
{
    int res = 0;

    if (seqIO._fileFormat == SeqIOFileFormat_::FILE_FORMAT_FASTA)
        res = seqIO._impl->writeRecord(id, seq, Fasta(), options);
    else if (seqIO._fileFormat == SeqIOFileFormat_::FILE_FORMAT_FASTQ)
        res = seqIO._impl->writeRecord(id, seq, qual, Fastq(), options);
    else
        res = 1;

    seqIO._isGood = seqIO._impl->_isGood;
    return res;
}

template <typename TId, typename TSequence, typename TQualities>
int writeRecord(SequenceStream & seqIO,
                TId const & id,
                TSequence const & seq,
                TQualities const & qual)
{
    return writeRecord(seqIO, id, seq, qual, seqIO.outputOptions);
}

template <typename TId, typename TSequence>
int writeRecord(SequenceStream & seqIO,
                TId const & id,
                TSequence const & seq,
                SequenceOutputOptions const & options)
{
    int res = 0;

    if (seqIO._fileFormat == SeqIOFileFormat_::FILE_FORMAT_FASTA)
        res = seqIO._impl->writeRecord(id, seq, Fasta(), options);
    else if (seqIO._fileFormat == SeqIOFileFormat_::FILE_FORMAT_FASTQ)
        res = seqIO._impl->writeRecord(id, seq, Fastq(), options);
    else
        res = 1;

    seqIO._isGood = seqIO._impl->_isGood;
    return res;
}

template <typename TId, typename TSequence>
int writeRecord(SequenceStream & seqIO,
                TId const & id,
                TSequence const & seq)
{
    return writeRecord(seqIO, id, seq, seqIO.outputOptions);
}

// ----------------------------------------------------------------------------
// Function writeAll()
// ----------------------------------------------------------------------------

/**
.Function.SequenceStream#writeAll
..class:Class.SequenceStream
..summary:Write sequence records from to a @Class.SequenceStream@ object.
..description:
The records are appended to the file if you have written out any previously.
When writing out @Spec.Dna5@, qualities are automatically taken from the sequence characters.
..signature:int writeAll(seqIO, ids, seqs[, options])
..signature:int writeAll(seqIO, ids, seqs, quals, [options])
..param.seqIO:The @Class.SequenceStream@ object to write to.
...type:Class.SequenceStream
..param.ids:Identifiers to write out.
...type:nolink:@Class.StringSet@ of @Shortcut.CharString@.
..param.seq:Sequences to write out.
...type:Class.StringSet
..param.quals:Qualities to write out. Optional.
...type:nolink:@Class.StringSet@ of @Shortcut.CharString@.
...remarks:Qualities are ignored if the file format does not suppor them.  If none are given for FASTQ, score 40 is written out for all.
..param.options:The configuration for writing FASTA and FASTQ files.
...type:Class.SequenceOutputOptions
..returns:An integer, $0$ on success, $1$ on errors.
...type:nolink:$int$
..example:Write out all sequences.
..example.code:
int main()
{
    seqan::SequenceStream seqIO("in.fasta", seqan::SequenceStream::WRITE);
    seqan::StringSet<seqan::CharString> ids;
    appendValue(ids, "seq1");
    appendValue(ids, "seq2");
    seqan::StringSet<seqan::Dna5String> seqs;
    appendValue(seqs, "CGAT");
    appendValue(seqs, "TTTT");

    int res = writeAll(seqIO, ids, seqs);
    if (res != 0)
    {
        std::cerr << "ERROR: Could not write records!\n";
        return 1;
    }

    return 0;
}
..include:seqan/seq_io.h
*/

template <typename TId, typename TIdSpec, typename TSequence, typename TSeqSpec, typename TQualities,
          typename TQualSpec>
int writeAll(SequenceStream & seqIO,
             StringSet<TId, TIdSpec> const & ids,
             StringSet<TSequence, TSeqSpec> const & seqs,
             StringSet<TQualities, TQualSpec> const & quals,
             SequenceOutputOptions const & options)
{
    int res = 0;

    if (seqIO._fileFormat == SeqIOFileFormat_::FILE_FORMAT_FASTA)
        res = seqIO._impl->writeAll(ids, seqs, Fasta(), options);
    else if (seqIO._fileFormat == SeqIOFileFormat_::FILE_FORMAT_FASTQ)
        res = seqIO._impl->writeAll(ids, seqs, quals, Fastq(), options);
    else
        res = 1;

    seqIO._isGood = seqIO._impl->_isGood;
    return res;
}

template <typename TId, typename TIdSpec, typename TSequence, typename TSeqSpec, typename TQualities,
          typename TQualSpec>
int writeAll(SequenceStream & seqIO,
             StringSet<TId, TIdSpec> const & ids,
             StringSet<TSequence, TSeqSpec> const & seqs,
             StringSet<TQualities, TQualSpec> const & quals)
{
    return writeAll(seqIO, ids, seqs, quals, seqIO.outputOptions);
}

template <typename TId, typename TIdSpec, typename TSequence, typename TSeqSpec>
int writeAll(SequenceStream & seqIO,
             StringSet<TId, TIdSpec> const & ids,
             StringSet<TSequence, TSeqSpec> const & seqs,
             SequenceOutputOptions const & options)
{
    int res = 0;

    if (seqIO._fileFormat == SeqIOFileFormat_::FILE_FORMAT_FASTA)
        res = seqIO._impl->writeAll(ids, seqs, Fasta(), options);
    else if (seqIO._fileFormat == SeqIOFileFormat_::FILE_FORMAT_FASTQ)
        res = seqIO._impl->writeAll(ids, seqs, Fastq(), options);
    else
        res = 1;

    seqIO._isGood = seqIO._impl->_isGood;
    return res;
}

template <typename TId, typename TIdSpec, typename TSequence, typename TSeqSpec>
int writeAll(SequenceStream & seqIO,
             StringSet<TId, TIdSpec> const & ids,
             StringSet<TSequence, TSeqSpec> const & seqs)
{
    return writeAll(seqIO, ids, seqs, seqIO.outputOptions);
}

}  // namespace seqan

#endif  // CORE_INCLUDE_SEQAN_SEQ_IO_SEQUENCE_SEQ_IO_H_
