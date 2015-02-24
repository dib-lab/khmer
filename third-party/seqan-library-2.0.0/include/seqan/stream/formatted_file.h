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
// Base class for high level access to formatted files. It supports different
// file formats, compression, auto detection from file extensions and content
// as well as reading/writing from/to stdin/stdout.
// ==========================================================================

#ifndef SEQAN_STREAM_SMART_FILE_H_
#define SEQAN_STREAM_SMART_FILE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

/*!
 * @defgroup FileFormats File Formats
 * @brief Tags for identifying file formats.
 */

/*!
 * @mfn FormattedFile#FileFormat
 * @brief Metafunction for retrieving the file format type of a formatted file.
 *
 * @signature FileFormat<TFormattedFile>::Type;
 *
 * @tparam TFormattedFile The formatted file type to query for its file format type.
 * @return Type           The resulting @link FileFormats @endlink file formats type.
 */

template <typename TFormattedFile>
struct FileFormat;

template <typename TFormattedFile, typename TStorageSpec>
struct FormattedFileContext;

template <typename TObject, typename TStorageSpec>
struct StorageSwitch;

// ============================================================================
// Concepts
// ============================================================================

// ----------------------------------------------------------------------------
// Concept FormattedFileHeaderConcept
// ----------------------------------------------------------------------------

/*!
 * @concept FormattedFileHeaderConcept
 * @extends DefaultConstructibleConcept
 * @extends CopyConstructibleConcept
 * @extends AssignableConcept
 * @headerfile <seqan/stream.h>
 * @brief Concept for header of formatted files.
 * @signature concept FormattedFileHeaderConcept;
 */

// ----------------------------------------------------------------------------
// Concept FormattedFileRecordConcept
// ----------------------------------------------------------------------------

/*!
 * @concept FormattedFileRecordConcept
 * @extends DefaultConstructibleConcept
 * @extends CopyConstructibleConcept
 * @extends AssignableConcept
 * @headerfile <seqan/stream.h>
 * @brief Concept for record of formatted files.
 * @signature concept FormattedFileRecordConcept;
 */

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class FormattedFile
// ----------------------------------------------------------------------------

/*!
 * @class FormattedFile
 * @headerfile <seqan/stream.h>
 * @brief Base class for formatted file I/O.
 *
 * @signature template <typename TFileFormat, typename TDirection[, typename TSpec]>
 *            struct FormattedFile;
 *
 * @tparam TFileFormat  A type specifying the file format.
 * @tparam TDirection The direction of the file, one of @link DirectionTags#Input Input
 *                    @endlink or @link DirectionTags#Output @endlink.
 * @tparam TSpec      A tag for the specialization, defauls to <tt>void</tt>.
 *
 * FormattedFile provides the following basic I/O operations on formatted files:
 *
 * <ul>
 * <li>Open a file given its filename or attach to an existing stream like stdin/stdout.</li>
 * <li>Guess the file format from the file content or filename extension.</li>
 * <li>Set the file format manually.</li>
 * <li>Access compressed or uncompressed files transparently.</li>
 * </ul>
 *
 * FormattedFile encapsulates a VirtualStream and provides access to its @link StreamConcept#DirectionIterator @endlink.
 * Each instance of FormattedFile keeps a @link FormattedFile#FormattedFileContext @endlink while reading or writing the formatted file.
 */

/*!
 * @class FormattedFileIn
 * @headerfile <seqan/stream.h>
 * @brief Base class for reading formatted files.
 * @signature typedef FormattedFile<TFileFormat, Input, TSpec> FormattedFileIn;
 * @extends FormattedFile
 *
 * A formatted input file can write @link FormattedFileHeaderConcept Header @endlink
 * and @link FormattedFileRecordConcept Records @endlink.
 */

/*!
 * @fn FormattedFileIn#readHeader
 * @brief Read one @link FormattedFileHeaderConcept @endlink from a @link FormattedFileIn @endlink object.
 *
 * @signature void readHeader(fileIn, header);
 *
 * @param[in,out] fileIn    The @link FormattedFileIn @endlink object to read from.
 * @param[in]     header    The @link FormattedFileHeaderConcept @endlink to read.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 */

/*!
 * @fn FormattedFileIn#readRecord
 * @brief Read one @link FormattedFileRecordConcept @endlink from a @link FormattedFileIn @endlink object.
 *
 * @signature void readRecord(fileIn, record);
 *
 * @param[in,out] fileIn    The @link FormattedFileIn @endlink object to read from.
 * @param[in]     record    The @link FormattedFileRecordConcept @endlink to read.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 */

/*!
 * @class FormattedFileOut
 * @headerfile <seqan/stream.h>
 * @brief Base class for writing formatted files.
 * @signature typedef FormattedFile<TFileFormat, Output, TSpec> FormattedFileOut;
 * @extends FormattedFile
 *
 * A formatted output file can write @link FormattedFileHeaderConcept Header @endlink
 * and @link FormattedFileRecordConcept Records @endlink.
 */

/*!
 * @fn FormattedFileOut#writeHeader
 * @brief Write one @link FormattedFileHeaderConcept @endlink to a @link FormattedFileOut @endlink object.
 *
 * @signature void writeHeader(fileOut, header);
 *
 * @param[in,out] fileOut   The @link FormattedFileOut @endlink object to write into.
 * @param[in]     header    The @link FormattedFileHeaderConcept @endlink to write.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 */

/*!
 * @fn FormattedFileOut#writeRecord
 * @brief Write one @link FormattedFileRecordConcept @endlink to a @link FormattedFileOut @endlink object.
 *
 * @signature void writeRecord(fileOut, record);
 *
 * @param[in,out] fileOut   The @link FormattedFileOut @endlink object to write into.
 * @param[in]     record    The @link FormattedFileRecordConcept @endlink to write.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 */

template <typename TFileFormat, typename TDirection, typename TSpec = void>
struct FormattedFile
{
    typedef VirtualStream<char, TDirection>                                     TStream;
    typedef typename Iterator<TStream, TDirection>::Type                        TIter;
    typedef typename FileFormat<FormattedFile>::Type                            TFileFormats;
    typedef typename FormattedFileContext<FormattedFile, Owner<> >::Type        TOwnerContext;
    typedef typename FormattedFileContext<FormattedFile, Dependent<> >::Type    TDependentContext;

    TStream             stream;
    TIter               iter;
    TFileFormats        format;
    TOwnerContext       data_context;
    TDependentContext   context;

    /*!
     * @fn FormattedFile#FormattedFile
     * @brief Provides default construction.
     *
     * @signature FormattedFile::FormattedFile();
     * @signature FormattedFile::FormattedFile(fileName[, openMode]);
     * @signature FormattedFile::FormattedFile(stream);
     * @signature FormattedFile::FormattedFile(other);
     * @signature FormattedFile::FormattedFile(otherContext);
     * @signature FormattedFile::FormattedFile(otherContext, fileName[, openMode]);
     * @signature FormattedFile::FormattedFile(otherContext, stream);
     *
     * @param[in] fileName     Path to file to open, <tt>char const *</tt>.
     * @param[in] openMode     Optionally, the file open mode, default obtained from <tt>TDirection</tt>.
     * @param[in] stream       A <tt>std::basic_istream&lt;&gt;</tt> to read from or <tt>std::basic_ostream&lt;&gt;</tt>
     *                         to write to, depending on <tt>TDirection</tt>.
     * @param[in] other        A second FormattedFile, this FormattedFile's dependent context will depend on <tt>other</tt>'s
     *                         dependent context.
     * @param[in] otherContext The dependent context of another FormattedFile, this FormattedFile's dependent context will depend on <tt>otherContext</tt>.
     *
     * @throw IOError The variants that accept the <tt>fileName</tt> parameter throw an exception of this
     *                type in case opening the file fails.
     */
    FormattedFile() : context(data_context)
    {}

    // filename based c'tors
    explicit FormattedFile(const char * fileName, int openMode = DefaultOpenMode<FormattedFile>::VALUE) :
        context(data_context)
    {
        _open(*this, fileName, openMode, True());
    }

    // stream based c'tors
    template <typename TValue>
    explicit
    FormattedFile(std::basic_istream<TValue> &istream,
              SEQAN_CTOR_ENABLE_IF(And<IsSameType<TDirection, Input>, IsSameType<TValue, char> >)) :
        context(data_context)
    {
        _open(*this, istream, _mapFileFormatToCompressionFormat(format), True());
        ignoreUnusedVariableWarning(dummy);
    }

    template <typename TValue, typename TFormat>
    FormattedFile(std::basic_ostream<TValue> &ostream,
              TFormat const &format,
              SEQAN_CTOR_ENABLE_IF(And<IsSameType<TDirection, Output>, IsSameType<TValue, char> >)) :
        context(data_context)
    {
        bool success = open(*this, ostream, format);
        ignoreUnusedVariableWarning(dummy);
        ignoreUnusedVariableWarning(success);
        SEQAN_ASSERT(success);
    }

    // now everything given another context
    explicit
    FormattedFile(TDependentContext &otherCtx) :
        context(otherCtx)
    {}

    explicit
    FormattedFile(FormattedFile &other) :
        context(other.context)
    {}

    FormattedFile(TDependentContext &otherCtx, const char *fileName, int openMode = DefaultOpenMode<FormattedFile>::VALUE) :
        context(otherCtx)
    {
        _open(*this, fileName, openMode, True());
    }

    template <typename TValue>
    explicit
    FormattedFile(TDependentContext &otherCtx,
              std::basic_istream<TValue> &istream,
              SEQAN_CTOR_ENABLE_IF(And<IsSameType<TDirection, Input>, IsSameType<TValue, char> >)) :
        context(otherCtx)
    {
        _open(*this, istream, _mapFileFormatToCompressionFormat(format), True());
        ignoreUnusedVariableWarning(dummy);
    }

    template <typename TValue, typename TFormat>
    FormattedFile(TDependentContext &otherCtx,
              std::basic_ostream<TValue> &ostream,
              TFormat const &format,
              SEQAN_CTOR_ENABLE_IF(And<IsSameType<TDirection, Output>, IsSameType<TValue, char> >)) :
        context(otherCtx)
    {
        bool success = open(*this, ostream, format);
        ignoreUnusedVariableWarning(dummy);
        ignoreUnusedVariableWarning(success);
        SEQAN_ASSERT(success);
    }

    ~FormattedFile()
    {
        close(*this);
    }

    /*!
     * @fn FormattedFile::operatorTDependentContext FormattedFile::operator TDependentContext
     * @brief Allows conversion to a dependent context for the FormattedFile
     * @signature TDependentContext & FormattedFile::operator TDependentContext();
     *
     * @return TDependentContext The dependent context of this FormattedFile.
     */

    operator TDependentContext & ()
    {
        return context;
    }

    /*!
     * @fn FormattedFile::getFileExtensions
     * @brief Static function that returns a list of allowed file format extension.
     *
     * @signature TExtensionVector getFileExtensions()
     *
     * @return TExtensionVector A <tt>std::vector&lt;std::string&gt;</tt> with the allowed file extensions.
     */
    static std::vector<std::string>
    getFileExtensions()
    {
        std::vector<std::string> extensions;

        _getCompressionExtensions(extensions,
                                  TFileFormats(),
                                  typename FileFormat<TStream>::Type(),
//                                  true);
                                  IsSameType<TDirection, Output>::VALUE);
        return extensions;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction DirectionIterator
// ----------------------------------------------------------------------------

template <typename TFileFormat, typename TDirection, typename TSpec>
struct DirectionIterator<FormattedFile<TFileFormat, TDirection, TSpec>, TDirection>
{
    typedef typename FormattedFile<TFileFormat, TDirection, TSpec>::TIter Type;
};

// ----------------------------------------------------------------------------
// Metafunction FormattedFileContext
// ----------------------------------------------------------------------------
// FormattedFile holds a Context either as owner or dependent (see StorageSwitch below).
// Note that the FormattedFile class contains twice the context: both as owner and dependent.
// NOTE(esiragusa): A statically-typed Holder should abstract this pattern.

/*!
 * @mfn FormattedFile#FormattedFileContext
 * @brief Returns the context type for a FormattedFile.
 *
 * @signature FormattedFileContext<TFormattedFile, TStorageSpec>::Type
 *
 * @tparam TFormattedFile   The FormattedFile to query.
 * @tparam TStorageSpec     The storage specification, passed as specialization to any @link StringSet @endlink
 *                          contained in the context.
 * @tparam Type             The resulting FormattedFile context type.
 */

template <typename TFormattedFile, typename TStorageSpec>
struct FormattedFileContext
{
    typedef Nothing Type;
};

// ----------------------------------------------------------------------------
// Metafunction StorageSwitch
// ----------------------------------------------------------------------------
// This metafunction is used to switch the relationship between the FormattedFile
// and its Context among aggregation and composition.
// NOTE(esiragusa): there was a more generic metafunction Member<> to do this.

// Composition (owner).
template <typename TObject, typename TStorageSpec>
struct StorageSwitch
{
    typedef TObject Type;
};

// Aggregation (dependent).
template <typename TObject, typename TSpec>
struct StorageSwitch<TObject, Dependent<TSpec> >
{
    typedef TObject * Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TFileFormat, typename TDirection, typename TSpec>
struct Value<FormattedFile<TFileFormat, TDirection, TSpec> > :
    Value<typename FormattedFile<TFileFormat, TDirection, TSpec>::TStream> {};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <typename TFileFormat, typename TDirection, typename TSpec>
struct Position<FormattedFile<TFileFormat, TDirection, TSpec> > :
    Position<typename FormattedFile<TFileFormat, TDirection, TSpec>::TStream> {};

// ----------------------------------------------------------------------------
// Metafunction DefaultOpenMode
// ----------------------------------------------------------------------------

template <typename TFileFormat, typename TDirection, typename TSpec, typename TDummy>
struct DefaultOpenMode<FormattedFile<TFileFormat, TDirection, TSpec>, TDummy> :
    DefaultOpenMode<typename FormattedFile<TFileFormat, TDirection, TSpec>::TStream, TDummy> {};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Helper Function _throwIf().
// ----------------------------------------------------------------------------
// Helper functions that reduce number of "uncaught exception" false positives in static analysis tools.

template <typename TException> void _throwIf(TException const & e, True const & /*tag*/) { SEQAN_THROW(e); }
template <typename TException> void _throwIf(TException const & /*e*/, False const & /*tag*/) { /*nop*/ }

// ----------------------------------------------------------------------------
// Function directionIterator()
// ----------------------------------------------------------------------------

template <typename TFileFormat, typename TDirection, typename TSpec>
inline typename FormattedFile<TFileFormat, TDirection, TSpec>::TIter &
directionIterator(FormattedFile<TFileFormat, TDirection, TSpec> & file, TDirection const &)
{
    return file.iter;
}

// ----------------------------------------------------------------------------
// Function format()
// ----------------------------------------------------------------------------

/*!
 * @fn FormattedFile#format
 * @brief Return the format of a FormattedFile.
 *
 * @signature TFormat format(file);
 *
 * @param[in] file The FormattedFile to check.
 * @return TFormat The type as returned from @link FormattedFile#FileFormat @endlink.
 */

template <typename TFileFormat, typename TDirection, typename TSpec>
inline typename FileFormat<FormattedFile<TFileFormat, TDirection, TSpec> >::Type &
format(FormattedFile<TFileFormat, TDirection, TSpec> & file)
{
    return file.format;
}

template <typename TFileFormat, typename TDirection, typename TSpec>
inline typename FileFormat<FormattedFile<TFileFormat, TDirection, TSpec> >::Type const &
format(FormattedFile<TFileFormat, TDirection, TSpec> const & file)
{
    return file.format;
}

// ----------------------------------------------------------------------------
// Function setFormat()
// ----------------------------------------------------------------------------

/*!
 * @fn FormattedFile#setFormat
 * @brief Set the format of a FormattedFile.
 *
 * @signature void setFormat(file, format);
 *
 * @param[in,out] file The FormattedFile to change.
 * @param[in]     format The @link FormattedFile#FileFormat @endlink to set.
 */

template <typename TFileFormat, typename TDirection, typename TSpec, typename TFormat>
inline void
setFormat(FormattedFile<TFileFormat, TDirection, TSpec> & file, TFormat format)
{
    assign(file.format, format);
}

// ----------------------------------------------------------------------------
// Function guessFormat()
// ----------------------------------------------------------------------------

/*!
 * @fn FormattedFile#guessFormat
 * @brief Guess the format of an open FormattedFile.
 *
 * @signature bool guessFormat(file);
 *
 * @param[in,out] file The open FormattedFile for which the format is to be guessed.
 */

template <typename TFileFormat, typename TSpec>
inline bool guessFormat(FormattedFile<TFileFormat, Input, TSpec> & file)
{
    return guessFormatFromStream(file.stream, file.format);
}

template <typename TFileFormat, typename TSpec>
inline bool guessFormat(FormattedFile<TFileFormat, Output, TSpec> &)
{
    return true;
}

// ----------------------------------------------------------------------------
// _mapFileFormatToCompressionFormat
// ----------------------------------------------------------------------------

template <typename TFormat>
inline Nothing
_mapFileFormatToCompressionFormat(TFormat)
{
    return Nothing();
}

template <typename TCompressedFileTypes>
inline void
_mapFileFormatToCompressionFormat(TagSelector<TCompressedFileTypes> &,
                                  TagSelector<void> const &)
{}

template <typename TCompressedFileTypes, typename TTagList>
inline void
_mapFileFormatToCompressionFormat(TagSelector<TCompressedFileTypes> & result,
                                  TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormat;
    typedef typename TagSelector<TTagList>::Base TBase;

    if (isEqual(format, TFormat()))
        assign(result, _mapFileFormatToCompressionFormat(TFormat()));
    else
        _mapFileFormatToCompressionFormat(result, static_cast<TBase const &>(format));
}

template <typename TFileFormatList>
inline TagSelector<CompressedFileTypes>
_mapFileFormatToCompressionFormat(TagSelector<TFileFormatList> format)
{
    TagSelector<CompressedFileTypes> compressionType;
    _mapFileFormatToCompressionFormat(compressionType, format);
    return compressionType;
}


template <typename TFormattedFile, typename TFormat>
inline void
_checkThatStreamOutputFormatIsSet(TFormattedFile const &, TFormat const &)
{}

template <typename TFileFormat, typename TSpec, typename TFileFormatList>
inline void
_checkThatStreamOutputFormatIsSet(FormattedFile<TFileFormat, Output, TSpec> const &, TagSelector<TFileFormatList> const &format)
{
    if (value(format) < 0)
        SEQAN_FAIL("FormattedFile: File format not specified, use setFormat().");
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

/*!
 * @fn FormattedFile#open
 * @brief Open a FormattedFile.
 *
 * @signature bool open(file, fileName);
 *
 * @param[in,out] file The FormattedFile to open.
 * @param[in]     fileName The name of the file open.
 * @return bool <tt>true</tt> in the case of success, <tt>false</tt> otherwise.
 */

template <typename TFileFormat, typename TDirection, typename TSpec,
          typename TStream, typename TCompressionFormat, typename TThrowExceptions>
inline bool _open(FormattedFile<TFileFormat, TDirection, TSpec> & file,
                  TStream &stream,
                  TCompressionFormat const &compressionFormat,
                  TThrowExceptions = True())
{
    if (!open(file.stream, stream, compressionFormat))
    {
        _throwIf(UnknownFileFormat(), TThrowExceptions());
        return false;
    }

    if (!guessFormat(file))
    {
        _throwIf(UnknownFileFormat(), TThrowExceptions());
        return false;
    }

    file.iter = directionIterator(file.stream, TDirection());

    return true;
}

template <typename TFileFormat, typename TDirection, typename TSpec, typename TStream>
inline bool open(FormattedFile<TFileFormat, TDirection, TSpec> & file,
                 TStream &stream)
{
    return _open(file, stream, _mapFileFormatToCompressionFormat(file.format), False());
}

template <typename TFileFormat, typename TSpec, typename TStream, typename TFormat_>
inline bool open(FormattedFile<TFileFormat, Output, TSpec> & file,
                 TStream &stream,
                 Tag<TFormat_> format)
{
    setFormat(file, format);
    return _open(file, stream, _mapFileFormatToCompressionFormat(file.format), False());
}

template <typename TFileFormat, typename TSpec, typename TStream, typename TFormats>
inline bool open(FormattedFile<TFileFormat, Output, TSpec> & file,
                 TStream &stream,
                 TagSelector<TFormats> format)
{
    setFormat(file, format);
    return _open(file, stream, _mapFileFormatToCompressionFormat(file.format), False());
}

// ----------------------------------------------------------------------------
// Function open(fileName)
// ----------------------------------------------------------------------------

template <typename TFileFormat, typename TDirection, typename TSpec, typename TThrowExceptions>
inline bool _open(FormattedFile<TFileFormat, TDirection, TSpec> & file,
                  const char *fileName,
                  int openMode = DefaultOpenMode<FormattedFile<TFileFormat, TDirection, TSpec> >::VALUE,
                  TThrowExceptions = True())
{
    if (!open(file.stream, fileName, openMode))
    {
        _throwIf(FileOpenError(fileName), TThrowExceptions());
        return false;
    }

    if (IsSameType<TDirection, Input>::VALUE && _isPipe(fileName))
    {
        if (!guessFormat(file))
        {
            // read from a pipe (without file extension)
            _throwIf(UnknownFileFormat(), TThrowExceptions());
            return false;
        }
    }
    else
    {
        Prefix<const char *>::Type basename = _getUncompressedBasename(fileName, format(file.stream));
        if (!guessFormatFromFilename(basename, file.format))    // read/write from/to a file (with extension)
        {
            close(file.stream);
            _throwIf(UnknownExtensionError(fileName), TThrowExceptions());
            return false;
        }
    }

    file.iter = directionIterator(file.stream, TDirection());
    return true;
}

template <typename TFileFormat, typename TDirection, typename TSpec>
inline bool open(FormattedFile<TFileFormat, TDirection, TSpec> & file,
                 const char *fileName,
                 int openMode = DefaultOpenMode<FormattedFile<TFileFormat, TDirection, TSpec> >::VALUE)
{
    return _open(file, fileName, openMode, False());
}

// ----------------------------------------------------------------------------
// Function close()
// ----------------------------------------------------------------------------

/*!
 * @fn FormattedFile#close
 * @brief Close a FormattedFile.
 *
 * @signature bool close(file);
 *
 * @param[in,out] file The FormattedFile to close.
 * @return bool <tt>true</tt> in the case of success, <tt>false</tt> otherwise.
 */

template <typename TFileFormat, typename TDirection, typename TSpec>
inline bool close(FormattedFile<TFileFormat, TDirection, TSpec> & file)
{
    setFormat(file, typename FileFormat<FormattedFile<TFileFormat, TDirection, TSpec> >::Type());
    file.iter = typename DirectionIterator<FormattedFile<TFileFormat, TDirection, TSpec>, TDirection>::Type();
    return close(file.stream);
}

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

/*!
 * @fn FormattedFile#atEnd
 * @brief Determines whether a FormattedFile is at the end.
 *
 * @signature bool atEnd(file);
 *
 * @param[in,out] file The FormattedFile to check.
 * @return bool <tt>true</tt> in the case of success, <tt>false</tt> otherwise.
 */

template <typename TFileFormat, typename TDirection, typename TSpec>
inline SEQAN_FUNC_ENABLE_IF(Is<InputStreamConcept<typename FormattedFile<TFileFormat, TDirection, TSpec>::TStream> >, bool)
atEnd(FormattedFile<TFileFormat, TDirection, TSpec> const & file)
{
    return atEnd(file.iter);
}

// ----------------------------------------------------------------------------
// Function position()
// ----------------------------------------------------------------------------

template <typename TFileFormat, typename TSpec>
inline typename Position<FormattedFile<TFileFormat, Output, TSpec> >::Type
position(FormattedFile<TFileFormat, Output, TSpec> & file)
{
    return file.stream.tellp();
}

template <typename TFileFormat, typename TSpec>
inline typename Position<FormattedFile<TFileFormat, Input, TSpec> >::Type
position(FormattedFile<TFileFormat, Input, TSpec> & file)
{
    return file.stream.tellg();
}

// ----------------------------------------------------------------------------
// Function setPosition()
// ----------------------------------------------------------------------------

template <typename TFileFormat, typename TSpec, typename TPosition>
inline bool
setPosition(FormattedFile<TFileFormat, Output, TSpec> & file, TPosition pos)
{
    return (TPosition)file.stream.rdbuf()->pubseekpos(pos, std::ios_base::out) == pos;
}

template <typename TFileFormat, typename TSpec, typename TPosition>
inline bool
setPosition(FormattedFile<TFileFormat, Input, TSpec> & file, TPosition pos)
{
    return (TPosition)file.stream.rdbuf()->pubseekpos(pos, std::ios_base::in) == pos;
}

// ----------------------------------------------------------------------------
// Function context()
// ----------------------------------------------------------------------------

/*!
 * @fn FormattedFile#context
 * @brief Return the FormattedFile's dependent context object.
 *
 * @signature TDependentContext & context(file);
 *
 * @param[in,out] file The FormattedFile to query for its context.
 *
 * @return TDependentContext The dependent context, type as returned from @link FormattedFile#FormattedFileContext @endlink.
 */

template <typename TFileFormat, typename TDirection, typename TSpec>
inline typename FormattedFileContext<FormattedFile<TFileFormat, TDirection, TSpec>, Dependent<> >::Type &
context(FormattedFile<TFileFormat, TDirection, TSpec> & file)
{
    return file.context;
}

template <typename TFileFormat, typename TDirection, typename TSpec>
inline typename FormattedFileContext<FormattedFile<TFileFormat, TDirection, TSpec>, Dependent<> >::Type const &
context(FormattedFile<TFileFormat, TDirection, TSpec> const & file)
{
    return file.context;
}

// ----------------------------------------------------------------------------
// Function _getCompressionFormatExtensions()
// ----------------------------------------------------------------------------

template <typename TStringSet, typename TFormat_, typename TCompressionFormats>
inline void
_getCompressionExtensions(
    TStringSet &stringSet,
    Tag<TFormat_> const & /*formatTag*/,
    TCompressionFormats const & compress,
    bool primaryExtensionOnly,
    Nothing)
{
    typedef Tag<TFormat_> TFormat;

    std::vector<std::string> compressionExtensions;
    _getFileExtensions(compressionExtensions, compress, primaryExtensionOnly);

    unsigned len = (primaryExtensionOnly)? 1 : sizeof(FileExtensions<TFormat>::VALUE) / sizeof(char*);
    for (unsigned i = 0; i < len; ++i)
        for (unsigned j = 0; j < compressionExtensions.size(); ++j)
        {
            size_t jj = (j == 0)? compressionExtensions.size() - 1 : j - 1;    // swap first and last compression extension
            appendValue(stringSet, (std::string)FileExtensions<TFormat>::VALUE[i] + compressionExtensions[jj]);
        }
}

template <typename TStringSet, typename TFormat_, typename TCompressionFormats, typename TCompression_>
inline void
_getCompressionExtensions(
    TStringSet &stringSet,
    Tag<TFormat_> const & formatTag,
    TCompressionFormats const &,
    bool primaryExtensionOnly,
    Tag<TCompression_>)
{
    _getFileExtensions(stringSet, formatTag, primaryExtensionOnly);
}

template <typename TStringSet, typename TFormat_, typename TCompressionFormats>
inline void
_getCompressionExtensions(
    TStringSet &stringSet,
    Tag<TFormat_> const & formatTag,
    TCompressionFormats const & compress,
    bool primaryExtensionOnly)
{
    _getCompressionExtensions(stringSet, formatTag, compress, primaryExtensionOnly, _mapFileFormatToCompressionFormat(formatTag));
}

template <typename TStringSet, typename TTag, typename TCompressionFormats>
inline void
_getCompressionExtensions(
    TStringSet &stringSet,
    TagList<TTag, void> const & /*formatTag*/,
    TCompressionFormats const & compress,
    bool primaryExtensionOnly = false)
{
    _getCompressionExtensions(stringSet, TTag(), compress, primaryExtensionOnly);
}

template <typename TStringSet, typename TTag, typename TSubList, typename TCompressionFormats>
inline void
_getCompressionExtensions(
    TStringSet &stringSet,
    TagList<TTag, TSubList> const & /*formatTag*/,
    TCompressionFormats const & compress,
    bool primaryExtensionOnly = false)
{
    _getCompressionExtensions(stringSet, TTag(), compress, primaryExtensionOnly);
    _getCompressionExtensions(stringSet, TSubList(), compress, primaryExtensionOnly);
}

template <typename TStringSet, typename TTagList, typename TCompressionFormats>
inline void
_getCompressionExtensions(
    TStringSet &stringSet,
    TagSelector<TTagList> const & /*formatTag*/,
    TCompressionFormats const & compress,
    bool primaryExtensionOnly = false)
{
    _getCompressionExtensions(stringSet, TTagList(), compress, primaryExtensionOnly);
}

// ----------------------------------------------------------------------------
// Function getFileExtensions()
// ----------------------------------------------------------------------------

/*!
 * @fn FormattedFile#getFileExtensions
 * @brief Static function that returns a list of allowed file format extension.
 *
 * @signature TExtensionVector getFileExtensions(file)
 *
 * @param[in] file The FormattedFile to query.
 * @return TExtensionVector A <tt>std::vector&lt;std::string&gt;</tt> with the allowed file extensions.
 *
 * This is a shortcut to @link FormattedFile#getFileExtensions @endlink.
 */

template <typename TFileFormat, typename TDirection, typename TSpec>
static std::vector<std::string>
getFileExtensions(FormattedFile<TFileFormat, TDirection, TSpec> const & file)
{
    return file.getFileExtensions();
}

}  // namespace seqan

#endif // SEQAN_STREAM_SMART_FILE_H_
