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
// Defines basic tags and flags for file access.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_FILE_INTERFACE_H_
#define SEQAN_INCLUDE_SEQAN_FILE_INTERFACE_H_

namespace seqan {

/*!
 * @class SyncFile
 * @extends File
 * @headerfile <seqan/file.h>
 * @brief File structure supporting synchronous input/output access.
 *
 * @signature template <[typename TSpec]>
 *            class File<Sync<TSpec> >;
 *
 * @tparam TSpec Further specializing type.  Default: <tt>void</tt>.
 *
 * This class supports pseudo-asynchronous access methods, i.e. the method to initiate an I/O request blocks until the
 * request completion.
 */

template <typename TSpec = void>
struct Sync;

/*!
 * @class AsyncFile
 * @extends File
 * @headerfile <seqan/file.h>
 * @brief File structure supporting asynchronous input/output access.
 *
 * @signature template <[typename TSpec]>
 *            class File<Async<TSpec> >;
 *
 * @tparam TSpec Further specializing type.  Default: <tt>void</tt>.
 */

template <typename TSpec = void>
struct Async;

/*!
 * @class File
 * @headerfile <seqan/file.h>
 * @brief Represents a file.
 *
 * @signature template <[typename TSpec]>
 *            class File<TSpec>;
 *
 * @tparam TSpec Specializing type.  Default: <tt>Async&lt;&gt;</tt>.
 */

template <typename TSpec = Async<> >
class File;

/*!
 * @enum FileOpenMode
 * @headerfile <seqan/file.h>
 * @brief Flags to select th eopen mode of a @link File @endlink or external string.
 *
 * These flags can be combined via the <tt>|</tt> operator (bitwise OR).  The defualt open mode is <tt>OPEN_RDWR |
 * OPEN_CREATE | OPEN_APPEND</tt>.
 *
 * If you omit the <tt>OPEN_APPEND</tt> flag in write mode, the file will be truncated to size 0 when opened.
 *
 * @section Examples
 *
 * Code example to test for read-only mode:
 *
 * @code{.cpp}
 * if (openMode & OPEN_MASK == OPEN_READ)
 *     // do something if opened in read-only mode
 * @endcode
 *
 * @val FileOpenMode OPEN_RDONLY
 * @brief Open in read-only mode.
 *
 * @val FileOpenMode OPEN_WRONLY
 * @brief Open in write-only mode.
 *
 * @val FileOpenMode OPEN_RDWR
 * @brief Open for reading and writing.
 *
 * @val FileOpenMode OPEN_CREATE
 * @brief Create the file if it does not yet exists.
 *
 * @val FileOpenMode OPEN_APPEND
 * @brief Keep the existing data.  If this flag is not given then the file is cleared in write mode.
 *
 * @val FileOpenMode OPEN_QUIET
 * @brief Don't print any warning message if the file could not be opened.
 *
 * @val FileOpenMode OPEN_MASK
 * @brief (Internal) Bitmask to extract the read/write open mode.
 *
 * @val FileOpenMode OPEN_ASYNC
 * @brief (Internal) Open the file for asynchronous file access.  For asynchronous file access, use the @link AsyncFile
 *                   @endlink.
 *
 * @val FileOpenMode OPEN_TEMPORARY
 * @brief (Internal) Open automatically delete the file after close.  Use the <tt>openTemp</tt> methods to open
 *        temporary files.
 */

// --------------------------------------------------------------------------
// Enum FileOpenMode
// --------------------------------------------------------------------------

enum FileOpenMode {
    OPEN_RDONLY     = 1,
    OPEN_WRONLY     = 2,
    OPEN_RDWR       = 3,
    OPEN_MASK       = 3,
    OPEN_CREATE     = 4,
    OPEN_APPEND     = 8,
    OPEN_ASYNC      = 16,
    OPEN_TEMPORARY    = 32,
    OPEN_QUIET        = 128
}; //IOREV is it intended that two labels share the same value? What is OPEN_MASK anyway?

// --------------------------------------------------------------------------
// Direction Tags
// --------------------------------------------------------------------------

struct Input_;
typedef Tag<Input_> Input;

struct Output_;
typedef Tag<Output_> Output;

struct Bidirectional_;
typedef Tag<Bidirectional_> Bidirectional;

// --------------------------------------------------------------------------
// Metafunction DefaultOpenMode
// --------------------------------------------------------------------------

// helper metafunction to avoid ambigous partial specializations of the 1st/2nd argument
template <typename TDirection>
struct DefaultFileOpenMode_
{
    enum { VALUE = OPEN_RDWR | OPEN_CREATE | OPEN_APPEND };
};

template <>
struct DefaultFileOpenMode_<Input>
{
    enum { VALUE = OPEN_RDONLY };
};

template <>
struct DefaultFileOpenMode_<Output>
{
    enum { VALUE = OPEN_WRONLY | OPEN_CREATE };
};

template <typename T, typename TDirection = Bidirectional>
struct DefaultOpenMode:
    DefaultFileOpenMode_<TDirection> {};

// --------------------------------------------------------------------------
// Metafunction DefaultOpenTempMode
// --------------------------------------------------------------------------

template <typename T>
struct DefaultOpenTempMode
{
    enum { VALUE = OPEN_RDWR | OPEN_CREATE };
};

enum FileSeekMode {
    SEEK_BEGIN   = 0,
    SEEK_CURRENT = 1
#ifndef SEEK_END
  , SEEK_END     = 2
#endif
}; //IOREV why not use constants SEEK_SET, SEEK_CUR, SEEK_END from cstdio?

//////////////////////////////////////////////////////////////////////////////
// result type of asynch. functions
// you have to call release(AsyncRequest<T>) after a finished *event based* transfer
struct AsyncDummyRequest
{
    AsyncDummyRequest()
    {}

    AsyncDummyRequest(AsyncDummyRequest &, Move)
    {}

#ifdef SEQAN_CXX11_STANDARD
    AsyncDummyRequest(AsyncDummyRequest &&)
    {}
#endif

private:
    AsyncDummyRequest(AsyncDummyRequest const &)
    {}
};

template <>
struct HasMoveConstructor<AsyncDummyRequest> : True {};



/*!
 * @class AsyncRequest
 * @headerfile <seqan/file.h>
 * @brief Associated with an asynchronous I/O request.
 *
 * @signature template <typenam TFile>
 *            struct AsyncRequest;
 *
 * @tparam TFile The file type.
 *
 * This structure is used to identify asynchronous requests after their initiation.
 */

    template < typename T >
    struct AsyncRequest
    {
//IOREV _stub_ this seems not to be implemented at all, most functions are commented
        typedef AsyncDummyRequest Type;
    };


// ============================================================================
// Exceptions
// ============================================================================

// ----------------------------------------------------------------------------
// Exception IOError
// ----------------------------------------------------------------------------

/*!
 * @class IOError
 * @headerfile <seqan/basic.h>
 * @brief Input/Output error exception.
 *
 * @signature typedef std::ios_base::failure IOError;
 *
 *
 * @fn IOError::IOError
 * @brief Constructor
 *
 * @signature IOError::IOError(msg, errorCode);
 *
 * @param[in] msg       Message as <tt>std::string</tt>.
 * @param[in] errorCode The error code as an <tt>int</tt>.
 */

typedef std::ios_base::failure  IOError;

// ----------------------------------------------------------------------------
// Exception FileOpenError
// ----------------------------------------------------------------------------

struct FileOpenError : IOError
{
public:
    FileOpenError(const char *fileName):
        IOError((std::string)"Could not open file " + fileName)
    {}

protected:
    FileOpenError(std::string const &msg):
        IOError(msg)
    {}
};

// ----------------------------------------------------------------------------
// Exception UnknownFileFormat
// ----------------------------------------------------------------------------

struct UnknownFileFormat : FileOpenError
{
    UnknownFileFormat():
        FileOpenError(std::string("Could not detect file format"))
    {}

    UnknownFileFormat(const char *fileName):
        FileOpenError(std::string("Could not detect file format of ") + fileName)
    {}
};

// ----------------------------------------------------------------------------
// Exception UnknownExtensionError
// ----------------------------------------------------------------------------

struct UnknownExtensionError : FileOpenError
{
    UnknownExtensionError(const char *fileName):
        FileOpenError(std::string("Unknown file extension of ") + fileName)
    {}
};

}  // namespace seqan;

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_FILE_INTERFACE_H_
