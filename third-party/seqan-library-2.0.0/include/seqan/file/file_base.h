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
// Defines basic file access functions.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_FILE_BASE_H_
#define SEQAN_INCLUDE_SEQAN_FILE_BASE_H_

/* IOREV
 * _doc_
 *
 * base class with SPecs and Tags
 * also contains standard calls for IO as wrappers around members
 * in system/file_sync.h and system/file_ssync.h
 * (hese files are built around c++ fstream IO)
 * well documented (in comparison to other files)
 *
 * SEQAN_DIRECTIO Macro mentioned here but not documented or tested
 */

namespace seqan {

// Manual Forward.
template < typename TSpec, typename TPos >
inline typename Position< File<TSpec> >::Type seek(File<TSpec> &me, TPos const fileOfs, int origin);
template < typename TSpec, typename TPos >
inline typename Position< File<TSpec> >::Type seek(File<TSpec> &me, TPos const fileOfs);

    //////////////////////////////////////////////////////////////////////////////
    // generic open/close interface

/*!
 * @fn File#open
 * @brief Opens a file, stream, or persistent string.
 *
 * @signature bool open(file, fileName, openMode);
 *
 * @param[in,out] file     The File to open.
 * @param[in]     fileName A <tt>char const *</tt> string containing the file name.
 * @param[in]     openMode Combination of flags defining how the file should be opened.  See @link FileOpenMode
 *                         @endlink for more details.  Type: <tt>int</tt>.  If you omit the <tt>OPEN_APPEND</tt> flag in
 *                         write mode, the file will be cleared when opened.  Default: <tt>OPEN_RDWR | OPEN_CREATE |
 *                         OPEN_APPEND</tt>.
 *
 * @return bool <tt>true</tt> on success, <tt>false</tt> on failure.
 */

    template < typename TSpec >
    inline bool open(File<TSpec> &me, const char *fileName, int openMode)
    {
//IOREV resolves to member in file_(a)sync.h which resolves to fstream
        return me.open(fileName, openMode);
    }

    template < typename TSpec >
    inline bool open(File<TSpec> &me, const char *fileName)
    {
//IOREV
        return open(me, fileName, DefaultOpenMode<File<TSpec> >::VALUE);
    }

/*!
 * @fn File#openTemp
 * @brief Opens a temporary file.
 *
 * @signature bool openTemp(file);
 *
 * @param[in,out] file The File object to open the temporary file.
 *
 * @return bool <tt>true</tt> on success, <tt>false</tt> on failure.
 *
 * @section Remarks
 *
 * After closing this file will be deleted automatically.  The openmode (see @link File#open @endlink) is <tt>OPEN_RDWR
 * | OPEN_CREATE</tt>.
 *
 * @return bool <tt>true</tt> on success <tt>false</tt> on failure.
 */

    template < typename TSpec >
    inline bool openTemp(File<TSpec> &me)
    {
//IOREV
        return me.openTemp();
    }

    template < typename TSpec >
    inline bool openTemp(File<TSpec> &me, int openMode)
    {
//IOREV
        return me.openTemp(openMode);
    }

    template < typename File >
    inline void reopen(File &, int)
    {
//IOREV _stub_ _nodoc_ This is currently a no-op. Is that intended?
    }

/*!
 * @fn File#close
 * @brief Close a file.
 *
 * @signature bool close(file);
 *
 * @param[in,out] file The File object to close.
 *
 * @return bool <tt>true</tt> on success, <tt>false</tt> on failure.
 */

    template < typename TSpec >
    inline bool close(File<TSpec> & me)
    {
//IOREV
        return me.close();
    }

    template < typename TSpec >
    inline unsigned sectorSize(File<TSpec> const & /*me*/)
    {
//IOREV _duplicate_ _nodoc_ duplicate or identical spec. in file_cstyle.h should'nt this be variable
        return 4096;
    }


    //////////////////////////////////////////////////////////////////////////////
    // generic read(At)/write(At) interface

/*!
 * @fn File#read
 * @brief Loads record from a file.
 *
 * @signature bool read(file, memPtr, count);
 *
 * @param[in,out] file   The File object.
 * @param[out]    memPtr A pointer to the first destination record in memory.
 * @param[in]     count  The amount of records to be read.
 *
 * @return bool <tt>true</tt> on success, <tt>false</tt> on failure.
 *
 * @section Remarks
 *
 * The records are read from the position pointed to by the current file pointer (see @link File#seek @endlink).
 */

    template < typename TSpec, typename TValue, typename TSize >
    inline bool read(File<TSpec> & me, TValue *memPtr, TSize const count)
    {
//IOREV
        typedef typename Size<File<TSpec> >::Type TFileSize;
        TFileSize nbytes = (TFileSize)count * (TFileSize)sizeof(TValue);
        return me.read(memPtr, nbytes) == nbytes;
    }

/*!
 * @fn File#write
 * @brief Saves records to a file.
 *
 * @signature bool write(file, memPtr, count);
 *
 * @param[in,out] file   The File object.
 * @param[in]     memPtr Pointer to the source for the data to write.
 * @param[in]     count  The number of records to write.
 *
 * @return bool <tt>true</tt> on success, <tt>false</tt> on failure.
 *
 * @section Remarks
 *
 * The records are written at the position pointed to by th ecurrent file pointer (see @link File#seek @endlink).
 */

    template < typename TSpec, typename TValue, typename TSize >
    inline bool write(File<TSpec> & me, TValue const *memPtr, TSize const count)
    {
//IOREV
        typedef typename Size<File<TSpec> >::Type TFileSize;
        TFileSize nbytes = (TFileSize)count * (TFileSize)sizeof(TValue);
        return me.write(memPtr, nbytes) == nbytes;
    }

/*!
 * @fn File#readAt
 * @brief Loads records from a specific position in a file.
 *
 * @signature bool readAt(file, memPtr, count, fileOfs);
 *
 * @param[in,out] file    The File object to read from.
 * @param[out]    memPtr  A pointer to the first destination record in memory.
 * @param[in]     count   The amount of records to be read.
 * @param[in]     fileOfs The absolute file position in bytes measured from the beginning.
 *
 * @return bool <tt>true</tt> on success and <tt>false</tt> on failure.
 */

    template < typename TFile, typename TValue, typename TSize, typename TPos >
    inline bool readAt(TFile & me, TValue *memPtr, TSize const count, TPos const fileOfs)
    {
//IOREV
        typedef typename Position<TFile>::Type pos_t;
        seek(me, (pos_t)fileOfs * (pos_t)sizeof(TValue));
        return read(me, memPtr, count);
    }

/*!
 * @fn File#writeAt
 * @brief Saves records to a specific position in a file.
 *
 * @signature bool writeAt(file, memPtr, count, fileOfs);
 *
 * @param[in,out] file    The File object to write to.
 * @param[in]     memPtr  Pointer to the memory to write.
 * @param[in]     count   The amount of records to be written.
 * @param[in]     fileOfs The absolute file position in bytes measured from the beginning.
 *
 * @return bool <tt>true</tt> on success, <tt>false</tt> on failure.
 */

    template < typename TFile, typename TValue, typename TSize, typename TPos >
    inline bool writeAt(TFile & me, TValue const *memPtr, TSize const count, TPos const fileOfs)
    {
//IOREV
        typedef typename Position<TFile>::Type pos_t;
        seek(me, (pos_t)fileOfs * (pos_t)sizeof(TValue));
        return write(me, memPtr, count);
    }



    //////////////////////////////////////////////////////////////////////////////
    // generic seek/tell/size/resize interface

/*!
 * @fn File#seek
 * @brief Changes the current file pointer.
 *
 * @signature TPosition seek(file, fileOfs[, origin]);
 *
 * @param[in,out] file    The File object to seek in.
 * @param[in]     fileOfs A file offset measured in bytes relative to <tt>origin</tt>.
 * @param[in]     origin  Selects the origin from where to calculate the new position.  One of <tt>SEEK_BEGIN</tt>,
 *                        <tt>SEEK_CURRENT</tt>, and <tt>SEEK_END<tt> (origin is beginning, current pointer, end of
 *                        the file).  Default: <tt>SEEK_BEGIN</tt>.
 *
 * @return TPosition The new file position measured in bytes from the beginning.
 */

    template < typename TSpec, typename TPos >
    inline typename Position< File<TSpec> >::Type seek(File<TSpec> &me, TPos const fileOfs, int origin)
    {
//IOREV
        typedef typename Position< File<TSpec> >::Type TFilePos;
        TFilePos newOfs = me.seek(fileOfs, origin);
        #if SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING
            if (origin == SEEK_BEGIN && newOfs != (TFilePos)fileOfs) {
                std::cerr << "seek returned " << std::hex << newOfs << " instead of " << fileOfs << std::dec << std::endl;
            }
        #endif
        return newOfs;
    }

    template < typename TSpec, typename TPos >
    inline typename Position< File<TSpec> >::Type seek(File<TSpec> &me, TPos const fileOfs)
    {
//IOREV
        return seek(me, fileOfs, SEEK_BEGIN);
    }

/*!
 * @fn File#tell
 * @brief Gets the current file pointer.
 *
 * @signature TPosition tell(file);
 *
 * @param[in] file The File object to query for the current position.
 *
 * @return TPosition The current position in the file.
 */

    template < typename TSpec >
    inline typename Position< File<TSpec> >::Type tell(File<TSpec> &me)
    {
//IOREV
        return me.tell();
    }

/*!
 * @fn File#rewind
 * @brief Sets the current file pointer to the beginning of a file.
 *
 * @signature void rewind(file);
 *
 * @param[in,out] file The file to reset the file pointer of.
 */

    template < typename File >
    inline void rewind(File &me)
    {
//IOREV
        seek(me, 0);
    }

/*!
 * @fn File#length
 * @brief Return the file size.
 *
 * @signature TSize length(file);
 *
 * @param[in] file The File object to query for its size.
 *
 * @return TSize The file size measured in bytes.
 */

    template < typename TSpec >
    inline typename Size<File<TSpec> >::Type length(File<TSpec> &me)
    {
//IOREV
        typename Size<File<TSpec> >::Type old_pos = tell(me);
        typename Size<File<TSpec> >::Type result = seek(me, 0, SEEK_END);
        seek(me, old_pos, SEEK_BEGIN);
        return result;
    }

/*!
 * @fn File#resize
 * @brief A file object.
 *
 * @signature void resize(file, newLength);
 *
 * @param[in,out] file      The File object to resize.
 * @param[in]     newLength The file size in bytes to resize to in bytes.
 */

    template < typename TSpec, typename TSize >
    inline void resize(File<TSpec> &me, TSize new_length)
    {
//IOREV possibly not standard-conformant, see resize() in file_cstyle.h
        typename Size<File<TSpec> >::Type old_pos = tell(me);
        seek(me, new_length, SEEK_BEGIN);
        setEof(me);
        seek(me, old_pos, SEEK_BEGIN);
    }

/*!
 * @fn File#setEof
 * @brief Sets the file end to the current pointer.
 *
 * @signature bool setEof(file);
 *
 * @param[in,out] file The File object to set the end of.
 *
 * @return bool <tt>true</tt> on success, <tt>false</tt> on failure.
 */

    template < typename TSpec >
    inline bool setEof(File<TSpec> &/*me*/)
    {
//IOREV _noop_ specialized for async file access but not for sync
        return true;
    }


    //////////////////////////////////////////////////////////////////////
    // Pseudo asynchronous Methods
    //////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////
    // callback based read/write
/*
    template < typename File, typename TValue, typename TSize,
               typename aCallback, typename aHint >
    inline typename AsyncRequest<File>::Type
    asyncRead(File & me, TValue *memPtr, TSize const count,
        aCallback* cb, aHint* hint)
    {
        result = read(me, memPtr, count);
        cb(hint);
        return NULL;
    }

    template < typename File, typename TValue, typename TSize,
               typename aCallback, typename aHint >
    inline typename AsyncRequest<File>::Type
    asyncWrite(File & me, TValue const *memPtr, TSize const count,
        aCallback* cb, aHint* hint)
    {
        write(me, memPtr, count);
        cb(hint);
        return NULL;
    }

    template < typename File, typename TValue, typename TSize, typename TPos,
               typename aCallback, typename aHint >
    inline typename AsyncRequest<File>::Type
    asyncReadAt(File & me, TValue *memPtr, TSize const count, TPos const fileOfs,
        aCallback* cb, aHint* hint)
    {
        readAt(me, memPtr, count, fileOfs);
        cb(hint);
        return NULL;
    }

    template < typename File, typename TValue, typename TSize, typename TPos,
               typename aCallback, typename aHint >
    inline typename AsyncRequest<File>::Type
    asyncWriteAt(File & me, TValue const *memPtr, TSize const count, TPos const fileOfs,
        aCallback* cb, aHint* hint)
    {
        result = writeAt(me, memPtr, count, fileOfs);
        cb(hint);
        return NULL;
    }


    //////////////////////////////////////////////////////////////////////
    // event based read/write

    template < typename File, typename TValue, typename TSize,
               typename aEvent >
    inline typename AsyncRequest<File>::Type
    asyncRead(File & me, TValue *memPtr, TSize const count,
        aEvent &event)
    {
        read(me, memPtr, count);
        event.signal();
        return NULL;
    }

    template < typename File, typename TValue, typename TSize,
               typename aEvent >
    inline typename AsyncRequest<File>::Type
    asyncWrite(File & me, TValue const *memPtr, TSize const count,
        aEvent &event)
    {
        write(me, memPtr, count);
        event.signal();
        return NULL;
    }

    template < typename File, typename TValue, typename TSize, typename TPos,
               typename aEvent >
    inline typename AsyncRequest<File>::Type
    asyncReadAt(File & me, TValue *memPtr, TSize const count, TPos const fileOfs,
        aEvent &event)
    {
        readAt(me, memPtr, count, fileOfs);
        event.signal();
        return NULL;
    }

    template < typename File, typename TValue, typename TSize, typename TPos,
               typename aEvent >
    inline typename AsyncRequest<File>::Type
    asyncWriteAt(File & me, TValue const *memPtr, TSize const count, TPos const fileOfs,
        aEvent &event)
    {
        writeAt(me, memPtr, count, fileOfs);
        event.signal();
        return NULL;
    }
*/

    //////////////////////////////////////////////////////////////////////
    // queue-less request based pseudo asychronous read/write

/*!
 * @fn File#asyncReadAt
 * @brief Asynchronously loads records from a specific position in a file.
 *
 * @signature bool asyncReadAt(file, memPtr, count, fileOfs, request);
 *
 * @param[in,out] file    The File object to read from.
 * @param[out]    memPtr  A pointer to the first destination record in memory.
 * @param[in]     count   The amount of records to be read.
 * @param[in]     fileOfs The absolute file position in bytes measured from the beginning.
 * @param[in]     request Reference to a structure that will be associated with this asynchronous request.  Type:
 *                        @link AsyncRequest @endlink.
 *
 * @return bool <tt>true</tt> on success, <tt>false</tt> on failure.
 */

    template < typename File, typename TValue, typename TSize, typename TPos,
               typename AsyncRequest >
    inline bool
    asyncReadAt(File & me, TValue *memPtr, TSize const count, TPos const fileOfs,
        AsyncRequest &)
    {
//IOREV _stub_ see general discussion about AsynRequest
        return readAt(me, memPtr, count, fileOfs);
    }

/*!
 * @fn File#asyncWriteAt
 * @brief Asynchronously writes records to a specific position in a file.
 *
 * @signature bool asyncWriteAt(file, memPtr, count, fileOfs, request);
 *
 * @param[in,out] file    The File object.
 * @param[in]     memPtr  A pointer to the first source record in memory.
 * @param[in]     count   The amount of records to be written.
 * @param[in]     fileOfs The absolute file position in bytes measured form the beginning.
 * @param[in]     request Reference to a structure that will be associated with this asynchronous request.
 *
 * @return bool <tt>true</tt> on succcess, <tt>false</tt> on failure.
 */

    template < typename File, typename TValue, typename TSize, typename TPos,
               typename TAsyncRequest >
    inline bool
    asyncWriteAt(File & me, TValue const *memPtr, TSize const count, TPos const fileOfs,
        TAsyncRequest &)
    {
//IOREV _stub_ see general discussion about AsynRequest
        return writeAt(me, memPtr, count, fileOfs);
    }


    //////////////////////////////////////////////////////////////////////
    // pseudo queue specific functions

/*!
 * @fn File#flush
 * @brief Waits for all open requests to complete.
 *
 * @signature void flush(file);
 *
 * @param[in,out] file The File object to flush.
 */

    template < typename TSpec >
    inline void flush(File<TSpec> &)
    {
//IOREV _noop_ specialized for async file access but not for sync
    }

/*!
 * @fn File#waitFor
 * @brief Waits for an asychronous request to complete.
 *
 * @signature bool waitFor(request[, timeout]);
 *
 * @param[in,out] request Reference to an AsyncRequest.
 * @param[in]     timeout A timeout value in milliseconds.  A value of 0 can be used to test for completion without
 *                        waiting.  Default: 0.
 *
 * @return bool <tt>true</tt> on completion, <tt>false</tt> on timeout.
 *
 * @section Remarks
 *
 * <tt>waitfor</tt> block sand suspends the calling thread process until <tt>request</tt> is completed or after
 * <tt>timeout</tt> milliseconds.
 */

    inline bool waitFor(AsyncDummyRequest &)
    {
//IOREV _noop_ see general discussion about AsynRequest
        return true;
    }

    template < typename TTime >
    inline bool waitFor(AsyncDummyRequest &, TTime, bool &inProgress)
    {
//IOREV _noop_ see general discussion about AsynRequest
        inProgress = false;
        return true;
    }

    // deprecated
    template < typename TSpec, typename AsyncRequest >
    inline void release(File<TSpec> &, AsyncRequest &)
    {
//IOREV _noop_ see general discussion about AsynRequest
    }

/*!
 * @fn File#cancel
 * @brief Cancels an asynchronous request.
 *
 * @signature bool cancel(file, request);
 *
 * @param[in,out] file    The File to cancel the request for.
 * @param[in]     request Reference to an AsyncRequest object.  Type: @link AsyncRequest @endlink.
 *
 * @return bool <tt>true</tt> on success, <tt>false</tt> on failure.
 */

    template < typename TSpec, typename AsyncRequest >
    inline bool cancel(File<TSpec> &, AsyncRequest &)
    {
//IOREV _noop_ see general discussion about AsynRequest
        return true;
    }


    // little helpers

    template <typename T1, typename T2> inline
    T1 enclosingBlocks(T1 _size, T2 _blockSize)
    {
//IOREV not sure what this does, but is used in several places
        return (_size + (T1)_blockSize - (T1)1) / (T1)_blockSize;
    }

    template <typename T1, typename T2> inline
    T1 alignSize(T1 _size, T2 _aligning)
    {
//IOREV not sure what this does, but is used in several places
        if (_size < (T1)_aligning)
            return _aligning;
        else
            return (_size / (T1)_aligning) * (T1)_aligning;
    }
}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_FILE_BASE_H_
