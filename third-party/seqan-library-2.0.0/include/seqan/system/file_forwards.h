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

#ifndef SEQAN_HEADER_FILE_MANUAL_FORWARDS_H
#define SEQAN_HEADER_FILE_MANUAL_FORWARDS_H

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

//////////////////////////////////////////////////////////////////////////////
// CLASSES
//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN {

//____________________________________________________________________________
// TagAllocateAligned_

struct TagAllocateAligned_;           // "include/seqan/file/file_async.h"(283)
struct AiocbWrapper;

//////////////////////////////////////////////////////////////////////////////
// TYPEDEFS
//////////////////////////////////////////////////////////////////////////////

//____________________________________________________________________________
// (*sighandler_t)(int)

typedef void (*sighandler_t)(int);           // "include/seqan/file/file_async.h"(258)

//____________________________________________________________________________
// TagAllocateAligned

typedef Tag<TagAllocateAligned_> const TagAllocateAligned;           // "include/seqan/file/file_async.h"(284)


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//____________________________________________________________________________
// allocate

template <typename T, typename TValue, typename TSize> inline void allocate(T const & me, TValue * & data, TSize count, TagAllocateAligned const);           // "include/seqan/file/file_async.h"(292)
template <typename TSpec, typename TValue, typename TSize> inline void allocate( File<Async<TSpec> > const & me, TValue * & data, TSize count);           // "include/seqan/file/file_async.h"(351)

//____________________________________________________________________________
// asyncReadAt

template <typename TSpec, typename TValue, typename TSize, typename TPos > bool asyncReadAt(File<Async<TSpec> > & me, TValue *memPtr, TSize const count, TPos const fileOfs, AiocbWrapper &request);           // "include/seqan/file/file_async.h"(136)

//____________________________________________________________________________
// asyncWriteAt

template <typename TSpec, typename TValue, typename TSize, typename TPos > bool asyncWriteAt(File<Async<TSpec> > & me, const TValue *memPtr, TSize const count, TPos const fileOfs, AiocbWrapper &request);           // "include/seqan/file/file_async.h"(160)

//____________________________________________________________________________
// cancel

template <typename TSpec> inline bool cancel(File<Async<TSpec> > & me, AiocbWrapper &request);           // "include/seqan/file/file_async.h"(242)

//____________________________________________________________________________
// deallocate

template <typename T, typename TValue, typename TSize> inline void deallocate( T const & me, TValue * data, TSize count, TagAllocateAligned const);           // "include/seqan/file/file_async.h"(310)
template <typename TSpec, typename TValue, typename TSize> inline void deallocate( File<Async<TSpec> > const & me, TValue * data, TSize count);           // "include/seqan/file/file_async.h"(360)

//____________________________________________________________________________
// error

inline int error(AiocbWrapper const &request);           // "include/seqan/file/file_async.h"(246)

//____________________________________________________________________________
// fileExists

inline bool fileExists(const char *fileName);           // "include/seqan/file/file_sync.h"(189)

//____________________________________________________________________________
// fileUnlink

inline bool fileUnlink(const char *fileName);           // "include/seqan/file/file_sync.h"(194)

//____________________________________________________________________________
// flush

template <typename TSpec> inline bool flush(File<Async<TSpec> > & me);           // "include/seqan/file/file_async.h"(182)

//____________________________________________________________________________
// printRequest

inline void printRequest(AiocbWrapper &request, const char *_hint);           // "include/seqan/file/file_async.h"(125)
inline void printRequest(AiocbWrapper &request);           // "include/seqan/file/file_async.h"(125)

//____________________________________________________________________________
// read

template <typename TSpec, typename TValue, typename TSize > inline bool read(File<Sync<TSpec> > & me, TValue *memPtr, TSize const count);           // "include/seqan/file/file_sync.h"(226)

//____________________________________________________________________________
// release

template <typename TSpec> inline void release(File<Async<TSpec> > & me, AiocbWrapper const &request);           // "include/seqan/file/file_async.h"(255)

//____________________________________________________________________________
// _returnValue

inline int _returnValue(AiocbWrapper &request);           // "include/seqan/file/file_async.h"(250)

//____________________________________________________________________________
// waitFor

inline bool waitFor(AiocbWrapper &request);           // "include/seqan/file/file_async.h"(189)
inline bool waitFor(AiocbWrapper &request, long timeoutMilliSec, bool &inProgress);           // "include/seqan/file/file_async.h"(204)

//____________________________________________________________________________
// waitForAny

template <typename TSize > inline TSize waitForAny(AiocbWrapper const * const contexts[], TSize count);           // "include/seqan/file/file_async.h"(223)
template <typename TSize > inline TSize waitForAny(AiocbWrapper const * const contexts[], TSize count, long timeoutMilliSec);           // "include/seqan/file/file_async.h"(231)

//____________________________________________________________________________
// write

template <typename TSpec, typename TValue, typename TSize > inline bool write(File<Sync<TSpec> > & me, TValue const *memPtr, TSize const count);           // "include/seqan/file/file_sync.h"(231)

} //namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////

#endif

