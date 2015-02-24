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

#ifndef SEQAN_HEADER_SYSTEM_MANUAL_FORWARDS_H
#define SEQAN_HEADER_SYSTEM_MANUAL_FORWARDS_H

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

//////////////////////////////////////////////////////////////////////////////
// CLASSES
//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN {

//____________________________________________________________________________
// Event

struct Event;           // "include/seqan/system/system_event.h"(18)

//____________________________________________________________________________
// Mutex

struct Mutex;           // "include/seqan/system/system_mutex.h"(16)

//____________________________________________________________________________
// Thread

template <typename Worker> struct Thread;           // "include/seqan/system/system_thread.h"(18)
//IOREV _todo_


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//____________________________________________________________________________
// close

inline bool close(Event &e);           // "include/seqan/system/system_event.h"(109)
inline bool close(Mutex &m);           // "include/seqan/system/system_mutex.h"(79)
template <typename TWorker> inline bool close(Thread<TWorker> &m);           // "include/seqan/system/system_thread.h"(97)

//____________________________________________________________________________
// kill

template <typename TWorker> inline bool kill(Thread<TWorker> &m);           // "include/seqan/system/system_thread.h"(102)

//____________________________________________________________________________
// lock

inline bool lock(Mutex &m);           // "include/seqan/system/system_mutex.h"(83)

//____________________________________________________________________________
// open

inline bool open(Event &e, bool initial);           // "include/seqan/system/system_event.h"(101)
inline bool open(Event &e);           // "include/seqan/system/system_event.h"(105)
inline bool open(Mutex &m, bool initial);           // "include/seqan/system/system_mutex.h"(71)
inline bool open(Mutex &m);           // "include/seqan/system/system_mutex.h"(75)
template <typename TWorker> inline bool open(Thread<TWorker> &m);           // "include/seqan/system/system_thread.h"(87)

//____________________________________________________________________________
// run

template <typename TWorker> inline bool run(Thread<TWorker> &m);           // "include/seqan/system/system_thread.h"(92)

//____________________________________________________________________________
// signal

inline bool signal(Event &e);           // "include/seqan/system/system_event.h"(131)

//____________________________________________________________________________
// unlock

inline bool unlock(Mutex &m);           // "include/seqan/system/system_mutex.h"(87)

//____________________________________________________________________________
// waitFor

inline bool waitFor(Event &e);           // "include/seqan/system/system_event.h"(113)
template <typename TTime > inline bool waitFor(Event &e, TTime timeoutMilliSec, bool &inProgress);           // "include/seqan/system/system_event.h"(118)
template <typename TWorker> inline bool waitFor(Thread<TWorker> &m);           // "include/seqan/system/system_thread.h"(107)

} //namespace SEQAN_NAMESPACE_MAIN


#endif

