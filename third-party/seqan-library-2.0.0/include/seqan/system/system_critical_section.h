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
// Critical Section class. In conjunction with a condition object it allows
// to suspend a thread until another wakes it up.
// ==========================================================================

#ifndef SEQAN_HEADER_SYSTEM_CRITICAL_SECTION_H_
#define SEQAN_HEADER_SYSTEM_CRITICAL_SECTION_H_

namespace seqan {

#ifdef PLATFORM_WINDOWS

struct CriticalSection
{
    CRITICAL_SECTION data_cs;

    CriticalSection()
    {
        InitializeCriticalSection(&data_cs);
    }

    ~CriticalSection()
    {
        DeleteCriticalSection(&data_cs);
    }
};

#else

struct CriticalSection
{
    pthread_mutex_t data_cs;

    CriticalSection()
    {
        int result = pthread_mutex_init(&data_cs, NULL);
        ignoreUnusedVariableWarning(result);
        SEQAN_ASSERT_EQ(result, 0);
    }

    ~CriticalSection()
    {
        int result = pthread_mutex_destroy(&data_cs);
        ignoreUnusedVariableWarning(result);
        SEQAN_ASSERT_EQ(result, 0);
    }
};

#endif

inline void
lock(CriticalSection &cs)
{
#ifdef PLATFORM_WINDOWS
    EnterCriticalSection(&cs.data_cs);
#else
    int result = pthread_mutex_lock(&cs.data_cs);
    ignoreUnusedVariableWarning(result);
    SEQAN_ASSERT_EQ(result, 0);
#endif
}

inline void
unlock(CriticalSection &cs)
{
#ifdef PLATFORM_WINDOWS
    LeaveCriticalSection(&cs.data_cs);
#else
    int result = pthread_mutex_unlock(&cs.data_cs);
    ignoreUnusedVariableWarning(result);
    SEQAN_ASSERT_EQ(result, 0);
#endif
}

}

#endif  // SEQAN_HEADER_SYSTEM_CRITICAL_SECTION_H_
