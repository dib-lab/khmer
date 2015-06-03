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
// Condition class. In conjunction with a critical section (given in the
// c'tor) this class allows to suspend a thread until another wakes it up via
// signal().
// ==========================================================================

#ifndef SEQAN_HEADER_SYSTEM_CONDITION_H_
#define SEQAN_HEADER_SYSTEM_CONDITION_H_

namespace seqan {

#ifdef PLATFORM_WINDOWS

struct Condition
{
    enum { Infinite = INFINITE };

    CONDITION_VARIABLE  data_cond;
    CriticalSection     *csPtr;

    explicit
    Condition(CriticalSection &cs) :
        csPtr(&cs)
    {
        InitializeConditionVariable(&data_cond);
    }

    ~Condition()
    {
        WakeAllConditionVariable(&data_cond);
    }
};

#else

struct Condition
{
    enum { Infinite = LONG_MAX };

    pthread_cond_t      data_cond;
    CriticalSection     *csPtr;

    explicit
    Condition(CriticalSection &cs) :
        csPtr(&cs)
    {
        int result = pthread_cond_init(&data_cond, NULL);
        ignoreUnusedVariableWarning(result);
        SEQAN_ASSERT_EQ(result, 0);
    }

    ~Condition()
    {
        int result = pthread_cond_destroy(&data_cond);
        ignoreUnusedVariableWarning(result);
        SEQAN_ASSERT_EQ(result, 0);
    }
};

#endif

inline void
waitFor(Condition &cond)
{
#ifdef PLATFORM_WINDOWS
    BOOL success = SleepConditionVariableCS(&cond.data_cond, &cond.csPtr->data_cs, INFINITE);
    ignoreUnusedVariableWarning(success);
    SEQAN_ASSERT(success);
#else
    int result = pthread_cond_wait(&cond.data_cond, &cond.csPtr->data_cs);
    ignoreUnusedVariableWarning(result);
    SEQAN_ASSERT_EQ(result, 0);
#endif
}

inline void
waitFor(Condition &cond, long timeoutMilliSec, bool & inProgress)
{
#ifdef PLATFORM_WINDOWS
    inProgress = (SleepConditionVariableCS(&cond.data_cond, &cond.csPtr->data_cs, timeoutMilliSec) == 0);
    if (inProgress)
        SEQAN_ASSERT_EQ(GetLastError(), ERROR_TIMEOUT);
#else
    if (timeoutMilliSec != Condition::Infinite)
    {
        timespec ts;
        ts.tv_sec = timeoutMilliSec / 1000;
        ts.tv_nsec = (timeoutMilliSec % 1000) * 1000;
        int result = pthread_cond_timedwait(&cond.data_cond, &cond.csPtr->data_cs, &ts);
        inProgress = (result == ETIMEDOUT);
        SEQAN_ASSERT(result == 0 || inProgress);
    }
    else
    {
        inProgress = false;
        waitFor(cond);
    }
#endif
}

inline void
signal(Condition &cond)
{
#ifdef PLATFORM_WINDOWS
    WakeAllConditionVariable(&cond.data_cond);
#else
    int result = pthread_cond_broadcast(&cond.data_cond);
    ignoreUnusedVariableWarning(result);
    SEQAN_ASSERT_EQ(result, 0);
#endif
}

}

#endif  // SEQAN_HEADER_SYSTEM_CONDITION_H_
