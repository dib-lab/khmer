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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#ifndef SEQAN_HEADER_SYSTEM_SEMAPHORE_H
#define SEQAN_HEADER_SYSTEM_SEMAPHORE_H

namespace SEQAN_NAMESPACE_MAIN
{

#ifdef PLATFORM_WINDOWS

    static SECURITY_ATTRIBUTES SemaphoreDefaultAttributes = {
        sizeof(SECURITY_ATTRIBUTES),
        NULL,
        true
    };

    struct Semaphore
    {
        typedef LONG Type;
        typedef HANDLE Handle;
        enum { MAX_VALUE = MAXLONG };
        
        Handle hSemaphore;

        Semaphore(Type init = 0, Type max = MAX_VALUE)
        {
            hSemaphore = CreateSemaphore(&SemaphoreDefaultAttributes, init, max, NULL);
            SEQAN_ASSERT_MSG(hSemaphore != NULL, "Could not create Semaphore!");
        }

        ~Semaphore()
        {
            int res = CloseHandle(hSemaphore);
            (void)res;  // Only used in assertion.
            SEQAN_ASSERT_MSG(res != 0, "Could not destroy Semaphore!");
        }

        bool lock(DWORD timeoutMilliSec = INFINITE) {
            return WaitForSingleObject(hSemaphore, timeoutMilliSec) != WAIT_TIMEOUT;
        }

        void unlock()
        {
            int res = ReleaseSemaphore(hSemaphore, 1, NULL);
            (void)res;  // Only used in assertion.
            SEQAN_ASSERT_MSG(res != 0, "Could not unlock Semaphore!");
        }

    private:

        Semaphore(Semaphore const &) {
            // resource copying is not yet supported (performance reason)
            // it needs a reference counting technique
        }

    };

#else

    struct Semaphore
    {
        typedef unsigned int Type;
        typedef sem_t* Handle;

        sem_t data, *hSemaphore;

        Semaphore(Type init = 0):
            hSemaphore(&data)
        {
            bool res = !sem_init(hSemaphore, 0, init);
            (void)res;  // Only used in assertion.
            SEQAN_ASSERT(res);
        }

        ~Semaphore() {
            bool res = sem_destroy(hSemaphore);
            (void)res;  // Only used in assertion.
            SEQAN_ASSERT(res);
        }

        void lock() {
            bool res = !sem_wait(hSemaphore);
            (void)res;  // Only used in assertion.
            SEQAN_ASSERT(res);
        }

        void unlock() {
            bool res = !sem_post(hSemaphore);
            (void)res;  // Only used in assertion.
            SEQAN_ASSERT(res);
        }

    private:

        Semaphore(Semaphore const &) {
            // resource copying is not yet supported (performance reason)
            // it needs a reference counting technique
        }

    };


#endif

}

#endif
