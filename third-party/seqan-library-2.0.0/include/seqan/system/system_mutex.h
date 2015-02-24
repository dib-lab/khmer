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

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#ifndef SEQAN_HEADER_SYSTEM_MUTEX_H
#define SEQAN_HEADER_SYSTEM_MUTEX_H

namespace SEQAN_NAMESPACE_MAIN
{

#ifdef PLATFORM_WINDOWS

    static SECURITY_ATTRIBUTES MutexDefaultAttributes = {
        sizeof(SECURITY_ATTRIBUTES),
        NULL,
        true
    };

    struct Mutex
    {
        typedef HANDLE Handle;

        Handle hMutex;

        Mutex():
            hMutex(NULL) {}

        Mutex(BOOL initial) {
            SEQAN_DO_SYS2(open(initial), "Could not create Mutex");
        }

        // Move constructors
        Mutex(Mutex & other, Move) :
            hMutex(other.hMutex)
        {
            other.hMutex = NULL;
        }

#ifdef SEQAN_CXX11_STANDARD
        Mutex(Mutex && other) :
            hMutex(other.hMutex)
        {
            other.hMutex = NULL;
        }
#endif

        ~Mutex() {
            if (*this)
                SEQAN_DO_SYS2(close(), "Could not destroy Mutex");
        }

        inline bool open(BOOL initial = false) {
            return (hMutex = CreateMutex(&MutexDefaultAttributes, initial, NULL)) != NULL;
        }

        inline bool close() {
            bool success = CloseHandle(hMutex);
            hMutex = NULL;
            return success;
        }

        inline bool lock(DWORD timeoutMilliSec = INFINITE) {
            return WaitForSingleObject(hMutex, timeoutMilliSec) != WAIT_TIMEOUT;
        }

        inline bool unlock() {
            return ReleaseMutex(hMutex) != 0;
        }

        inline operator bool() const {
            return hMutex != NULL;
        }

    private:

        Mutex(Mutex const &) :
            hMutex(NULL)
        {
            // we only support move construction (no copy-construction)
        }
    };

#else

    struct Mutex
    {
        typedef pthread_mutex_t* Handle;

        pthread_mutex_t data, *hMutex;

        Mutex():
            hMutex(NULL)
        {}

        Mutex(bool initial)
        {
            SEQAN_DO_SYS(open(initial));
        }

        // Move constructors
        Mutex(Mutex & other, Move) :
            hMutex(other.hMutex)
        {
            other.hMutex = NULL;
        }

#ifdef SEQAN_CXX11_STANDARD
        Mutex(Mutex && other) :
            hMutex(other.hMutex)
        {
            other.hMutex = NULL;
        }
#endif

        ~Mutex() {
            if (*this)
                SEQAN_DO_SYS(close());
        }

        inline bool open(bool initial = false)
        {
            if (!pthread_mutex_init(&data, NULL) && (hMutex = &data)) {
                if (initial) return lock();
                return true;
            } else
                return false;
        }

        inline bool close() {
            bool success = (pthread_mutex_destroy(hMutex) == 0);
            hMutex = NULL;
            return success;
        }

        inline bool lock() {
            return !pthread_mutex_lock(hMutex);
        }

        inline bool unlock() {
            return !pthread_mutex_unlock(hMutex);
        }

        inline operator bool() const {
            return hMutex != NULL;
        }

    private:

        Mutex(Mutex const &) :
            hMutex(NULL)
        {
            // we only support move construction (no copy-construction)
        }

    };

#endif

    template <>
    struct HasMoveConstructor<Mutex> : True {};

    //////////////////////////////////////////////////////////////////////////////
    // global mutex functions

    inline bool open(Mutex &m, bool initial) {
        return m.open(initial);
    }

    inline bool open(Mutex &m) {
        return open(m, false);
    }

    inline bool close(Mutex &m) {
        return m.close();
    }

    inline bool lock(Mutex &m) {
        return m.lock();
    }

    inline bool unlock(Mutex &m) {
        return m.unlock();
    }

}

#endif
