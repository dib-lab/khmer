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

#ifndef SEQAN_HEADER_SYSTEM_THREAD_H
#define SEQAN_HEADER_SYSTEM_THREAD_H

namespace SEQAN_NAMESPACE_MAIN
{

#ifdef PLATFORM_WINDOWS

    static SECURITY_ATTRIBUTES ThreadDefaultAttributes = {
        sizeof(SECURITY_ATTRIBUTES),
        NULL,
        true
    };

    template <typename Worker>
    struct Thread
    {
//IOREV _notio_
        typedef HANDLE Handle;

        Handle hThread;
        DWORD  hThreadID;
        Worker worker;

        Thread() {}

        template <typename TArg>
        Thread(TArg &arg):
            worker(arg) {}

        ~Thread() {
            if (*this) {
                cancel();
                wait();
            }
        }

        inline bool open(BOOL initital = false) {
            return hThread = CreateThread(
                &ThreadDefaultAttributes,    // default security attributes 
                0,                           // use default stack size  
                &_start,                     // thread function 
                this,                        // argument to thread function 
                0,                           // use default creation flags 
                &hThreadID);                 // returns the thread identifier 
        }

        inline bool close() {
            if (CloseHandle(hThread)) return true;
			hThread = NULL;
			return false;
        }

        inline bool cancel(DWORD exitCode = 0) {
            return !TerminateThread(hThread, exitCode);
        }

        inline bool wait(DWORD timeoutMilliSec = INFINITE) {
            return WaitForSingleObject(hThread, timeoutMilliSec) != WAIT_TIMEOUT;
        }

        inline operator bool() const {
            return hThread != NULL;
        }

    private:

        Thread(Thread const &) {
            // resource copying is not yet supported (performance reason)
            // it needs a reference counting technique
        }

        static DWORD WINAPI _start(LPVOID _this) {
            reinterpret_cast<Thread*>(_this)->worker.run(&reinterpret_cast<Thread*>(_this));
			return 0;	// return value should indicate success/failure
        }
    };
    
#else

    template <typename Worker>
    struct Thread
    {
//IOREV _notio_
        typedef pthread_t* Handle;

        pthread_t data, *hThread;
        Worker worker;

        Thread() {}

        template <typename TArg>
        Thread(TArg &arg):
            worker(arg) {}

        ~Thread() {
            if (*this) {
                cancel();
                wait();
            }
        }

        inline bool open()
        {
            if (!pthread_create(&data, NULL, _start, this) && (hThread = &data)) {
                return true;
            } else
                return false;
        }

        inline bool close() {
            return cancel() && wait() && !(hThread == NULL);
        }

        inline bool cancel() {
            return !(pthread_cancel(data));
        }

        inline bool wait() {
            return !(pthread_join(data, NULL));
        }

        inline bool wait(void* &retVal) {
            return !(pthread_join(data, &retVal));
        }

        inline bool detach() {
            return !(pthread_detach(data));
        }

        inline operator bool() const {
            return hThread != NULL;
        }

    private:

        Thread(Thread const &) {
            // resource copying is not yet supported (performance reason)
            // it needs a reference counting technique
        }

        static void* _start(void* _this) {
            reinterpret_cast<Thread*>(_this)->worker.run(&reinterpret_cast<Thread*>(_this));
			return 0;
        }
    };
    
#endif


	//////////////////////////////////////////////////////////////////////////////
	// global thread functions

	template <typename TWorker>
	inline bool open(Thread<TWorker> &m) {
		return m.open();
	}

	template <typename TWorker>
	inline bool run(Thread<TWorker> &m) {
		return m.open();
	}

	template <typename TWorker>
	inline bool close(Thread<TWorker> &m) {
		return m.close();
	}

	template <typename TWorker>
	inline bool kill(Thread<TWorker> &m) {
		return m.close();
	}

	template <typename TWorker>
	inline bool waitFor(Thread<TWorker> &m) {
		return m.wait();
	}

}

#endif
