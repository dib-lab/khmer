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
// Code for profiling.
// ==========================================================================

// TODO(holtgrew): This could use some cleanup.

#include <ctime>

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_PROFILING_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_PROFILING_H_

namespace seqan
{

// todo: substitute defines with inlines
#ifndef SEQAN_PROFILE

    #define SEQAN_PROSET(i,v)           do {} while (false)
    #define SEQAN_PROADD(i,v)           do {} while (false)
    #define SEQAN_PROSUB(i,v)           do {} while (false)
    #define SEQAN_PROVAL(i)             0
    #define SEQAN_PROEXTRAS(i)          do {} while (false)
    #define SEQAN_PROMARK(m)            do {} while (false)
    #define SEQAN_PROENDMARK(m)         do {} while (false)
    #define SEQAN_PRORESET              do {} while (false)
    #define SEQAN_PROGETTIME            0
    #define SEQAN_PROTIMESTART(a)       do {} while (false)
    #define SEQAN_PROTIMEDIFF(a)        0
    #define SEQAN_PROTIMEUPDATE(a)      0
    // replace malloc and free in external tools
    // with SEQAN_PROMALLOC and SEQAN_PROFREE to profile
    // their memory usage
    #define SEQAN_PROMALLOC(s)          malloc(s)
    #define SEQAN_PROFREE(p)            free(p)

#else

    #define SEQAN_PROSET(i,v)           _profileSet(i,v)
    #define SEQAN_PROADD(i,v)           _profileAdd(i,v)
    #define SEQAN_PROSUB(i,v)           _profileSub(i,v)
    #define SEQAN_PROVAL(i)             (ProfileData_<>::_proValue[i])
    #define SEQAN_PROEXTRAS(i)          {ProfileData_<>::_proExtraCount = i;}
    #define SEQAN_PROMARK(m)            _profileMark(m)
    #define SEQAN_PROENDMARK(m)         _profileEndMark(m)
    #define SEQAN_PRORESET              _profileReset()
    #define SEQAN_PROGETTIME            sysTime()
    #define SEQAN_PROTIMESTART(a)       _proFloat a = sysTime()
    #define SEQAN_PROTIMEDIFF(a)        (sysTime() - a)
    #define SEQAN_PROTIMEUPDATE(a)      (_profileUpdate(a))
    #define SEQAN_PROMALLOC(s)          _profileMalloc(s)
    #define SEQAN_PROFREE(p)            _profileFree(p)

#endif

#ifdef PLATFORM_WINDOWS
    typedef __int64   ProfileInt_; //IOREV _notio_
#else
    typedef int64_t ProfileInt_; //IOREV _notio_
#endif

    typedef double    _proFloat;


    typedef _proFloat ProfileTimeValue_; //IOREV _notio_

    enum ProfileConstants_ {
        SEQAN_PROPAGESIZE         = 4096, // B in byte
        SEQAN_PROFLOAT            = 0,
        SEQAN_PROINT              = 1,
        SEQAN_PROTIME             = 2,
        SEQAN_PROTYPEMASK         = 3,
        SEQAN_PROSTATE            = 4
    };

    enum ProfileValueIndex_ {
        SEQAN_PROSYSTIME          = 0,
        SEQAN_PROCPUTIME          = 1,
        SEQAN_PROMEMORY           = 2,    // current memory usage (state value)
        SEQAN_PROIO               = 3,    // IOs done (measured in Blocks of size B)
        SEQAN_PROIORANDOM         = 4,    // IOs calls done (read/write calls done)
        SEQAN_PROIOVOLUME         = 5,    // current disk usage (state value)
        SEQAN_PRODEPTH            = 6,    // algorithmic rec. depth or loop count
        SEQAN_PROOPENFILES        = 7,    // currently opened files
        SEQAN_PROIWAIT            = 8,    // waiting time (initiating)
        SEQAN_PROCWAIT            = 9,    // waiting time (completing)
        SEQAN_PROEXTRA1           = 10,
        SEQAN_PROEXTRA2           = 11,
        SEQAN_PROEXTRA3           = 12,
        SEQAN_PROINDEXCOUNT       = 13,
        SEQAN_PROEXTRACOUNT       = 3
    };

    const char ProfileValueType_[] = {
        SEQAN_PROTIME, 
        SEQAN_PROTIME, 
        SEQAN_PROINT + SEQAN_PROSTATE, 
        SEQAN_PROINT,
        SEQAN_PROINT,
        SEQAN_PROINT + SEQAN_PROSTATE, 
        SEQAN_PROINT + SEQAN_PROSTATE, 
        SEQAN_PROINT + SEQAN_PROSTATE, 
        SEQAN_PROFLOAT,
        SEQAN_PROFLOAT,
        SEQAN_PROFLOAT + SEQAN_PROSTATE,
        SEQAN_PROFLOAT + SEQAN_PROSTATE,
        SEQAN_PROFLOAT + SEQAN_PROSTATE
    };

    typedef ProfileTimeValue_ ProfileTStates_[SEQAN_PROINDEXCOUNT]; //IOREV _notio_
    typedef _proFloat  ProfileTTimes[SEQAN_PROINDEXCOUNT]; //IOREV _notio_



    struct ProfileFile_;
//IOREV

    template <typename T = void>
    struct ProfileData_
    {
//IOREV _notio_
        static ProfileTStates_  _proValue;
        static ProfileTTimes    _proLastUpdate;
        static int          _proExtraCount;
        
        static clock_t      _proCpuTimeLast;            // clock_t wraps around every 72mins
        static ProfileInt_      _proCpuTimeOffset;          // we have to work around this

        static ProfileFile_*    _proPFile;
        static ProfileFile_*    _proPFileStream;
    };

    template <typename T> ProfileTStates_   ProfileData_<T>::_proValue = {};
    template <typename T> ProfileTStates_   ProfileData_<T>::_proLastUpdate = {};
    template <typename T> int           ProfileData_<T>::_proExtraCount = 0;
    template <typename T> clock_t       ProfileData_<T>::_proCpuTimeLast = 0;
    template <typename T> ProfileInt_       ProfileData_<T>::_proCpuTimeOffset = 0;
    template <typename T> ProfileFile_*     ProfileData_<T>::_proPFile = NULL;
    template <typename T> ProfileFile_*     ProfileData_<T>::_proPFileStream = NULL;


    inline ProfileFile_* & _proPFile()          { return ProfileData_<>::_proPFile; }
//IOREV
    inline ProfileFile_* & _proPFileStream()    { return ProfileData_<>::_proPFileStream; }
//IOREV

/**
.Function.cpuTime
..cat:Miscellaneous
..summary:Returns the cpu time in seconds.
..signature:cpuTime()
..returns:A $double$, cpu time stamp in seconds.
...type:nolink:double
..remarks:
Calls $clock$ to retrieve the processor time used by the running thread.
This implies that the thread's processor time does not tick if the thread is suspended.
While this has its advantages, benchmarks should generally focus on wall clock time, not processor time.
Wall clock time is returned by @Function.sysTime@.
..see:Function.sysTime
..include:seqan/basic.h
*/


// HINT: The unit of all time functions is second.
    inline _proFloat cpuTime() {
        clock_t now = clock();
        if (ProfileData_<>::_proCpuTimeLast > now) {        // test for time wrap
            ProfileData_<>::_proCpuTimeOffset += (~0u);     // got one
            ProfileData_<>::_proCpuTimeOffset ++;
//          printf("\n!!WRAP!! old:%d, now:%d    ofs:%d\n",ProfileData_<>::_proCpuTimeLast,now,ProfileData_<>::_proCpuTimeOffset);
        }
        ProfileData_<>::_proCpuTimeLast = now;
        return (ProfileData_<>::_proCpuTimeOffset + now) / (_proFloat)CLOCKS_PER_SEC;
    }


/**
.Function.sysTime
..cat:Miscellaneous
..summary:Returns the system time in seconds.
..signature:sysTime()
..returns:A $double$, system time stamp in seconds.
...type:nolink:double
..remarks:In contrast to @Function.cpuTime@, the system time corresponds to the wall clock time under Linux and Mac OS X.
Under Windows @Function.sysTime@ returns the result of @Function.cpuTime@.
..remarks:Use this for benchmarking uner Linux and Mac Os X.
..remarks:Calls $clock_gettime$ under Linux and $gettimeofday$ under Mac OS X.
..see:Function.cpuTime
..example.text:
We can use @Function.sysTime@ to instrument our code for profiling/timing information quite robustly.
The following demonstrates how the Function.sysTime is used in many SeqAn apps for collecting timing information.
..example.code:
bool printTiming = true;

// ...

double startTime = sysTime();
// Do some complex calculation.
if (printTiming)
    std::cerr << "Some complex calculation too " << sysTime() - startTime << " s." << std::endl;
..include:seqan/basic.h
*/

    #ifdef PLATFORM_WINDOWS
//        inline _proFloat sysTime() { return GetTickCount() * 1e-3; }
        inline _proFloat sysTime() { return ( (_proFloat) clock() ) / CLOCKS_PER_SEC; }
    #else

        #include <unistd.h>
        #if _POSIX_TIMERS > 0
            #ifndef SEQAN_USE_CLOCKGETTIME
            #define SEQAN_USE_CLOCKGETTIME
            #endif
        #endif
        
        #ifndef SEQAN_USE_CLOCKGETTIME
        /* some systems e.g. darwin have no clock_gettime */
        
            #include <sys/time.h>
            
            inline _proFloat sysTime() {
                struct timeval tp;
                gettimeofday(&tp, NULL);
                return tp.tv_sec + tp.tv_usec * 1e-6;
            }

        #else

            inline _proFloat sysTime() {
                struct timespec tp;
                clock_gettime(CLOCK_MONOTONIC, &tp);
                return tp.tv_sec + tp.tv_nsec * 1e-9;
            }

        #endif

    #endif

    
    struct ProfileFile_ {
//IOREV not generic, uses FILE* instead of File() and custom IO

        FILE   *out;
        bool   running;

        _proFloat dumpStep;            // 0 .. manual dump mode, >0 .. live stream
        _proFloat dumpNext;        

        ProfileTStates_ all, last;
        ::std::string mark;
        unsigned    lines;

        ProfileFile_() {
            running = false;
        }

        ProfileFile_(char const *fname, _proFloat _dumpStep = 300.0) { // five minutes default dump interval
            running = false;
            start(fname, _dumpStep);
        }

        ~ProfileFile_() {
            if (running) stop();
        }

        inline void start(char const *fname, _proFloat _dumpStep = 300.0, bool append = false) {
            if (append)
                out = fopen(fname, "a");
            else {
                out = fopen(fname, "w");
                dumpHeader();
            }

            if (!out) printf("WARNING: proFile could not be opened.\n");

            setTime(ProfileData_<>::_proValue);
            syncAll(all);
            syncAll(last);
            running      = true;
            lines        = 0;
            dumpStep     = _dumpStep;
            dumpNext     = sysTime();
            dump(last);
        }

        inline void stop() {
            dump(last);
            maximize(all, last);
            if (dumpStep == 0) {
                mark = "Zusammenfassung";
                dump(all);
            }
            fclose(out);
            running = false;
        }

        inline void syncTime(ProfileTStates_ &dst) {
            ::std::memcpy(dst, ProfileData_<>::_proValue, 2 * sizeof(ProfileTimeValue_));
        }

        inline void sync(ProfileTStates_ &dst) {
            ::std::memcpy(&(dst[2]), &(ProfileData_<>::_proValue[2]), sizeof(ProfileTStates_) - 2 * sizeof(ProfileTimeValue_));
        }

        inline void syncAll(ProfileTStates_ &dst) {
            ::std::memcpy(dst, ProfileData_<>::_proValue, sizeof(ProfileTStates_));
        }

        inline static void setTime(ProfileTStates_ &dst) {
            dst[0] = sysTime();
            dst[1] = cpuTime();
        }

        inline void maximize(ProfileTStates_ &dst, ProfileTStates_ const &src) {
            for(int i = 0; i < SEQAN_PROINDEXCOUNT; ++i)
                if (((ProfileValueType_[i] & SEQAN_PROSTATE) != 0))
                    if (dst[i] < src[i])
                        dst[i] = src[i];
        }

        inline void dumpTab() {
            if (!bol)
                fprintf(out, " \t");
            bol = false;
        }

        inline void dumpEndl() { fprintf(out, "\n"); }

        inline void dumpHeader() {
            fprintf(out, "\"Echtzeit\"\t\"CPU-Zeit\"\t\"Speicher\"\t\"I/O-Zugriffe\"\t\"wahlfreie I/Os\"\t\"I/O-Volumen\"\t\"Rekursionstiefe\"\t\"Offene Dateien\"\t\"Idle-Zeit vor I/O\"\t\"Idle-Zeit nach I/O\"\n");
        }

        inline void dumpTime(_proFloat seconds) {
            if (seconds < 0) {
                fputc('-', out);
                seconds = -seconds;
            }
            int secs    = (int)seconds;
            int mins    = secs/60;  secs -= 60*mins;
            int hours   = mins/60;  mins -= 60*hours;
            fprintf(out, "%d:%02d:%02d", hours, mins, secs);
        }

        inline void dumpTimeEx(_proFloat seconds) {
            int milli   = (int)(seconds * 1000.0);
            int secs    = (int)seconds;
            int mins    = secs/60;  secs -= 60*mins;
            int hours   = mins/60;  mins -= 60*hours;
            fprintf(out, "%d:%02d:%02d.%03d", hours, mins, secs, milli);
        }

        inline void dumpValue(ProfileTStates_ &stat, int valNum) {
            _proFloat f = stat[valNum];
            if ((ProfileValueType_[valNum] & SEQAN_PROSTATE) == 0)
                f = ProfileData_<>::_proValue[valNum] - f;

            switch (ProfileValueType_[valNum] & SEQAN_PROTYPEMASK) {
                case SEQAN_PROINT:                                      // state value -> print last seen maximum
                    fprintf(out, "%.0f", f);
                    break;

                case SEQAN_PROFLOAT:
                    fprintf(out, "%f", f);
                    break;

                case SEQAN_PROTIME:
                    dumpTimeEx(f);
            }
        }

        inline void dumpSysValues(ProfileTStates_ &stat) {
            for(int i = 0; i < SEQAN_PROINDEXCOUNT - SEQAN_PROEXTRACOUNT; ++i) {
                dumpTab();
                dumpValue(stat, i);
            }
        }

        inline void dumpExtraValues(ProfileTStates_ &stat) {
            for(int i = 0; i < ProfileData_<>::_proExtraCount; ++i) {
                dumpTab();
                dumpValue(stat, SEQAN_PROINDEXCOUNT - SEQAN_PROEXTRACOUNT + i);
            }
    }
    
        inline void dumpMark() {
            if (!mark.empty()) {
                dumpTab();
                fprintf(out, "\"%s\"", mark.c_str());
                mark.erase();
            }
        }

        inline void dump(ProfileTStates_ &stat) {
            setTime(ProfileData_<>::_proValue);
            dumpNext += dumpStep;
            bol = true;
            bool _flush = ((dumpStep == 0.0)) || ((lines & 16) == 0);

            dumpSysValues(stat);
            dumpExtraValues(stat);
            dumpMark();
            dumpEndl();
            if (_flush) fflush(out);
            ++lines;
        }

        inline void signalDumpTest(_proFloat now) {
            if (dumpStep > 0 && now > dumpNext && running) {
                dump(last);
                maximize(all, last);
                sync(last);
            }
        }

        inline void signalNewMax(int valNum) {
            if (running)
                if (last[valNum] < ProfileData_<>::_proValue[valNum])
                    last[valNum] = ProfileData_<>::_proValue[valNum];
        }

        inline void setMark(const char *text) {
            if (running) {
                mark = text;
                if (dumpStep == 0.0) {
                    dump(last);                 // manual dump;
                    maximize(all, last);
                    sync(last);
                }
            }
        }
        
        inline void reset() {
            syncTime(last);
        }

        inline void setEndMark(const char *text) {
            if (running) {
                setMark(text);
                reset();
            }
        }

    private:
        
        bool bol;   // begin of line
    };



/*
    inline void _profileSignalDumpTest(_proFloat now);
    inline void _profileSignalNewMax(int valNum);
    inline void _profileMark(const char *text);
    inline void _profileEndMark(const char *text);
    inline void _profileReset();

    inline void _profileSet(int valNum, _proFloat value);
    inline void _profileAdd(int valNum, _proFloat value);
    inline void _profileSub(int valNum, _proFloat value);
    
    // simple interface for external programs
    inline void *_profileMalloc(size_t size);
    inline void _profileFree(void *_ptr);
*/

    inline void _profileSignalDumpTest(_proFloat now) {
//IOREV _notio_
        if (ProfileData_<>::_proPFileStream) ProfileData_<>::_proPFileStream->signalDumpTest(now);
    }

    inline void _profileSignalNewMax(int valNum) {
//IOREV _notio_
        if (((ProfileValueType_[valNum] & SEQAN_PROSTATE) != 0)) {
            if (ProfileData_<>::_proPFileStream) ProfileData_<>::_proPFileStream->signalNewMax(valNum);
            if (ProfileData_<>::_proPFile)       ProfileData_<>::_proPFile->signalNewMax(valNum);
        }
    }

    inline void _profileMark(const char *text) {
//IOREV _notio_
        if (ProfileData_<>::_proPFileStream) ProfileData_<>::_proPFileStream->setMark(text);
        if (ProfileData_<>::_proPFile)       ProfileData_<>::_proPFile->setMark(text);
    }

    inline void _profileEndMark(const char *text) {
//IOREV _notio_
        if (ProfileData_<>::_proPFileStream) { ProfileData_<>::_proPFileStream->setEndMark(text); }
        if (ProfileData_<>::_proPFile)       { ProfileData_<>::_proPFile->setEndMark(text); }
    }

    inline void _profileReset() {
//IOREV _notio_
        if (ProfileData_<>::_proPFileStream) { ProfileData_<>::_proPFileStream->reset(); }
        if (ProfileData_<>::_proPFile)       { ProfileData_<>::_proPFile->reset(); }
    }




    template <typename TValue>
    inline void _profileSet(ProfileValueIndex_ valNum, TValue value) {
//IOREV _notio_
        _proFloat now = sysTime();
        ProfileData_<>::_proLastUpdate[valNum] = now;
        if (ProfileData_<>::_proValue[valNum] < value) {
            ProfileData_<>::_proValue[valNum] = value;
            _profileSignalNewMax(valNum);
        } else
            ProfileData_<>::_proValue[valNum] = value;
        _profileSignalDumpTest(now);
    }

    template <typename TValue>
    inline void _profileAdd(ProfileValueIndex_ valNum, TValue value) {
//IOREV _notio_
        _proFloat now = sysTime();
        ProfileData_<>::_proValue[valNum] += value;
        ProfileData_<>::_proLastUpdate[valNum] = now;
        if (valNum == SEQAN_PROIO) _profileAdd(SEQAN_PROIORANDOM, 1);
        _profileSignalNewMax(valNum);
        _profileSignalDumpTest(now);
    }

    template <typename TValue>
    inline void _profileSub(ProfileValueIndex_ valNum, TValue value) {
//IOREV _notio_
        _proFloat now = sysTime();
        ProfileData_<>::_proValue[valNum] -= value;
        ProfileData_<>::_proLastUpdate[valNum] = now;
        _profileSignalDumpTest(now);
    }
    
    // simple interface for external programs
    inline void *_profileMalloc(size_t size) {
//IOREV _notio_
        size_t *ptr = reinterpret_cast<size_t*>(malloc(size + sizeof(size_t)));
        if (ptr) {
            _profileAdd(SEQAN_PROMEMORY, (_proFloat)(*ptr = size));
//          printf("_profileMalloc %x size %d\n", ptr, size);
            ++ptr;
        }
        return ptr;
    }

    inline void _profileFree(void *_ptr) {
//IOREV _notio_
        size_t *ptr = reinterpret_cast<size_t*>(_ptr);
        if (ptr) {
            --ptr;
//          printf("_profileFree   %x size %d\n", _ptr, *ptr);
            _profileSub(SEQAN_PROMEMORY, (_proFloat)*ptr);
        }
        free(ptr);
    }

    inline _proFloat _profileUpdate(_proFloat& a) {
//IOREV _notio_
        _proFloat x = sysTime() - a;
        a += x;
        return x;
    }
}

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_PROFILING_H_
