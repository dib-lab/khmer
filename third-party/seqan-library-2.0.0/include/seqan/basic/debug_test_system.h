// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
// Copyright (c) 2013 NVIDIA Corporation
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// The SeqAn testing infrastructure.  Based on ideas from the OpenMS
// "ClassTest.h".
// ==========================================================================

// TODO(holtgrew): This could use some cleanup.

// SEQAN_NO_GENERATED_FORWARDS

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_DEBUG_TEST_SYSTEM_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_DEBUG_TEST_SYSTEM_H_

#include <iostream>  // stdout, stderr
#include <iomanip>
#include <cstring>   // strrpos
#include <cstdlib>   // exit()
#include <cstdio>
#include <cstdarg>   // va_start, va_list, va_end
#include <algorithm> // min()
#include <set>
#include <vector>
#include <string>
#include <typeinfo>

#ifdef PLATFORM_WINDOWS
#include <Windows.h>    // DeleteFile()
#else  // #ifdef PLATFORM_WINDOWS
#include <unistd.h>     // unlink()
#include <sys/stat.h>   // mkdir()
#include <dirent.h>     // DIR
#if SEQAN_HAS_EXECINFO
#include <execinfo.h>   // backtrace(), backtrace_symbols()
#endif  // #if SEQAN_HAS_EXECINFO
#include <cxxabi.h>     // __cxa_demangle()
#include <signal.h>
#endif  // #ifdef PLATFORM_WINDOWS

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Demangler
// ----------------------------------------------------------------------------
// Holds the name of a given C++ type T.
// NOTE(esiragusa): this class could become a subclass of CStyle String...

namespace seqan {

template <typename T>
struct Demangler
{
#ifdef PLATFORM_GCC
    char *data_begin;
#else
    const char *data_begin;
#endif

    Demangler()
    {
        T t;
        _demangle(*this, t);
    }

    Demangler(T const & t)
    {
        _demangle(*this, t);
    }

    ~Demangler()
    {
#ifdef PLATFORM_GCC
        free(data_begin);
#endif
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _demangle(Demangler)
// ----------------------------------------------------------------------------

template <typename T>
inline void _demangle(Demangler<T> & me, T const & t)
{
#ifdef PLATFORM_GCC
    int status;
    me.data_begin = abi::__cxa_demangle(typeid(t).name(), NULL, NULL, &status);
#else
    me.data_begin = typeid(t).name();
#endif
}

// ----------------------------------------------------------------------------
// Function toCString(Demangler)
// ----------------------------------------------------------------------------

template <typename T>
inline const char * toCString(Demangler<T> const & me)
{

    return me.data_begin;
}

}

/*!
 * @defgroup AssertMacros Assertion and Check Macros
 * @brief The assertion and check macros provided by SeqAn.
 *
 * Assertions are checks performed at runtime when debugging is enabled.  Debugging is enabled by defining the
 * preprocessor symbol <tt>SEQAN_ENABLE_DEBUG</tt> as <tt>1</tt> (the default is to set it to <tt>0</tt> if the common C
 * macro <tt>NDEBUG</tt> is defined and to set it to <tt>1</tt> otherwise.  When using the SeqAn build system or the
 * CMake FindSeqAn.cmake module, this is automatically set appropriately.
 *
 * The SEQAN_CHECK and SEQAN_FAIL macro always lead to an exit of the program with a non-0 return value.
 */

/*!
 * @macro AssertMacros#SEQAN_FAIL
 * @headerfile <seqan/basic.h>
 * @brief Force abortion of program, regardless of debugging settings.
 *
 * @signature SEQAN_FAIL(msg[, args]);
 *
 * @param[in] msg  A format string.
 * @param[in] args An optional list of arguments that are used for filling msg.
 *
 * @section Remarks
 *
 * Use this if something really unexpected happens inside your functions and there is no way to report this through the
 * API.  A good example would be logic errors, e.g. invalid values.
 *
 * @section Examples
 *
 * In the following example, the <tt>SEQAN_FAIL</tt> is there if a possible value is added to <tt>MyEnum</tt> but the
 * function <tt>foo</tt> is not updated accordingly.
 *
 * @code{.cpp}
 * enum MyEnum
 * {
 *   VALUE_ONE,
 *   VALUE_TWO
 * };
 *
 * bool foo(MyEnum x)
 * {
 *     switch (x)
 *     {
 *     case VALUE_ONE:
 *         // do something
 *         return true;
 *     case VALUE_TWO:
 *         // do something
 *         return true;
 *     }
 *
 *     SEQAN_FAIL("Logic error. Should never reach here. x == %d.", x);
 *     return false;
 * }
 * @endcode
 */

#define SEQAN_FAIL(...)                                                 \
    do {                                                                \
        ::seqan::ClassTest::forceFail(__FILE__, __LINE__,               \
                                      __VA_ARGS__);                     \
        ::seqan::ClassTest::fail();                                     \
    } while (false)

/*!
 * @macro AssertMacros#SEQAN_CHECK
 * @headerfile <seqan/basic.h>
 * @brief Force abortion of program if a condition is not met, regardless of debugging settings.
 *
 * @signature SEQAN_CHECK(condition, msg[, args]);
 *
 * @param[in] condition An expression that is checked.
 * @param[in] msg       A format string.
 * @param[in] args      An optional list of arguments.
 *
 * @section Remarks
 *
 * Use this if something really unexpected happens inside your functions and there is no way to report this through the
 * API.  A good example would be logic errors, e.g. invalid values.
 *
 * @section Examples
 *
 * In the following example, the <tt>SEQAN_CHECK</tt> stops program execution if a value is added to <tt>MyEnum</tt> but
 * the function <tt>foo</tt> is not updated accordingly.
 *
 * @code{.cpp}
 * enum MyEnum
 * {
 *   VALUE_ONE,
 *   VALUE_TWO
 * };
 *
 * bool foo(MyEnum x)
 * {
 *     SEQAN_CHECK((x == VALUE_ONE || x == VALUE_TWO), "Invalid value for x == %d.", x);
 *
 *     switch (x)
 *     {
 *     case VALUE_ONE:
 *         // do something
 *         return true;
 *     case VALUE_TWO:
 *         // do something
 *         return true;
 *     }
 *
 *     return false;  // Should never reach here, checked above with SEQAN_CHECK.
 * }
 * @endcode
 */

#define SEQAN_CHECK(_arg1, ...)                                         \
    do {                                                                \
        if (!::seqan::ClassTest::testTrue(__FILE__, __LINE__,           \
                                          (_arg1), # _arg1,              \
                                          __VA_ARGS__)) {               \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)

// SeqAn's has three global debug/testing levels: testing, debug and
// release.  Depending on the level, the SEQAN_ASSERT_* and
// SEQAN_CHECKPOINT macros will be enabled.
//
// Note that this is independent of the <cassert> assertions and
// NDEBUG being defined.
//
// The levels are enabled by the values of the macros
// SEQAN_ENABLE_TESTING and SEQAN_ENABLE_DEBUG.  By setting a macro to
// 0, one disables the level and by setting the macro to 1, one
// enables a level.  Enabling testing also enables debug, overriding a
// value of 0 for SEQAN_ENABLE_DEBUG.
//
// If the level is release (both the macros for debug and testing are
// 0), the assertions will be disabled.  If the level is debug then
// the assertions will be enabled.  If the level is testing then the
// checkpoint macros will also be enabled.
//
// The default is to enable debugging but disable testing.
//
// You can print the current level using the function seqan::printDebugLevel().

/*!
 * @macro TestSystemMacros#SEQAN_ENABLE_TESTING
 * @headerfile <seqan/basic.h>
 * @brief Indicates whether testing is enabled.
 *
 * @signature SEQAN_ENABLE_TESTING
 *
 * When set to 1, testing is enabled.  If it is undefined or set to 0, testing is disabled.  This means the macros for
 * the tests (SEQAN_BEGIN_TESTSUITE, SEQAN_DEFINE_TEST, SEQAN_CALL_TEST, and SEQAN_END_TESTSUITE) will be enabled.  This
 * makes failing assertions raise exceptions instead of calling <tt>abort()</tt> (which terminates the program).
 *
 * By default, this is set to 0.
 *
 * If you want to change this value in your C++ program code you have to define this value before including any SeqAn header!
 *
 * If set to 1 then @link TestSystemMacros#SEQAN_ENABLE_DEBUG @endlink is forced to 1 as well.
 *
 * @see TestSystemMacros#SEQAN_ENABLE_DEBUG
 */

// Set default for SEQAN_ENABLE_TESTING.
#ifndef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 0
#endif  // #ifndef SEQAN_ENABLE_TESTING

/*!
 * @macro TestSystemMacros#SEQAN_ENABLE_DEBUG
 * @headerfile <seqan/basic.h>
 * @brief Indicates whether debugging is enabled.
 *
 * @signature SEQAN_ENABLE_DEBUG
 *
 * When enabled (set to 1) then debugging is enabled.  This means the assertion macros are expanded to actual test code.
 * If debugging (and testing) is disabled then the SeqAn assertion macros expand to no instructions.
 *
 * By default, thi sis set to 0 if <tt>NDEBUG</tt> is defined and set to 1 if <tt>NDEBUG</tt> is not defined.
 *
 * If you want to change this value then you have to define this value before including any SeqAn header.
 *
 * Force-enabled if SEQAN_ENABLE_TESTING is set to 1.
 *
 * @see TestSystemMacros#SEQAN_ENABLE_TESTING
 */

// Set default for SEQAN_ENABLE_DEBUG.
#ifndef SEQAN_ENABLE_DEBUG
  #ifdef NDEBUG
    #define SEQAN_ENABLE_DEBUG 0
  #else  // #ifdef NDEBUG
    #define SEQAN_ENABLE_DEBUG 1
  #endif  // #ifdef NDEBUG
#endif  // #ifndef SEQAN_ENABLE_DEBUG

// Force-enable debugging if testing is enabled.
#if SEQAN_ENABLE_TESTING
#undef SEQAN_ENABLE_DEBUG
#define SEQAN_ENABLE_DEBUG 1
#endif  // #if SEQAN_ENABLE_TESTING

// Allow disabling checkpoints independent of testing.
#ifndef SEQAN_ENABLE_CHECKPOINTS
#define SEQAN_ENABLE_CHECKPOINTS 0 // SEQAN_ENABLE_TESTING
#endif  // #ifndef SEQAN_ENABLE_CHECKPOINTS

/*!
 * @macro TestSystemMacros#SEQAN_TYPEDEF_FOR_DEBUG
 * @headerfile <seqan/basic.h>
 * @brief When using typedefs that are only used in debug mode then they have to be marked with macro.
 *
 * @signature SEQAN_TYPEDE_FOR_DEBUG
 *
 * @section Examples
 *
 * @code{.cpp}
 * typedef int TInt SEQAN_TYPEDEF_FOR_DEBUG;
 * @endcode
 */

#if !SEQAN_ENABLE_DEBUG
#  if defined(__GNUC__) && ((__GNUC__ > 4) || ((__GNUC__ == 4) && (__GNUC_MINOR__ >= 7)))
#    define SEQAN_TYPEDEF_FOR_DEBUG __attribute__((unused))
#  else
#    define SEQAN_TYPEDEF_FOR_DEBUG
#  endif
#else
#  define SEQAN_TYPEDEF_FOR_DEBUG
#endif

// TODO(holtgrew): This one is for profiling and in tests.
#if defined(__GNUC__) && ((__GNUC__ > 4) || ((__GNUC__ == 4) && (__GNUC_MINOR__ >= 7)))
#  define SEQAN_UNUSED_TYPEDEF __attribute__((unused))
#else
#  define SEQAN_UNUSED_TYPEDEF
#endif

namespace seqan {

// SEQAN_CXX_FLAGS_ contains the compiler flags, SEQAN_CXX_FLAGS is a string
// literal with this value.
#if !defined(SEQAN_CXX_FLAGS_)
#define SEQAN_CXX_FLAGS_ SEQAN_CXX_FLAGS_NOT_SET
#endif //  !defined(SEQAN_CXX_FLAGS__)
#define SEQAN_MKSTRING_(str) # str
#define SEQAN_MKSTRING(str) SEQAN_MKSTRING_(str)
#define SEQAN_CXX_FLAGS SEQAN_MKSTRING(SEQAN_CXX_FLAGS_)
//#undef SEQAN_MKSTRING
//#undef SEQAN_MKSTRING_

/*!
 * @fn printDebugLevel
 * @headerfile <seqan/basic.h>
 * @brief Print the current SeqAn debug level and the compiler flags to the given stream.
 *
 * @signature void printDebugLevel(stream);
 *
 * @param[in,out] stream A std::ostream where the information about the levels are streamed to.
 */

template <typename TStream>
void printDebugLevel(TStream & stream)
{
    stream << "SEQAN_ENABLE_DEBUG == " << SEQAN_ENABLE_DEBUG << std::endl;
    stream << "SEQAN_ENABLE_TESTING == " << SEQAN_ENABLE_TESTING << std::endl;
    stream << "SEQAN_ENABLE_CHECKPOINTS == " << SEQAN_ENABLE_CHECKPOINTS << std::endl;
    stream << "SEQAN_CXX_FLAGS == \"" << SEQAN_CXX_FLAGS << "\"" << std::endl;
}

#if defined(PLATFORM_WINDOWS) || !SEQAN_HAS_EXECINFO

template <typename TSize>
void printStackTrace(TSize /*maxFrames*/)
{}

#else

// print a demangled stack backtrace of the caller function
// TODO(esiragusa): use Demangler.
template <typename TSize>
void printStackTrace(TSize maxFrames)
{
    void * addrlist[256];
    char temp[4096];
    char addr[20];
    char offset[20];

    size_t size;
    int status;
    char * symname;
    char * demangled;

    std::cerr << std::endl << "stack trace:" << std::endl;

    int addrlist_len = backtrace(addrlist, maxFrames);
    char ** symbollist = backtrace_symbols(addrlist, addrlist_len);
    for (int i = 1; i < addrlist_len; ++i)
    {
        offset[0] = 0;
        addr[0] = 0;
        demangled = NULL;

        // LINUX FORMAT:
        //          ./sam2svg [0x473b8c]
        //          /lib/libc.so.6 [0x7f40d2526f60]
        //          ./sam2svg(_Z2f3v+0x10) [0x47200c]
        //          ./sam2svg(_Z2f2v+0xd) [0x472021]
        //          ./sam2svg(main+0x1367) [0x4735fc]
        //          /lib/libc.so.6(__libc_start_main+0xe6) [0x7f40d25131a6]
        //

        if (3 == sscanf(symbollist[i], "%*[^(](%4095[^+]+%[^)]) %s", temp, offset, addr))
        {
            symname = temp;
            if (NULL != (demangled = abi::__cxa_demangle(temp, NULL, &size, &status)))
            {
                symname = demangled;
            }
        }
        // MAC OS X FORMAT:
        //          1   sam2svg                             0x0000000100003a39 _ZN5seqanL28signalHandlerPrintStackTraceEi + 21
        //          2   libSystem.B.dylib                   0x00007fff87a6d67a _sigtramp + 26
        //          3   libSystem.B.dylib                   0x00007fff87a76df7 tiny_free_do_recirc_to_depot + 980
        //          4   sam2svg                             0x00000001000021b9 _Z2f2v + 9
        //          5   sam2svg                             0x00000001000034b1 main + 4546
        //          6   sam2svg                             0x0000000100002190 start + 52
        else if (3 == sscanf(symbollist[i], "%*d %*s %s %s %*s %s", addr, temp, offset))
        {
            symname = temp;
            if (NULL != (demangled = abi::__cxa_demangle(temp, NULL, &size, &status)))
            {
                symname = demangled;
            }
        }
        // LINUX FORMAT:
        //          ./sam2svg [0x473b8c]
        //          /lib/libc.so.6 [0x7f40d2526f60]
        else if (2 == sscanf(symbollist[i], "%s %s", temp, addr))
        {
            symname = temp;
        }
        // DEFAULT:
        else
        {
            symname = symbollist[i];
        }

        std::cerr << std::setw(3) << i - 1;
        std::cerr << std::setw(20) << addr;
        std::cerr << "  " << symname;
        if (offset[0] != 0)
            std::cerr << " + " << offset;
        std::cerr << std::endl;

        free(demangled);
    }
    std::cerr << std::endl;
    // Only the array must be freed according to man page, not the contents.
    free(symbollist);
}

static void signalHandlerPrintStackTrace(int signum)
{
    std::cerr << std::endl;
    printStackTrace(20);
    signal(signum, SIG_DFL);
    kill(getpid(), signum);
}

inline int _deploySignalHandlers()
{
    signal(SIGSEGV, signalHandlerPrintStackTrace);      // segfault
    signal(SIGFPE, signalHandlerPrintStackTrace);       // divide by zero
    // ...
    return 0;
}

#if SEQAN_ENABLE_DEBUG

// automatically deploy signal handlers that output the stack trace on a trap (in debug mode)

template <typename T>
struct SignalHandlersDummy_
{
    static const int i;
};

template <typename T>
const int SignalHandlersDummy_<T>::i = _deploySignalHandlers();

namespace {
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#endif  // ifdef __clang__
volatile int signalHandlersDummy_ = SignalHandlersDummy_<void>::i;
#ifdef __clang__
#pragma clang diagnostic pop
#endif  // ifdef __clang__
}

#endif  // #if SEQAN_ENABLE_DEBUG
#endif  // defined(PLATFORM_WINDOWS) || !SEQAN_HAS_EXECINFO


// Namespace for the testing infrastructure.
//
// This namespace contains the variables and functions that are used
// in the macros below to perform the tests.
namespace ClassTest {
// Raised when an assertion fails in test mode.
struct AssertionFailedException {};

// Container for static global data for the tests.
struct StaticData
{
    // Number of tests that were run.
    static int & testCount()
    {
        static int result = 0;
        return result;
    }

    // Number of errors that occurred.
    static int & errorCount()
    {
        static int result = 0;
        return result;
    }

    // Number of skipped tests.
    static int & skippedCount()
    {
        static int result = 0;
        return result;
    }

    // Flag whether there was an error in this test.
    static bool & thisTestOk()
    {
        static bool result = 0;
        return result;
    }

    // Flag whether this test was skipped.
    static bool & thisTestSkipped()
    {
        static bool result = 0;
        return result;
    }

    // Name of the current test.
    static const char * & currentTestName()
    {
        const char * defaultValue = "";
        static const char * result = const_cast<char *>(defaultValue);
        return result;
    }

    // Base path to the binary.  Extrapolated from __FILE__.
    static char * & basePath()
    {
        const char * defaultValue = ".";
        static char * result = const_cast<char *>(defaultValue);
        return result;
    }

    static char const * _computePathToRoot()
    {
        // Get path to include.
        const char * file = __FILE__;
        int pos = -1;
        for (size_t i = 0; i < strlen(file) - strlen("include"); ++i)
        {
            if (strncmp(file + i, "include", strlen("include")) == 0)
            {
                pos = i;
            }
        }
        for (; pos > 0 && *(file + pos - 1) != '/' &&  *(file + pos - 1) != '\\'; --pos)
            continue;
        if (pos == -1)
        {
            std::cerr << "Could not extrapolate path to repository from __FILE__ == \""
                      << __FILE__ << "\"" << std::endl;
            exit(1);
        }

        static char buffer[1024];
        strncpy(&buffer[0], file, pos);
        buffer[pos - 1] = '\0';
        return &buffer[0];
    }

    // Base path to the directory containing "core" and "extras."
    // Extrapolated from __FILE__.
    static char const * pathToRoot()
    {
        const char * result = 0;
        if (!result)
            result = _computePathToRoot();
        return result;
    }

    // Total number of checkpoints in header file.
    static int & totalCheckPointCount()
    {
        static int result = 0;
        return result;
    }

    // Total number of checkpoints found in binary files.
    static int & foundCheckPointCount()
    {
        static int result = 0;
        return result;
    }

    // Names of temporary files as returned by tempFileName.  This
    // global state is used to remove any existing such files
    // after completing the testsuite.
    static::std::vector<std::string> & tempFileNames()
    {
        static::std::vector<std::string> filenames;
        return filenames;
    }
};

// Open a temporary file, unlink it, return posix handle.  Note: This has not been tested yet.
// TODO(holtgrew): Not used yet and Windows code does not work.
/*
inline
int openTempFile() {
#ifdef PLATFORM_WINDOWS
    char * fileName = _tempnam(NULL, "SQN");
    if (!fileName) {
        ::std::cerr << "Cannot create a unique temporary filename" << ::std::endl;
        exit(1);
    }
    int result = open(fileName, _O_RDWR | OPEN_TEMPORARY);
    free(fileName);
    return result;
#else  // A Unix...
    char filenameBuffer[100];
    strcpy(filenameBuffer, "/tmp/SEQANXXXXXXXXXX");
    int result = mkstemp(filenameBuffer);
    unlink(filenameBuffer);
    return result;
#endif  // ifdef PLATFORM_WINDOWS
}
*/

// Return the path to a temporary file, in a static buffer in this
// function.  This is not thread safe!
inline
const char * tempFileName()
{
    static char fileNameBuffer[1000];
#ifdef PLATFORM_WINDOWS
    static char filePathBuffer[1000];
    //  Gets the temp path env string (no guarantee it's a valid path).
    DWORD dwRetVal = 0;
    dwRetVal = GetTempPath(1000,            // length of the buffer
                           filePathBuffer); // buffer for path
    if (dwRetVal > 1000 || (dwRetVal == 0))
    {
        std::cerr << "GetTempPath failed" << std::endl;
        exit(1);
    }

    UINT uRetVal   = 0;
    uRetVal = GetTempFileName(filePathBuffer,   // directory for tmp files
                              TEXT("SEQAN."),   // temp file name prefix
                              0,                // create unique name
                              fileNameBuffer);  // buffer for name

    if (uRetVal == 0)
    {
        std::cerr << "GetTempFileName failed" << std::endl;
        exit(1);
    }

    DeleteFile(fileNameBuffer);
    CreateDirectoryA(fileNameBuffer, NULL);
    StaticData::tempFileNames().push_back(fileNameBuffer);
    strcat(fileNameBuffer, "\\test_file");
    return fileNameBuffer;

#else  // ifdef PLATFORM_WINDOWS_VS
    strcpy(fileNameBuffer, "/tmp/SEQAN.XXXXXXXXXXXXXXXXXXXX");
    mode_t cur_umask = umask(S_IRWXO | S_IRWXG);  // to silence Coverity warning
    int _tmp = mkstemp(fileNameBuffer);
    (void) _tmp;
    umask(cur_umask);
    unlink(fileNameBuffer);
    mkdir(fileNameBuffer, 0777);

    StaticData::tempFileNames().push_back(fileNameBuffer);

    strcat(fileNameBuffer, "/test_file");
    return fileNameBuffer;

#endif  // ifdef PLATFORM_WINDOWS
}

// Initialize the testing infrastructure.
//
// Used through SEQAN_BEGIN_TESTSUITE(test_name)
inline
void beginTestSuite(const char * testSuiteName, const char * argv0)
{
    // First things first: Print test suite name and current debug level.
    std::cout << "TEST SUITE " << testSuiteName << std::endl;
    printDebugLevel(std::cout);
    (void)testSuiteName;
    StaticData::testCount() = 0;
    StaticData::skippedCount() = 0;
    StaticData::errorCount() = 0;
    StaticData::totalCheckPointCount() = 0;
    StaticData::foundCheckPointCount() = 0;
    // Get path to argv0.
    const char * end = argv0;
    const char * ptr = std::min(strchr(argv0, '\\'), strchr(argv0, '/'));     // On Windows, we can have both \ and /.
    for (; ptr != 0; ptr = std::min(strchr(ptr + 1, '\\'), strchr(ptr + 1, '/')))
        end = ptr;
    int rpos = end - argv0;
    if (rpos <= 0)
    {
        StaticData::basePath() = new char[2];
        strcpy(StaticData::basePath(), ".");
    }
    else
    {
        int len = rpos;
        StaticData::basePath() = new char[len];
        strncpy(StaticData::basePath(), argv0, len);
    }

#ifdef PLATFORM_WINDOWS_VS
    // Set CRT reporting such that everything goes to stderr and there are
    // no popups causing timeouts.
    _set_error_mode(_OUT_TO_STDERR);
    _CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
    _CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDERR);
    _CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_FILE);
    _CrtSetReportFile(_CRT_ERROR, _CRTDBG_FILE_STDERR);
    _CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);
    _CrtSetReportFile(_CRT_ASSERT, _CRTDBG_FILE_STDERR);
#endif  // PLATFORM_WINDOWS_VS
}

// Run test suite finalization.
//
// Used through SEQAN_END_TESTSUITE
//
// Prints a bottom banner with the error count and returns the
// program's return code.
inline
int endTestSuite()
{
    delete[] StaticData::basePath();

    std::cout << "**************************************" << std::endl;
    std::cout << " Total Check Points : " << StaticData::totalCheckPointCount() << std::endl;
    std::cout << " Found Check Points : " << StaticData::foundCheckPointCount() << std::endl;
    std::cout << " Lost Check Points  : " << StaticData::totalCheckPointCount() - StaticData::foundCheckPointCount() << std::endl;
    std::cout << "--------------------------------------" << std::endl;
    std::cout << " Total Tests: " << StaticData::testCount() << std::endl;
    std::cout << " Skipped:     " << StaticData::skippedCount() << std::endl;
    std::cout << " Errors:      " << StaticData::errorCount() << std::endl;
    std::cout << "**************************************" << std::endl;
    // TODO(holtgrew): Re-enable that all check points have to be found for the test to return 1;
    /*
    if (StaticData::totalCheckPointCount() != StaticData::foundCheckPointCount())
        return 1;
    */
    // Delete all temporary files that still exist.
    for (unsigned i = 0; i < StaticData::tempFileNames().size(); ++i)
    {
#ifdef PLATFORM_WINDOWS
        HANDLE hFind;
        WIN32_FIND_DATA data;

        std::string temp = StaticData::tempFileNames()[i].c_str() + std::string("\\*");
        hFind = FindFirstFile(temp.c_str(), &data);
        if (hFind != INVALID_HANDLE_VALUE)
        {
            do
            {
                std::string tempp = StaticData::tempFileNames()[i].c_str() + std::string("\\") + data.cFileName;
                if (strcmp(data.cFileName, ".") == 0 || strcmp(data.cFileName, "..") == 0)
                    continue;  // Skip these.
                if (!DeleteFile(tempp.c_str()))
                    std::cerr << "WARNING: Could not delete file " << tempp << "\n";
            }
            while (FindNextFile(hFind, &data));
            FindClose(hFind);
        }

        if (!RemoveDirectory(StaticData::tempFileNames()[i].c_str()))
            std::cerr << "WARNING: Could not delete directory " << StaticData::tempFileNames()[i] << "\n";
#else  // #ifdef PLATFORM_WINDOWS
        DIR * dpdf;
        struct dirent * epdf;

        dpdf = opendir(StaticData::tempFileNames()[i].c_str());
        if (dpdf != NULL)
        {
            while ((epdf = readdir(dpdf)) != NULL)
            {
                std::string temp = StaticData::tempFileNames()[i].c_str() + std::string("/") + std::string(epdf->d_name);
                unlink(temp.c_str());
            }
        }

        rmdir(StaticData::tempFileNames()[i].c_str());
        if (closedir(dpdf) != 0)
            std::cerr << "WARNING: Could not delete directory " << StaticData::tempFileNames()[i] << "\n";
#endif  // #ifdef PLATFORM_WINDOWS
    }

    if (StaticData::errorCount() != 0)
        return 1;

    return 0;
}

// Run test initialization.
inline
void beginTest(const char * testName)
{
    StaticData::currentTestName() = testName;
    StaticData::thisTestOk() = true;
    StaticData::thisTestSkipped() = false;
    StaticData::testCount() += 1;
}

// Run test finalization.
inline
void endTest()
{
    if (StaticData::thisTestSkipped())
    {
        std::cout << StaticData::currentTestName() << " SKIPPED" << std::endl;
    }
    else if (StaticData::thisTestOk())
    {
        std::cout << StaticData::currentTestName() << " OK" << std::endl;
    }
    else
    {
        std::cerr << StaticData::currentTestName() << " FAILED" << std::endl;
    }
}

// Marks the current test as "skipped".
inline
void skipCurrentTest()
{
    StaticData::thisTestSkipped() = true;
    StaticData::skippedCount() += 1;
}

// Called by the macro SEQAN_ASSERT_FAIL.
inline void forceFail(const char * file, int line,
                      const char * comment, ...)
{
    StaticData::errorCount() += 1;
    std::cerr << file << ":" << line << " FAILED! ";
    if (comment)
    {
        std::cerr << " (";
        va_list args;
        va_start(args, comment);
        vfprintf(stderr, comment, args);
        va_end(args);
        std::cerr << ")";
    }
    std::cerr << std::endl;
}

// Similar to forceFail above, but accepting a va_list parameter.
inline void vforceFail(const char * file, int line,
                       const char * comment, va_list argp)
{
    StaticData::errorCount() += 1;
    std::cerr << file << ":" << line << " FAILED! ";
    if (comment)
    {
        std::cerr << " (";
        vfprintf(stderr, comment, argp);
        std::cerr << ")";
    }
    std::cerr << std::endl;
}

// Same as forceFail above, but with comment set to 0.
inline void forceFail(const char * file, int line)
{
    forceFail(file, line, 0);
}

// Called by the macro SEQAN_ASSERT_EQ.
//
// Tests that the given two value are equal.  Returns true iff the
// two values are equal.
template <typename T1, typename T2>
bool testEqual(char const * file, int line,
               T1 const & value1, char const * expression1,
               T2 const & value2, char const * expression2,
               char const * comment, ...)
{
    if (!(value1 == value2))
    {
        // Increase global error count.
        StaticData::thisTestOk() = false;
        StaticData::errorCount() += 1;
        // Print assertion failure text, with comment if any is given.
        std::cerr << file << ":" << line << " Assertion failed : "
                  << expression1 << " == " << expression2 << " was: " << value1
                  << " != " << value2;
        if (comment)
        {
            std::cerr << " (";
            va_list args;
            va_start(args, comment);
            vfprintf(stderr, comment, args);
            va_end(args);
            std::cerr << ")";
        }
        std::cerr << std::endl;
        return false;
    }
    return true;
}

// Similar to testEqual above, but accepts a va_list instead of variadic
// parameters.
template <typename T1, typename T2>
bool vtestEqual(const char * file, int line,
                const T1 & value1, const char * expression1,
                const T2 & value2, const char * expression2,
                const char * comment, va_list argp)
{
    if (!(value1 == value2))
    {
        // Increase global error count.
        StaticData::thisTestOk() = false;
        StaticData::errorCount() += 1;
        // Print assertion failure text, with comment if any is given.
        std::cerr << file << ":" << line << " Assertion failed : "
                  << expression1 << " == " << expression2 << " was: " << value1
                  << " != " << value2;
        if (comment)
        {
            std::cerr << " (";
            vfprintf(stderr, comment, argp);
            std::cerr << ")";
        }
        std::cerr << std::endl;
        return false;
    }
    return true;
}

// Same as testEqual above, but with comment set to 0.
template <typename T1, typename T2>
bool testEqual(const char * file, int line,
               const T1 & value1, const char * expression1,
               const T2 & value2, const char * expression2)
{
    return testEqual(file, line, value1, expression1, value2, expression2, 0);
}

// Called by the macro SEQAN_ASSERT_IN_DELTA.
//
// Tests that the given two value are equal.  Returns true iff the
// two values are equal.
template <typename T1, typename T2, typename T3>
bool testInDelta(const char * file, int line,
                 const T1 & value1, const char * expression1,
                 const T2 & value2, const char * expression2,
                 const T3 & value3, const char * expression3,
                 const char * comment, ...)
{
    if (!(value1 >= value2 - value3 && value1 <= value2 + value3))
    {
        // Increase global error count.
        StaticData::thisTestOk() = false;
        StaticData::errorCount() += 1;
        // Print assertion failure text, with comment if any is given.
        std::cerr << file << ":" << line << " Assertion failed : "
                  << expression1 << " in [" << expression2 << " - " << expression3
                  << ", " << expression2 << " + " << expression3 << "] was: " << value1
                  << " not in [" << value2 - value3 << ", " << value2 + value3 << "]";
        if (comment)
        {
            std::cerr << " (";
            va_list args;
            va_start(args, comment);
            vfprintf(stderr, comment, args);
            va_end(args);
            std::cerr << ")";
        }
        std::cerr << std::endl;
        return false;
    }
    return true;
}

// Similar to testInDelta above, but accepts a va_list instead of variadic
// parameters.
template <typename T1, typename T2, typename T3>
bool vtestInDelta(const char * file, int line,
                  const T1 & value1, const char * expression1,
                  const T2 & value2, const char * expression2,
                  const T3 & value3, const char * expression3,
                  const char * comment, va_list argp)
{
    if (!(value1 >= value2 - value3 && value1 <= value2 + value3))
    {
        // Increase global error count.
        StaticData::thisTestOk() = false;
        StaticData::errorCount() += 1;
        // Print assertion failure text, with comment if any is given.
        std::cerr << file << ":" << line << " Assertion failed : "
                  << expression1 << " in [" << expression2 << " - " << expression3
                  << ", " << expression2 << " + " << expression3 << "] was: " << value1
                  << " not in [" << value2 - value3 << ", " << value2 + value3 << "]";
        if (comment)
        {
            std::cerr << " (";
            vfprintf(stderr, comment, argp);
            std::cerr << ")";
        }
        std::cerr << std::endl;
        return false;
    }
    return true;
}

// Same as testInDelta above, but with comment set to 0.
template <typename T1, typename T2, typename T3>
bool testInDelta(const char * file, int line,
                 const T1 & value1, const char * expression1,
                 const T2 & value2, const char * expression2,
                 const T3 & value3, const char * expression3)
{
    return testInDelta(file, line, value1, expression1, value2, expression2, value3, expression3, 0);
}

// Called by the macro SEQAN_ASSERT_NEQ.
//
// Tests that the given two value are not equal.  Returns true iff
// the two values are equal.
template <typename T1, typename T2>
bool testNotEqual(const char * file, int line,
                  const T1 & value1, const char * expression1,
                  const T2 & value2, const char * expression2,
                  const char * comment, ...)
{
    if (!(value1 != value2))
    {
        // Increase global error count.
        StaticData::thisTestOk() = false;
        StaticData::errorCount() += 1;
        // Print assertion failure text, with comment if any is given.
        std::cerr << file << ":" << line << " Assertion failed : "
                  << expression1 << " != " << expression2 << " was: " << value1
                  << " == " << value2;
        if (comment)
        {
            std::cerr << " (";
            va_list args;
            va_start(args, comment);
            vfprintf(stderr, comment, args);
            va_end(args);
            std::cerr << ")";
        }
        std::cerr << std::endl;
        return false;
    }
    return true;
}

// Similar to testNotEqual above, but accepts a va_list instead of variadic
// parameters.
template <typename T1, typename T2>
bool vtestNotEqual(const char * file, int line,
                   const T1 & value1, const char * expression1,
                   const T2 & value2, const char * expression2,
                   const char * comment, va_list argp)
{
    if (!(value1 != value2))
    {
        // Increase global error count.
        StaticData::thisTestOk() = false;
        StaticData::errorCount() += 1;
        // Print assertion failure text, with comment if any is given.
        std::cerr << file << ":" << line << " Assertion failed : "
                  << expression1 << " != " << expression2 << " was: " << value1
                  << " == " << value2;
        if (comment)
        {
            std::cerr << " (";
            vfprintf(stderr, comment, argp);
            std::cerr << ")";
        }
        std::cerr << std::endl;
        return false;
    }
    return true;
}

// Same as testNotEqual above, but with comment set to 0.
template <typename T1, typename T2>
bool testNotEqual(const char * file, int line,
                  const T1 & value1, const char * expression1,
                  const T2 & value2, const char * expression2)
{
    return testNotEqual(file, line, value1, expression1, value2, expression2, 0);
}

// Called by the macro SEQAN_ASSERT_GEQ.
//
// Tests that the first value is greater than or equal to the
// second one.  Returns true iff the test yields true.
template <typename T1, typename T2>
bool testGeq(const char * file, int line,
             const T1 & value1, const char * expression1,
             const T2 & value2, const char * expression2,
             const char * comment, ...)
{
    if (!(value1 >= value2))
    {
        // Increase global error count.
        StaticData::thisTestOk() = false;
        StaticData::errorCount() += 1;
        // Print assertion failure text, with comment if any is given.
        std::cerr << file << ":" << line << " Assertion failed : "
                  << expression1 << " >= " << expression2 << " was: " << value1
                  << " < " << value2;
        if (comment)
        {
            std::cerr << " (";
            va_list args;
            va_start(args, comment);
            vfprintf(stderr, comment, args);
            va_end(args);
            std::cerr << ")";
        }
        std::cerr << std::endl;
        return false;
    }
    return true;
}

// Similar to testGeq above, but accepts a va_list instead of variadic
// parameters.
template <typename T1, typename T2>
bool vtestGeq(const char * file, int line,
              const T1 & value1, const char * expression1,
              const T2 & value2, const char * expression2,
              const char * comment, va_list argp)
{
    if (!(value1 >= value2))
    {
        // Increase global error count.
        StaticData::thisTestOk() = false;
        StaticData::errorCount() += 1;
        // Print assertion failure text, with comment if any is given.
        std::cerr << file << ":" << line << " Assertion failed : "
                  << expression1 << " >= " << expression2 << " was: " << value1
                  << " < " << value2;
        if (comment)
        {
            std::cerr << " (";
            vfprintf(stderr, comment, argp);
            std::cerr << ")";
        }
        std::cerr << std::endl;
        return false;
    }
    return true;
}

// Same as testGeq above, but with comment set to 0.
template <typename T1, typename T2>
bool testGeq(const char * file, int line,
             const T1 & value1, const char * expression1,
             const T2 & value2, const char * expression2)
{
    return testGeq(file, line, value1, expression1, value2, expression2, 0);
}

// Called by the macro SEQAN_ASSERT_GT.
//
// Tests that the first value is greater than the second one.
// Returns true iff the test yields true.
template <typename T1, typename T2>
bool testGt(const char * file, int line,
            const T1 & value1, const char * expression1,
            const T2 & value2, const char * expression2,
            const char * comment, ...)
{
    if (!(value1 > value2))
    {
        // Increase global error count.
        StaticData::thisTestOk() = false;
        StaticData::errorCount() += 1;
        // Print assertion failure text, with comment if any is given.
        std::cerr << file << ":" << line << " Assertion failed : "
                  << expression1 << " > " << expression2 << " was: " << value1
                  << " <= " << value2;
        if (comment)
        {
            std::cerr << " (";
            va_list args;
            va_start(args, comment);
            vfprintf(stderr, comment, args);
            va_end(args);
            std::cerr << ")";
        }
        std::cerr << std::endl;
        return false;
    }
    return true;
}

// Similar to testGt above, but accepts a va_list instead of variadic
// parameters.
template <typename T1, typename T2>
bool vtestGt(const char * file, int line,
             const T1 & value1, const char * expression1,
             const T2 & value2, const char * expression2,
             const char * comment, va_list argp)
{
    if (!(value1 > value2))
    {
        // Increase global error count.
        StaticData::thisTestOk() = false;
        StaticData::errorCount() += 1;
        // Print assertion failure text, with comment if any is given.
        std::cerr << file << ":" << line << " Assertion failed : "
                  << expression1 << " > " << expression2 << " was: " << value1
                  << " <= " << value2;
        if (comment)
        {
            std::cerr << " (";
            vfprintf(stderr, comment, argp);
            std::cerr << ")";
        }
        std::cerr << std::endl;
        return false;
    }
    return true;
}

// Same as testGt above, but with comment set to 0.
template <typename T1, typename T2>
bool testGt(const char * file, int line,
            const T1 & value1, const char * expression1,
            const T2 & value2, const char * expression2)
{
    return testGt(file, line, value1, expression1, value2, expression2, 0);
}

// Called by the macro SEQAN_ASSERT_LEQ.
//
// Tests that the first value is less than or equal to the second
// one.  Returns true iff the test yields true.
template <typename T1, typename T2>
bool testLeq(const char * file, int line,
             const T1 & value1, const char * expression1,
             const T2 & value2, const char * expression2,
             const char * comment, ...)
{
    if (!(value1 <= value2))
    {
        // Increase global error count.
        StaticData::thisTestOk() = false;
        StaticData::errorCount() += 1;
        // Print assertion failure text, with comment if any is given.
        std::cerr << file << ":" << line << " Assertion failed : "
                  << expression1 << " <= " << expression2 << " was: " << value1
                  << " > " << value2;
        if (comment)
        {
            std::cerr << " (";
            va_list args;
            va_start(args, comment);
            vfprintf(stderr, comment, args);
            va_end(args);
            std::cerr << ")";
        }
        std::cerr << std::endl;
        return false;
    }
    return true;
}

// Similar to testLeq above, but accepts a va_list instead of variadic
// parameters.
template <typename T1, typename T2>
bool vtestLeq(const char * file, int line,
              const T1 & value1, const char * expression1,
              const T2 & value2, const char * expression2,
              const char * comment, va_list argp)
{
    if (!(value1 <= value2))
    {
        // Increase global error count.
        StaticData::thisTestOk() = false;
        StaticData::errorCount() += 1;
        // Print assertion failure text, with comment if any is given.
        std::cerr << file << ":" << line << " Assertion failed : "
                  << expression1 << " <= " << expression2 << " was: " << value1
                  << " > " << value2;
        if (comment)
        {
            std::cerr << " (";
            vfprintf(stderr, comment, argp);
            std::cerr << ")";
        }
        std::cerr << std::endl;
        return false;
    }
    return true;
}

// Same as testLeq above, but with comment set to 0.
template <typename T1, typename T2>
bool testLeq(const char * file, int line,
             const T1 & value1, const char * expression1,
             const T2 & value2, const char * expression2)
{
    return testLeq(file, line, value1, expression1, value2, expression2, 0);
}

// Called by the macro SEQAN_ASSERT_LT.
//
// Tests that the first value is greater than the second one.
// Returns true iff the test yields true.
template <typename T1, typename T2>
bool testLt(const char * file, int line,
            const T1 & value1, const char * expression1,
            const T2 & value2, const char * expression2,
            const char * comment, ...)
{
    if (!(value1 < value2))
    {
        // Increase global error count.
        StaticData::thisTestOk() = false;
        StaticData::errorCount() += 1;
        // Print assertion failure text, with comment if any is given.
        std::cerr << file << ":" << line << " Assertion failed : "
                  << expression1 << " < " << expression2 << " was: " << value1
                  << " >= " << value2;
        if (comment)
        {
            std::cerr << " (";
            va_list args;
            va_start(args, comment);
            vfprintf(stderr, comment, args);
            va_end(args);
            std::cerr << ")";
        }
        std::cerr << std::endl;
        return false;
    }
    return true;
}

// Similar to testLt above, but accepts a va_list instead of variadic
// parameters.
template <typename T1, typename T2>
bool vtestLt(const char * file, int line,
             const T1 & value1, const char * expression1,
             const T2 & value2, const char * expression2,
             const char * comment, va_list argp)
{
    if (!(value1 < value2))
    {
        // Increase global error count.
        StaticData::thisTestOk() = false;
        StaticData::errorCount() += 1;
        // Print assertion failure text, with comment if any is given.
        std::cerr << file << ":" << line << " Assertion failed : "
                  << expression1 << " < " << expression2 << " was: " << value1
                  << " >= " << value2;
        if (comment)
        {
            std::cerr << " (";
            vfprintf(stderr, comment, argp);
            std::cerr << ")";
        }
        std::cerr << std::endl;
        return false;
    }
    return true;
}

// Same as testLt above, but comment is 0.
template <typename T1, typename T2>
bool testLt(const char * file, int line,
            const T1 & value1, const char * expression1,
            const T2 & value2, const char * expression2)
{
    return testLt(file, line, value1, expression1, value2, expression2, 0);
}

// Called by the macro SEQAN_ASSERT.
//
// Test that the given argument evaluates to true.
template <typename T>
bool testTrue(const char * file, int line,
              const T & value_, const char * expression_,
              const char * comment, ...)
{
    if (!(value_))
    {
        // Increase global error count.
        StaticData::thisTestOk() = false;
        StaticData::errorCount() += 1;
        // Print assertion failure text, with comment if any is given.
        std::cerr << file << ":" << line << " Assertion failed : "
                  << expression_ << " should be true but was " << (value_);
        if (comment)
        {
            std::cerr << " (";
            va_list args;
            va_start(args, comment);
            vfprintf(stderr, comment, args);
            va_end(args);
            std::cerr << ")";
        }
        std::cerr << std::endl;
        return false;
    }
    return true;
}

// Similar to testTrue above, but accepts a va_list instead of variadic
// parameters.
template <typename T>
bool vtestTrue(const char * file, int line,
               const T & value_, const char * expression_,
               const char * comment, va_list argp)
{
    if (!(value_))
    {
        // Increase global error count.
        StaticData::thisTestOk() = false;
        StaticData::errorCount() += 1;
        // Print assertion failure text, with comment if any is given.
        std::cerr << file << ":" << line << " Assertion failed : "
                  << expression_ << " should be true but was " << (value_);
        if (comment)
        {
            std::cerr << " (";
            vfprintf(stderr, comment, argp);
            std::cerr << ")";
        }
        std::cerr << std::endl;
        return false;
    }
    return true;
}

// Same as testTrue above, but comment will automatically be set to 0.
template <typename T>
bool testTrue(const char * file, int line,
              const T & value_, const char * expression_)
{
    return testTrue(file, line, value_, expression_, 0);
}

// Called by the macro SEQAN_ASSERT.
//
// Test that the given argument evaluates to false.
template <typename T>
bool testFalse(const char * file, int line,
               const T & value_, const char * expression_,
               const char * comment, ...)
{
    if (value_)
    {
        // Increase global error count.
        StaticData::thisTestOk() = false;
        StaticData::errorCount() += 1;
        // Print assertion failure text, with comment if any is given.
        std::cerr << file << ":" << line << " Assertion failed : "
                  << expression_ << " should be false but was " << (value_);
        if (comment)
        {
            std::cerr << " (";
            va_list args;
            va_start(args, comment);
            vfprintf(stderr, comment, args);
            va_end(args);
            std::cerr << ")";
        }
        std::cerr << std::endl;
        return false;
    }
    return true;
}

// Similar to testFalse above, but accepts a va_list instead of variadic
// parameters.
template <typename T>
bool vtestFalse(const char * file, int line,
                const T & value_, const char * expression_,
                const char * comment, va_list argp)
{
    if (value_)
    {
        // Increase global error count.
        StaticData::thisTestOk() = false;
        StaticData::errorCount() += 1;
        // Print assertion failure text, with comment if any is given.
        std::cerr << file << ":" << line << " Assertion failed : "
                  << expression_ << " should be false but was " << (value_);
        if (comment)
        {
            std::cerr << " (";
            vfprintf(stderr, comment, argp);
            std::cerr << ")";
        }
        std::cerr << std::endl;
        return false;
    }
    return true;
}

// Same as testFalse above, but comment will automatically be set to 0.
template <typename T>
bool testFalse(const char * file, int line,
               const T & value_, const char * expression_)
{
    return testFalse(file, line, value_, expression_, 0);
}

// Represents a check point in a file.
struct CheckPoint
{
    // Path to the file.
    const char * file;
    // Line in the file.
    unsigned int line;

    // Less-than comparator for check points.
    bool operator<(const CheckPoint & other) const
    {
        int c = strcmp(file, other.file);
        if (c < 0)
            return true;

        if (c == 0 && line < other.line)
            return true;

        return false;
    }

};

// Wrapper for a set of check points.
// TODO(holtgrew): Simply store the set?
struct CheckPointStore
{
    static::std::set<CheckPoint> & data()
    {
        static::std::set<CheckPoint> result;
        return result;
    }
};

// Puts the given check point into the CheckPointStore's data.
inline bool
registerCheckPoint(unsigned int line, const char * file)
{
    const char * file_name = strrchr(file, '/');
    const char * file_name_2 = strrchr(file, '\\');
    if (file_name_2 > file_name)
        file_name = file_name_2;
    if (!file_name)
        file_name = file;
    else
        ++file_name;

    CheckPoint cp = {file_name, line};
        #ifdef _OMP
        #pragma omp critical
        #endif  // #ifdef _OMP
    CheckPointStore::data().insert(cp);
    return true;
}

// Test whether the given check point exists in the check point
// store.
inline void
testCheckPoint(const char * file, unsigned int line)
{
    StaticData::totalCheckPointCount() += 1;
    CheckPoint cp = {file, line};
    if (CheckPointStore::data().find(cp) == CheckPointStore::data().end())
    {
        std::cerr << file << ":" << line << "  -- Check point lost."
                  << std::endl;
        return;
    }
    StaticData::foundCheckPointCount() += 1;
}

// Verify the check points for the given file.
inline void
verifyCheckPoints(const char * file)
{
    char const * file_name = strrchr(file, '/');
    char const * file_name_2 = strrchr(file, '\\');
    if (file_name_2 > file_name)
        file_name = file_name_2;
    if (!file_name)
        file_name = file;
    else
        ++file_name;



    int len = strlen(StaticData::pathToRoot()) +
              strlen("/") + strlen(file) + 1;
    char * absolutePath = new char[len];
    absolutePath[0] = '\0';
    strcat(absolutePath, StaticData::pathToRoot());
    strcat(absolutePath, "/");
    strcat(absolutePath, file);

    FILE * fl = ::std::fopen(absolutePath, "r");
    delete[] absolutePath;
    if (!fl)
    {
        std::cerr << file << " -- verifyCheckPoints could not find this file." << std::endl;
    }
    unsigned int line_number = 1;
    char buf[1 << 16];

    while (::std::fgets(buf, sizeof(buf), fl))
    {
        if (::std::strstr(buf, "SEQAN_CHECKPOINT"))
        {
            testCheckPoint(file_name, line_number);
        }
        ++line_number;
    }

    ::std::fclose(fl);
}

#if SEQAN_ENABLE_TESTING
// If in testing mode then raise an AssertionFailedException.
inline void fail()
{
    StaticData::thisTestOk() = false;
    printStackTrace(20);
    throw AssertionFailedException();
}

#else
// If not in testing mode then quit with an abort.
inline void fail()
{
    printStackTrace(20);
    abort();
}

#endif  // #if SEQAN_ENABLE_TESTING

}  // namespace ClassTest

/*!
 * @macro TestSystemMacros#SEQAN_DEFINE_TEST
 * @headerfile <seqan/basic.h>
 * @brief Expand to test definition.
 *
 * @signature SEQAN_DEFINE_TEST(test_name)
 *
 * This macro expands to the definition of a $void$ function with <tt>SEQAN_TEST_ + test_name</tt> as its name.
 *
 * @section Example
 *
 * @code{.cpp}
 * SEQAN_DEFINE_TEST(test_name)
 * {
 *     SEQAN_ASSERT_LT(0, 3);
 * }
 * @endcode
 */

// This macro expands to function header for one test.
#define SEQAN_DEFINE_TEST(test_name)                    \
    template <bool speed_up_dummy_to_prevent_compilation_of_unused_tests_> \
    void SEQAN_TEST_ ## test_name()

/*!
 * @defgroup TestSystemMacros Test System Macros
 * @brief Macros for the test system.
 */

/*!
 * @macro TestSystemMacros#SEQAN_BEGIN_TESTSUITE
 * @headerfile <seqan/basic.h>
 * @brief Expand to a test suite beginning.
 *
 * @signature SEQAN_BEGIN_TESTSUITE(name)
 *
 * @param[in] name The name of the test suite.
 *
 * This macro expands to a <tt>main()</tt> function and some initialization code that sets up the test system.
 *
 * @section Examples
 *
 * @code{.cpp}
 * #include <seqan/basic.h>
 *
 * SEQAN_BEGIN_TESTSUITE(test_foo)
 * {
 *    SEQAN_CALL_TEST(test_foo_my_test);
 * }
 * SEQAN_END_TESTSUITE
 * @endcode
 */

#if SEQAN_ENABLE_TESTING
// This macro expands to startup code for a test file.
#define SEQAN_BEGIN_TESTSUITE(suite_name)                       \
    int main(int argc, char ** argv) {                           \
        (void) argc;                                                \
        ::seqan::ClassTest::beginTestSuite(# suite_name, argv[0]);

/*!
 * @macro TestSystemMacros#SEQAN_END_TESTSUITE
 * @headerfile <seqan/basic.h>
 * @brief Expand to test suite ending.
 *
 * @signature SEQAN_END_TESTSUITE
 *
 * This macro expands to finalization code for a test suite.
 *
 * @section Examples
 *
 * @code{.cpp}
 * #include <seqan/basic.h>
 *
 * SEQAN_BEGIN_TESTSUITE(test_foo)
 * {
 *     SEQAN_CALL_TEST(test_foo_my_test);
 * }
 * SEQAN_END_TESTSUITE
 * @endcode
 */

// This macro expands to shutdown code for a test file.
#define SEQAN_END_TESTSUITE                     \
    return ::seqan::ClassTest::endTestSuite();  \
    }

/*!
 * @macro TestSystemMacros#SEQAN_CALL_TEST
 * @headerfile <seqan/basic.h>
 * @brief Expand to calling a test.
 *
 * @signature SEQAN_CALL_TEST(test_name);
 *
 * This expects the test to be defined with SEQAN_DEFINE_TEST.  This macro will expand to code that calls the code
 * inside a try/catch block. Use this macro within a test suite, only.
 *
 * @section Examples
 *
 * @code{.cpp}
 * // Within a test suite.
 * SEQAN_CALL_TEST(test_name);
 * @endcode
 */

// This macro expands to code to call a given test.
#define SEQAN_CALL_TEST(test_name)                                      \
    do {                                                                \
        seqan::ClassTest::beginTest(# test_name);                       \
        try {                                                           \
            SEQAN_TEST_ ## test_name<true>();                           \
        } catch (seqan::ClassTest::AssertionFailedException e) {        \
            /* Swallow exception, go on with next test. */              \
            (void) e;  /* Get rid of unused variable warning. */        \
        } catch (std::exception const & e) {                            \
            std::cerr << "Unexpected exception of type "                \
                      << toCString(seqan::Demangler<std::exception>(e)) \
                      << "; message: " << e.what() << "\n";             \
            seqan::ClassTest::StaticData::thisTestOk() = false;         \
            seqan::ClassTest::StaticData::errorCount() += 1;            \
        } catch (...) {                                                 \
            std::cerr << "Unexpected exception of unknown type\n";      \
            seqan::ClassTest::StaticData::thisTestOk() = false;         \
            seqan::ClassTest::StaticData::errorCount() += 1;            \
        }                                                               \
        seqan::ClassTest::endTest();                                    \
    } while (false)

/*!
 * @macro TestSystemMacros#SEQAN_SKIP_TEST
 * @headerfile <seqan/basic.h>
 * @brief Force the test to return without failing and mark it as skipped.
 *
 * @signature SEQAN_SKIP_TEST;
 *
 * @section Examples
 *
 * @code{.cpp}
 * SEQAN_DEFINE_TEST(test_skipped)
 * {
 *     SEQAN_SKIP_TEST;
 * }
 * @endcode
 */

// This macro returns from the current function and logs a "skipped"
// event for the current test.
#define SEQAN_SKIP_TEST                         \
    do {                                        \
        ::seqan::ClassTest::skipCurrentTest();  \
        return;                                 \
    } while (false)
#endif  // #if SEQAN_ENABLE_TESTING

// variadic macros are not supported by VS 2003 and before
#if !defined(_MSC_VER) || (_MSC_VER >= 1400)

#if SEQAN_ENABLE_DEBUG && !defined(__CUDA_ARCH__)

/*!
 * @macro AssertMacros#SEQAN_ASSERT
 * @headerfile <seqan/basic.h>
 * @brief Test that the given expression can be coerced to <tt>true</tt>.
 *
 * @signature SEQAN_ASSERT(expression);
 * @signature SEQAN_ASSERT_MSG(expression, message[, parameters]);
 *
 * @param[in] expression An expression to check for being true.
 * @param[in] message    A format string.
 * @param[in] parameters An optional list of parameters.
 *
 * @section Remarks
 *
 * The main advantage of this macro is that it prints the values of its argument on failures.  Note that the
 * <tt>operator&lt;&lt;</tt> to the type of <tt>std::cerr</tt> has to be defined for the type of both expression
 * parameters. Otherwise, simply use the equivalent SEQAN_ASSERT @call.
 *
 * See SEQAN_CHECK and SEQAN_FAIL for (conditionally) aborting your program regardless of debug settings.
 *
 * @section Examples
 *
 * @code{.cpp}
 * SEQAN_ASSERT(0);  // will fail
 * SEQAN_ASSERT(1);  // will run through
 * SEQAN_ASSERT_MSG(0, "message %d", 2);  // Will fail with message.
 * @endcode
 */

/*!
 * @macro AssertMacros#SEQAN_ASSERT_NOT
 * @headerfile <seqan/basic.h>
 * @brief Test that the given expression can be coerced to <tt>false</tt>.
 *
 * @signature SEQAN_ASSERT_NOT(expression)
 * @signature SEQAN_ASSERT_NOT_MSG(expression, message[, parameters])
 *
 * @param[in] expression An expression to check for being false.
 * @param[in] message    A format string.
 * @param[in] parameters An optional list of parameters.
 *
 * @section Remarks
 *
 * The main advantage of this macro is that it prints the values of its argument on failures.  Note that the
 * <tt>operator&lt;&lt;</tt> to the type of <tt>std::cerr</tt> has to be defined for the type of both expression
 * parameters.  Otherwise, simply use the equivalent SEQAN_ASSERT call.
 *
 * See SEQAN_CHECK and SEQAN_FAIL for (conditionally) aborting your program regardless of debug settings.
 *
 * @section Examples
 *
 * @code{.cpp}
 * SEQAN_ASSERT_NOT(0);  // will run through
 * SEQAN_ASSERT_NOT(1);  // will fail
 * SEQAN_ASSERT_NOT_MSG(0, "msg %s", "test");  // will fail with message
 * @endcode
 */

/*!
 * @macro AssertMacros#SEQAN_ASSERT_EQ
 * @headerfile <seqan/basic.h>
 * @brief Test that two given expressions are equal, as defined by the matching call to the <tt>operator=(,)</tt>.

 * @signature SEQAN_ASSERT_EQ(expression1, expression2);
 * @signature SEQAN_ASSERT_EQ_MSG(expression1, expression2, comment[, parameters]);
 *
 * @param[in] expression1 The first expression.
 * @param[in] expression2 The second expression.
 * @param[in] comment     A C-string (<tt>char const *</tt>) to use as a format string for printing a message
 *                        on failure.
 * @param[in] parameters  An optional parameter that is put into <tt>printf()</tt> with format string
 *                        <tt>comment</tt>.
 *
 * The main advantage of this macro is that it prints the values of its argument on failures.  Note that the
 * <tt>operator&lt;&lt;</tt> to the type of <tt>std::cerr</tt> has to be defined for the type of both expression
 * parameters.  Otherwise, simply use the equivalent SEQAN_ASSERT call.
 *
 * See SEQAN_CHECK and SEQAN_FAIL for (conditionally) aborting your program regardless of debug settings.
 *
 * @section Examples
 *
 * @code{.cpp}
 * SEQAN_ASSERT_EQ(0, false);  // will run through
 * SEQAN_ASSERT_EQ(1, false);  // will fail
 * SEQAN_ASSERT_EQ(1, "foo");  // will not compile
 * SEQAN_ASSERT_EQ_MSG(1, false, "msg");  // will fail with message
 * @endcode
 */

/*!
 * @macro AssertMacros#SEQAN_ASSERT_NEQ
 * @headerfile <seqan/basic.h>
 * @brief Test that two given expressions are not equal, as defined by the matching call to the <tt>operator!=(,)</tt>.
 *
 * @signature SEQAN_ASSERT_NEQ(expression1, expression2);
 * @signature SEQAN_ASSERT_NEQ_MSG(expression1, expression2, comment[, parameters]);
 *
 * @param[in] expression1 The first expression.
 * @param[in] expression2 The second expression.
 * @param[in] comment     A C-string (<tt>char const *</tt>) to use as a format string for printing a message
 *                        on failure.
 * @param[in] parameters  An optional parameter that is put into <tt>printf()</tt> with format string
 *                        <tt>comment</tt>.
 *
 * The main advantage of this macro is that it prints the values of its argument on failures.  Note that the
 * <tt>operator&lt;&lt;</tt> to the type of <tt>std::cerr</tt> has to be defined for the type of both expression
 * parameters.  Otherwise, simply use the equivalent SEQAN_ASSERT call.
 *
 * See SEQAN_CHECK and SEQAN_FAIL for (conditionally) aborting your program regardless of debug settings.
 *
 * @section Examples
 *
 * @code{.cpp}
 * SEQAN_ASSERT_NEQ(0, false);  // will fail
 * SEQAN_ASSERT_NEQ(1, false);  // will run through
 * SEQAN_ASSERT_NEQ(1, "foo");  // will not compile
 * SEQAN_ASSERT_NEQ_MSG(1, false, "msg");  // will fail with message
 * @endcode
 */

/*!
 * @macro AssertMacros#SEQAN_ASSERT_LT
 * @headerfile <seqan/basic.h>
 * @brief Test that the two given expressions are in the less-than relation as defined by the matching call to
 *        operator<(,).
 *
 * @signature SEQAN_ASSERT_LT(expression1, expression2);
 * @signature SEQAN_ASSERT_LT(expression1, expression2, comment[, parameters]);
 *
 * @param[in] expression1 The first expression.
 * @param[in] expression2 The second expression.
 * @param[in] comment     A C-string (<tt>char const *</tt>) to use as a format string for printing a message
 *                        on failure.
 * @param[in] parameters  An optional parameter that is put into <tt>printf()</tt> with format string
 *                        <tt>comment</tt>.
 *
 * The main advantage of this macro is that it prints the values of its argument on failures.  Note that the
 * <tt>operator&lt;&lt;</tt> to the type of <tt>std::cerr</tt> has to be defined for the type of both expression
 * parameters.  Otherwise, simply use the equivalent SEQAN_ASSERT call.
 *
 * See SEQAN_CHECK and SEQAN_FAIL for (conditionally) aborting your program regardless of debug settings.
 *
 * @section Examples
 *
 * @code{.cpp}
 * SEQAN_ASSERT_LT(0, 1);  // will run through
 * SEQAN_ASSERT_LT(1, 1);  // will not run through
 * SEQAN_ASSERT_LT_MSG(1, 1, "msg");  // will fail with message
 * @endcode
 */

/*!
 * @macro AssertMacros#SEQAN_ASSERT_LEQ
 *
 * @brief Test that the two given expressions are in the less-than-or-equal
 *        relation as defined by the matching call to operator<=(,).
 *
 * @signature SEQAN_ASSERT_LEQ(expression1, expression2)
 * @signature SEQAN_ASSERT_LEQ_MSG(expression1, expression2, comment[,
 *            parameters])
 *
 * @param[in] expression1 The first expression.
 * @param[in] expression2 The second expression.
 * @param[in] comment     A C-string (<tt>char const *</tt>) to use as a format string for printing a message
 *                        on failure.
 * @param[in] parameters  An optional parameter that is put into <tt>printf()</tt> with format string
 *                        <tt>comment</tt>.
 *
 * The main advantage of this macro is that it prints the values of its argument
 * on failures. Note that the <tt>operator&lt;&lt;</tt> to the type of
 * <tt>std::cerr</tt> has to be defined for the type of both expression
 * parameters. Otherwise, simply use the equivalent  SEQAN_ASSERT
 * call.
 *
 * See  SEQAN_CHECK  and  SEQAN_FAIL  for
 * (conditionally) aborting your program regardless of debug settings.
 *
 * @section Examples
 *
 * @code{.cpp}
 * SEQAN_ASSERT_LEQ(1, 1);  // will run through
 * SEQAN_ASSERT_LEQ(1, 2);  // will not run through
 * SEQAN_ASSERT_LEQ_MSG(1, 2, "msg");  // will fail with message
 * @endcode
 */

/*!
 * @macro AssertMacros#SEQAN_ASSERT_GT
 *
 * @brief Test that the two given expressions are in the greather-than relation
 *        as defined by the matching call to operator>(,).
 *
 * @signature SEQAN_ASSERT_GT(expression1, expression2);
 * @signature SEQAN_ASSERT_GT_MSG(expression1, expression2, comment[, parameters]);
 *
 * @param[in] expression1 The first expression.
 * @param[in] expression2 The second expression.
 * @param[in] comment     A C-string (<tt>char const *</tt>) to use as a format string for printing a message
 *                        on failure.
 * @param[in] parameters  An optional parameter that is put into <tt>printf()</tt> with format string
 *                        <tt>comment</tt>.
 *
 * The main advantage of this macro is that it prints the values of its argument
 * on failures. Note that the <tt>operator&lt;&lt;</tt> to the type of
 * <tt>std::cerr</tt> has to be defined for the type of both expression
 * parameters. Otherwise, simply use the equivalent  SEQAN_ASSERT
 * call.
 *
 * See  SEQAN_CHECK  and  SEQAN_FAIL  for
 * (conditionally) aborting your program regardless of debug settings.
 *
 * @section Examples
 *
 * @code{.cpp}
 * SEQAN_ASSERT_GT(2, 1);  // will run through
 * SEQAN_ASSERT_GT(1, 1);  // will not run through
 * SEQAN_ASSERT_GT_MSG(1, 1, "msg");  // will fail with message
 * @endcode
 */

/*!
 * @macro AssertMacros#SEQAN_ASSERT_GEQ
 *
 * @brief Test that the two given expressions are in the greater-than-or-equal
 *        relation as defined by the matching call to operator>=(,).
 *
 * @signature SEQAN_ASSERT_GEQ(expression1, expression2);
 * @signature SEQAN_ASSERT_GEQ_MSG(expression1, expression2, comment[, parameters]);
 *
 * @param[in] expression1 The first expression.
 * @param[in] expression2 The second expression.
 * @param[in] comment     A C-string (<tt>char const *</tt>) to use as a format string for printing a message
 *                        on failure.
 * @param[in] parameters  An optional parameter that is put into <tt>printf()</tt> with format string
 *                        <tt>comment</tt>.
 *
 * The main advantage of this macro is that it prints the values of its argument on failures.  Note that the
 * <tt>operator&lt;&lt;</tt> to the type of <tt>std::cerr</tt> has to be defined for the type of both expression
 * parameters.  Otherwise, simply use the equivalent SEQAN_ASSERT call.
 *
 * See SEQAN_CHECK and SEQAN_FAIL for (conditionally) aborting your program regardless of debug settings.
 *
 * @section Examples
 *
 * @code{.cpp}
 * SEQAN_ASSERT_GEQ(1, 1);  // will run through
 * SEQAN_ASSERT_GEQ(0, 1);  // will not run through
 * SEQAN_ASSERT_GEQ_MSG(0, 1, "msg");  // will fail with message
 * @endcode
 */

/*!
 * @macro AssertMacros#SEQAN_ASSERT_IN_DELTA
 *
 * @brief Test that a value <tt>y</tt> lies within an <tt>delta</tt> environment of a value <tt>x</tt>.
 *
 * @signature SEQAN_ASSERT_IN_DELTA(x, y, delta);
 * @signature SEQAN_ASSERT_IN_DELTA_MSG(x, y, delta, comment[, parameters]);
 *
 * @param[in] x           The value to center the environment in.
 * @param[in] y           The value to check whether it falls within the environment.
 * @param[in] delta       The environment size.
 * @param[in] comment     A C-string (<tt>char const *</tt>) to use as a format string for printing a message
 *                        on failure.
 * @param[in] parameters  An optional parameter that is put into <tt>printf()</tt> with format string
 *                        <tt>comment</tt>.
 *
 * The main advantage of this macro is that it prints the values of its argument on failures.  Note that the
 * <tt>operator&lt;&lt;</tt> to the type of <tt>std::cerr</tt> has to be defined for the type of both expression
 * parameters.  Otherwise, simply use the equivalent SEQAN_ASSERT call.
 *
 * See SEQAN_CHECK and SEQAN_FAIL for (conditionally) aborting your program regardless of debug settings.
 *
 * @section Examples
 *
 * @code{.cpp}
 * SEQAN_ASSERT_IN_DELTA(0, 0, 0.1);  // will run through
 * SEQAN_ASSERT_IN_DELTA(1, -2, 1);  // will fail
 * SEQAN_ASSERT_IN_DELTA(1, "foo");  // will not compile
 * SEQAN_ASSERT_IN_DELTA_MSG(1, 0, 0.1, "msg");  // will fail with message
 * @endcode
 */

// Force a test failure.
//
// Usage:  SEQAN_ASSERT_FAIL("Failed at position %d", pos);
#define SEQAN_ASSERT_FAIL(...)                                          \
    do {                                                                \
        ::seqan::ClassTest::forceFail(__FILE__, __LINE__,               \
                                      __VA_ARGS__);                     \
        ::seqan::ClassTest::fail();                                     \
    } while (false)


// Equality assertion without a comment.
//
// Usage:  SEQAN_ASSERT_EQ(4, 4);
#define SEQAN_ASSERT_EQ(_arg1, _arg2)                                   \
    do {                                                                \
        if (!::seqan::ClassTest::testEqual(__FILE__, __LINE__,          \
                                           (_arg1), # _arg1,             \
                                           (_arg2), # _arg2)) {          \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// Equality assertion with a comment.
//
// Usage:  SEQAN_ASSERT_EQ(4, 4);
#define SEQAN_ASSERT_EQ_MSG(_arg1, _arg2, ...)                          \
    do {                                                                \
        if (!::seqan::ClassTest::testEqual(__FILE__, __LINE__,          \
                                           (_arg1), # _arg1,             \
                                           (_arg2), # _arg2,             \
                                           __VA_ARGS__)) {              \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// In-delta-environment assertion without a comment.
//
// Usage:  SEQAN_ASSERT_IN_DELTA(4.1, 4, 0.1);
#define SEQAN_ASSERT_IN_DELTA(_arg1, _arg2, _arg3)                      \
    do {                                                                \
        if (!::seqan::ClassTest::testInDelta(__FILE__, __LINE__,        \
                                             (_arg1), # _arg1,           \
                                             (_arg2), # _arg2,           \
                                             (_arg3), # _arg3)) {        \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// In-delta-environment assertion witha comment.
//
// Usage:  SEQAN_ASSERT_IN_DELTA_MSG(4.1, 4, 0.1, "3.9 <= 4.1 <= 4.1");
#define SEQAN_ASSERT_IN_DELTA_MSG(_arg1, _arg2, _arg3, ...)             \
    do {                                                                \
        if (!::seqan::ClassTest::testInDelta(__FILE__, __LINE__,        \
                                             (_arg1), # _arg1,           \
                                             (_arg2), # _arg2,           \
                                             (_arg3), # _arg3,           \
                                             __VA_ARGS__)) {            \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// Inequality assertion without a comment.
//
// Usage:  SEQAN_ASSERT_NEQ(4, 5);
#define SEQAN_ASSERT_NEQ(_arg1, _arg2)                                  \
    do {                                                                \
        if (!::seqan::ClassTest::testNotEqual(__FILE__, __LINE__,       \
                                              (_arg1), # _arg1,          \
                                              (_arg2), # _arg2)) {       \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// Inequality assertion with a comment.
//
// Usage:  SEQAN_ASSERT_NEQ(4, 5);
#define SEQAN_ASSERT_NEQ_MSG(_arg1, _arg2, ...)                         \
    do {                                                                \
        if (!::seqan::ClassTest::testNotEqual(__FILE__, __LINE__,       \
                                              (_arg1), # _arg1,          \
                                              (_arg2), # _arg2,          \
                                              __VA_ARGS__)) {           \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// Less-than-or-equal assertion without a comment.
#define SEQAN_ASSERT_LEQ(_arg1, _arg2)                                  \
    do {                                                                \
        if (!::seqan::ClassTest::testLeq(__FILE__, __LINE__,            \
                                         (_arg1), # _arg1,               \
                                         (_arg2), # _arg2)) {            \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// Less-than-or-equal assertion with a comment.
#define SEQAN_ASSERT_LEQ_MSG(_arg1, _arg2, ...)                         \
    do {                                                                \
        if (!::seqan::ClassTest::testLeq(__FILE__, __LINE__,            \
                                         (_arg1), # _arg1,               \
                                         (_arg2), # _arg2,               \
                                         __VA_ARGS__)) {                \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// Less-than assertion without a comment.
#define SEQAN_ASSERT_LT(_arg1, _arg2)                                   \
    do {                                                                \
        if (!::seqan::ClassTest::testLt(__FILE__, __LINE__,             \
                                        (_arg1), # _arg1,                \
                                        (_arg2), # _arg2)) {             \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// Less-than assertion with a comment.
#define SEQAN_ASSERT_LT_MSG(_arg1, _arg2, ...)                          \
    do {                                                                \
        if (!::seqan::ClassTest::testLt(__FILE__, __LINE__,             \
                                        (_arg1), # _arg1,                \
                                        (_arg2), # _arg2,                \
                                        __VA_ARGS__)) {                 \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// Greater-than-or-equal assertion without a comment.
#define SEQAN_ASSERT_GEQ(_arg1, _arg2)                                  \
    do {                                                                \
        if (!::seqan::ClassTest::testGeq(__FILE__, __LINE__,            \
                                         (_arg1), # _arg1,               \
                                         (_arg2), # _arg2)) {            \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// Greater-than-or-equal assertion with a comment.
#define SEQAN_ASSERT_GEQ_MSG(_arg1, _arg2, ...)                         \
    do {                                                                \
        if (!::seqan::ClassTest::testGeq(__FILE__, __LINE__,            \
                                         (_arg1), # _arg1,               \
                                         (_arg2), # _arg2,               \
                                         __VA_ARGS__)) {                \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// Greater-than assertion without a comment.
#define SEQAN_ASSERT_GT(_arg1, _arg2)                                   \
    do {                                                                \
        if (!::seqan::ClassTest::testGt(__FILE__, __LINE__,             \
                                        (_arg1), # _arg1,                \
                                        (_arg2), # _arg2)) {             \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// Greater-than assertion with a comment.
#define SEQAN_ASSERT_GT_MSG(_arg1, _arg2, ...)                          \
    do {                                                                \
        if (!::seqan::ClassTest::testGt(__FILE__, __LINE__,             \
                                        (_arg1), # _arg1,                \
                                        (_arg2), # _arg2,                \
                                        __VA_ARGS__)) {                 \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// TODO(holtgrew): Rename to SEQAN_ASSERT once that name is free.;
// Trueness assertion with a comment.
//
// Usage:  SEQAN_ASSERT(false);
#define SEQAN_ASSERT(_arg1)                                        \
    do {                                                                \
        if (!::seqan::ClassTest::testTrue(__FILE__, __LINE__,           \
                                          (_arg1), # _arg1)) {           \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// TODO(holtgrew): Rename to SEQAN_ASSERT once that name is free.;
// Trueness assertion with a comment.
#define SEQAN_ASSERT_MSG(_arg1, ...)                               \
    do {                                                                \
        if (!::seqan::ClassTest::testTrue(__FILE__, __LINE__,           \
                                          (_arg1), # _arg1,              \
                                          __VA_ARGS__)) {             \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// Falseness assertion without a comment.
//
// Usage:  SEQAN_ASSERT_NOT(false);
#define SEQAN_ASSERT_NOT(_arg1)                                       \
    do {                                                              \
        if (!::seqan::ClassTest::testFalse(__FILE__, __LINE__,        \
                                           (_arg1), # _arg1)) {        \
            ::seqan::ClassTest::fail();                               \
        }                                                             \
    } while (false)


// Falseness assertion with a comment.
#define SEQAN_ASSERT_NOT_MSG(_arg1, ...)                              \
    do {                                                              \
        if (!::seqan::ClassTest::testFalse(__FILE__, __LINE__,        \
                                           (_arg1), # _arg1,           \
                                           __VA_ARGS__)) {          \
            ::seqan::ClassTest::fail();                               \
        }                                                             \
    } while (false)


#elif SEQAN_ENABLE_DEBUG && defined(__CUDA_ARCH__)

#define SEQAN_ASSERT_EQ(_arg1, _arg2) do { assert(_arg1 == _arg2); } while (false)
#define SEQAN_ASSERT_EQ_MSG(_arg1, _arg2, ...) do { assert(_arg1 == _arg2); } while (false)
#define SEQAN_ASSERT_NEQ(_arg1, _arg2) do { assert(_arg1 != _arg2); } while (false)
#define SEQAN_ASSERT_NEQ_MSG(_arg1, _arg2, ...) do { assert(_arg1 != _arg2); } while (false)
#define SEQAN_ASSERT_LEQ(_arg1, _arg2) do { assert(_arg1 <= _arg2); } while (false)
#define SEQAN_ASSERT_LEQ_MSG(_arg1, _arg2, ...) do { assert(_arg1 <= _arg2); } while (false)
#define SEQAN_ASSERT_LT(_arg1, _arg2) do { assert(_arg1 < _arg2); } while (false)
#define SEQAN_ASSERT_LT_MSG(_arg1, _arg2, ...) do { assert(_arg1 < _arg2); } while (false)
#define SEQAN_ASSERT_GEQ(_arg1, _arg2) do { assert(_arg1 >= _arg2); } while (false)
#define SEQAN_ASSERT_GEQ_MSG(_arg1, _arg2, ...) do { assert(_arg1 >= _arg2); } while (false)
#define SEQAN_ASSERT_GT(_arg1, _arg2) do { assert(_arg1 > _arg2); } while (false)
#define SEQAN_ASSERT_GT_MSG(_arg1, _arg2, ...) do { assert(_arg1 > _arg2); } while (false)
#define SEQAN_ASSERT(_arg1) do { assert(_arg1); } while (false)
#define SEQAN_ASSERT_MSG(_arg1, ...) do { assert(_arg1); } while (false)
#define SEQAN_ASSERT_NOT(_arg1) do { assert(!_arg1); } while (false)
#define SEQAN_ASSERT_NOT_MSG(_arg1, ...) do { assert(!_arg1); } while (false)
#define SEQAN_ASSERT_FAIL(...) do { assert(false); } while (false)

#else

#define SEQAN_ASSERT_EQ(_arg1, _arg2) do {} while (false)
#define SEQAN_ASSERT_EQ_MSG(_arg1, _arg2, ...) do {} while (false)
#define SEQAN_ASSERT_NEQ(_arg1, _arg2) do {} while (false)
#define SEQAN_ASSERT_NEQ_MSG(_arg1, _arg2, ...) do {} while (false)
#define SEQAN_ASSERT_LEQ(_arg1, _arg2) do {} while (false)
#define SEQAN_ASSERT_LEQ_MSG(_arg1, _arg2, ...) do {} while (false)
#define SEQAN_ASSERT_LT(_arg1, _arg2) do {} while (false)
#define SEQAN_ASSERT_LT_MSG(_arg1, _arg2, ...) do {} while (false)
#define SEQAN_ASSERT_GEQ(_arg1, _arg2) do {} while (false)
#define SEQAN_ASSERT_GEQ_MSG(_arg1, _arg2, ...) do {} while (false)
#define SEQAN_ASSERT_GT(_arg1, _arg2) do {} while (false)
#define SEQAN_ASSERT_GT_MSG(_arg1, _arg2, ...) do {} while (false)
#define SEQAN_ASSERT(_arg1) do {} while (false)
#define SEQAN_ASSERT_MSG(_arg1, ...) do {} while (false)
#define SEQAN_ASSERT_NOT(_arg1) do {} while (false)
#define SEQAN_ASSERT_NOT_MSG(_arg1, ...) do {} while (false)
#define SEQAN_ASSERT_FAIL(...) do {} while (false)

#endif  // #if defined(SEQAN_ENABLE_DEBUG) && !defined(__CUDA_ARCH__)

#else // no variadic macros

#if SEQAN_ENABLE_DEBUG
inline void SEQAN_ASSERT_FAIL(const char * comment, ...)
{
    va_list args;
    va_start(args, comment);
    ::seqan::ClassTest::vforceFail("", 0, comment, args);
    ::seqan::ClassTest::fail();
    va_end(args);
}

template <typename T1, typename T2, typename T3>
void SEQAN_ASSERT_IN_DELTA(T1 const & _arg1, T2 const & _arg2, T3 const & _arg3)
{
    if (!::seqan::ClassTest::testInDelta("", 0, _arg1, "", _arg2, "", _arg3, ""))
        ::seqan::ClassTest::fail();
}

template <typename T1, typename T2, typename T3>
void SEQAN_ASSERT_IN_DELTA_MSG(T1 const & _arg1, T2 const & _arg2, T3 const & _arg3, const char * comment, ...)
{
    va_list args;
    va_start(args, comment);
    if (!::seqan::ClassTest::vtestInDelta("", 0, _arg1, "", _arg2, "", _arg3, "", comment, args))
        ::seqan::ClassTest::fail();
    va_end(args);
}

template <typename T1, typename T2>
void SEQAN_ASSERT_EQ(T1 const & _arg1, T2 const & _arg2)
{
    if (!::seqan::ClassTest::testEqual("", 0, _arg1, "", _arg2, ""))
        ::seqan::ClassTest::fail();
}

template <typename T1, typename T2>
void SEQAN_ASSERT_EQ_MSG(T1 const & _arg1, T2 const & _arg2, const char * comment, ...)
{
    va_list args;
    va_start(args, comment);
    if (!::seqan::ClassTest::vtestEqual("", 0, _arg1, "", _arg2, "", comment, args))
        ::seqan::ClassTest::fail();
    va_end(args);
}

template <typename T1, typename T2>
void SEQAN_ASSERT_NEQ(T1 const & _arg1, T2 const & _arg2)
{
    if (!::seqan::ClassTest::testNotEqual("", _arg1, "", _arg2, ""))
        ::seqan::ClassTest::fail();
}

template <typename T1, typename T2>
void SEQAN_ASSERT_NEQ_MSG(T1 const & _arg1, T2 const & _arg2, const char * comment, ...)
{
    va_list args;
    va_start(args, comment);
    if (!::seqan::ClassTest::vtestNotEqual("", _arg1, "", _arg2, "", comment, args))
        ::seqan::ClassTest::fail();
    va_end(args);
}

template <typename T1, typename T2>
void SEQAN_ASSERT_LEQ(T1 const & _arg1, T2 const & _arg2)
{
    if (!::seqan::ClassTest::testLeq("", 0, _arg1, "", _arg2, ""))
        ::seqan::ClassTest::fail();
}

template <typename T1, typename T2>
void SEQAN_ASSERT_LEQ_MSG(T1 const & _arg1, T2 const & _arg2, const char * comment, ...)
{
    va_list args;
    va_start(args, comment);
    if (!::seqan::ClassTest::vtestLeq("", 0, _arg1, "", _arg2, "", comment, args))
        ::seqan::ClassTest::fail();
    va_end(args);
}

template <typename T1, typename T2>
void SEQAN_ASSERT_LT(T1 const & _arg1, T2 const & _arg2)
{
    if (!::seqan::ClassTest::testLt("", 0, _arg1, "", _arg2, ""))
        ::seqan::ClassTest::fail();
}

template <typename T1, typename T2>
void SEQAN_ASSERT_LT_MSG(T1 const & _arg1, T2 const & _arg2, const char * comment, ...)
{
    va_list args;
    va_start(args, comment);
    if (!::seqan::ClassTest::vtestLt("", 0, _arg1, "", _arg2, "", comment, args))
        ::seqan::ClassTest::fail();
    va_end(args);
}

template <typename T1, typename T2>
void SEQAN_ASSERT_GEQ(T1 const & _arg1, T2 const & _arg2)
{
    if (!::seqan::ClassTest::testGeq("", 0, _arg1, "", _arg2, ""))
        ::seqan::ClassTest::fail();
}

template <typename T1, typename T2>
void SEQAN_ASSERT_GEQ_MSG(T1 const & _arg1, T2 const & _arg2, const char * comment, ...)
{
    va_list args;
    va_start(args, comment);
    if (!::seqan::ClassTest::vtestGeq("", 0, _arg1, "", _arg2, "", comment, args))
        ::seqan::ClassTest::fail();
    va_end(args);
}

template <typename T1, typename T2>
void SEQAN_ASSERT_GT(T1 const & _arg1, T2 const & _arg2)
{
    if (!::seqan::ClassTest::testGt("", 0, _arg1, "", _arg2, ""))
        ::seqan::ClassTest::fail();
}

template <typename T1, typename T2>
void SEQAN_ASSERT_GT_MSG(T1 const & _arg1, T2 const & _arg2, const char * comment, ...)
{
    va_list args;
    va_start(args, comment);
    if (!::seqan::ClassTest::vtestGt("", 0, _arg1, "", _arg2, "", comment, args))
        ::seqan::ClassTest::fail();
    va_end(args);
}

template <typename T1>
void SEQAN_ASSERT(T1 const & _arg1)
{
    if (!::seqan::ClassTest::testTrue("", 0, _arg1, ""))
        ::seqan::ClassTest::fail();
}

template <typename T1>
void SEQAN_ASSERT_MSG(T1 const & _arg1, const char * comment, ...)
{
    va_list args;
    va_start(args, comment);
    if (!::seqan::ClassTest::vtestTrue("", 0, _arg1, "", comment, args))
        ::seqan::ClassTest::fail();
    va_end(args);
}

template <typename T1>
void SEQAN_ASSERT_NOT(T1 const & _arg1)
{
    if (!::seqan::ClassTest::testFalse("", 0, _arg1, ""))
        ::seqan::ClassTest::fail();
}

template <typename T1>
void SEQAN_ASSERT_NOT_MSG(T1 const & _arg1, const char * comment, ...)
{
    va_list args;
    va_start(args, comment);
    if (!::seqan::ClassTest::vtestFalse("", 0, _arg1, "", comment, args))
        ::seqan::ClassTest::fail();
    va_end(args);
}

#else // #if SEQAN_ENABLE_DEBUG

inline void SEQAN_ASSERT_FAIL(const char * comment, ...) {}
template <typename T1, typename T2, typename T3>
void SEQAN_ASSERT_IN_DELTA(T1 const & _arg1, T2 const & _arg2, T3 const & _arg3) {}
template <typename T1, typename T2, typename T3>
void SEQAN_ASSERT_IN_DELTA_MSG(T1 const & _arg1, T2 const & _arg2, T3 const & _arg3, const char * comment, ...) {}
template <typename T1, typename T2>
void SEQAN_ASSERT_EQ(T1 const & _arg1, T2 const & _arg2) {}
template <typename T1, typename T2>
void SEQAN_ASSERT_EQ_MSG(T1 const & _arg1, T2 const & _arg2, const char * comment, ...) {}
template <typename T1, typename T2>
void SEQAN_ASSERT_NEQ(T1 const & _arg1, T2 const & _arg2) {}
template <typename T1, typename T2>
void SEQAN_ASSERT_NEQ_MSG(T1 const & _arg1, T2 const & _arg2, const char * comment, ...) {}
template <typename T1, typename T2>
void SEQAN_ASSERT_LEQ(T1 const & _arg1, T2 const & _arg2) {}
template <typename T1, typename T2>
void SEQAN_ASSERT_LEQ_MSG(T1 const & _arg1, T2 const & _arg2, const char * comment, ...) {}
template <typename T1, typename T2>
void SEQAN_ASSERT_LT(T1 const & _arg1, T2 const & _arg2) {}
template <typename T1, typename T2>
void SEQAN_ASSERT_LT_MSG(T1 const & _arg1, T2 const & _arg2, const char * comment, ...) {}
template <typename T1, typename T2>
void SEQAN_ASSERT_GEQ(T1 const & _arg1, T2 const & _arg2) {}
template <typename T1, typename T2>
void SEQAN_ASSERT_GEQ_MSG(T1 const & _arg1, T2 const & _arg2, const char * comment, ...) {}
template <typename T1, typename T2>
void SEQAN_ASSERT_GT(T1 const & _arg1, T2 const & _arg2) {}
template <typename T1, typename T2>
void SEQAN_ASSERT_GT_MSG(T1 const & _arg1, T2 const & _arg2, const char * comment, ...) {}
template <typename T1>
void SEQAN_ASSERT(T1 const & _arg1) {}
template <typename T1>
void SEQAN_ASSERT_MSG(T1 const & _arg1, const char * comment, ...) {}
template <typename T1>
void SEQAN_ASSERT_NOT(T1 const & _arg1) {}
template <typename T1>
void SEQAN_ASSERT_NOT_MSG(T1 const & _arg1, const char * comment, ...) {}

#endif // #if SEQAN_ENABLE_DEBUG

#endif // no variadic macros

// Returns a string (of type char*) with the path to the called binary.
//
// Use this to locate files relative to the test binary.
#define SEQAN_PROGRAM_PATH                      \
    ::seqan::ClassTest::StaticData::basePath()

/*!
 * @macro SEQAN_PATH_TO_ROOT
 * @headerfile <seqan/basic.h>
 * @brief Return path to the checkout root directory.
 *
 * @signature TCharPtr SEQAN_PATH_TO_ROOT()
 *
 * @return TCharPtr <tt>char const *</tt>, string with the path to the parent directory of the tests directory.
 *
 * This only works when using the SeqAn SVN checkout!
 *
 * The pointed to string is initialized on program startup by the code generated by SEQAN_BEGIN_TESTSUITE.
 *
 * @section Examples
 *
 * @code{.cpp}
 * CharString buffer = SEQAN_PATH_TO_ROOT();
 * append(buffer, "/tests/files/example.txt");
 *
 * FILE *f = fopen(toCString(buffer), "w");
 * fprintf(f, "Test Data");
 * fclose(f);
 * @endcode
 *
 * @see SEQAN_TEMP_FILENAME
 */

// TODO(holtgrew): Subject to change wiht restructuring.
// Returns a const char * string with the path to the projects directory.
#define SEQAN_PATH_TO_ROOT()                      \
    ::seqan::ClassTest::StaticData::pathToRoot()


// Returns the POSIX int file handle to an open file.
// TODO(holtgrewe): Uncomment if openTempFile has been implemented.
// #define SEQAN_OPEN_TEMP_FILE() (::seqan::ClassTest::openTempFile())

/*!
 * @macro SEQAN_TEMP_FILENAME
 * @headerfile <seqan/basic.h>
 * @brief Generates the name to a temporary file.
 *
 * @signature TCharType SEQAN_TEMP_FILENAME();
 *
 * @return TCharType <tt>char const *</tt>, string with the path to a temporary file.
 *
 * @section Remarks
 *
 * The pointed to string is stored in a buffer and is overwritten by the next call to this macro. Copy it out if you
 * need it.
 *
 * @section Examples
 *
 * @code{.cpp}
 * const char *p = SEQAN_TEMP_FILENAME();
 * buffer char tempFilename[1000];
 * strcpy(tempFilename, p);
 * FILE *f = fopen(tempFilename, "w");
 * fprintf(f, "Test Data");
 * fclose(f);
 * @endcode
 * @see SEQAN_PATH_TO_ROOT
 */

// Returns a temporary filename.
#define SEQAN_TEMP_FILENAME() (::seqan::ClassTest::tempFileName())


#if SEQAN_ENABLE_CHECKPOINTS

// Create a check point at the point where the macro is placed.
// TODO(holtgrew): Should be called SEQAN_CHECK_POINT to be consistent.
#define SEQAN_CHECKPOINT                                        \
    ::seqan::ClassTest::registerCheckPoint(__LINE__, __FILE__);

// Call the check point verification code for the given file.
#define SEQAN_VERIFY_CHECKPOINTS(filename)          \
    ::seqan::ClassTest::verifyCheckPoints(filename)

#else  // #if SEQAN_ENABLE_CHECKPOINTS

#define SEQAN_CHECKPOINT

// If checkpoints are to be verified if testing is disabled then print
// a warning.
#define SEQAN_VERIFY_CHECKPOINTS(filename)                              \
    do {                                                                \
        fprintf(stderr, ("WARNING: Check point verification is "        \
                         "disabled. Trying to verify %s from %s:%d.\n"), \
                filename, __FILE__, __LINE__);                          \
    } while (false)

#endif  // #if SEQAN_ENABLE_CHECKPOINTS

#if !SEQAN_ENABLE_TESTING

#define SEQAN_BEGIN_TESTSUITE(suite_name)                               \
    int main(int argc, char ** argv) {                                   \
        (void) argc;                                                        \
        (void) argv;                                                        \
        fprintf(stderr, "Warning: SEQAN_ENABLE_TESTING is wrong and you used the macro SEQAN_BEGIN_TESTSUITE!\n");
#define SEQAN_END_TESTSUITE \
    return 0;                                   \
    }
#define SEQAN_CALL_TEST(test_name) do { SEQAN_TEST_ ## test_name(); } while (false)
#define SEQAN_SKIP_TEST do {} while (false)

#endif  // #if !SEQAN_ENABLE_TESTING

}  // namespace seqan

#endif  // SEQAN_INCLUDE_SEQAN_BASIC_DEBUG_TEST_SYSTEM_H_
