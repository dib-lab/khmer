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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// Definition of basic exceptions.
// ==========================================================================

#ifndef SEQAN_BASIC_BASIC_EXCEPTION_H_
#define SEQAN_BASIC_BASIC_EXCEPTION_H_

// ============================================================================
// Prerequisites
// ============================================================================

#include <typeinfo>
#include <exception>
#include <stdexcept>

#ifdef PLATFORM_GCC
#include <cxxabi.h>
#endif

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename T>
struct Tag;

//struct Nothing_;
//typedef Tag<Nothing_> Nothing;

// ============================================================================
// Macros
// ============================================================================

/*!
 * @defgroup ExceptionHandling SeqAn Exception Handling
 * @brief Macros supporting exception handling on various platforms.
 */

/*!
 * @macro ExceptionHandling#SEQAN_EXCEPTIONS
 * @headerfile <seqan/basic.h>
 * @brief Determines whether exceptions are enabled or not.
 *
 * @signature #define SEQAN_EXCEPTIONS
 *
 * @see ExceptionHandling#SEQAN_TRY
 * @see ExceptionHandling#SEQAN_CATCH
 * @see ExceptionHandling#SEQAN_THROW
 * @see Exception
 */

#define SEQAN_EXCEPTIONS    __EXCEPTIONS

/*!
 * @macro ExceptionHandling#SEQAN_TRY
 * @headerfile <seqan/basic.h>
 * @brief Replaces the C++ try keyword.
 *
 * @signature SEQAN_TRY {} SEQAN_CATCH() {}
 *
 * When exceptions are disabled, i.e. SEQAN_EXCEPTIONS is set to false, the code inside the try block is always executed.
 *
 * @see ExceptionHandling#SEQAN_CATCH
 * @see ExceptionHandling#SEQAN_THROW
 * @see Exception
 *
 * @section Examples
 *
 * @code{.cpp}
 *
 * SEQAN_TRY
 * {
 *     SEQAN_THROW(Exception)
 * }
 * SEQAN_CATCH(Exception const & e)
 * {
 *     std::cerr << e.what() << std::endl;
 * }
 *
 * @endcode
 */

/*!
 * @macro ExceptionHandling#SEQAN_CATCH
 * @headerfile <seqan/basic.h>
 * @brief Replaces the C++ catch keyword.
 *
 * @signature SEQAN_TRY {} SEQAN_CATCH() {}
 *
 * When exceptions are disabled, i.e. SEQAN_EXCEPTIONS is set to false, the code inside the catch block is never executed.
 *
 * @see ExceptionHandling#SEQAN_TRY
 * @see ExceptionHandling#SEQAN_THROW
 * @see Exception
 *
 * @section Examples
 *
 * See @link ExceptionHandling#SEQAN_TRY @endlink for a full example.
 */

/*!
 * @macro ExceptionHandling#SEQAN_THROW
 * @headerfile <seqan/basic.h>
 * @brief Replaces the C++ throw keyword.
 *
 * @signature SEQAN_THROW(Exception);
 *
 * When exceptions are disabled, i.e. AssertMacros#SEQAN_EXCEPTIONS is set to false, the macro turns into SEQAN_FAIL.
 *
 * @see ExceptionHandling#SEQAN_TRY
 * @see ExceptionHandling#SEQAN_CATCH
 * @see AssertMacros#SEQAN_FAIL
 * @see Exception
 *
 * @section Examples
 *
 * See @link ExceptionHandling#SEQAN_TRY @endlink for a full example.
 */

#ifdef SEQAN_EXCEPTIONS

#define SEQAN_TRY           try
#define SEQAN_CATCH(E)      catch(E)
#define SEQAN_THROW(E)      throw E
#define SEQAN_RETHROW       throw

#else

#define SEQAN_TRY           if (true)
#define SEQAN_CATCH(E)      if (false)
//#define SEQAN_CATCH(E)      for (E ; false; )
#define SEQAN_THROW(E)      SEQAN_FAIL(#E)
#define SEQAN_RETHROW

#endif // #ifdef SEQAN_EXCEPTIONS

// ============================================================================
// Exceptions
// ============================================================================

// ----------------------------------------------------------------------------
// Basic Exception
// ----------------------------------------------------------------------------

/*!
 * @class Exception
 * @headerfile <seqan/basic.h>
 * @brief Generic SeqAn exception.
 * @signature typedef std::exception Exception;
 *
 * @fn Exception::Exception
 * @brief Constructor.
 *
 * @signature Exception::Exception(msg);
 * @param[in] msg The message as a <tt>std::string</tt>.
 */

typedef std::exception          Exception;

// ----------------------------------------------------------------------------
// Exception BadAlloc
// ----------------------------------------------------------------------------

/*!
 * @class BadAlloc
 * @headerfile <seqan/basic.h>
 * @brief Generic SeqAn exception.
 * @signature typedef std::bad_alloc BadAlloc;
 *
 * @fn BadAlloc::BadAlloc
 * @brief Constructor.
 *
 * @signature BadAlloc::BadAlloc(msg);
 * @param[in] msg The message as a <tt>std::string</tt>.
 */

typedef std::bad_alloc          BadAlloc;

// ----------------------------------------------------------------------------
// Exception BadCast
// ----------------------------------------------------------------------------

/*!
 * @class BadCast
 * @headerfile <seqan/basic.h>
 * @brief Generic SeqAn exception.
 * @signature typedef std::bad_cast BadCast;
 *
 * @fn BadCast::BadCast
 * @brief Constructor.
 *
 * @signature BadCast::BadCast(msg);
 * @param[in] msg The message as a <tt>std::string</tt>.
 */

typedef std::bad_cast           BadCast;

// ----------------------------------------------------------------------------
// Exceptions Bad*
// ----------------------------------------------------------------------------
// NOTE(esiragusa): These exceptions can be introduced as long as we need them.

//typedef std::bad_exception      BadException;
//typedef std::bad_typeid         BadTypeId;
//typedef std::bad_function_call  BadFunctionCall;
//typedef std::bad_weak_ptr       BadWeakPtr;

// ----------------------------------------------------------------------------
// Exception RuntimeError
// ----------------------------------------------------------------------------

/*!
 * @class RuntimeError
 * @headerfile <seqan/basic.h>
 * @brief Runtime error exception.
 * @signature typedef std::runtime_error RuntimeError;
 *
 *
 * @fn RuntimeError::RuntimeError
 * @brief Constructor.
 *
 * @signature RuntimeError::RuntimeError(msg);
 * @param[in] msg The message as a <tt>std::string</tt>.
 */

typedef std::runtime_error      RuntimeError;

// ----------------------------------------------------------------------------
// Exception LogicError
// ----------------------------------------------------------------------------
// NOTE(esiragusa): Always prefer SEQAN_ASSERT to logic error exceptions.

//typedef std::logic_error        LogicError;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction ExceptionMessage
// ----------------------------------------------------------------------------

template <typename T, typename TSpec = void>
struct ExceptionMessage
{
    static const std::string VALUE;
};

template <typename T, typename TSpec>
const std::string ExceptionMessage<T, TSpec>::VALUE;

// ----------------------------------------------------------------------------
// Function getExceptionMessage()
// ----------------------------------------------------------------------------

template <typename TFunctor, typename TContext>
inline std::string const &
getExceptionMessage(TFunctor const &, TContext const &)
{
    return ExceptionMessage<TFunctor, TContext>::VALUE;
}

// ============================================================================
// Functors
// ============================================================================

// ----------------------------------------------------------------------------
// Functor AssertFunctor
// ----------------------------------------------------------------------------

template <typename TFunctor, typename TException, typename TContext = Nothing, bool RETURN_VALUE = false>
struct AssertFunctor
{
    TFunctor func;

    AssertFunctor() {}

    AssertFunctor(TFunctor & func) :
        func(func)
    {}

    std::string escapeChar(unsigned char val)
    {
        if (val <= '\r')
        {
            static const char * const escapeCodes[14] = {
                "\\0",  "\\1",  "\\2",  "\\3",  "\\4",  "\\5",  "\\6",  "\\a",
                "\\b",  "\\t",  "\\n",  "\\v",  "\\f",  "\\r" };
            return std::string(escapeCodes[val]);
        }
        else if (' ' <= val && val < 128u)
            return std::string() + (char)val;
        else
        {
            char buffer[6]; // 5 + 1, e.g. "\0xff" + trailing zero
            sprintf(buffer, "\\%#2x", (unsigned)val);
            return std::string(buffer);
        }
    }

    template <typename TValue>
    bool operator() (TValue const & val)
    {
        if (SEQAN_UNLIKELY(!func(val)))
            throw TException(std::string("Unexpected character '") + escapeChar(val) + "' found. " +
                             getExceptionMessage(func, TContext()));
        return RETURN_VALUE;
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function globalExceptionHandler()
// ----------------------------------------------------------------------------

#if defined(SEQAN_EXCEPTIONS) && !defined(SEQAN_NO_GLOBAL_EXCEPTION_HANDLER)
// Declare global exception handler.
static void globalExceptionHandler();

// Install global exception handler.
static const std::terminate_handler _globalExceptionHandler = std::set_terminate(globalExceptionHandler);

inline static void globalExceptionHandler()
{
    SEQAN_TRY
    {
        SEQAN_RETHROW;
    }
    SEQAN_CATCH(Exception & e)
    {
        SEQAN_FAIL("Uncaught exception of type %s: %s", toCString(Demangler<Exception>(e)), e.what());
    }
    SEQAN_CATCH(...)
    {
        SEQAN_FAIL("Uncaught exception of unknown type.\n");
    }
}
#endif  // #if defined(SEQAN_EXCEPTIONS) && !defined(SEQAN_NO_GLOBAL_EXCEPTION_HANDLER)

}  // namespace seqan

#endif  // #ifndef SEQAN_BASIC_BASIC_EXCEPTION_H_
