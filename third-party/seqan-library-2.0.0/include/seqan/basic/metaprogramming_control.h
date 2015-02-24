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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Metaprogramming control structures.
//
// This header defines metaprogramming control structure such as conditionals
// and loops.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_METAPROGRAMMING_CONTROL_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_METAPROGRAMMING_CONTROL_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Switch;  Supporting Tags Case, NilCase.
// ----------------------------------------------------------------------------

/*!
 * @defgroup MetafunctionSwitch Switch Metafunction Types
 * @brief Tags for the metafunction Switch and the metafunction itself.
 *
 * @section Example
 *
 * The following shows a complete example of using the Switch statement.
 *
 * @snippet demos/basic/metaprogramming_switch.cpp switch demo
 */

// TODO(holtgrew): Use Tag<>?
/*!
 * @tag MetafunctionSwitch#NilCase
 * @brief Tag for terminating the case in Switch statement.
 *
 * @signature struct NilCase {};
 */

/*!
 * @mfn MetafunctionSwitch#Case
 * @brief Tag for one case.
 *
 * @signature template <int TAG, typename TResult, typename TNext>
 *            struct Case;
 *
 * @tparam TAG     The <tt>int</tt> tag number to use.
 * @tparam TResult The type to return in <tt>Case<...>::Type</tt> if matches.
 * @tparam TNext   The next <tt>Case</tt>.
 */

/*!
 * @mfn MetafunctionSwitch#Switch
 * @brief Switch statement for metaprogramming.
 *
 * @signature Switch<TAG, TCase>::Type
 *
 * @tparam TAG   <tt>int</tt> with the current value.
 * @tparam TCase First <tt>Case</tt> statement.
 *
 * @return Type The selected type.
 */

const int DEFAULT_ = ~(~0u >> 1); // initialize with the smallest int

struct NilCase {};

template <int TAG, typename Type_, typename Next_ = NilCase>
struct Case
{
    enum { TAG_ = TAG };
    typedef Type_ Type;
    typedef Next_ Next;
};

template <int TAG, typename Case_>
struct Switch
{
    typedef typename Case_::Next NextCase_;
    enum
    {
        CASE_TAG_ = Case_::TAG_,
        FOUND_    = (CASE_TAG_ == TAG || CASE_TAG_ == DEFAULT_)
    };

    typedef typename IfC<FOUND_,
                         typename Case_::Type,
                         typename Switch<TAG, NextCase_>::Type
                         >::Type Type;
};

template <int TAG>
struct Switch<TAG, NilCase>
{
    typedef NilCase Type;
};

// ----------------------------------------------------------------------------
// Metafunction Loop
// ----------------------------------------------------------------------------

/*!
 * @class Loop
 * @headerfile <seqan/basic.h>
 * @brief Helper for loops.
 *
 * @signature template <typename TWorker, unsigned COUNT>
 *            struct Loop;
 *
 * @tparam TWorker A struct with a static inline void function called <tt>body</tt>.  <tt>body</tt> should have two
 *                 parameters, one for passing in values and state from the outside and the second is an int.  The
 *                 function will be called <tt>COUNT</tt> times with the same reference for the first one and the values
 *                 <tt>COUNT</tt>, <tt>COUNT - 1</tt>, ..., <tt>1</tt> for the second parameter.
 * @tparam COUNT   An <tt>int</tt> constant.
 *
 * @section Example
 *
 * We define the following worker to print an integer.  The first argument is of <tt>Nothing</tt> as a dummy.  Note that
 * the parameter is not const.
 *
 * @snippet demos/basic/metaprogramming_control.cpp print worker
 *
 * The following shows an example calling <tt>PrintWorker::body()</tt> through Loop.  We have to create a local variable
 * since the first parameter is not const.  The reason for this is that the parameter can also be used for a mutable
 * object that holds some state.
 *
 * @snippet demos/basic/metaprogramming_control.cpp print worker call loop reverse
 *
 * @see LoopReverse
 *
 * @fn Loop::run
 * @brief Run the loop body.
 *
 * @signature Loop::run(arg, i);
 *
 * @param[in,out] arg The argument to pass to the worker's <tt>body()</tt> function.
 * @param[in]     i   The <tt>int</tt> passed to the <tt>body()</tt> function.
 */

// Example of a loop Worker class.  Could be removed, serves no
// purpose.
struct WorkerNothing
{
    template <typename Arg>
    static inline void body(Arg & /*arg*/, int /*I*/)
    {}
};

template <typename Worker, int I>
class Loop {
public:
    template <typename Arg>
    static inline void run(Arg & arg)
    {
        Loop<Worker, I - 1>::run(arg);
        Worker::body(arg, I);
    }
};

template <typename Worker>
class Loop<Worker, 0>
{
public:
    // end of loop
    template <typename Arg>
    static inline void run(Arg &)
    {}
};

// ----------------------------------------------------------------------------
// Metafunction LoopReverse
// ----------------------------------------------------------------------------

/*!
 * @class LoopReverse
 * @brief Helper for reverse loops.
 *
 * @headerfile <seqan/basic.h>
 *
 * @signature template <typename TWorker, unsigned COUNT>
 *            struct LoopReverse;
 *
 * @tparam TWorker A struct with a static inline void function called <tt>body</tt>.  <tt>body</tt> should have two
 *                 parameters, one for passing in values and state from the outside and the second is an int.  The
 *                 function will be called <tt>COUNT</tt> times with the same reference for the first one and the values
 *                 <tt>COUNT</tt>, <tt>COUNT - 1</tt>, ..., <tt>1</tt> for the second parameter.
 * @tparam COUNT   An <tt>int</tt> constant.
 *
 * @section Example
 *
 * We define the following worker to print an integer.  The first argument is of <tt>Nothing</tt> as a dummy.  Note that
 * the parameter is not const.
 *
 * @snippet demos/basic/metaprogramming_control.cpp print worker
 *
 * The following shows an example calling <tt>PrintWorker::body()</tt> through LoopReverse.  We have to create a local
 * variable since the first parameter is not const.  The reason for this is that the parameter can also be used for a
 * mutable object that holds some state.
 *
 * @snippet demos/basic/metaprogramming_control.cpp print worker call loop
 *
 * @see Loop
 *
 * @fn LoopReverse::run
 * @brief Run the loop body.
 *
 * @signature LoopReverse::run(arg, i);
 *
 * @param[in,out] arg The argument to pass to the worker's <tt>body()</tt> function.
 * @param[in]     i   The <tt>int</tt> passed to the <tt>body()</tt> function.
 */

template <typename Worker, int I>
class LoopReverse
{
public:
    template <typename Arg>
    static inline void run(Arg & arg)
    {
        Worker::body(arg, I);
        LoopReverse<Worker, I - 1>::run(arg);
    }
};

template <typename Worker>
class LoopReverse<Worker, 0>
{
  public:
    // end of loop
    template <typename Arg>
    static inline void run(Arg &) {}
};

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_METAPROGRAMMING_CONTROL_H_
