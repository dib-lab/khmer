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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Metaprogramming control structures.
//
// This header defines metaprogramming control structure such as conditionals
// and loops.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_METAPROGRAMMING_CONTROL_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_METAPROGRAMMING_CONTROL_H_

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
 * @snippet core/demos/basic/metaprogramming_switch.cpp switch demo
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

/**
.Tag.NilCase
..cat:Metaprogramming
..summary:Metaprogramming $default:$ case expression.
..signature:NilCase
..remarks:$NilCase$ is returned by metafunction @Metafunction.Switch@ if no case matched.
..include:seqan/basic.h
..remarks:The documentation of @Metafunction.Switch@ gives an example.
..see:Tag.Case
..see:Metafunction.Switch

.Tag.Case
..cat:Metaprogramming
..summary:Metaprogramming $case$ expression.
..signature:Case<TAG[, NextCase]>
..param.TAG:
..param.NextCase:Optional next case of type @Tag.Case@, defaults to @Tag.NilCase@.
...type:Tag.NilCase
...type:Tag.Case
..include:seqan/basic.h
..remarks:The documentation of @Metafunction.Switch@ gives an example.
..see:Tag.NilCase
..see:Metafunction.Switch

.Metafunction.Switch
..cat:Metaprogramming
..summary:Metaprogramming $switch$ expression.
..signature:Switch<TAG, Case>::Type
..param.TAG:The case label.
...type:nolink:$int$
..param.Case:Cascade of $Case$ tags.
...type:Tag.Case
...type:Tag.NilCase
..returns:The selected type from the @Tag.Case@ cascade or @Tag.NilCase@.
..include:seqan/basic.h
..see:Tag.NilCase
..see:Tag.Case
..example.code:
int switchTest(Nothing const &) { return -1; }
int switchTest(False const &) { return 0; }
int switchTest(True const &) { return 1; }
int switchTest(NilCase const &) { return 2; }

template <int X>
struct SwitchTest
{
    typedef typename Switch<
        X,
        Case<-1, Nothing,
        Case<0, False,
        Case<1, True
        > > > >::Type Type;
};

SEQAN_DEFINE_TEST(test_metaprogramming_switch)
{
    typedef typename SwitchTest<-1>::Type T1;
    typedef typename SwitchTest< 0>::Type T2;
    typedef typename SwitchTest< 1>::Type T3;
    typedef typename SwitchTest< 2>::Type T4;

    std::cout << switchTest(T1()) << std::endl;  // => -1
    std::cout << switchTest(T2()) << std::endl;  // =>  0
    std::cout << switchTest(T3()) << std::endl;  // =>  1
    std::cout << switchTest(T4()) << std::endl;  // =>  2
 }
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

// TODO(holtgrew): This is more of a class than a function or metafunction.

/*!
 * @fn Loop
 * @headerfile <seqan/basic.h>
 *
 * @signature void Loop<TWorker, COUNT>::run(args);
 *
 * @tparam TWorker A struct with a static inline void function called <tt>body</tt>.  <tt>body</tt> should have two
 *                 parameters, one for passing in values and state from the outside and the second is an int.  The
 *                 function will be called <tt>COUNT</tt> times with the same reference for the first one and the values
 *                 <tt>COUNT</tt>, <tt>COUNT - 1</tt>, ..., <tt>1</tt> for the second parameter.
 * @param COUNT    An <tt>int</tt> constant.
 *
 * @section Example
 *
 * We define the following worker to print an integer.  The first argument is of <tt>Nothing</tt> as a dummy.  Note that
 * the parameter is not const.
 *
 * @snippet core/demos/basic/metaprogramming_control.cpp print worker
 *
 * The following shows an example calling <tt>PrintWorker::body()</tt> through Loop.  We have to create a local variable
 * since the first parameter is not const.  The reason for this is that the parameter can also be used for a mutable
 * object that holds some state.
 *
 * @snippet core/demos/basic/metaprogramming_control.cpp print worker call loop reverse
 *
 * @see LoopReverse
 */

/**
.Metafunction.Loop:
..cat:Metaprogramming
..summary:Metafunction returning a function that iterates over a static integer range.
..signature:Loop<Worker, I>::run(Arg & arg)
..param.Worker:A worker $struct$. It has to implement the static (and preferably finline/inline) function $body$ that accepts two parameters. The first one will be a reference to $arg$, as given to $run()$.  The second will be the current value of the iterating variable.
..param.I:The upper limit for the iteration.
..param.arg:The argument to be passed into the workers' $body()$ function.
..remarks:The loop will go from 1 up to and including $I$.
..see:Metafunction.LoopReverse
..include:seqan/basic.h
..example.text:Print the values 1, 2, ..., $I-1$, $I$.
..example.code:
struct PrintWorker
{
    static inline void body(Nothing & arg, int I)
    {
        (void)arg;  // ignored
        printf("%d\n", I);
    }
};

Loop<PrintWorker, 10>::run(Nothing());
// This will print the numbers 1, 2, ..., 9, 10.
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

// TODO(holtgrew): This is more of a class than a function or metafunction.

/*!
 * @fn LoopReverse
 * @headerfile <seqan/basic.h>
 *
 * @signature void LoopReverse<TWorker, COUNT>::run(args);
 *
 * @tparam TWorker A struct with a static inline void function called <tt>body</tt>.  <tt>body</tt> should have two
 *                 parameters, one for passing in values and state from the outside and the second is an int.  The
 *                 function will be called <tt>COUNT</tt> times with the same reference for the first one and the values
 *                 <tt>COUNT</tt>, <tt>COUNT - 1</tt>, ..., <tt>1</tt> for the second parameter.
 * @param COUNT    An <tt>int</tt> constant.
 *
 * @section Example
 *
 * We define the following worker to print an integer.  The first argument is of <tt>Nothing</tt> as a dummy.  Note that
 * the parameter is not const.
 *
 * @snippet core/demos/basic/metaprogramming_control.cpp print worker
 *
 * The following shows an example calling <tt>PrintWorker::body()</tt> through LoopReverse.  We have to create a local
 * variable since the first parameter is not const.  The reason for this is that the parameter can also be used for a
 * mutable object that holds some state.
 *
 * @snippet core/demos/basic/metaprogramming_control.cpp print worker call loop
 *
 * @see Loop
 */

/**
.Metafunction.LoopReverse:
..cat:Metaprogramming
..summary:Metafunction returning a function that iterates over a static integer range in reverse order.
..signature:LoopReverse<Worker, I>::run(Arg & arg)
..param.Worker:A worker $struct$. It has to implement the static (and preferably finline/inline) function $body$ that accepts two parameters. The first one will be a reference to $arg$, as given to $run()$.  The second will be the current value of the iterating variable.
..param.I:The upper limit for the iteration.
..param.arg:The argument to be passed into the workers' $body()$ function.
..remarks:The loop will go from $I$ down to and including 1.
..include:seqan/basic.h
..see:Metafunction.Loop
..example.text:Print the values $I$, $I - 1$, ..., 2, 1.
..example.code:
struct PrintWorker
{
    static inline body(Nothing & arg, int I)
    {
        (void)arg;  // ignored
        printf("%d\n", I);
    }
};

Loop<PrintWorker, 10>::run(Nothing());
// This will print the numbers 1, 2, ..., 9, 10.
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

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_METAPROGRAMMING_CONTROL_H_
