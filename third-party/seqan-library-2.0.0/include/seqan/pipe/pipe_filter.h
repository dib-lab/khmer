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

#ifndef SEQAN_HEADER_PIPE_FILTER_H
#define SEQAN_HEADER_PIPE_FILTER_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

    template <typename TValue, typename TResult = typename Value<TValue, 1>::Type>
    struct filterI1 : public std::unary_function<TValue, TResult>
    {
        inline TResult operator() (const TValue & x) const
        {
            return x.i1;
        }
    };

    template <typename TValue, typename TResult = typename Value<TValue, 2>::Type>
    struct filterI2 : public std::unary_function<TValue, TResult>
    {
        inline TResult operator() (const TValue & x) const
        {
            return x.i2;
        }
    };

    template <typename TValue, typename TResult = typename Value<TValue, 3>::Type>
    struct filterI3 : public std::unary_function<TValue, TResult>
    {
        inline TResult operator() (const TValue & x) const
        {
            return x.i3;
        }
    };


    template < typename TFunctor >
    struct Filter;

    template < typename TInput, typename TFunctor >
    struct Value< Pipe< TInput, Filter<TFunctor> > >
    {
        typedef typename TFunctor::result_type Type;
    };

/*!
 * @class Filter
 * @extends Pipe
 * @headerfile <seqan/pipe.h>
 * @brief Applies a specific function to the input stream.
 *
 * @signature template <typename TInput, typename TFunctor>
 *            struct Pipe<TInput, Filter<TFunctor> >;
 *
 * @tparam TFunctor A unary function (see STL's <tt>unary_function</tt>).  The argument type of <tt>TFunctor</tt>
 *                  must be <tt>VALUE&lt;TInput&gt;::Type</tt>.
 * @tparam TInput   The type of the pipeline module this module reads from.
 *
 * The output type of this pipe is the result type of <tt>TFunctor</tt>.
 */

    //////////////////////////////////////////////////////////////////////////////
    // filter class
    template <typename TInput, typename TFunctor >
    struct Pipe< TInput, Filter<TFunctor> >
    {
        TInput      &in;
        TFunctor    F;

/*!
 * @fn Filter::Pipe
 * @brief Constructor
 *
 * @signature Pipe::Pipe(in[, func]);
 *
 * @param[in] in   Reference to an input pipe.
 * @param[in] func A <tt>TFunctor</tt> object.
 */

        Pipe(TInput& _in):
            in(_in) {}

        Pipe(TInput& _in, const TFunctor& F_) :
            in(_in),
            F(F_) {}

        inline typename Value<Pipe>::Type const operator*() const
        {
            return F(*in);
        }

        Pipe & operator++()
        {
            ++in;
            return *this;
        }

    };

//}

}

#endif
