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

#ifndef SEQAN_HEADER_PIPE_CASTER_H
#define SEQAN_HEADER_PIPE_CASTER_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

    struct CasterReinterpret;
    struct CasterConvert;

    template < typename TValue, typename TSpec = CasterReinterpret >
    struct Caster;

    template < typename TInput, typename TValue, typename TSpec >
    struct Value< Pipe< TInput, Caster<TValue, TSpec> > > {
        typedef TValue Type;
    };

/*!
 * @class Caster
 * @extends Pipe
 * @headerfile <seqan/pipe.h>
 * @brief Casts the input type in a specific output type.
 *
 * @signature template <typename TInput, typename TValue[, typename TSpec]>
 *            class Pipe<TInput, Caster<TValue, TSpec> >;
 *
 * @tparam TInput The type of the pipeline module this module reads from.
 * @tparam TValue The new output type.
 * @tparam TSpec  <tt>CasterReinterpret</tt> (default) or <tt>CasterConvert</tt>.
 *
 * The input stream is casted using <tt>reinterpret_cast&lt;TValue&gt;</tt>.
 */

    //////////////////////////////////////////////////////////////////////////////
    // caster pipe
    template <typename TInput, typename TValue >
    struct Pipe< TInput, Caster<TValue, CasterReinterpret> >
    {
        TInput      &in;

        Pipe(TInput& _in):
            in(_in) {}

        inline TValue const & operator*() const {
            return reinterpret_cast<TValue const &>(*in);
        }

        Pipe& operator++() {
            ++in;
            return *this;
        }
    };

    template <typename TInput, typename TValue >
    struct Pipe< TInput, Caster<TValue, CasterConvert> >
    {
        TInput      &in;

        Pipe(TInput& _in):
            in(_in) {}

        inline TValue operator*() const {
            return TValue(*in);
        }

        Pipe& operator++() {
            ++in;
            return *this;
        }
    };
//}

}

#endif
