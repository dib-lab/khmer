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

#ifndef SEQAN_HEADER_PIPE_ITERATOR_H
#define SEQAN_HEADER_PIPE_ITERATOR_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

    //////////////////////////////////////////////////////////////////////////////
    // input pipe iterator
    // interface between pipe modules and algorithms working with iterators

    template <typename TInput>
    struct IPipeIterator
    {
        //////////////////////////////////////////////////////////////////////////////
        // public iterator interface

        typedef IPipeIterator                    iterator;
        typedef std::forward_iterator_tag        iterator_category;
        typedef typename Value<TInput>::Type    value_type;
        typedef const value_type                const_value_type;
        typedef typename Size<TInput>::Type        difference_type;
        typedef value_type*                        pointer;
        typedef value_type&                        reference;
        typedef const value_type&                const_reference;

        TInput*            in;
        difference_type    rest;

        IPipeIterator():
            in(NULL),
            rest(0) {}

        IPipeIterator(TInput &_in):
            in(&_in),
            rest(length(_in))
        {
            beginRead(*in);
        }

        IPipeIterator(const iterator &I):
            in(I.in),
            rest(I.rest) {}

        const_value_type operator* () const {
            return **in;
        }

        //const_reference operator* () const {
        //    return **in;
        //}

        iterator& operator++ () {
            ++(*in);
            if (!--rest) endRead(*in);
            return *this;
        }

        iterator operator++ (int) {
            iterator before = *this;
            ++*this;
            return before;
        }

        iterator operator+ (difference_type delta) {
            for(; delta != 0; --delta)
                ++*this;
            return *this;
        }

        difference_type operator- (const iterator &I) const {
            return I.rest - rest;
        }

        bool Equal_(const iterator &I) const {
            return rest == I.rest;
        }
    };

    template <typename TInput>
    bool operator==(const IPipeIterator<TInput>& _Left, const IPipeIterator<TInput>& Right_)
    {
        return _Left.Equal_(Right_);
    }

    template <typename TInput>
    bool operator!=(const IPipeIterator<TInput>& _Left, const IPipeIterator<TInput>& Right_)
    {
        return !(_Left == Right_);
    }



    //////////////////////////////////////////////////////////////////////////////
    // output pipe iterator
    // interface between algorithms working with iterators and pipe modules

    template <typename TOutput>
    struct OPipeIterator
    {
        //////////////////////////////////////////////////////////////////////////////
        // public iterator interface

        typedef OPipeIterator                    iterator;
        typedef std::forward_iterator_tag        iterator_category;
        typedef typename Value<TOutput>::Type    value_type;
        typedef typename Size<TOutput>::Type    difference_type;
        typedef iterator*                        pointer;
        typedef iterator&                        reference;
        typedef value_type const &                const_reference;

        TOutput*        out;
        difference_type    rest;

        OPipeIterator():
            out(NULL),
            rest(0) {}

        OPipeIterator(TOutput &_out):
            out(&_out),
            rest(length(_out))
        {
            beginWrite(*out);
        }

        OPipeIterator(const iterator &I):
            out(I.out),
            rest(I.rest) {}

        reference operator* () {
            return *this;
        }

        iterator operator= (const_reference _v) {
            out->push(_v);
            return *this;
        }

        iterator& operator++ () {
            ++(*out);
            if (!--rest) endWrite(*out);
            return *this;
        }

        iterator operator++ (int) {
            iterator before = *this;
            ++*this;
            return before;
        }

        difference_type operator- (const iterator &I) const {
            return I.rest - rest;
        }

        bool Equal_(const iterator &I) const {
            return rest == I.rest;
        }
    };

    template <typename TOutput>
    bool operator==(const OPipeIterator<TOutput>& _Left, const OPipeIterator<TOutput>& Right_)
    {
        return _Left.Equal_(Right_);
    }

    template <typename TOutput>
    bool operator!=(const OPipeIterator<TOutput>& _Left, const OPipeIterator<TOutput>& Right_)
    {
        return !(_Left == Right_);
    }


    template < typename TInput, typename TSpec, typename TIteratorSpec >
    struct Iterator< Pipe< TInput, TSpec >, TIteratorSpec> {
        typedef IPipeIterator< Pipe< TInput, TSpec > > Type;
    };

    template < typename TInput >
    struct Value< IPipeIterator< TInput > > {
        typedef typename Value<TInput>::Type Type;
    };

    template < typename TInput >
    struct Size< IPipeIterator< TInput > > {
        typedef typename Size<TInput>::Type Type;
    };

    template < typename TOutput >
    struct Value< OPipeIterator< TOutput > > {
        typedef typename Value<TOutput>::Type Type;
    };

    template < typename TOutput >
    struct Size< OPipeIterator< TOutput > > {
        typedef typename Size<TOutput>::Type Type;
    };


    template < typename TInput, typename TSpec >
    IPipeIterator< Pipe< TInput, TSpec > >
    begin(Pipe< TInput, TSpec > &pipe) {
        return IPipeIterator< Pipe< TInput, TSpec > >(pipe);
    }

    template < typename TInput, typename TSpec >
    IPipeIterator< Pipe< TInput, TSpec > >
    end(Pipe< TInput, TSpec > &/*pipe*/) {
        return IPipeIterator< Pipe< TInput, TSpec > >();
    }


    template < typename TInput >
    inline typename Difference<TInput>::Type
    difference(IPipeIterator<TInput> first, IPipeIterator<TInput> last) {
        return last - first;
    }

    template < typename TOutput >
    inline typename Difference<TOutput>::Type
    difference(OPipeIterator<TOutput> first, OPipeIterator<TOutput> last) {
        return last - first;
    }

//}

}

#endif
