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
// Thread-safe / lock-free sequence operations.
// ==========================================================================

#ifndef SEQAN_PARALLEL_PARALLEL_SEQUENCE_H_
#define SEQAN_PARALLEL_PARALLEL_SEQUENCE_H_

namespace seqan {

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class ConcurrentAppender
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec = void>
struct ConcurrentAppender
{
    TString *       data_host;
    ReadWriteLock   lock;

    ConcurrentAppender() :
        data_host()
    {}

    ConcurrentAppender(TString & string)
    {
        setHost(*this, string);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec>
struct Value<ConcurrentAppender<TString, TSpec> > : Value<TString> {};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function host()
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec>
inline TString &
host(ConcurrentAppender<TString, TSpec> & me)
{
    return *me.data_host;
}

template <typename TString, typename TSpec>
inline TString const &
host(ConcurrentAppender<TString, TSpec> const & me)
{
    return *me.data_host;
}

// ----------------------------------------------------------------------------
// Function setHost()
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec>
inline void
setHost(ConcurrentAppender<TString, TSpec> & me, TString & string)
{
    me.data_host = &string;
}

// ----------------------------------------------------------------------------
// Function _incLength()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TParallel>
inline typename Size<String<TValue, Alloc<TSpec> > >::Type
_incLength(String<TValue, Alloc<TSpec> > & me, Tag<TParallel> const & tag)
{
    return atomicInc(me.data_end, tag) - begin(me, Standard());
}

// ----------------------------------------------------------------------------
// Function _decLength()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TParallel>
inline typename Size<String<TValue, Alloc<TSpec> > >::Type
_decLength(String<TValue, Alloc<TSpec> > & me, Tag<TParallel> const & tag)
{
    return atomicDec(me.data_end, tag) - begin(me, Standard());
}

// ----------------------------------------------------------------------------
// Function appendValue()
// ----------------------------------------------------------------------------

template <typename TTargetValue, typename TTargetSpec, typename TValue>
inline void
appendValue(String<TTargetValue, TTargetSpec> & me, TValue const & _value, Insist, Parallel)
{
    valueConstruct(begin(me, Standard()) + _incLength(me, Parallel()) - 1, _value);
}

template <typename TString, typename TSpec, typename TValue, typename TExpand, typename TParallel>
inline void
appendValue(ConcurrentAppender<TString, TSpec> & me,
            TValue const & val,
            Tag<TExpand> const & expandTag,
            Tag<TParallel> const & parallelTag)
{
    appendValue(host(me), val, expandTag, parallelTag);
}

template <typename TString, typename TSpec, typename TValue, typename TExpand>
inline void
appendValue(ConcurrentAppender<TString, TSpec> & me,
            TValue const & val,
            Tag<TExpand> const & expandTag,
            Parallel)
{
    typedef typename Size<TString>::Type TSize;

    TString & string = host(me);

    while (true)
    {
        // try to append the value
        {
            ScopedReadLock<> readLock(me.lock);
            TSize newLen = _incLength(string, Parallel());
            if (newLen <= capacity(string))
            {
                valueConstruct(begin(string, Standard()) + newLen - 1, val);
                break;
            }
            _decLength(string, Parallel());
        }

        // try to extend capacity
        {
            ScopedWriteLock<> writeLock(me.lock);
            TSize cap = capacity(string);
            if (cap == length(string))
                reserve(string, cap + 1, expandTag);
        }
    }
}

}  // namespace seqan

#endif  // #ifndef SEQAN_PARALLEL_PARALLEL_SEQUENCE_H_
