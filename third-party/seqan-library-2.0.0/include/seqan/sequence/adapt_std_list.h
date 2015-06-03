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
// Adaptions for STL vectors to SeqAn sequences.
// ==========================================================================

// SEQAN_NO_GENERATED_FORWARDS: No forwards are generated for this file.

#ifndef SEQAN_SEQUENCE_ADAPT_STD_LIST_H_
#define SEQAN_SEQUENCE_ADAPT_STD_LIST_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// ===========================================================================
// Metafunctions
// ===========================================================================

template <typename TValue, typename TAlloc>
struct IsContiguous< std::list<TValue, TAlloc> >
{
    enum { VALUE = false };
};

template <typename  TValue, typename TAlloc>
struct IsContiguous< std::list<TValue, TAlloc> const>
        : IsContiguous< std::list<TValue, TAlloc> > {};

template <typename TValue, typename TAlloc>
struct Value< std::list<TValue, TAlloc> >
{
    typedef typename std::list<TValue, TAlloc>::value_type Type;
};

template <typename TValue, typename TAlloc>
struct Value< std::list<TValue, TAlloc> const>
        : Value< std::list<TValue, TAlloc> > {};

// TODO(holtgrew): GetValue is a reference?! I thought the reverse was true in respect to Value<>.
template <typename TValue, typename TAlloc>
struct GetValue< std::list<TValue, TAlloc> >
{
    typedef typename std::list<TValue, TAlloc>::reference Type;
};

template <typename TValue, typename TAlloc>
struct GetValue< std::list<TValue, TAlloc> const>
{
    typedef typename std::list<TValue, TAlloc>::const_reference Type;
};

template <typename TValue, typename TAlloc>
struct Reference< std::list<TValue, TAlloc> >
{
    typedef typename std::list<TValue, TAlloc>::reference Type;
};

template <typename TValue,  typename TAlloc>
struct Reference< std::list<TValue, TAlloc> const>
{
    typedef typename std::list<TValue,  TAlloc>::const_reference Type;
};

template <typename TValue, typename TAlloc>
struct Iterator< std::list<TValue, TAlloc>, Rooted>
{
    typedef std::list<TValue, TAlloc> TString_;
    typedef Iter<TString_, StdIteratorAdaptor> TIterator_;
    typedef Iter<TString_, AdaptorIterator<TIterator_> > Type;
};

template <typename TValue, typename TAlloc>
struct Iterator< std::list<TValue, TAlloc> const, Rooted>
{
    typedef std::list<TValue, TAlloc> const TString_;
    typedef Iter<TString_, StdIteratorAdaptor> TIterator_;
    typedef Iter<TString_, AdaptorIterator<TIterator_> > Type;
};

template <typename TValue, typename TAlloc>
struct Iterator< std::list<TValue, TAlloc>, Standard>
{
    typedef Iter< std::list<TValue, TAlloc>, StdIteratorAdaptor> Type;
};

template <typename TValue, typename TAlloc>
struct Iterator< std::list<TValue, TAlloc> const, Standard>
{
    typedef Iter< std::list<TValue, TAlloc> const, StdIteratorAdaptor> Type;
};

template <typename TValue, typename TAlloc>
struct Position< std::list<TValue, TAlloc> >
{
    typedef typename std::list<TValue, TAlloc>::size_type Type;
};

template <typename TValue, typename TAlloc>
struct Position< std::list<TValue, TAlloc> const>
        : Position< std::list<TValue, TAlloc> > {};

template <typename TValue, typename TAlloc>
struct Size< std::list<TValue, TAlloc> >
{
    typedef typename std::list<TValue, TAlloc>::size_type Type;
};

template <typename TValue, typename TAlloc>
struct Size< std::list<TValue, TAlloc> const>
        : Size< std::list<TValue, TAlloc> > {};

template <typename TValue, typename TAlloc>
struct StdContainerIterator< std::list<TValue, TAlloc> >
{
    typedef std::list<TValue, TAlloc> TContainer_;
    typedef typename TContainer_::iterator Type;
};

template <typename TValue, typename TAlloc>
struct StdContainerIterator< std::list<TValue, TAlloc> const>
{
    typedef std::list<TValue, TAlloc> TContainer_;
    typedef typename TContainer_::const_iterator Type;
};

// ===========================================================================
// Functions
// ===========================================================================

template<typename TValue>
inline
typename Iterator<std::list<TValue>, Standard>::Type
begin(std::list<TValue> & list,
      Standard const &)
{
    SEQAN_CHECKPOINT;
    return list.begin();
}

template <typename TValue>
inline
typename Iterator<std::list<TValue> const, Standard>::Type
begin(std::list<TValue> const & list,
      Standard const &)
{
    SEQAN_CHECKPOINT;
    return list.begin();
}

template<typename TValue>
inline
typename Iterator<std::list<TValue>, Standard>::Type
end(std::list<TValue> & list,
      Standard const &)
{
    SEQAN_CHECKPOINT;
    return list.end();
}

template <typename TValue>
inline
typename Iterator<std::list<TValue> const, Standard>::Type
end(std::list<TValue> const & list,
      Standard const &)
{
    SEQAN_CHECKPOINT;
    return list.end();
}

template <typename TValue>
inline typename Reference<std::list<TValue> >::Type
front(std::list<TValue> & list)
{
    SEQAN_CHECKPOINT;
    return list.front();
}

template <typename TValue>
inline typename Reference<std::list<TValue> const>::Type
front(std::list<TValue> const & list)
{
    SEQAN_CHECKPOINT;
    return list.front();
}

template <typename TValue>
inline typename Reference<std::list<TValue> >::Type
back(std::list<TValue> & list)
{
    SEQAN_CHECKPOINT;
    return list.back();
}

template <typename TValue>
inline typename Reference<std::list<TValue> const>::Type
back(std::list<TValue> const & list)
{
    SEQAN_CHECKPOINT;
    return list.back();
}

template <typename TValue>
inline typename Size<std::list<TValue> >::Type
length(std::list<TValue> & list)
{
    SEQAN_CHECKPOINT;
    return list.size();
}

template <typename TValue>
inline typename Size<std::list<TValue> const>::Type
length(std::list<TValue> const & list)
{
    SEQAN_CHECKPOINT;
    return list.size();
}

template <typename T, typename T2>
inline void
prependValue(std::list<T> & list,
             T2 value)
{
    SEQAN_CHECKPOINT;
    list.push_front(value);
}

template <typename T, typename T2>
inline void
appendValue(std::list<T> & list,
            T2 value)
{
    SEQAN_CHECKPOINT;
    list.push_back(value);
}

template <typename T>
inline void
clear(std::list<T> & list)
{
    SEQAN_CHECKPOINT;
    list.clear();
}

template <typename T, typename TTag>
inline void
reserve(std::list<T> & /*list*/, TTag const & /*tag*/)
{
    SEQAN_CHECKPOINT;
}

template <typename T>
inline typename Size<std::list<T> >::Type
capacity(std::list<T> const & list)
{
    SEQAN_CHECKPOINT;
    return length(list);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQUENCE_ADAPT_STD_LIST_H_
