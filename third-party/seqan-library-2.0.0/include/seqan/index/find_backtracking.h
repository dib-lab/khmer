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
// Approximate string matching via backtracking on VSTrees
// ==========================================================================

#ifndef SEQAN_FIND_BACKTRACKING_H_
#define SEQAN_FIND_BACKTRACKING_H_

//#define SEQAN_DEBUG

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TDistance, typename TSpec>
struct Backtracking;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TPrefix, typename TDistance>
struct BacktrackingState_ {};

// ============================================================================

template <typename TSuffix, typename TDistance>
struct SuffixAligner_
{};

template <typename TSuffix>
struct SuffixAligner_<TSuffix, HammingDistance>
{
    typedef typename Size<TSuffix>::Type    TSize;

    TSize               suffix_length;
    TSize               suffix_position;

    SuffixAligner_() :
        suffix_length(0),
        suffix_position(0)
    {}
};

// ============================================================================

template <typename TPrefix, typename TDistance>
struct PrefixAligner_
{};

template <typename TPrefix>
struct PrefixAligner_<TPrefix, HammingDistance>
{
    typedef typename Size<TPrefix>::Type    TSize;

    TSize               prefix_length;
    TSize               prefix_position;
    TSize               global_position;
    unsigned            errors;

    PrefixAligner_() :
        prefix_length(0),
        prefix_position(0),
        global_position(0),
        errors(0)
    {}
};

// ============================================================================

template <typename TText, typename TSpec, typename TDistance, typename TBacktrackingSpec>
class Finder<Index<TText, TSpec>, Backtracking<TDistance, TBacktrackingSpec> >
{
protected:
    typedef Index<TText, TSpec>                                         TIndex;
    typedef typename Iterator<TIndex, TopDown<> >::Type                 TIndexIterator;
    typedef String<TIndexIterator>                                      TParentStack;

    typedef typename Fibre<TIndex, FibreSA>::Type                       TSA;
    typedef typename Iterator<TSA const, Standard>::Type                TIterator;
    typedef typename Size<TIndex>::Type                                 TSize;

    typedef typename EdgeLabel<TIndexIterator>::Type                    TSuffix;
    typedef SuffixAligner_<TSuffix, TDistance>                           TSuffixAligner;

public:
    bool                index_iterator_at_root;

    Holder<TIndex>      index;
    TIndexIterator      index_iterator;
    TParentStack        index_parents;
    Pair<TSize>         index_range;

    TIterator           data_iterator;
    TSize               data_length;
    Pair<TIterator>     range;

    TSuffixAligner      suffix_aligner;

    Finder() {}

    Finder(TIndex & index) :
        index(index),
        index_iterator(index)
    {
        clear(*this);
    }

    Finder(TIndex const & index) :
        index(index),
        index_iterator(index)
    {
        clear(*this);
    }

    ~Finder()
    {
#ifdef SEQAN_DEBUG
        std::cout << "Finder Parents Height: " << length(this->index_parents) << std::endl;
#endif
    }

};

// ============================================================================

template <typename TNeedle, typename TDistance, typename TBacktrackingSpec>
class Pattern<TNeedle, Backtracking<TDistance, TBacktrackingSpec> >
{
protected:
    typedef TNeedle                                     TPrefix;
    typedef PrefixAligner_<TPrefix, TDistance>           TPrefixAligner;

    typedef typename BacktrackingState_<TPrefix, TDistance>::Type    TState;
    typedef String<TState>                              TStateStack;

public:
    Holder<TNeedle>         data_host;
    TPrefixAligner          prefix_aligner;
    TStateStack             state;

    Pattern() {}

    Pattern(TNeedle const & needle)
    {
        setHost(*this, needle);
    }

    ~Pattern()
    {
        // Empty stack
        clear(this->state);
    }

};

template <typename TNeedle, typename TSpec, typename TDistance, typename TBacktrackingSpec>
class Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> >
{
protected:
    typedef Index<TNeedle, TSpec>                                       TIndex;
    typedef typename Iterator<TIndex, TopDown<> >::Type                 TIndexIterator;
    typedef String<TIndexIterator>                                      TParentStack;

    typedef typename Fibre<TIndex, FibreSA>::Type                       TSA;
    typedef typename Iterator<TSA const, Standard>::Type                TIterator;
    typedef typename Size<TIndex>::Type                                 TSize;

    typedef typename EdgeLabel<TIndexIterator>::Type                    TPrefix;
    typedef PrefixAligner_<TPrefix, TDistance>                           TPrefixAligner;

    typedef typename BacktrackingState_<TPrefix, TDistance>::Type                    TState;
    typedef String<TState>                                              TStateStack;

    typedef String<bool>                                                TEndStack;

public:
    bool                index_iterator_at_root;

    Holder<TIndex>      data_host;
    TIndexIterator      index_iterator;
    TParentStack        index_parents;
    Pair<TSize>         index_range;

    TIterator           data_iterator;
    TSize               data_length;
    Pair<TIterator>     range;

    TPrefixAligner      prefix_aligner;
    TStateStack         state;
    TEndStack           atEnd;

    unsigned            exact;
    bool                search;

    // TODO(esiragusa): Remove depth, isLeaf(it) should return true when repLength(it) >= depth
    TSize               depth;

    Pattern() {}

    Pattern(TIndex & index, TSize depth) :
        data_host(index),
        index_iterator(index),
        depth(depth)
    {
        setHost(*this, index);
    }

    Pattern(TIndex const & index, TSize depth) :
        data_host(index),
        index_iterator(index),
        depth(depth)
    {
        setHost(*this, index);
    }

    ~Pattern()
    {
#ifdef SEQAN_DEBUG
        std::cout << "BacktrackingState_ Height: " << length(this->state) << std::endl;
        std::cout << "atEnd Height: " << length(this->atEnd) << std::endl;
        std::cout << "Pattern Parents Height: " << length(this->index_parents) << std::endl;
#endif

        // Empty stacks
        clear(this->index_parents);
        clear(this->state);
        clear(this->atEnd);
    }

};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TPrefix>
struct BacktrackingState_<TPrefix, HammingDistance>
{
    typedef typename Size<TPrefix>::Type            TSize;
    typedef Pair<TSize, unsigned>                   Type;
};

// ============================================================================

template <typename TNeedle, typename TSpec, typename TDistance, typename TBacktrackingSpec>
struct Position<Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > >:
    SAValue<Index<TNeedle, TSpec> >
{};

// ============================================================================
// Functions
// ============================================================================

template <typename TPrefix>
inline unsigned
getScore(PrefixAligner_<TPrefix, HammingDistance> const & me, bool atEnd)
{
    return atEnd ? me.errors : MaxValue<unsigned>::VALUE;
}

template <typename TPrefix>
inline unsigned
getMinScore(PrefixAligner_<TPrefix, HammingDistance> const & me)
{
    return me.errors;
}

template <typename TPrefix>
inline typename BacktrackingState_<TPrefix, HammingDistance>::Type
getInitialState(PrefixAligner_<TPrefix, HammingDistance> &)
{
    typedef typename BacktrackingState_<TPrefix, HammingDistance>::Type  TState;
    return TState(0, 0);
}

template <typename TPrefix>
inline typename BacktrackingState_<TPrefix, HammingDistance>::Type
getState(PrefixAligner_<TPrefix, HammingDistance> & me)
{
    typedef typename BacktrackingState_<TPrefix, HammingDistance>::Type  TState;
    return TState(me.global_position, me.errors);
}

template <typename TPrefix, typename TState>
inline void
setState(PrefixAligner_<TPrefix, HammingDistance> & me, TState const & state)
{
    me.global_position = state.i1;
    me.errors = state.i2;
#ifdef SEQAN_DEBUG
    std::cout << "Old BacktrackingState_:      " << "(" << state.i1 << ", " << state.i2 << ")" << std::endl;
#endif
}

template <typename TSuffix, typename TState>
inline void
setState(SuffixAligner_<TSuffix, HammingDistance> &, TState const &)
{
//    me.global_position = state.i1;
}

template <typename TSuffix, typename TSize>
inline void setLength(SuffixAligner_<TSuffix, HammingDistance> & me, TSize suffix_length)
{
    me.suffix_length = suffix_length;
}

template <typename TPrefix, typename TSize>
inline void setLength(PrefixAligner_<TPrefix, HammingDistance> & me, TSize prefix_length)
{
    me.prefix_length = prefix_length;
}

template <typename TSuffix, typename TSize>
inline void setPosition(SuffixAligner_<TSuffix, HammingDistance> & me, TSize suffix_position)
{
    me.suffix_position = suffix_position;
}

template <typename TPrefix, typename TSize>
inline void setPosition(PrefixAligner_<TPrefix, HammingDistance> & me, TSize prefix_position)
{
    me.prefix_position = prefix_position;
}

// ============================================================================

template <typename TSuffix, typename TPrefix, typename TErrors>
inline bool align(SuffixAligner_<TSuffix, HammingDistance> & suffix_aligner,
                  PrefixAligner_<TPrefix, HammingDistance> & prefix_aligner,
                  TSuffix & suffix,
                  TPrefix & prefix,
                  TErrors errors)
{
    return _align(suffix_aligner, prefix_aligner, suffix, prefix, errors,
                  typename IsSequence<TSuffix>::Type(), typename IsSequence<TPrefix>::Type());
}

template <typename TSuffix, typename TPrefix, typename TErrors>
inline bool _align(SuffixAligner_<TSuffix, HammingDistance> & suffix_aligner,
                   PrefixAligner_<TPrefix, HammingDistance> & prefix_aligner,
                   TSuffix & suffix,
                   TPrefix & prefix,
                   TErrors errors,
                   False const & /* tag */,
                   False const & /* tag */)
{
    typedef typename Size<TPrefix>::Type            TSize;

    TSize extension = _min(suffix_aligner.suffix_length - suffix_aligner.suffix_position,
                           prefix_aligner.prefix_length - prefix_aligner.prefix_position);

    if (extension > 0)
    {
        if (suffix != prefix)
            if (++prefix_aligner.errors > errors)
                return false;

        suffix_aligner.suffix_position++;
        prefix_aligner.prefix_position++;
        prefix_aligner.global_position++;
    }

    return true;
}

template <typename TSuffix, typename TPrefix, typename TErrors>
inline bool _align(SuffixAligner_<TSuffix, HammingDistance> & suffix_aligner,
                   PrefixAligner_<TPrefix, HammingDistance> & prefix_aligner,
                   TSuffix & suffix,
                   TPrefix & prefix,
                   TErrors errors,
                   True const & /* tag */,
                   False const & /* tag */)
{
    typedef typename Size<TPrefix>::Type            TSize;

    TSize extension = _min(suffix_aligner.suffix_length - suffix_aligner.suffix_position,
                           prefix_aligner.prefix_length - prefix_aligner.prefix_position);

    if (extension > 0)
    {
        if (suffix[suffix_aligner.suffix_position] != prefix)
            if (++prefix_aligner.errors > errors)
                return false;

        suffix_aligner.suffix_position++;
        prefix_aligner.prefix_position++;
        prefix_aligner.global_position++;
    }

    return true;
}

template <typename TSuffix, typename TPrefix, typename TErrors>
inline bool _align(SuffixAligner_<TSuffix, HammingDistance> & suffix_aligner,
                   PrefixAligner_<TPrefix, HammingDistance> & prefix_aligner,
                   TSuffix & suffix,
                   TPrefix & prefix,
                   TErrors errors,
                   False const & /* tag */,
                   True const & /* tag */)
{
    typedef typename Size<TPrefix>::Type            TSize;

    TSize extension = _min(suffix_aligner.suffix_length - suffix_aligner.suffix_position,
                           prefix_aligner.prefix_length - prefix_aligner.prefix_position);

    if (extension > 0)
    {
        if (suffix != prefix[prefix_aligner.prefix_position])
            if (++prefix_aligner.errors > errors)
                return false;

        suffix_aligner.suffix_position++;
        prefix_aligner.prefix_position++;
        prefix_aligner.global_position++;
    }

    return true;
}

template <typename TSuffix, typename TPrefix, typename TErrors>
inline bool _align(SuffixAligner_<TSuffix, HammingDistance> & suffix_aligner,
                   PrefixAligner_<TPrefix, HammingDistance> & prefix_aligner,
                   TSuffix & suffix,
                   TPrefix & prefix,
                   TErrors errors,
                   True const & /* tag */,
                   True const & /* tag */)
{
    typedef typename Size<TPrefix>::Type            TSize;

    TSize extension = _min(suffix_aligner.suffix_length - suffix_aligner.suffix_position,
                           prefix_aligner.prefix_length - prefix_aligner.prefix_position);

    TSize endPosition = prefix_aligner.global_position + extension;

    while (prefix_aligner.global_position < endPosition)
    {
        if (suffix[suffix_aligner.suffix_position] != prefix[prefix_aligner.prefix_position])
            if (++prefix_aligner.errors > errors)
                return false;

        suffix_aligner.suffix_position++;
        prefix_aligner.prefix_position++;
        prefix_aligner.global_position++;
    }

    return true;
}

// ============================================================================

template <typename TSuffix, typename TPrefix>
inline bool match(SuffixAligner_<TSuffix, HammingDistance> & suffix_aligner,
                  PrefixAligner_<TPrefix, HammingDistance> & prefix_aligner,
                  TSuffix const & suffix,
                  TPrefix const & prefix)
{
    return _match(suffix_aligner, prefix_aligner, suffix, prefix,
                  typename IsSequence<TSuffix>::Type(), typename IsSequence<TPrefix>::Type());
}

template <typename TSuffix, typename TPrefix>
inline bool _match(SuffixAligner_<TSuffix, HammingDistance> & suffix_aligner,
                   PrefixAligner_<TPrefix, HammingDistance> & prefix_aligner,
                   TSuffix const & suffix,
                   TPrefix const & prefix,
                   False const & /* tag */,
                   False const & /* tag */)
{
    typedef typename Size<TPrefix>::Type            TSize;

    TSize extension = _min(suffix_aligner.suffix_length - suffix_aligner.suffix_position,
                           prefix_aligner.prefix_length - prefix_aligner.prefix_position);

    if (extension > 0)
    {
        if (suffix != prefix)
            return false;

        suffix_aligner.suffix_position++;
        prefix_aligner.prefix_position++;
        prefix_aligner.global_position++;
    }

    return true;
}

template <typename TSuffix, typename TPrefix>
inline bool _match(SuffixAligner_<TSuffix, HammingDistance> & suffix_aligner,
                   PrefixAligner_<TPrefix, HammingDistance> & prefix_aligner,
                   TSuffix const & suffix,
                   TPrefix const & prefix,
                   True const & /* tag */,
                   False const & /* tag */)
{
    typedef typename Size<TPrefix>::Type            TSize;

    TSize extension = _min(suffix_aligner.suffix_length - suffix_aligner.suffix_position,
                           prefix_aligner.prefix_length - prefix_aligner.prefix_position);

    if (extension > 0)
    {
        if (suffix[suffix_aligner.suffix_position] != prefix)
            return false;

        suffix_aligner.suffix_position++;
        prefix_aligner.prefix_position++;
        prefix_aligner.global_position++;
    }

    return true;
}

template <typename TSuffix, typename TPrefix>
inline bool _match(SuffixAligner_<TSuffix, HammingDistance> & suffix_aligner,
                   PrefixAligner_<TPrefix, HammingDistance> & prefix_aligner,
                   TSuffix const & suffix,
                   TPrefix const & prefix,
                   False const & /* tag */,
                   True const & /* tag */)
{
    typedef typename Size<TPrefix>::Type            TSize;

    TSize extension = _min(suffix_aligner.suffix_length - suffix_aligner.suffix_position,
                           prefix_aligner.prefix_length - prefix_aligner.prefix_position);

    if (extension > 0)
    {
        if (suffix != prefix[prefix_aligner.prefix_position])
            return false;

        suffix_aligner.suffix_position++;
        prefix_aligner.prefix_position++;
        prefix_aligner.global_position++;
    }

    return true;
}

template <typename TSuffix, typename TPrefix>
inline bool _match(SuffixAligner_<TSuffix, HammingDistance> & suffix_aligner,
                   PrefixAligner_<TPrefix, HammingDistance> & prefix_aligner,
                   TSuffix const & suffix,
                   TPrefix const & prefix,
                   True const & /* tag */,
                   True const & /* tag */)
{
    typedef typename Size<TPrefix>::Type            TSize;

    TSize extension = _min(suffix_aligner.suffix_length - suffix_aligner.suffix_position,
                           prefix_aligner.prefix_length - prefix_aligner.prefix_position);

    TSize endPosition = prefix_aligner.global_position + extension;

    while (prefix_aligner.global_position < endPosition)
    {
        if (suffix[suffix_aligner.suffix_position] != prefix[prefix_aligner.prefix_position])
            return false;

        suffix_aligner.suffix_position++;
        prefix_aligner.prefix_position++;
        prefix_aligner.global_position++;
    }

    return true;
}

// ============================================================================

template <typename TPrefix>
inline bool _atEnd(PrefixAligner_<TPrefix, HammingDistance> const & me)
{
    return me.prefix_position == me.prefix_length;
}

template <typename TSuffix>
inline bool _atEnd(SuffixAligner_<TSuffix, HammingDistance> const & me)
{
    return me.suffix_position == me.suffix_length;
}

// ============================================================================

template <typename TText, typename TSpec, typename TDistance, typename TBacktrackingSpec>
inline void
clear(Finder<Index<TText, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & me)
{
    typedef Index<TText, TSpec>                                         TIndex;
    typedef typename Iterator<TIndex, TopDown<> >::Type                 TIndexIterator;

    typedef typename Fibre<TIndex, FibreSA>::Type                       TSA;
    typedef typename Iterator<TSA const, Standard>::Type                TIterator;

    typedef typename EdgeLabel<TIndexIterator>::Type                    TSuffix;
    typedef SuffixAligner_<TSuffix, TDistance>                           TSuffixAligner;

    // Init backtracking on root node
    me.index_iterator = TIndexIterator(host(me));
    me.index_iterator_at_root = true;

    // Empty parent stack
    clear(me.index_parents);

    // Init data iterator on empty range
    hostIterator(me) = begin(indexSA(host(me)), Standard());
    me.range.i1 = me.range.i2 = TIterator();
    me.data_length = 0;

    // Call SuffixAligner_ constructor
    me.suffix_aligner = TSuffixAligner();
}

// ============================================================================

template <typename TNeedle, typename TDistance, typename TBacktrackingSpec, typename TNewNeedle>
void setHost(Pattern<TNeedle, Backtracking<TDistance, TBacktrackingSpec> > & me, TNewNeedle const & needle)
{
    typedef TNeedle                                     TPrefix;
    typedef PrefixAligner_<TPrefix, TDistance>           TPrefixAligner;
    //typedef typename BacktrackingState_<TPrefix, TDistance>::Type    TState;

    SEQAN_ASSERT_NOT(empty(needle));
    setValue(me.data_host, needle);

    // Call PrefixAligner_ constructor
    me.prefix_aligner = TPrefixAligner();

    // Init stack on empty word
    clear(me.state);
    appendValue(me.state, getInitialState(me.prefix_aligner));
}

template <typename TNeedle, typename TSpec, typename TDistance, typename TBacktrackingSpec>
void _moveIteratorAtRoot(Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & me)
{
    typedef Index<TNeedle, TSpec>                                       TIndex;
    typedef typename Iterator<TIndex, TopDown<> >::Type                 TIndexIterator;
    me.index_iterator = TIndexIterator(host(me));
}

template <typename TNeedle, typename TSpec, typename TDistance, typename TBacktrackingSpec, typename TNewNeedle, typename TNewSpec>
void setHost(Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & me,
             Index<TNewNeedle, TNewSpec> const & index)
{
    typedef Index<TNeedle, TSpec>                                       TIndex;
    typedef typename Iterator< TIndex, TopDown<> >::Type                TIndexIterator;
    typedef typename Fibre<TIndex, FibreSA>::Type                       TSA;
    typedef typename Iterator<TSA const, Standard>::Type                TIterator;

    typedef typename EdgeLabel<TIndexIterator>::Type                    TPrefix;
    typedef PrefixAligner_<TPrefix, TDistance>                           TPrefixAligner;
    //typedef typename BacktrackingState_<TPrefix, TDistance>::Type                    TState;

    // TODO(esiragusa): Update index holder
    setValue(me.data_host, index);

    // Init backtracking on root node
    _moveIteratorAtRoot(me);
//    me.index_iterator = TIndexIterator(host(me));
    me.index_iterator_at_root = true;

    // Empty parent stack
    clear(me.index_parents);

    // Init data iterator on empty range
    hostIterator(me) = begin(indexSA(host(me)), Standard());
    me.range.i1 = me.range.i2 = TIterator();
    me.data_length = 0;

    // Init stack on empty word
    clear(me.state);
    appendValue(me.state, getInitialState(me.prefix_aligner));

    // Empty atEnd stack
    clear(me.atEnd);
//    appendValue(me.atEnd, true);

    // Empty exact search stack
    me.exact = 0;
    me.search = false;

    // Call PrefixAligner_ constructor
    me.prefix_aligner = TPrefixAligner();
}

template <typename TNeedle, typename TDistance, typename TBacktrackingSpec, typename TNewNeedle>
void setHost(Pattern<TNeedle, Backtracking<TDistance, TBacktrackingSpec> > & me, TNewNeedle & needle)
{
    setHost(me, reinterpret_cast<TNeedle const &>(needle));
}

// ============================================================================

template < typename TNeedle, typename TSpec, typename TDistance, typename TBacktrackingSpec >
inline typename Iterator< typename Fibre<Index<TNeedle, TSpec>, FibreSA>::Type const, Standard>::Type const &
hostIterator(Pattern< Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > const & me)
{
    return me.data_iterator;
}

template < typename TNeedle, typename TSpec, typename TDistance, typename TBacktrackingSpec >
inline typename Iterator< typename Fibre<Index<TNeedle, TSpec>, FibreSA>::Type const, Standard >::Type &
hostIterator(Pattern< Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & me)
{
    return me.data_iterator;
}

// ============================================================================

template <typename TNeedle, typename TSpec, typename TDistance, typename TBacktrackingSpec>
inline bool
empty(Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & me)
{
    return me.range.i1 == me.range.i2;
}

template <typename TNeedle, typename TSpec, typename TDistance, typename TBacktrackingSpec>
inline bool
atBegin(Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & me)
{
    return (empty(me) || hostIterator(me) == me.range.i1);
}

template <typename TNeedle, typename TSpec, typename TDistance, typename TBacktrackingSpec>
inline bool
atEnd(Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & me)
{
    return (empty(me) || hostIterator(me) == me.range.i2);
}

template <typename TNeedle, typename TSpec, typename TDistance, typename TBacktrackingSpec>
inline void
goBegin(Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & me)
{
    hostIterator(me) = me.range.i1;
}

template <typename TNeedle, typename TSpec, typename TDistance, typename TBacktrackingSpec>
inline void
goEnd(Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & me)
{
    hostIterator(me) = me.range.i2;
}

// ============================================================================

template <typename TNeedle, typename TSpec, typename TDistance, typename TBacktrackingSpec>
inline typename Position<Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > >::Type
beginPosition(Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & me)
{
    SEQAN_ASSERT_NOT(empty(me));
    return *me.data_iterator;
}

template <typename TNeedle, typename TSpec, typename TDistance, typename TBacktrackingSpec>
inline typename Position<Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > >::Type
endPosition(Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & me)
{
    return posAdd(beginPosition(me), me.data_length);
}

template <typename TNeedle, typename TSpec, typename TDistance, typename TBacktrackingSpec>
inline typename Position<Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > >::Type
position(Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & me)
{
    return beginPosition(me);
}

// ============================================================================

template <typename TText, typename TSpec, typename TDistance, typename TBacktrackingSpec, typename TNeedle>
inline bool
_resume(Finder<Index<TText, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & finder,
        Pattern<TNeedle, Backtracking<TDistance, TBacktrackingSpec> > & pattern)
{
    // Resume positively from root only the first time
    if (finder.index_iterator_at_root)
    {
        finder.index_iterator_at_root = false;
        return true;
    }

    return _cut(finder, pattern);
}

template <typename TText, typename TSpec, typename TDistance, typename TBacktrackingSpec, typename TNeedle, typename TErrors>
inline bool
_backtrack(Finder<Index<TText, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & finder,
           Pattern<TNeedle, Backtracking<TDistance, TBacktrackingSpec> > & pattern,
           TErrors errors)
{
    typedef Index<TText, TSpec>                             TIndex;
    typedef typename Iterator<TIndex, TopDown<> >::Type     TIndexIterator;
    typedef typename EdgeLabel<TIndexIterator>::Type        TSuffix;
    typedef typename Size<TIndex>::Type                     TSuffixSize;

    typedef typename BacktrackingState_<TNeedle, TDistance>::Type        TState;

    setLength(pattern.prefix_aligner, length(needle(pattern)));
    setPosition(pattern.prefix_aligner, 0);

    do
    {
        SEQAN_ASSERT_NOT(empty(pattern.state));

#ifdef SEQAN_DEBUG
        std::cout << "Stack Height:   " << length(pattern.state) << std::endl;
        std::cout << "Suffix:         " <<
        prefix(representative(finder.index_iterator),
               _min(repLength(finder.index_iterator), length(needle(pattern))))
        << std::endl;
#endif
        // Restore last state
        setState(finder.suffix_aligner, back(pattern.state));
        setState(pattern.prefix_aligner, back(pattern.state));

        // Update current prefix
        setPosition(pattern.prefix_aligner, pattern.prefix_aligner.global_position);

        // Update current suffix
        TSuffix suffixEdge = parentEdgeLabel(finder.index_iterator);
        TSuffixSize suffixLength = parentEdgeLength(finder.index_iterator);
        TSuffixSize suffixPos = pattern.prefix_aligner.global_position - parentRepLength(finder.index_iterator);
        setLength(finder.suffix_aligner, suffixLength);
        setPosition(finder.suffix_aligner, suffixPos);

        // Align suffix with pattern
        align(finder.suffix_aligner, pattern.prefix_aligner, suffixEdge, needle(pattern), errors);

        // A complete match was found
        if (getScore(pattern.prefix_aligner, _atEnd(pattern.prefix_aligner)) <= errors)
        {
            finder.index_range = range(finder.index_iterator);
            return true;
        }
        // Reduce to exact suffix search (speedup)
        else if (getMinScore(pattern.prefix_aligner) == errors)
        {
            TIndexIterator index_iterator(finder.index_iterator);

            // A complete match was found
            if (goDown(index_iterator, suffix(needle(pattern), pattern.prefix_aligner.global_position)))
            {
                // Move aligners to end of pattern
                pattern.prefix_aligner.global_position = length(needle(pattern));

                finder.index_range = range(index_iterator);
                return true;
            }
            else
                _cut(finder, pattern);
        }
        // Walk down text index only if an alignment is still possible
        else if (getMinScore(pattern.prefix_aligner) <= errors && !isLeaf(finder.index_iterator))
        {
            appendValue(finder.index_parents, finder.index_iterator);
            goDown(finder.index_iterator);
//            appendValue(pattern.state, getState(pattern.prefix_aligner));
            appendValue(pattern.state, TState(pattern.prefix_aligner.global_position, pattern.prefix_aligner.errors));
        }
        // Otherwise cut branch
        else
            _cut(finder, pattern);
    }
    while (!isRoot(finder.index_iterator));

    return false;
}

template <typename TText, typename TSpec, typename TDistance, typename TBacktrackingSpec, typename TNeedle>
inline bool
_cut(Finder<Index<TText, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & finder,
     Pattern<TNeedle, Backtracking<TDistance, TBacktrackingSpec> > & pattern)
{
    // Walk to next node
    if (!goRight(finder.index_iterator))
        while (!isRoot(finder.index_iterator))
        {
            SEQAN_ASSERT_NOT(empty(finder.index_parents));
            finder.index_iterator = back(finder.index_parents);
            eraseBack(finder.index_parents);

            SEQAN_ASSERT_NOT(empty(pattern.state));
            eraseBack(pattern.state);
            if (goRight(finder.index_iterator))
                break;
        }

    return !isRoot(finder.index_iterator);
}

// ============================================================================

template <typename TText, typename TTextSpec, typename TDistance, typename TBacktrackingSpec, typename TNeedle, typename TNeedleSpec, typename TErrors>
inline bool
_resume(Finder<Index<TText, TTextSpec>, Backtracking<TDistance, TBacktrackingSpec> > & finder,
        Pattern<Index<TNeedle, TNeedleSpec>, Backtracking<TDistance, TBacktrackingSpec> > & pattern,
        TErrors errors)
{
    // Resume positively from root only the first time
    if (finder.index_iterator_at_root)
    {
        finder.index_iterator_at_root = false;
        return _backtrack(finder, pattern, errors);
    }

    if (pattern.search)
    {
        if (_cut_exact(finder, pattern) && _search(finder, pattern))
            return true;

        SEQAN_ASSERT_NOT(pattern.exact);

        SEQAN_ASSERT_NOT(empty(finder.index_parents));
        finder.index_iterator = back(finder.index_parents);
        eraseBack(finder.index_parents);
        SEQAN_ASSERT_NOT(empty(pattern.index_parents));
        pattern.index_iterator = back(pattern.index_parents);
        eraseBack(pattern.index_parents);
        SEQAN_ASSERT_NOT(empty(pattern.state));
        eraseBack(pattern.state);

        pattern.search = false;
    }

    return _cut(finder, pattern) && _backtrack(finder, pattern, errors);
}

template <typename TText, typename TTextSpec, typename TDistance, typename TBacktrackingSpec, typename TNeedle, typename TNeedleSpec, typename TErrors>
inline bool
_backtrack(Finder<Index<TText, TTextSpec>, Backtracking<TDistance, TBacktrackingSpec> > & finder,
           Pattern<Index<TNeedle, TNeedleSpec>, Backtracking<TDistance, TBacktrackingSpec> > & pattern,
           TErrors errors)
{
    typedef Index<TText, TTextSpec>                                     TTextIndex;
    typedef typename Iterator<TTextIndex, TopDown<> >::Type             TTextIndexIterator;
    typedef typename EdgeLabel<TTextIndexIterator>::Type                TSuffix;
    typedef typename Size<TTextIndex>::Type                             TSuffixSize;

    typedef Index<TNeedle, TNeedleSpec>                                 TNeedleIndex;
    typedef typename Iterator<TNeedleIndex, TopDown<> >::Type           TNeedleIndexIterator;
    typedef typename EdgeLabel<TNeedleIndexIterator>::Type              TPrefix;
    //typedef typename Size<TNeedleIndex>::Type                           TPrefixSize;

    typedef typename BacktrackingState_<TNeedle, TDistance>::Type                    TState;

    SEQAN_ASSERT_NOT(pattern.exact);

    do
    {
        SEQAN_ASSERT_NOT(empty(pattern.state));

#ifdef SEQAN_DEBUG
        std::cout << "Stack Height:   " << length(pattern.state) << std::endl;
        std::cout << "Suffix:         "    << representative(finder.index_iterator) << std::endl;
        std::cout << "Prefix:         " << representative(pattern.index_iterator) << std::endl;
#endif

        // Lookup last state
        setState(finder.suffix_aligner, back(pattern.state));
        setState(pattern.prefix_aligner, back(pattern.state));

        // Update current suffix
        TSuffix suffixEdge = parentEdgeLabel(finder.index_iterator);
        TSuffixSize suffixLength = parentEdgeLength(finder.index_iterator);
        TSuffixSize suffixPos = pattern.prefix_aligner.global_position - parentRepLength(finder.index_iterator);
        setLength(finder.suffix_aligner, suffixLength);
        setPosition(finder.suffix_aligner, suffixPos);

        // Update current prefix
        TPrefix prefixEdge = parentEdgeLabel(pattern.index_iterator);
        TSuffixSize prefixLength = parentEdgeLength(pattern.index_iterator);
        TSuffixSize prefixPos = pattern.prefix_aligner.global_position - parentRepLength(pattern.index_iterator);
        setLength(pattern.prefix_aligner, _min(prefixLength, pattern.depth - parentRepLength(pattern.index_iterator)));
        setPosition(pattern.prefix_aligner, prefixPos);

        // Align suffix with prefix
        align(finder.suffix_aligner, pattern.prefix_aligner, suffixEdge, prefixEdge, errors);

        // A complete match was found
//        if (getScore(pattern.prefix_aligner, isLeaf(pattern.index_iterator)) <= errors)
        if (getScore(pattern.prefix_aligner, (pattern.prefix_aligner.global_position >= pattern.depth)) <= errors)
        {
            finder.index_range = range(finder.index_iterator);
            pattern.index_range = range(pattern.index_iterator);
            return true;
        }
        // Reduce to exact suffix search (speedup)
        else if (getMinScore(pattern.prefix_aligner) == errors)
        {
            pattern.search = true;

            appendValue(finder.index_parents, finder.index_iterator);
            appendValue(pattern.index_parents, pattern.index_iterator);
//            appendValue(pattern.state, getState(pattern.prefix_aligner));
            appendValue(pattern.state, TState(pattern.prefix_aligner.global_position, pattern.prefix_aligner.errors));

            if (_search(finder, pattern))
                return true;

            SEQAN_ASSERT_NOT(pattern.exact);

            SEQAN_ASSERT_NOT(empty(finder.index_parents));
            finder.index_iterator = back(finder.index_parents);
            eraseBack(finder.index_parents);
            SEQAN_ASSERT_NOT(empty(pattern.index_parents));
            pattern.index_iterator = back(pattern.index_parents);
            eraseBack(pattern.index_parents);
            SEQAN_ASSERT_NOT(empty(pattern.state));
            eraseBack(pattern.state);

            pattern.search = false;

            _cut(finder, pattern);
        }
        else if (getMinScore(pattern.prefix_aligner) <= errors)
        {
            if (_atEnd(finder.suffix_aligner))
            {
                if (!isLeaf(finder.index_iterator))
                {
                    appendValue(finder.index_parents, finder.index_iterator);
                    appendValue(pattern.index_parents, pattern.index_iterator);
//                    appendValue(pattern.state, getState(pattern.prefix_aligner));
                    appendValue(pattern.state, TState(pattern.prefix_aligner.global_position, pattern.prefix_aligner.errors));
                    appendValue(pattern.atEnd, false);

                    goDown(finder.index_iterator);
                }
                else
                    _cut(finder, pattern);
            }
            else if (_atEnd(pattern.prefix_aligner))
            {
                if (prefixLength < pattern.depth && !isLeaf(pattern.index_iterator))
                {
                    appendValue(finder.index_parents, finder.index_iterator);
                    appendValue(pattern.index_parents, pattern.index_iterator);
//                    appendValue(pattern.state, getState(pattern.prefix_aligner));
                    appendValue(pattern.state, TState(pattern.prefix_aligner.global_position, pattern.prefix_aligner.errors));
                    appendValue(pattern.atEnd, true);

                    goDown(pattern.index_iterator);
                }
                else
                    _cut(finder, pattern);
            }
        }
        else
            _cut(finder, pattern);
    }
    while (!(isRoot(pattern.index_iterator) && isRoot(finder.index_iterator)));

    return false;
}

template <typename TText, typename TTextSpec, typename TDistance, typename TBacktrackingSpec, typename TNeedle, typename TNeedleSpec>
inline bool
_cut(Finder<Index<TText, TTextSpec>, Backtracking<TDistance, TBacktrackingSpec> > & finder,
     Pattern<Index<TNeedle, TNeedleSpec>, Backtracking<TDistance, TBacktrackingSpec> > & pattern)
{
    SEQAN_ASSERT_NOT(pattern.exact);

    if (empty(pattern.atEnd))
        return false;

    do
    {
        SEQAN_ASSERT_NOT(empty(pattern.atEnd));

        if (!back(pattern.atEnd))
        {
            if (goRight(finder.index_iterator))
            {
                SEQAN_ASSERT_NOT(empty(pattern.index_parents));
                pattern.index_iterator = back(pattern.index_parents);
                break;
            }
        }
        else
        {
            if (goRight(pattern.index_iterator))
            {
                SEQAN_ASSERT_NOT(empty(finder.index_parents));
                finder.index_iterator = back(finder.index_parents);
                break;
            }
        }

        SEQAN_ASSERT_NOT(empty(finder.index_parents));
        finder.index_iterator = back(finder.index_parents);
        eraseBack(finder.index_parents);
        SEQAN_ASSERT_NOT(empty(pattern.index_parents));
        pattern.index_iterator = back(pattern.index_parents);
        eraseBack(pattern.index_parents);
        SEQAN_ASSERT_NOT(empty(pattern.state));
        eraseBack(pattern.state);
        SEQAN_ASSERT_NOT(empty(pattern.atEnd));
        eraseBack(pattern.atEnd);
    }
    while (!(isRoot(pattern.index_iterator) && isRoot(finder.index_iterator)));

    return !(isRoot(pattern.index_iterator) && isRoot(finder.index_iterator));
}

template <typename TText, typename TTextSpec, typename TDistance, typename TBacktrackingSpec, typename TNeedle, typename TNeedleSpec>
inline unsigned
_search(Finder<Index<TText, TTextSpec>, Backtracking<TDistance, TBacktrackingSpec> > & finder,
        Pattern<Index<TNeedle, TNeedleSpec>, Backtracking<TDistance, TBacktrackingSpec> > & pattern)
{
    typedef Index<TText, TTextSpec>                                     TTextIndex;
    typedef typename Iterator<TTextIndex, TopDown<> >::Type             TTextIndexIterator;
    typedef typename EdgeLabel<TTextIndexIterator>::Type                TSuffix;
    typedef typename Size<TTextIndex>::Type                             TSuffixSize;

    typedef Index<TNeedle, TNeedleSpec>                                 TNeedleIndex;
    typedef typename Iterator<TNeedleIndex, TopDown<> >::Type           TNeedleIndexIterator;
    typedef typename EdgeLabel<TNeedleIndexIterator>::Type              TPrefix;
    typedef typename Size<TNeedleIndex>::Type                           TPrefixSize;

    typedef typename BacktrackingState_<TNeedle, TDistance>::Type                    TState;

    do
    {
        SEQAN_ASSERT_NOT(empty(pattern.state));

#ifdef SEQAN_DEBUG
        std::cout << "Stack Height:   " << length(pattern.state) << std::endl;
        std::cout << "Exact Height:   " << pattern.exact << std::endl;
        std::cout << "Suffix:         "    << representative(finder.index_iterator) << std::endl;
        std::cout << "Prefix:         " << representative(pattern.index_iterator) << std::endl;
#endif

        // Lookup last state
        setState(finder.suffix_aligner, back(pattern.state));
        setState(pattern.prefix_aligner, back(pattern.state));

        // Update current suffix
        TSuffix suffixEdge = parentEdgeLabel(finder.index_iterator);
        TSuffixSize suffixLength = parentEdgeLength(finder.index_iterator);
        TSuffixSize suffixPos = pattern.prefix_aligner.global_position - parentRepLength(finder.index_iterator);
        setLength(finder.suffix_aligner, suffixLength);
        setPosition(finder.suffix_aligner, suffixPos);

        // Update current prefix
        TPrefix prefixEdge = parentEdgeLabel(pattern.index_iterator);
        TSuffixSize prefixLength = parentEdgeLength(pattern.index_iterator);
        TSuffixSize prefixPos = pattern.prefix_aligner.global_position - parentRepLength(pattern.index_iterator);
        setLength(pattern.prefix_aligner, _min(prefixLength, pattern.depth - parentRepLength(pattern.index_iterator)));
        setPosition(pattern.prefix_aligner, prefixPos);

        // Align exactly suffix with prefix
        if (match(finder.suffix_aligner, pattern.prefix_aligner, suffixEdge, prefixEdge))
        {
            if (!_atEnd(pattern.prefix_aligner) && _atEnd(finder.suffix_aligner))
            {
                TPrefixSize extension = pattern.prefix_aligner.prefix_length - pattern.prefix_aligner.prefix_position;

                TTextIndexIterator index_iterator(finder.index_iterator);
                // NOTE(esiragusa): This works if prefixEdge is a sequence, not if it is a simple type.
                if (!goDown(finder.index_iterator,
                            infixWithLength(prefixEdge, pattern.prefix_aligner.prefix_position, extension)))
                {
                    finder.index_iterator = index_iterator;
                    if (!_cut_exact(finder, pattern))
                        break;

                    continue;
                }

                pattern.prefix_aligner.global_position += extension;
                // NOTE(esiragusa): It is not necessary to update length and position as long as they are not used below.
//                suffixEdge = parentEdgeLabel(finder.index_iterator);
//                suffixLength = parentEdgeLength(finder.index_iterator);
//                suffixPos = pattern.prefix_aligner.global_position - parentRepLength(finder.index_iterator);
//                setLength(finder.suffix_aligner, suffixLength);
//                setPosition(finder.suffix_aligner, suffixPos);
            }

            // A complete match was found
            if (pattern.prefix_aligner.global_position >= pattern.depth)
            {
                finder.index_range = range(finder.index_iterator);
                pattern.index_range = range(pattern.index_iterator);
                return true;
            }
            else if (prefixLength < pattern.depth && !isLeaf(pattern.index_iterator))
            {
                ++pattern.exact;

                appendValue(finder.index_parents, finder.index_iterator);
                appendValue(pattern.index_parents, pattern.index_iterator);
//                appendValue(pattern.state, getState(pattern.prefix_aligner));
                appendValue(pattern.state, TState(pattern.prefix_aligner.global_position, pattern.prefix_aligner.errors));

                goDown(pattern.index_iterator);
            }
        }
        else if (!_cut_exact(finder, pattern))
            break;

    }
    while (!(isRoot(pattern.index_iterator) && isRoot(finder.index_iterator)));

    SEQAN_ASSERT_NOT(pattern.exact);

    return false;
}

template <typename TText, typename TTextSpec, typename TDistance, typename TBacktrackingSpec, typename TNeedle, typename TNeedleSpec>
inline bool
_cut_exact(Finder<Index<TText, TTextSpec>, Backtracking<TDistance, TBacktrackingSpec> > & finder,
           Pattern<Index<TNeedle, TNeedleSpec>, Backtracking<TDistance, TBacktrackingSpec> > & pattern)
{
    do
    {
        if (!pattern.exact)
            return false;

        if (goRight(pattern.index_iterator))
        {
            SEQAN_ASSERT_NOT(empty(finder.index_parents));
            finder.index_iterator = back(finder.index_parents);
            break;
        }

        SEQAN_ASSERT_NOT(empty(finder.index_parents));
        finder.index_iterator = back(finder.index_parents);
        eraseBack(finder.index_parents);
        SEQAN_ASSERT_NOT(empty(pattern.index_parents));
        pattern.index_iterator = back(pattern.index_parents);
        eraseBack(pattern.index_parents);
        SEQAN_ASSERT_NOT(empty(pattern.state));
        eraseBack(pattern.state);

        --pattern.exact;
    }
    while (!(isRoot(pattern.index_iterator) && isRoot(finder.index_iterator)));

    return !(isRoot(pattern.index_iterator) && isRoot(finder.index_iterator));
}

// ============================================================================

template <typename TText, typename TSpec, typename TDistance, typename TBacktrackingSpec, typename TNeedle, typename TErrors>
inline bool
find(Finder<Index<TText, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & finder,
     Pattern<TNeedle, Backtracking<TDistance, TBacktrackingSpec> > & pattern,
     TErrors errors)
{

    // Try to get another match from current node
    if (!atEnd(finder))
        goNext(hostIterator(finder));

    // Try to get more matches by backtracking some more
    if (atEnd(finder))
    {
        // Resume from last matching node and backtrack until next match
        if (_resume(finder, pattern) && _backtrack(finder, pattern, errors))
        {
            // Set data iterator range to the interval containing matches
            hostIterator(finder) = begin(indexSA(host(finder)), Standard());
            finder.range.i1 = hostIterator(finder) + finder.index_range.i1;
            finder.range.i2 = hostIterator(finder) + finder.index_range.i2;
            hostIterator(finder) = finder.range.i1;

            // Set match length
            _setFinderLength(finder, pattern.prefix_aligner.global_position);
        }
        // No more matches
        else
        {
            hostIterator(finder) = begin(indexSA(host(finder)), Standard());
            finder.range.i1 = hostIterator(finder);
            finder.range.i2 = hostIterator(finder);
        }
    }

    return !atEnd(finder);
}

template <typename TText, typename TTextSpec, typename TDistance, typename TBacktrackingSpec, typename TNeedle, typename TNeedleSpec, typename TErrors>
inline bool
find(Finder<Index<TText, TTextSpec>, Backtracking<TDistance, TBacktrackingSpec> > & finder,
     Pattern<Index<TNeedle, TNeedleSpec>, Backtracking<TDistance, TBacktrackingSpec> > & pattern,
     TErrors errors)
{

    // Try to get another match from current pattern node
    goNext(hostIterator(pattern));

    // Try to get another match from current text node
    if (atEnd(pattern))
    {
        goNext(hostIterator(finder));

        if (!atEnd(finder))
        {
            // Set data iterator range to the interval containing matches
            hostIterator(pattern) = begin(indexSA(host(pattern)), Standard());
            pattern.range.i1 = hostIterator(pattern) + pattern.index_range.i1;
            pattern.range.i2 = hostIterator(pattern) + pattern.index_range.i2;
            hostIterator(pattern) = pattern.range.i1;
        }
        // Try to get more matches by backtracking some more
        else
        {
            // Resume from last matching node and backtrack until next match
//                if (_resume(finder, pattern) && _backtrack(finder, pattern, errors))
            if (_resume(finder, pattern, errors))
            {
                // Set data iterator range to the interval containing matches
                hostIterator(finder) = begin(indexSA(host(finder)), Standard());
                finder.range.i1 = hostIterator(finder) + finder.index_range.i1;
                finder.range.i2 = hostIterator(finder) + finder.index_range.i2;
                hostIterator(finder) = finder.range.i1;

                hostIterator(pattern) = begin(indexSA(host(pattern)), Standard());
                pattern.range.i1 = hostIterator(pattern) + pattern.index_range.i1;
                pattern.range.i2 = hostIterator(pattern) + pattern.index_range.i2;
                hostIterator(pattern) = pattern.range.i1;

                // Set match length
                _setFinderLength(finder, pattern.prefix_aligner.global_position);
                pattern.data_length = pattern.prefix_aligner.global_position;
            }
            // No more matches
            else
            {
                hostIterator(finder) = begin(indexSA(host(finder)), Standard());
                finder.range.i1 = hostIterator(finder);
                finder.range.i2 = hostIterator(finder);

                hostIterator(pattern) = begin(indexSA(host(pattern)), Standard());
                pattern.range.i1 = hostIterator(pattern);
                pattern.range.i2 = hostIterator(pattern);
            }
        }

    }

    return !(atEnd(finder) && atEnd(pattern));
}

//template <typename TText, typename TSpec, typename TDistance, typename TBacktrackingSpec, typename TNeedle>
//inline bool
//find(Finder<Index<TText, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & finder,
//     Pattern<TNeedle, Backtracking<TDistance, TBacktrackingSpec> > & pattern)
//{
//    return find(finder, pattern, 0u);
//}

}

#endif  // #ifndef SEQAN_FIND_BACKTRACKING_H_
