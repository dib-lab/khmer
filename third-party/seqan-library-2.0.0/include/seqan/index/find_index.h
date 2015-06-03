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

#ifndef SEQAN_HEADER_INDEX_FIND_H
#define SEQAN_HEADER_INDEX_FIND_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

    template < typename TText, typename TSpec, typename TSpecFinder >
    struct Position< Finder< Index<TText, TSpec>, TSpecFinder > >:
        SAValue< Index<TText, TSpec> > {};


//////////////////////////////////////////////////////////////////////////////
// generic Finder class for all indices containing a suffix array or
// similar table, where a query result is an interval in this table
//
// your index must specialize the function _findFirstIndex and set
// finder.range to the interval containing the query hits. See:
//
//    template < typename TText, typename TSpec, typename TSpecFinder, typename TPattern >
//    inline void _findFirstIndex(
//        Finder< Index<TText, TSpec>, TSpecFinder > &finder,
//        TPattern const &pattern,
//        FinderMlr const)
//    {
//        Index<TText, TSpec> &index = haystack(finder);
//        indexRequire(index, EsaSA());
//        finder.range = equalRangeSAIterator(indexText(index), indexSA(index), pattern);
//    }

    template < typename TText, typename TSpec, typename TSpecFinder >
    class Finder< Index<TText, TSpec>, TSpecFinder >
    {
    protected:
        typedef Index<TText, TSpec>                                TIndex;
        typedef typename Fibre<TIndex, FibreSA>::Type            TSA;
        typedef typename Iterator<TSA const, Standard>::Type    TIterator;
        typedef typename Size<TIndex>::Type                        TSize;

    public:
        Holder<TIndex>    index;
        Pair<TIterator>    range;
        TIterator        data_iterator;
        TSize            data_length;

        Finder() : data_iterator(), data_length(0)
        {
            clear(*this);
        }
        Finder(TIndex &_index): index(_index), data_iterator(), data_length(0)
        {
            clear(*this);
        }
        Finder(TIndex const &_index): index(_index), data_iterator(), data_length(0)
        {
            clear(*this);
        }
    };

//____________________________________________________________________________

    template < typename TText, typename TSpec, typename TSpecFinder >
    inline typename Parameter_< Index<TText, TSpec> >::Type
    host(Finder< Index<TText, TSpec>, TSpecFinder > & me)
    {
SEQAN_CHECKPOINT
        return value(me.index);
    }

    template < typename TText, typename TSpec, typename TSpecFinder >
    inline typename Parameter_< Index<TText, TSpec> >::Type
    host(Finder< Index<TText, TSpec>, TSpecFinder > const & me)
    {
SEQAN_CHECKPOINT
        return value(me.index);
    }

    template < typename TText, typename TSpec, typename TSpecFinder >
    inline typename Parameter_< Index<TText, TSpec> >::Type
    container(Finder< Index<TText, TSpec>, TSpecFinder > & me)
    {
SEQAN_CHECKPOINT
        return value(me.index);
    }

    template < typename TText, typename TSpec, typename TSpecFinder >
    inline typename Parameter_< Index<TText, TSpec> >::Type
    container(Finder< Index<TText, TSpec>, TSpecFinder > const & me)
    {
SEQAN_CHECKPOINT
        return value(me.index);
    }

//____________________________________________________________________________

    template < typename TText, typename TSpec, typename TSpecFinder >
    inline void
    setHost(
        Finder< Index<TText, TSpec>, TSpecFinder > & me,
        typename Parameter_<Index<TText, TSpec> >::Type container_)
    {
        me.index = container_;
    }

    template < typename TText, typename TSpec, typename TSpecFinder >
    inline void
    setContainer(
        Finder< Index<TText, TSpec>, TSpecFinder > & me,
        typename Parameter_<Index<TText, TSpec> >::Type container_)
    {
        me.index = container_;
    }

//____________________________________________________________________________

    template < typename TText, typename TSpec, typename TSpecFinder >
    inline typename Iterator< typename Fibre<Index<TText, TSpec>, FibreSA>::Type const, Standard>::Type &
    hostIterator(Finder< Index<TText, TSpec>, TSpecFinder > & me)
    {
SEQAN_CHECKPOINT
        return me.data_iterator;
    }

    template < typename TText, typename TSpec, typename TSpecFinder >
    inline typename Iterator< typename Fibre<Index<TText, TSpec>, FibreSA>::Type const, Standard>::Type const &
    hostIterator(Finder< Index<TText, TSpec>, TSpecFinder > const & me)
    {
SEQAN_CHECKPOINT
        return me.data_iterator;
    }


//____________________________________________________________________________

    template < typename TText, typename TSpec, typename TSpecFinder >
    inline bool
    empty(Finder< Index<TText, TSpec>, TSpecFinder > & me)
    {
SEQAN_CHECKPOINT
        return me.range.i1 == me.range.i2;
    }

    template < typename TText, typename TSpec, typename TSpecFinder >
    inline void
    clear(Finder< Index<TText, TSpec>, TSpecFinder > & me)
    {
SEQAN_CHECKPOINT
        typedef Index<TText, TSpec>                                TIndex;
        typedef typename Fibre<TIndex, FibreSA>::Type            TSA;
        typedef typename Iterator<TSA const, Standard>::Type    TIterator;
        me.range.i1 = me.range.i2 = TIterator();
        me.data_length = 0;
    }

//____________________________________________________________________________

    template < typename TText, typename TSpec, typename TSpecFinder >
    inline bool
    atBegin(Finder< Index<TText, TSpec>, TSpecFinder > & me)
    {
SEQAN_CHECKPOINT
        return (empty(me) || hostIterator(me) == me.range.i1);
    }

    template < typename TText, typename TSpec, typename TSpecFinder >
    inline bool
    atEnd(Finder< Index<TText, TSpec>, TSpecFinder > & me)
    {
SEQAN_CHECKPOINT
        return (empty(me) || hostIterator(me) == me.range.i2);
    }

//____________________________________________________________________________

    template < typename TText, typename TSpec, typename TSpecFinder >
    inline void
    goBegin(Finder< Index<TText, TSpec>, TSpecFinder > & me)
    {
SEQAN_CHECKPOINT
        hostIterator(me) = me.range.i1;
    }

    template < typename TText, typename TSpec, typename TSpecFinder >
    inline void
    goEnd(Finder< Index<TText, TSpec>, TSpecFinder > & me)
    {
SEQAN_CHECKPOINT
        hostIterator(me) = me.range.i2;
    }

//____________________________________________________________________________

    template < typename TText, typename TSpec, typename TSpecFinder >
    inline typename Position< Finder< Index<TText, TSpec>, TSpecFinder > >::Type
    beginPosition(Finder< Index<TText, TSpec>, TSpecFinder > & me)
    {
        SEQAN_CHECKPOINT;
        SEQAN_ASSERT_NOT(empty(me));
        return *me.data_iterator;
    }

    template < typename TText, typename TSpec, typename TSpecFinder >
    inline typename Position< Finder< Index<TText, TSpec>, TSpecFinder > >::Type
    beginPosition(Finder< Index<TText, TSpec>, TSpecFinder > const & me)
    {
        SEQAN_CHECKPOINT;
        SEQAN_ASSERT_NOT(empty(me));
        return *me.data_iterator;
    }

//____________________________________________________________________________

    template < typename TText, typename TSpec, typename TSpecFinder >
    inline typename Position< Finder< Index<TText, TSpec>, TSpecFinder > >::Type
    endPosition(Finder< Index<TText, TSpec>, TSpecFinder > & me)
    {
SEQAN_CHECKPOINT
        return posAdd(beginPosition(me), me.data_length);
    }

    template < typename TText, typename TSpec, typename TSpecFinder >
    inline typename Position< Finder< Index<TText, TSpec>, TSpecFinder > >::Type
    endPosition(Finder< Index<TText, TSpec>, TSpecFinder > const & me)
    {
SEQAN_CHECKPOINT
        return posAdd(beginPosition(me), me.data_length);
    }

//____________________________________________________________________________

    template < typename TText, typename TSpec, typename TSpecFinder >
    inline typename Position< Finder< Index<TText, TSpec>, TSpecFinder > >::Type
    position(Finder< Index<TText, TSpec>, TSpecFinder > & me)
    {
SEQAN_CHECKPOINT
        return beginPosition(me);
    }

    template < typename TText, typename TSpec, typename TSpecFinder >
    inline typename Position< Finder< Index<TText, TSpec>, TSpecFinder > >::Type
    position(Finder< Index<TText, TSpec>, TSpecFinder > const & me)
    {
SEQAN_CHECKPOINT
        return beginPosition(me);
    }

//////////////////////////////////////////////////////////////////////////////
// find

    template < typename TText, typename TSpec, typename TSpecFinder, typename TPattern >
    inline bool find(
        Finder<Index<TText, TSpec>, TSpecFinder> &finder,
        TPattern const &pattern)
    {
        if (empty(finder))
        {
            _findFirstIndex(finder, needle(pattern), TSpecFinder());
            _setFinderLength(finder, length(needle(pattern)));
            hostIterator(finder) = finder.range.i1;
        } else
            ++hostIterator(finder);
        return !atEnd(finder);
    }

    template < typename TText, typename TSpec, typename TSpecFinder >
    inline bool find(Finder<Index<TText, TSpec>, TSpecFinder> &finder)
    {
        if (empty(finder)) return false;
        ++hostIterator(finder);
        return !atEnd(finder);
    }

}

#endif

