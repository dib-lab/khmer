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

#ifndef SEQAN_HEADER_INDEX_ESA_ALGS_H
#define SEQAN_HEADER_INDEX_ESA_ALGS_H

namespace SEQAN_NAMESPACE_MAIN
{

    //////////////////////////////////////////////////////////////////////////////
    // more sophisticated algorithms on enhanced suffix arrays of 1 sequence
    // (also called virtual suffix trees)
    //////////////////////////////////////////////////////////////////////////////


/*!
 * @class SuperMaxRepeatsIterator Super Max Repeats Iterator
 * @extends BottomUpIterator
 * @headerfile <seqan/index.h>
 *
 * @brief Iterator to search for all supermaximal repeats.
 *
 * @signature Iterator<TContainer, SuperMaxRepeats>::Type;
 * @signature template <typename TContainer>
 *            class Iter<TContainer, VSTree<BottomUp<SuperMaxRepeats> > >;
 *
 * @tparam TContainer Type of an index that can be iterated with a bottom-up
 *                    iterator. Types: @link IndexEsa @endlink
 *
 *
 * @note Note Instead of using the class Iter directly we recommend to use the result of the metafunction
 *               Iterator&lt;TContainer, SuperMaxRepeats&gt;::Type (which is Iter&lt;TContainer, VSTree&lt;BottomUp&lt;SuperMaxRepeats&gt; &g;t &gt;).
 *
 * @see DemoSupermaximalRepeats
 */

/*!
 * @fn SuperMaxRepeatsIterator::Iter
 * @brief The constructor
 *
 * @signature Iter::Iter(index[, minLength]);
 * @signature Iter::Iter(iterator);
 *
 * @param[in] index     The index to be used for the iteration. Types: @link IndexEsa @endlink
 * @param[in] minLength Minimum length of the supermaximal repeats, default value is 1.
 * @param[in] iterator  Another SuperMaxRepeats iterator. Types: @link SuperMaxRepeatsIterator @endlink
 */
    //////////////////////////////////////////////////////////////////////////////
    // super-maximal repeats - suffix tree version
    //////////////////////////////////////////////////////////////////////////////

    template < typename TSTree >
    struct GetVSTreeIteratorTraits< Iter< TSTree, VSTree< BottomUp<SuperMaxRepeats> > > > {
        typedef PostorderEmptyEdges    Type;
    };

    template < typename TSTree >
    class Iter< TSTree, VSTree< BottomUp<SuperMaxRepeats> > >:
        public Iter< TSTree, VSTree< BottomUp<> > >
    {
    public:
        typedef Iter< TSTree, VSTree< BottomUp<> > >    TBase;
        typedef typename Size<TSTree>::Type                TSize;
        typedef typename Value<TSTree>::Type            TValue;
//____________________________________________________________________________

        TSize                        minLength;
        typename Set<TValue>::Type    charSet;
//____________________________________________________________________________

        Iter() :
            TBase(),
            minLength(0)
        {}

        Iter(TSTree &_tree):
            TBase(_tree),
            minLength(1)
        {
            indexRequire(_tree, EsaChildtab());
            indexRequire(_tree, EsaBwt());
            goNext(*this);    // the iterator starts in a suffix, i.e. not a repeat
        }

        Iter(TSTree &_tree, MinimalCtor):
            TBase(_tree, MinimalCtor()) {}

        Iter(TSTree &_tree, TSize _minLength):
            TBase(_tree),
            minLength(_minLength)
        {
            indexRequire(_tree, EsaChildtab());
            indexRequire(_tree, EsaBwt());
            goNext(*this);    // the iterator starts in a suffix, i.e. not a repeat
        }

        Iter(Iter const &_origin):
            TBase((TBase const &)_origin),
            minLength(_origin.minLength),
            charSet(_origin.charSet) {}
    };

    template < typename TSTree >
    inline void goNext(Iter< TSTree, VSTree< BottomUp<SuperMaxRepeats> > > &it) {
        do {
            goNext(it, PostorderEmptyEdges());
        } while (!atEnd(it) &&
                 !(    childrenAreLeaves(it) &&
                    (repLength(it) >= it.minLength) &&
                    !isPartiallyLeftExtensible(it, it.charSet)) );
    }


/*!
 * @class SuperMaxRepeatsFastIterator Super Max Repeats Fast Iterator
 * @extends BottomUpIterator
 * @headerfile <seqan/index.h>
 *
 * @brief Iterator to search for all supermaximal repeats (for enh. suffix arrays only).
 *
 * @signature Iterator<TContainer, SuperMaxRepeatsFast>::Type;
 * @signature template <typename TContainer>
 *            class Iter<TContainer, VSTree<BottomUp<SuperMaxRepeatsFast> > >;
 *
 * @tparam TContainer Type of an index based on enhanced suffix array. Types:
 *                    @link IndexEsa @endlink
 *
 * @note Instead of using the class Iter directly we recommend to use the result of the metafunction
 *       Iterator&lt;TContainer, SuperMaxRepeatsFast&gt;::Type (which is Iter<TContainer, VSTree< BottomUp<SuperMaxRepeatsFast&gt; &gt; &gt;).
 */

/*!
 * @fn SuperMaxRepeatsFastIterator::Iter
 * @brief The constructor
 *
 * @signature Iter::Iter(index[, minLength]);
 * @signature Iter::Iter(iterator);
 *
 * @param[in] index     The index to be used for the iteration. Types: @link IndexEsa @endlink
 * @param[in] minLength Minimum length of the supermaximal repeats, default value is 1.
 * @param[in] iterator  Another SuperMaxRepeatsFast iterator. Types: @link SuperMaxRepeatsFastIterator @endlink
 */

    //////////////////////////////////////////////////////////////////////////////
    // supermaximal repeats - specialized for Enhanced Suffix Arrays
    //////////////////////////////////////////////////////////////////////////////

    template < typename TSTree >
    struct GetVSTreeIteratorTraits< Iter< TSTree, VSTree< BottomUp<SuperMaxRepeatsFast> > > > {
        typedef PostorderEmptyEdges    Type;
    };

    template < typename TText, typename TSpec >
    class Iter< Index<TText, IndexEsa<TSpec> >, VSTree< BottomUp<SuperMaxRepeatsFast> > >:
        public Iter< Index<TText, IndexEsa<TSpec> >, VSTree< BottomUp<> > >
    {
    public:
        typedef Index< TText, IndexEsa<TSpec> >        TIndex;
        typedef Iter< TIndex, VSTree< BottomUp<> > >    TBase;
        typedef typename Size<TIndex>::Type                TSize;
        typedef typename Value<TIndex>::Type            TValue;

        typedef typename Iterator<typename Fibre<TIndex, EsaLcp>::Type const>::Type    TLCPIter;
//____________________________________________________________________________

        TSize        minLength;
        TLCPIter    lIter, lEnd;    // lcp table iterators (optimization)
        TSize        lValueLast;        // current l-value of interval
        bool         rising;            // is the left interval border valid
        typename Set<TValue>::Type    charSet;
//____________________________________________________________________________

        Iter() :
            TBase(),
            minLength(0),
            lValueLast(0),
            rising(false)
        {}

        Iter(TIndex &_index):
            TBase(_index),
            minLength(1),
            lValueLast(0),
            rising(true)
        {
            this->vDesc.range = Pair<TSize>(0,0);
            indexRequire(_index, EsaBwt());
            lIter = begin(indexLcp(container(*this)));
            lEnd  = end(indexLcp(container(*this)));
            goNext(*this);
        }

        Iter(TIndex &_index, MinimalCtor):
            TBase(_index, MinimalCtor()) {}

        Iter(TIndex &_index, TSize _minLength):
            TBase(_index),
            minLength(_minLength),
            lValueLast(0),
            rising(true)
        {
            this->vDesc.i1 = Pair<TSize>(0,0);
            indexRequire(_index, EsaBwt());
            lIter = begin(indexLcp(container(*this)));
            lEnd  = end(indexLcp(container(*this)));
            goNext(*this);
        }

        Iter(Iter const &_origin):
            TBase((TBase const &)_origin),
            minLength(_origin.minLength),
            lIter(_origin.lIter),
            lEnd(_origin.lEnd),
            lValueLast(_origin.lValueLast),
            rising(_origin.rising) {}
    };

    template < typename TText, typename TSpec >
    inline void goNext(Iter< Index<TText, IndexEsa<TSpec> >, VSTree< BottomUp<SuperMaxRepeatsFast> > > &it)
    {
        typedef Index<TText, IndexEsa<TSpec> >        TIndex;
        typename Size<TIndex>::Type                    lcp;

        while (it.lIter != it.lEnd) {
            lcp = *it.lIter;

            if (lcp < it.lValueLast) {
                if (it.rising) {
                    if (it.lValueLast > it.minLength) {
                        _dfsLcp(it) = it.lValueLast;
                        ++_dfsRange(it).i2;
                        ++it.lIter;
                        it.lValueLast = lcp;
                        if (!isPartiallyLeftExtensible(it, it.charSet)) return;
                        continue;
                    }
                    it.rising = false;
                }
            } else
            if (lcp > it.lValueLast) {
                _dfsRange(it).i1 = _dfsRange(it).i2;
                it.rising = true;
            }

            ++_dfsRange(it).i2;
            ++it.lIter;
            it.lValueLast = lcp;
        }
        _dfsRange(it).i2 = 0;
        return;
    }


/*!
 * @class MaxRepeatsIterator Max Repeats Iterator
 * @extends BottomUpIterator
 * @headerfile <seqan/index.h>
 *
 * @brief Iterator to search for all maximal repeats.
 *
 * @signature Iterator<TContainer, MaxRepeats>::Type;
 * @signature template <typename TContainer>
 *            class Iter<TContainer, VSTree< BottomUp<MaxRepeats> > >;
 *
 * @tparam TContainer Type of an index that can be iterated with a bottom-up
 *                    iterator. Types: IndexEsa
 *
 * @note Instead of using the class Iter directly we recommend to use the result of the metafunction
 *       Iterator&lt;TContainer, MaxRepeats&gt;::Type (which is Iter<TContainer, VSTree< BottomUp<MaxRepeats&gt; &gt; &gt;).
 *
 * @see DemoMaximalRepeats
 */

/*!
 * @fn MaxRepeatsIterator::Iter
 * @brief The constructor
 *
 * @signature Iter::Iter(index[, minLength]);
 * @signature Iter::Iter(iterator);
 *
 * @param[in] index     The index to be used for the iteration. Types: @link IndexEsa @endlink
 * @param[in] minLength Minimum length of the supermaximal repeats, default value is 1.
 * @param[in] iterator  Another MaxRepeats iterator. Types: @link MaxRepeatsIterator @endlink
 */

    //////////////////////////////////////////////////////////////////////////////
    // maximal repeats - suffix tree version
    //////////////////////////////////////////////////////////////////////////////

    // contains a list of indices of the same bwt value (fraction)
    template <typename TSize>
    struct FractionHeader_ {
        TSize    begin, end;
        TSize    size;
        FractionHeader_() : begin(0), end(0), size(0) {}
        FractionHeader_(TSize _begin, TSize _end, TSize _size):
            begin(_begin), end(_end), size(_size) {}
    };

    // contains a set of fractions (one for each bwt value)
    // and a fraction for the undefined bwt value (for the virtual character at position -1)
    template <typename TValue, typename TSize>
    struct FractionCompound_ {
        typedef FractionHeader_<TSize>            TFractionHeader;
        typedef Pair<TValue, TFractionHeader>    TFraction;    // TFraction = (c,(begin,end))
        typedef typename Set<TFraction>::Type    TSet;        // c..char, begin/end indices in posList

        TSet            set;
        TFractionHeader    leftmost;

        FractionCompound_():
            leftmost(0,0,0) {}
    };

    // tests existence of left maximal repeats
    // right maximality is implicitly given
    template <typename TValue, typename TSize>
    int _haveMaximalRepeats(
        FractionCompound_<TValue, TSize> const &a,
        FractionCompound_<TValue, TSize> const &b)
    {
        TSize cs = length(a.set), ps = length(b.set);

        if (a.leftmost.size > 0) ++cs;
        if (b.leftmost.size > 0) ++ps;

        if (cs == 0 || ps == 0) return false;
        if (cs  > 1 || ps  > 1) return true;

        if (a.leftmost.size > 0 || b.leftmost.size > 0)
            return true;

        return (keyOf(begin(a.set)) != keyOf(begin(b.set)));
    }


    template <typename TValue, typename TSize>
    int _haveMaximalRepeats(
        FractionCompound_<TValue, TSize> const &a,
        FractionCompound_<TValue, TSize> const &b,
        TValue &equalKey)
    {
        TSize cs = length(a.set), ps = length(b.set);

        if (a.leftmost.size > 0) ++cs;
        if (b.leftmost.size > 0) ++ps;

        if (cs == 0 || ps == 0) return 0;
        if (cs  > 1 || ps  > 1) return 2;    // more than 2

        if (a.leftmost.size > 0 || b.leftmost.size > 0)
            return 2;

        if ((equalKey = keyOf(begin(a.set))) != keyOf(begin(b.set)))
            return 2;
        else
            return 1;
    }


    template < typename TSTree, typename TSpec >
    struct GetVSTreeIteratorTraits< Iter< TSTree, VSTree< BottomUp<MaxRepeats_<TSpec> > > > > {
        typedef PostorderEmptyEdges    Type;
    };

    template < typename TSTree, typename TSpec >
    class Iter< TSTree, VSTree< BottomUp<MaxRepeats_<TSpec> > > >:
        public Iter< TSTree, VSTree< BottomUp<> > >
    {
    public:
        typedef Iter< TSTree, VSTree< BottomUp<> > >    TBase;
        typedef typename Value<TSTree>::Type            TValue;
        typedef typename Size<TSTree>::Type                TSize;

        typedef FractionCompound_<TValue, TSize>        TFractionCompound;
        typedef String<TFractionCompound>               TSetStack;
        typedef String<TSize>                            TPositionList;

        typedef typename TFractionCompound::TSet        TSet;
        typedef typename Iterator<TSet>::Type            TSetIterator;
        typedef typename Iterator<TSet const>::Type        TConstSetIterator;

        typedef typename HistoryStackEntry_<TBase>::Type TStackEntry;

//____________________________________________________________________________

        TSize            minLength;
        TSetStack        setStack;
        TPositionList    posList;    // this list is indexed just as SA is and contains the next entry's index
        bool            canMerge;    // is false, if parent node appears after its first child on stack
//____________________________________________________________________________

        Iter() :
            TBase(),
            minLength(0),
            canMerge(false)
        {}

        Iter(TSTree &_index):
            TBase(_index, MinimalCtor()),
            minLength(1),
            canMerge(true)
        {
            indexRequire(_index, EsaSA());
            indexRequire(_index, EsaLcp());
            indexRequire(_index, EsaBwt());
            resize(posList, length(_index));

            if (!empty(indexSA(_index)))
            {
                TStackEntry e;
                e.range.i1 = 0;
                e.range.i2 = 0;
                _dfsOnPush(*this, e);
                goNext(*this);
            }
        }

        Iter(TSTree &_tree, MinimalCtor):
            TBase(_tree, MinimalCtor()) {}

        Iter(TSTree &_index, TSize _minLength):
            TBase(_index, MinimalCtor()),
            minLength(_minLength),
            canMerge(true)
        {
            indexRequire(_index, EsaSA());
            indexRequire(_index, EsaLcp());
            indexRequire(_index, EsaBwt());
            resize(posList, length(_index));

            if (!empty(indexSA(_index)))
            {
                TStackEntry e;
                e.range.i1 = 0;
                e.range.i2 = 0;
                _dfsOnPush(*this, e);
                goNext(*this);
            }
        }

        Iter(Iter const &_origin):
            TBase((TBase const &)_origin),
            minLength(_origin.minLength),
            setStack(_origin.setStack),
            posList(_origin.posList),
            canMerge(_origin.canMerge) {}

//____________________________________________________________________________

        inline bool hasRepeats()
        {
            if (length(setStack) < 2) return false;
            return _haveMaximalRepeats(back(setStack), setStack[length(setStack) - 2]) > 0;
        }

        inline TSize countRepeats() const
        {
            if (length(setStack) < 2) return 0;

            TFractionCompound const &child  = back(setStack);
            TFractionCompound const &parent = backPrev(setStack);

            TConstSetIterator childFraction    = begin(child.set);
            TConstSetIterator childEnd        = end(child.set);
            TConstSetIterator parentFraction    = begin(parent.set);
            TConstSetIterator parentEnd        = end(parent.set);

            TSize sum = 0;
            for(; childFraction != childEnd; ++childFraction) {
                for(; parentFraction != parentEnd; ++parentFraction) {
                    if (keyOf(childFraction) != keyOf(parentFraction))
                        sum += objectOf(childFraction).size * objectOf(parentFraction).size;

                    sum += child.leftmost.size * objectOf(parentFraction).size;
                }
                sum += objectOf(childFraction).size * parent.leftmost.size;
            }
            sum += child.leftmost.size * parent.leftmost.size;
            return sum;
        }
//____________________________________________________________________________

        inline void _dump() const {
            std::cerr << "SETSTACK of " << representative(*this) << ":" << std::endl;
            typename Iterator<TSetStack const>::Type it = begin(setStack), itEnd = end(setStack);
            while (it != itEnd) {
                TSet const &set = (*it).set;
                typename Iterator<TSet const>::Type sit = begin(set), sitEnd = end(set);

                while (sit != sitEnd) {
                    std::cerr << keyOf(sit) << "::";
                    typename TFractionCompound::TFractionHeader head = objectOf(sit);
                    TSize i = head.begin;
                    while (!_isSizeInval(i)) {
                        std::cerr << saAt(i,container(*this)) << "  ";
                        i = posList[i];
                    }
                    std::cerr << std::endl;
                    ++sit;
                }

                if ((*it).leftmost.size > 0)
                {
                    std::cerr << "\"\"::";
                    TSize i = (*it).leftmost.begin;
                    while (!_isSizeInval(i)) {
                        std::cerr << saAt(i,container(*this)) << "  ";
                        i = posList[i];
                    }
                    std::cerr << std::endl;
                }

                std::cerr << "_________________________" << std::endl;
                ++it;
            }
        }
    };
//____________________________________________________________________________

    // add bwt partitions of child to parent node
    template < typename TSTree, typename TSpec, typename TValue, typename TSize >
    inline void _fractionMerge(
        Iter<TSTree, VSTree< BottomUp<TSpec> > > &it,
        FractionCompound_<TValue, TSize> &parent,
        FractionCompound_<TValue, TSize> &child)
    {
        typedef FractionCompound_<TValue, TSize>    TCompound;
        typedef typename TCompound::TFraction        TFraction;
        typedef typename TCompound::TFractionHeader    TFractionHeader;
        typedef typename TCompound::TSet            TSet;
        typedef typename Iterator<TSet>::Type        TSetIterator;

        TSetIterator _end = end(child.set);
        for(TSetIterator i = begin(child.set); i != _end; ++i) {
            if (in(keyOf(i), parent.set)) {    // append child fraction to parent's fraction
                TFractionHeader &parent_header = objectOf(find(keyOf(i), parent.set));
                TFractionHeader const &child_header = objectOf(i);
                it.posList[parent_header.end] = child_header.begin;
                parent_header.end = child_header.end;
                parent_header.size += child_header.size;
            } else
                insert(TFraction(keyOf(i), objectOf(i)), parent.set);    // insert child fraction in parent's set
        }
        if (parent.leftmost.size > 0) {
            if (child.leftmost.size > 0) {
                it.posList[parent.leftmost.end] = child.leftmost.begin;
                parent.leftmost.end = child.leftmost.end;
                parent.leftmost.size += child.leftmost.size;
            }
        } else
            parent.leftmost = child.leftmost;
    }

    // maximal repeat push/leaf handlers of lcp-dfs-traversal
    template < typename TSTree, typename TElement, typename TSpec >
    inline void _dfsOnPush(Iter<TSTree, VSTree< BottomUp<MaxRepeats_<TSpec> > > > &it, TElement const &e)
    {
        typedef Iter<TSTree, VSTree< BottomUp<> > > TBase;
        _dfsOnPush((TBase&)it, e);

        if (it.canMerge)
            resize(it.setStack, length(it.setStack) + 1);
/*
        std::cerr << "PUSH ";
        _dumpHistoryStack(it);
        it._dump();
*/    }

    template < typename TSTree, typename TSpec >
    inline void _dfsOnLeaf(Iter<TSTree, VSTree< BottomUp<MaxRepeats_<TSpec> > > > &it)
    {
        typedef Iter<TSTree, VSTree< BottomUp<> > > TBase;
        _dfsOnLeaf((TBase&)it);

        typedef typename Value<TSTree>::Type    TValue;
        typedef typename Size<TSTree>::Type        TSize;
        typedef typename SAValue<TSTree>::Type    TSAValue;
        typedef FractionHeader_<TSize>            TFractionHeader;
        typedef Pair<TValue, TFractionHeader>    TFraction;

        resize(it.setStack, length(it.setStack) + 1);

        TSTree &index = container(it);

        TSize       i = _dfsRange(it).i1;
        TSAValue    lPos;
        posLocalize(lPos, getOccurrence(it), stringSetLimits(index));

        if (!posAtFirstLocal(lPos))
            insert(
                TFraction(
                    bwtAt(i, container(it)),
                    TFractionHeader(i, i, 1)),
                back(it.setStack).set);
        else
            back(it.setStack).leftmost = TFractionHeader(i, i, 1);

        _setSizeInval(it.posList[i]);
/*
        std::cerr << "LEAF ";
        _dumpHistoryStack(it);
        it._dump();
*/    }

//____________________________________________________________________________

    template < typename TSTree, typename TSpec >
    inline void goNext(Iter< TSTree, VSTree< BottomUp<MaxRepeats_<TSpec> > > > &it)
    {
        typedef Iter< TSTree, VSTree< BottomUp<MaxRepeats_<TSpec> > > > TIter;
        do {
            if (it.canMerge && length(it.setStack) >= 2)
            {
                typename Size<typename TIter::TSetStack>::Type len = length(it.setStack);
                _fractionMerge(it, it.setStack[len - 2], back(it.setStack));
                resize(it.setStack, len - 1);
            }
            goNext(it, PostorderEmptyEdges());
            if (empty(it.history))
                it.canMerge = false;
            else
                it.canMerge = !_dfsReversedOrder(it);
        } while (!eof(it) && !(it.canMerge && (repLength(it) >= it.minLength) && it.hasRepeats()));
    }
//____________________________________________________________________________

    template < typename TSTree, typename TSpec >
    inline typename VertexDescriptor<TSTree>::Type
    value(Iter< TSTree, VSTree< BottomUp<MaxRepeats_<TSpec> > > > const &it)
    {
        if (empty(it.history))
            return it.vDesc;
        typedef typename VertexDescriptor<TSTree>::Type TDesc;
        return TDesc(back(it.history).range.i1, it.vDesc.range.i2, 0);
    }
//____________________________________________________________________________

    template < typename TSTree, typename TSpec >
    inline typename Size<TSTree>::Type
    repLength(Iter< TSTree, VSTree< BottomUp<MaxRepeats_<TSpec> > > > const &it)
    {
        return back(it.history).range.i2;
    }
//____________________________________________________________________________

/*!
 * @fn MaxRepeatsIterator#length
 * @headerfile <seqan/index.h>
 * @brief Return the number of repeats.
 *
 * @signature TSize length(it);
 *
 * @param[in] it The MaxRepeatsIterator to query.
 *
 * @return TSize The number of found repeats.
 */


    template < typename TSTree, typename TSpec >
    inline typename Size<TSTree>::Type
    length(Iter< TSTree, VSTree< BottomUp<MaxRepeats_<TSpec> > > > const &it) {
        return it.countRepeats();
    }
//____________________________________________________________________________

    template < typename TSTree, class TSpec >
    inline typename Iterator< Iter< TSTree, VSTree< BottomUp<MaxRepeats_<TSpec> > > > >::Type
    begin(Iter< TSTree, VSTree< BottomUp<MaxRepeats_<TSpec> > > > &it)
    {
        typedef Iter< TSTree, VSTree< BottomUp<MaxRepeats_<TSpec> > > >    TIter;
        return typename Iterator<TIter>::Type(it);
    }
//____________________________________________________________________________

    template < typename TSTree, class TSpec >
    inline typename Iterator< Iter< TSTree, VSTree< BottomUp<MaxRepeats_<TSpec> > > > >::Type
    end(Iter< TSTree, VSTree< BottomUp<MaxRepeats_<TSpec> > > > &it)
    {
        typedef Iter< TSTree, VSTree< BottomUp<MaxRepeats_<TSpec> > > >    TIter;
        return typename Iterator<TIter>::Type(it, MinimalCtor());
    }




    //////////////////////////////////////////////////////////////////////////////
    // maximal repeat representation
    //////////////////////////////////////////////////////////////////////////////

    template <typename TSTree>
    struct MaxRepeat {
//        Iter< TSTree, VSTree<BottomUp<MaxRepeats> > > &it;
    };

    template <typename TSTree>
    struct Value< MaxRepeat<TSTree> > {
        typedef Pair< typename SAValue<TSTree>::Type > Type;
    };

    template <typename TSTree>
    struct Size< MaxRepeat<TSTree> > {
        typedef typename Size<TSTree>::Type Type;
    };


    template <typename TSTree>
    inline typename Size< MaxRepeat<TSTree> >::Type
    length(MaxRepeat<TSTree> const &repeat) {
        return repeat.it.countRepeats();
    }

/*
    template <typename TSTree>
    inline typename Iterator< MaxRepeat<TSTree> >::Type
    begin(MaxRepeat<TSTree> &repeat) {
        return Iterator< MaxRepeat<TSTree> >::Type(repeat.it);
    }

    template <typename TSTree>
    inline typename Iterator< MaxRepeat<TSTree> const >::Type
    begin(MaxRepeat<TSTree> const &repeat) {
        return Iterator< MaxRepeat<TSTree> >::Type(repeat.it);
    }
*/


    template <typename TSTree>
    class Iter< MaxRepeat<TSTree>, MaxRepeatOccurrences > {
    public:

        typedef typename Value<TSTree>::Type    TValue;
        typedef typename Size<TSTree>::Type        TSize;
        typedef    Pair<TSize>                        TPair;

        typedef FractionCompound_<TValue, TSize>    TFractionCompound;
        typedef typename TFractionCompound::TSet    TSet;
        typedef typename Iterator<TSet const>::Type    TSetIterator;

        TSize            childPtr, parentPtr;
        TSetIterator    childFraction, childEnd;
        TSetIterator    parentFraction, parentEnd;
        bool            _atEnd;
        TPair            tmp;
        bool            leftmostChild, leftmostParent;

        Iter<TSTree, VSTree<BottomUp<MaxRepeats> > > const *maxIt;

        inline Iter():
            childPtr(0), parentPtr(0),
            _atEnd(true),
            tmp(TPair(0, 0)),
            leftmostChild(false), leftmostParent(false) {}

        inline Iter(Iter<TSTree, VSTree<BottomUp<MaxRepeats> > > const &_maxIt):
            childPtr(0), parentPtr(0),
            _atEnd(true),
            tmp(TPair(0, 0)),
            leftmostChild(false), leftmostParent(false),
            maxIt(&_maxIt)
        {
            _init();
        }

        inline Iter(Iter<TSTree, VSTree<BottomUp<MaxRepeats> > > const &_maxIt, MinimalCtor):
            maxIt(&_maxIt),
            _atEnd(true) {}

        inline void _reinitParentFraction()
        {
            if (leftmostParent)
            {
                TFractionCompound const &parent = maxIt->setStack[length(maxIt->setStack) - 2];
                parentPtr = parent.leftmost.begin;
            } else
                parentPtr = objectOf(parentFraction).begin;
        }

        inline bool _innerStep() {
            if (_isSizeInval(childPtr = maxIt->posList[childPtr])) {
                if (_isSizeInval(parentPtr = maxIt->posList[parentPtr]))
                {
                    _reinitParentFraction();
                    return false;
                }
                childPtr = objectOf(childFraction).begin;
            }
            return true;
        }

        inline void _firstParentFraction() {
            TFractionCompound const &parent = maxIt->setStack[length(maxIt->setStack) - 2];

            parentFraction    = begin(parent.set);
            parentEnd        = end(parent.set);

            if (parentFraction != parentEnd) {
                leftmostParent = false;
                parentPtr = objectOf(parentFraction).begin;
            } else {
                leftmostParent = true;
                parentPtr = parent.leftmost.begin;
            }
        }

        inline void _firstChildFraction() {
            TFractionCompound const &child = back(maxIt->setStack);

            childFraction    = begin(child.set);
            childEnd        = end(child.set);

            if (childFraction != childEnd) {
                leftmostChild = false;
                childPtr = objectOf(childFraction).begin;
            } else {
                leftmostChild = true;
                childPtr = child.leftmost.begin;
            }
        }

        inline bool _nextParentFraction()
        {
            if (leftmostParent)
                return false;

            if (++parentFraction == parentEnd) {
                if (maxIt->setStack[length(maxIt->setStack) - 2].leftmost.size > 0) {
                    leftmostParent = true;
                    parentPtr = maxIt->setStack[length(maxIt->setStack) - 2].leftmost.begin;
                } else
                    return false;
            } else
                parentPtr = objectOf(parentFraction).begin;

            return true;
        }

        inline bool _nextChildFraction() {
            if (leftmostChild)
                return false;

            if (++childFraction == childEnd) {
                if (back(maxIt->setStack).leftmost.size > 0) {
                    leftmostChild = true;
                    childPtr = back(maxIt->setStack).leftmost.begin;
                } else
                    return false;
            } else
                childPtr = objectOf(childFraction).begin;

            return true;
        }

        inline bool _outerStep() {
            do {
                if (!_nextChildFraction()) {
                    _firstChildFraction();
                    if (!_nextParentFraction()) {
                        _atEnd = true;
                        return false;
                    }
                }
                if (leftmostChild || leftmostParent) break;
            } while (keyOf(childFraction) == keyOf(parentFraction));        // ignore occurrences with equal bwt entries
            return true;
        }

        inline void _init()
        {
            if (length(maxIt->setStack) < 2) {
                _atEnd = true;
                return;
            }

            _firstChildFraction();
            _firstParentFraction();

            if (!leftmostChild && !leftmostParent &&
                (keyOf(childFraction) == keyOf(parentFraction)))
                _atEnd = !_outerStep();
            else
                _atEnd = false;

            if (!_atEnd) {
                tmp.i1 = saAt(parentPtr, container(*maxIt));
                tmp.i2 = saAt(childPtr, container(*maxIt));
            }
        }
    };

//____________________________________________________________________________

    template < typename TRepeat >
    inline typename Value< Iter<TRepeat, MaxRepeatOccurrences> >::Type &
    value(Iter<TRepeat, MaxRepeatOccurrences> const &it)  {
        return it.tmp;
    }

    template < typename TRepeat >
    inline typename Value< Iter<TRepeat, MaxRepeatOccurrences> >::Type &
    value(Iter<TRepeat, MaxRepeatOccurrences> &it)  {
        return it.tmp;
    }
//____________________________________________________________________________

    template < typename TRepeat >
    inline Iter<TRepeat, MaxRepeatOccurrences> &
    goNext(Iter<TRepeat, MaxRepeatOccurrences> &it)  {
        if (it._innerStep()) {
            it.tmp.i1 = saAt(it.parentPtr, container(*it.maxIt));
            it.tmp.i2 = saAt(it.childPtr, container(*it.maxIt));
            return it;
        }
        if (it._outerStep()) {
            it.tmp.i1 = saAt(it.parentPtr, container(*it.maxIt));
            it.tmp.i2 = saAt(it.childPtr, container(*it.maxIt));
        }
        return it;
    }
//____________________________________________________________________________

    template < typename TRepeat >
    inline void
    goBegin(Iter<TRepeat, MaxRepeatOccurrences> &it)
    {
        it._init;
    }

//____________________________________________________________________________

    template < typename TRepeat >
    inline void
    goEnd(Iter<TRepeat, MaxRepeatOccurrences> &it)
    {
        it._atEnd = true;
    }

    template < typename TRepeat >
    inline bool atEnd(Iter<TRepeat, MaxRepeatOccurrences> const &it) {
        return it._atEnd;
    }

    template < typename TRepeat >
    inline bool atEnd(Iter<TRepeat, MaxRepeatOccurrences> &it) {
        return it._atEnd;
    }

//____________________________________________________________________________

    template < typename TRepeat >
    inline bool
    operator == (
        Iter<TRepeat, MaxRepeatOccurrences> const &itA,
        Iter<TRepeat, MaxRepeatOccurrences> const &itB)
    {
        if (itA._atEnd && itB._atEnd) return true;
        if (itA._atEnd || itB._atEnd) return false;
        return (itA.childPtr == itB.childPtr) && (itA.parentPtr == itB.parentPtr);
    }

    template < typename TRepeat >
    inline bool
    operator != (
        Iter<TRepeat, MaxRepeatOccurrences> const &itA,
        Iter<TRepeat, MaxRepeatOccurrences> const &itB)
    {
        if (itA._atEnd && itB._atEnd) return false;
        if (itA._atEnd || itB._atEnd) return true;
        return (itA.childPtr != itB.childPtr) || (itA.parentPtr != itB.parentPtr);
    }
//____________________________________________________________________________


    template <typename TSTree, typename TSpec>
    struct Size< Iter<TSTree, VSTree< BottomUp<MaxRepeats_<TSpec> > > > > {
        typedef typename Size<TSTree>::Type Type;
    };

    template <typename TSTree>
    struct Size< Iter<MaxRepeat<TSTree>, MaxRepeatOccurrences> > {
        typedef typename Size<TSTree>::Type Type;
    };


    //////////////////////////////////////////////////////////////////////////////
    // Iterator wrappers
    //////////////////////////////////////////////////////////////////////////////

    // iterates over all supermaximal repeats
    template <typename TSTree>
    struct Iterator< TSTree, SuperMaxRepeats > {
        typedef Iter< TSTree, VSTree< BottomUp<SuperMaxRepeats> > > Type;
    };
//____________________________________________________________________________

    // iterates over all maximal unique matches
    template <typename TSTree>
    struct Iterator< TSTree, SuperMaxRepeatsFast > {
        typedef Iter< TSTree, VSTree< BottomUp<SuperMaxRepeatsFast> > > Type;
    };
//____________________________________________________________________________

    // iterates over all maximal repeat structures
    template <typename TSTree, typename TSpec>
    struct Iterator< TSTree, MaxRepeats_<TSpec> > {
        typedef Iter< TSTree, VSTree< BottomUp<MaxRepeats_<TSpec> > > > Type;
    };
//____________________________________________________________________________

    // iterates over all maximal repeat pairs of a repeat structure

    // Iterator of Iterator<TIndex, MaxRepeats>
    template <typename TSTree, typename TSpec>
    struct Iterator< Iter< TSTree, TSpec >, MaxRepeatOccurrences > {
        typedef Iter< MaxRepeat<TSTree>, MaxRepeatOccurrences > Type;
    };

    template <typename TSTree, typename TSpec>
    struct DefaultIteratorSpec< Iter< TSTree, VSTree< BottomUp<MaxRepeats_<TSpec> > > > > {
        typedef MaxRepeatOccurrences Type;
    };

    // alternative (use MaxRepeat<TIndex> as a wrapper in between)
    template <typename TSTree>
    struct Iterator< MaxRepeat<TSTree>, MaxRepeatOccurrences > {
        typedef Iter <MaxRepeat<TSTree>, MaxRepeatOccurrences > Type;
    };

    template <typename TSTree>
    struct DefaultIteratorSpec< MaxRepeat<TSTree> > {
        typedef MaxRepeatOccurrences Type;
    };

//}

}

#endif
