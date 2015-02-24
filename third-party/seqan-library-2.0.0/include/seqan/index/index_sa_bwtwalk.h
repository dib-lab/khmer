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
// Authors: Tobias Marschall <tobias.marschall@cs.uni-dortmund.de>
//          Marcel Martin <marcel.martin@tu-dortmund.de>
//          David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_INDEX_SA_BWTWALK_H
#define SEQAN_HEADER_INDEX_SA_BWTWALK_H

// TODO export lcp table by defining:
//    template < typename TSA, typename TLCP, typename TText, typename TAlgSpec >
//    inline void createSAAndLCP(TSA &sa, TLCP &lcp, TText const &s, TAlgSpec const &alg)
//
// TODO: The in-place permutation that uses an extra bit vector is not tested automatically,
//       since it is only used for strings longer than 2**31 or 2**63 characters.

// definitions:
// lexnextpos[p] = pos[rank[p]+1]
// lexprevpos[p] = pos[rank[p]-1]
// lexxorpos = lexprevpos ^ lexnextpos

namespace SEQAN_NAMESPACE_MAIN
{
    // Public tags
    struct BwtWalkFast_ {};
    typedef Tag<BwtWalkFast_> const BwtWalkFast;

    struct BwtWalkInPlace_ {};
    typedef Tag<BwtWalkInPlace_> const BwtWalkInPlace;

    template <typename TSpec = BwtWalkFast >
    struct BwtWalk {};

    // Public function
    // NOTE: This function uses the supremum value of Value<TSA>::Type as NIL symbol
    //       for temporary calculations. Therefore, the caller must ensure that
    //       length(s) < MaxValue<TValue>::VALUE. Otherwise, behaviour is undefined!
    template < typename TSA, typename TText, typename TSpec >
    inline void createSuffixArray(
        TSA &SA,
        TText const &s,
        BwtWalk<TSpec> const &)
    {
        typedef typename AllowsFastRandomAccess<TSA>::Type TAllowsFastRandomAccess;
        _createSuffixArray(SA, s, BwtWalk<TSpec>(), TAllowsFastRandomAccess());
    }

    // List functions, for internal use (that's why there's a '_' in front)

    // Inserts p into the suffix list between the given predecessor and successor
    // NOTE: Intentionally not using references for predecessor and successor!
    template < typename TArray, typename TValue >
    inline void _insertBetween(
        TArray &lexprevpos,
        TArray &lexnextpos,
        TValue p,
        TValue predecessor,
        TValue successor)
    {
        value(lexprevpos, p) = predecessor;
        value(lexnextpos, p) = successor;
        const TValue NIL = MaxValue<TValue>::VALUE;
        if (predecessor != NIL) value(lexnextpos, predecessor) = p;
        if (successor != NIL) value(lexprevpos, successor) = p;
    }

    // Inserts p into the suffix list before the given successor
    template < typename TArray, typename TValue >
    inline void _insertBefore(
        TArray &lexprevpos,
        TArray &lexnextpos,
        TValue p,
        TValue successor)
    {
        _insertBetween(lexprevpos, lexnextpos, p, getValue(lexprevpos, successor), successor);
    }

    // Inserts p into the suffix list after the given predecessor
    template < typename TArray, typename TValue >
    inline void _insertAfter(
        TArray &lexprevpos,
        TArray &lexnextpos,
        TValue p,
        TValue predecessor)
    {
        _insertBetween(lexprevpos, lexnextpos, p, predecessor, getValue(lexnextpos, predecessor));
    }

    // Inserts p into the suffix list between the given predecessor and successor
    template < typename TArray, typename TValue >
    inline void _insertXorBetween(
        TArray &lexxorpos,
        TValue p,
        TValue predecessor,
        TValue successor)
    {
        const TValue NIL = MaxValue<TValue>::VALUE;
        value(lexxorpos, p) = predecessor ^ successor;
        if (predecessor != NIL) value(lexxorpos, predecessor) = getValue(lexxorpos, predecessor) ^ successor ^ p;
        if (successor != NIL) value(lexxorpos, successor) = getValue(lexxorpos, successor) ^ predecessor ^ p;
    }

    // Creates a suffix array. This is the version for fast random access.
    template < typename TSA, typename TText >
    inline void _createSuffixArray(
        TSA &SA,
        TText const &s,
        BwtWalk<BwtWalkFast> const &,
        True const &)
    {
        typedef typename Value<TSA>::Type TValue;
        typedef typename Iterator<TSA, Standard>::Type TSaIter;
        const TValue NIL = MaxValue<TValue>::VALUE;

        if (empty(s)) return;

        TSA lexnextpos;
        // Since SA is capable of fast random access, we can temporarily use it as lexprevpos array.
        TSA &lexprevpos = SA;

        resize(lexnextpos, length(s), NIL, Exact());
        std::fill(begin(lexprevpos, Standard()), end(lexprevpos, Standard()), NIL);

        // finally create suffix array
        //   1) first suffix is returned...
        TValue p = _createSuffixList(lexprevpos, lexnextpos, s);
        //   2) sequentially write out suffix array
        TSaIter saIt = begin(SA, Standard());
        while (p != NIL) {
            *saIt = p;
            p = getValue(lexnextpos, p);
            ++saIt;
        }
    }

    // Creates a suffix array. This is the version for external memory.
    // (SA is written sequentially from start to finish.)
    template < typename TSA, typename TText >
    inline void _createSuffixArray(
        TSA &SA,
        TText const& s,
        BwtWalk<BwtWalkFast> const &,
        False const &)
    {
        typedef typename Value<TSA>::Type TValue;
        typedef typename Iterator<TSA, Standard>::Type TSaIter;
        const TValue NIL = MaxValue<TValue>::VALUE;
        typedef String<TValue> TArray;

        if (empty(s)) return;

        TArray lexnextpos;
        TArray lexprevpos;

        resize(lexnextpos, length(s), NIL, Exact());
        resize(lexprevpos, length(s), NIL, Exact());

        // finally create suffix array
        //   1) first suffix is returned...
        TValue p = _createSuffixList(lexprevpos, lexnextpos, s);
        //   2) sequentially write out suffix array
        TSaIter saIt = begin(SA, Standard());
        while (p != NIL) {
            *saIt = p;
            p = getValue(lexnextpos, p);
            ++saIt;
        }
    }

    // Note: Writing to SA is done sequential only, i.e. well-suited to be directly
    //       written to external storage.
    template < typename TSA, typename TText >
    inline void _createSuffixArray(
        TSA &SA,
        TText const &s,
        BwtWalk<BwtWalkInPlace> const &,
        False const &)
    {
        typedef typename Value<TSA>::Type TValue;
        typedef typename Iterator<TSA, Standard>::Type TSaIter;
        const TValue NIL = MaxValue<TValue>::VALUE;
        typedef String<TValue> TArray;

        if (empty(s)) return;

        TArray lexxorpos;
        // invalid ^ invalid == 0
        resize(lexxorpos, length(s), 0, Exact());

        // finally create suffix array
        //   1) first suffix is returned...
        TValue p = _createXoredSuffixList(lexxorpos, s);
        //   2) sequentially write out suffix array
        TValue previous = NIL;
        TSaIter saIt = begin(SA, Standard());
        while (p != NIL) {
            *saIt = p;
            TValue tmp = p;
            p = getValue(lexxorpos, p) ^ previous;
            previous = tmp;
            ++saIt;
        }
    }

    // Variant for fast random access SA
    // This version computes the suffix array in place, but only if the text s is
    // short enough such that the most significant bit (MSB) of the entries in the
    // SA can be used as a flag bit, which is needed to invert a permutation in place.
    // Otherwise, an additional length(s) bits of memory are used.

    // This is a helper function for the case "length(s) >= (NIL && NOTFLAGBIT)".
    template < typename TSA, typename TText, typename TValue, typename TArray >
    inline void _createSuffixArray(
        TSA &SA,
        TArray & lexxorpos,
        TValue & p,
        TText const &s,
        TValue const & NIL,
        TValue const &,
        TValue const &,
        BwtWalk<BwtWalkInPlace> const &,
        True const &,
        True const &)
    {
        // If the string is too long, we cannot use the MSB of the SA entries as a flag.
        // Instead, we allocate a separate bitvector
        String<bool, Packed<> > visited;
        resize(visited, length(s), false);

        // 2) convert into inverse suffix array
        TValue n = 0;
        TValue previous = NIL;
        while (p != NIL) {
            TValue tmp = p;
            p = getValue(lexxorpos, p) ^ previous;
            previous = tmp;
            value(lexxorpos, previous) = n++;
        }

        // 3) invert array in place to obtain suffix array
        for (n = 0; n < length(s); ++n) {
            TValue i = getValue(SA, n);
            // check if current position has already been processed
            if (!getValue(visited, i)) {
                // if not, invert current cycle
                TValue previous = n;
                while (i != n) {
                    TValue next = getValue(SA, i);
                    value(SA, i) = previous;
                    value(visited, i) = true;
                    previous = i;
                    i = next;
                }
                value(visited, n) = true;
                value(SA, i) = previous;
            }
        }
    }

    // This is a helper function for the case "length(s) < (NIL && NOTFLAGBIT)".
    template < typename TSA, typename TText, typename TValue, typename TArray >
    inline void _createSuffixArray(
        TSA &SA,
        TArray & lexxorpos,
        TValue & p,
        TText const &s,
        TValue const & NIL,
        TValue const & FLAGBIT,
        TValue const & NOTFLAGBIT,
        BwtWalk<BwtWalkInPlace> const &,
        True const &,
        False const &)
    {
        // The string is short enough, we can use the MSB of the SA entries as a flag.

        // 2) convert into inverse suffix array and set every flag
        TValue n = 0;
        TValue previous = NIL;
        while (p != NIL) {
            TValue tmp = p;
            p = getValue(lexxorpos, p) ^ previous;
            previous = tmp;
            value(lexxorpos, previous) = FLAGBIT | n++;
        }

        // 3) invert array in place to obtain suffix array
        for (n = 0; n < length(s); ++n) {
            TValue i = getValue(SA, n);
            // check if current position has already been processed
            if ((i & FLAGBIT) != 0) {
                // if not, invert current cycle
                i &= NOTFLAGBIT;
                TValue previous = n;
                while (i != n) {
                    TValue next = getValue(SA, i) & NOTFLAGBIT;
                    value(SA, i) = previous;
                    previous = i;
                    i = next;
                }
                value(SA, i) = previous;
            }
        }
    }

    // This takes care of the case where the suffix array's entry type is
    // longer than the string length type to circumvent some warnings.
    template < typename TSA, typename TText, typename TValue, typename TArray >
    inline void _createSuffixArrayHandleTooLongSuffixArrayType(
        TSA &SA,
        TArray & lexxorpos,
        TValue & p,
        TText const &s,
        TValue const & NIL,
        TValue const & FLAGBIT,
        TValue const & NOTFLAGBIT,
        BwtWalk<BwtWalkInPlace> const &,
        True const &,
        True const &)
    {
        _createSuffixArray(SA, lexxorpos, p, s, NIL, FLAGBIT, NOTFLAGBIT, BwtWalk<BwtWalkInPlace>(), True(), False());
    }

    // This takes care of the case where the suffix array's entry type is NOT
    // longer than the string length type to circumvent some warnings.
    template < typename TSA, typename TText, typename TValue, typename TArray >
    inline void _createSuffixArrayHandleTooLongSuffixArrayType(
        TSA &SA,
        TArray & lexxorpos,
        TValue & p,
        TText const &s,
        TValue const & NIL,
        TValue const & FLAGBIT,
        TValue const & NOTFLAGBIT,
        BwtWalk<BwtWalkInPlace> const &,
        True const &,
        False const &)
    {
        if (length(s) >= (NIL & NOTFLAGBIT))
            _createSuffixArray(SA, lexxorpos, p, s, NIL, FLAGBIT, NOTFLAGBIT, BwtWalk<BwtWalkInPlace>(), True(), True());
        else
            _createSuffixArray(SA, lexxorpos, p, s, NIL, FLAGBIT, NOTFLAGBIT, BwtWalk<BwtWalkInPlace>(), True(), False());
    }

    // This is the main entry point.
    template < typename TSA, typename TText >
    inline void _createSuffixArray(
        TSA &SA,
        TText const &s,
        BwtWalk<BwtWalkInPlace> const &,
        True const &)
    {
        typedef typename Value<TSA>::Type TValue;
        const TValue NIL = MaxValue<TValue>::VALUE;

        const TValue FLAGBIT = (TValue)1 << (BitsPerValue<TValue>::VALUE-1);
        const TValue NOTFLAGBIT = ~FLAGBIT;

        if (empty(s)) return;

        // since SA is capable of fast random access, we can temporarily use it as the lexxorpos array
        TSA & lexxorpos = SA;
        std::fill(begin(lexxorpos, Standard()), end(lexxorpos, Standard()), 0);

        // finally create suffix array
        //   1) first suffix is returned ...
        TValue p = _createXoredSuffixList(lexxorpos, s);

        typedef typename Size<TText>::Type TSize;
        typedef typename Eval<!(sizeof(TSize) <= sizeof(TValue))>::Type TBool;
        _createSuffixArrayHandleTooLongSuffixArrayType(SA, lexxorpos, p, s, NIL, FLAGBIT, NOTFLAGBIT, BwtWalk<BwtWalkInPlace>(), True(), TBool());
    }

    template < typename TArray, typename TValue >
    inline void _findNewCharInsertionPosition(
        TArray &lexfirstpos,
        TArray &lexlastpos,
        TValue &pMinus,
        TValue &pPlus,
        unsigned int cOrd)
    {
        const TValue NIL = MaxValue<TValue>::VALUE;

        // current char does not exist
        // first find pMinus ...
        pMinus = NIL;
        unsigned int cMinus = cOrd;
        while (cMinus > 0) {
            --cMinus;
            if (getValue(lexlastpos, cMinus) != NIL) {
                pMinus = getValue(lexlastpos, cMinus);
                break;
            }
        }

        // ... and pPlus
        pPlus = NIL;
        for (unsigned int cPlus = cOrd+1; cPlus < length(lexfirstpos); ++cPlus) {
            if (getValue(lexfirstpos, cPlus) != NIL) {
                pPlus = getValue(lexfirstpos, cPlus);
                break;
            }
        }
    }

    // returns begin of suffix list
    template <typename TArray, typename TText>
    inline typename Value<TArray>::Type _createSuffixList(
        TArray &lexprevpos,
        TArray &lexnextpos,
        TText const &s)
    {
        typedef typename Iterator<const TText, Standard>::Type TTextIter;
        typedef typename Value<TArray>::Type TValue;

        typedef typename Value<const TText>::Type TChar;
        const unsigned int ALPHABETSIZE = ValueSize<TChar>::VALUE;
        const TValue NIL = MaxValue<TValue>::VALUE;
        typedef String<TValue, Array<ALPHABETSIZE> > TAlphabetArray;

        TAlphabetArray lexfirstpos, lexlastpos;

        resize(lexfirstpos, ALPHABETSIZE, NIL, Exact());
        resize(lexlastpos, ALPHABETSIZE, NIL, Exact());

        // main loop
        TTextIter it = end(s, Standard());
        TValue p = length(s);
        while (it != begin(s, Standard())) {
            --it;
            --p;
            TChar c = *it;
            // Instead of using c as an index, we use its ordValue.
            unsigned int cOrd = ordValue(c);
            if (getValue(lexfirstpos, cOrd) == NIL) {
                TValue pMinus, pPlus;
                _findNewCharInsertionPosition(lexfirstpos, lexlastpos, pMinus, pPlus, cOrd);
                _insertBetween(lexprevpos, lexnextpos, p, pMinus, pPlus);
                // this is our first suffix starting with c, remember in lexfirstpos/lexlastpos
                value(lexfirstpos, cOrd) = p;
                value(lexlastpos, cOrd) = p;
            } else {
                // current char already exists (case 2. in the paper)
                TValue pLeft = getValue(lexprevpos, p+1); // walking left
                TValue pRight = getValue(lexnextpos, p+1); // walking right
                // walk along Bwt in both directions
                while (true) {
                    // check if we have found the insertion position while walking to the left side
                    if (pLeft == NIL) {
                        if (back(s) == c) {
                            // this is the lexicographically second suffix starting with c, the first
                            // suffix starting with c consists only of one character (the last char in s)
                            // this case could be avoided by adding a $ symbol to the end of s.
                            _insertAfter(lexprevpos, lexnextpos, p, (TValue)(length(s)-1));
                            if (getValue(lexlastpos, cOrd) == length(s)-1) value(lexlastpos, cOrd) = p;
                        } else {
                            // this is the lexicographically first suffix starting with c,
                            // insert before lexfirstpos[c]
                            _insertBefore(lexprevpos, lexnextpos, p, getValue(lexfirstpos, cOrd));
                            value(lexfirstpos, cOrd) = p;
                        }
                        break;
                    } else {
                        // check if character has been found walking to the left
                        if (getValue(s, pLeft-1) == c) {
                            // insert in the middle of the list
                            _insertAfter(lexprevpos, lexnextpos, p, pLeft-1);
                            // check if just inserted suffix is now the last suffix starting with c
                            if (getValue(lexlastpos, cOrd) == pLeft-1) value(lexlastpos, cOrd) = p;
                            break;
                        }
                    }
                    pLeft = getValue(lexprevpos, pLeft);

                    // check if we have found the insertion position while walking to the right side
                    if (pRight == NIL) {
                        // this is the lexicographically last suffix starting with c,
                        // insert before lexfirstpos[c]
                        _insertAfter(lexprevpos, lexnextpos, p, getValue(lexlastpos, cOrd));
                        value(lexlastpos, cOrd) = p;
                        break;
                    } else {
                        // check if character has been found walking to the right
                        if (getValue(s, pRight-1) == c) {
                            // insert in the middle of the list
                            _insertBefore(lexprevpos, lexnextpos, p, pRight-1);
                            // check if just inserted suffix is now the first suffix starting with c
                            if (getValue(lexfirstpos, cOrd) == pRight-1) value(lexfirstpos, cOrd) = p;
                            break;
                        }
                    }
                    pRight = getValue(lexnextpos, pRight);
                }
            }
        }

        typename Iterator<TAlphabetArray, Standard>::Type fit = begin(lexfirstpos, Standard());
        while (*fit == NIL) ++fit;
        return *fit;
    }

    // returns beginning of suffix list
    template <typename TArray, typename TText>
    inline typename Value<TArray>::Type _createXoredSuffixList(
        TArray &lexxorpos,
        TText const &s)
    {
        typedef typename Iterator<const TText, Standard>::Type TTextIter;
        typedef typename Value<TArray>::Type TValue;

        typedef typename Value<const TText>::Type TChar;
        const unsigned int ALPHABETSIZE = ValueSize<TChar>::VALUE;
        const TValue NIL = MaxValue<TValue>::VALUE;
        typedef String<TValue, Array<ALPHABETSIZE> > TAlphabetArray;

        TAlphabetArray lexfirstpos;
        TAlphabetArray lexlastpos;
        resize(lexfirstpos, ALPHABETSIZE, NIL, Exact());
        resize(lexlastpos, ALPHABETSIZE, NIL, Exact());

        // main loop
        TTextIter it = end(s, Standard());
        TValue p = length(s);
        // store predecessor and successor from last round (p+1)
        TValue lastPredecessor = NIL;
        TValue lastSuccessor = NIL;
        while (it != begin(s, Standard())) {
            --it;
            --p;
            TChar c = *it;
            // Instead of using c as an index, we use its ordValue.
            unsigned int cOrd = ordValue(c);
            // to insert into xored list, we need to find predecessor _and_ successor
            TValue predecessor = NIL;
            TValue successor = NIL;
            if (getValue(lexfirstpos, cOrd) == NIL) {
                _findNewCharInsertionPosition(lexfirstpos, lexlastpos, predecessor, successor, cOrd);
                // this is our first suffix starting with c, remember in lexfirstpos/lexlastpos
                value(lexfirstpos, cOrd) = p;
                value(lexlastpos, cOrd) = p;
            } else {
                // current char already exists
                bool skipWalkRight = false;
                // ----- walk left ------------------ (search for predecessor) ------
                TValue pLast = p+1;
                TValue pLeft = getValue(lexxorpos, p+1) ^ lastSuccessor; // walking left
                while (true) {
                    // check if we have found the insertion position while walking to the left side
                    if (pLeft == NIL) {
                        if (back(s) == c) {
                            // this is the lexicographically second suffix starting with c, the first
                            // suffix starting with c consists only of one character (the last char in s)
                            // this case could be avoided by adding a $ symbol to the end of s.
                            predecessor = (TValue)(length(s)-1);
                            if (getValue(lexlastpos, cOrd) == length(s)-1) value(lexlastpos, cOrd) = p;
                        } else {
                            // this is the lexicographically first suffix starting with c,
                            // insert before lexfirstpos[c]
                            successor = getValue(lexfirstpos, cOrd);
                            // first find pMinus ...
                            unsigned int cMinus = cOrd;
                            while (cMinus > 0) {
                                --cMinus;
                                if (getValue(lexlastpos, cMinus) != NIL) {
                                    predecessor = getValue(lexlastpos, cMinus);
                                    break;
                                }
                            }
                            skipWalkRight = true;
                            value(lexfirstpos, cOrd)=p;
                        }
                        break;
                    } else {
                        // check if character has been found walking to the left
                        if (getValue(s, pLeft-1) == c) {
                            // insert in the middle of the list
                            predecessor = pLeft-1;
                            // check if just inserted suffix is now the last suffix starting with c
                            if (getValue(lexlastpos, cOrd) == pLeft-1) value(lexlastpos, cOrd) = p;
                            break;
                        }
                    }
                    TValue tmp = pLeft;
                    pLeft = getValue(lexxorpos, pLeft) ^ pLast;
                    pLast = tmp;
                }

                // ----- walk right ------------------ (search for successor) ------
                pLast = p+1;
                TValue pRight = getValue(lexxorpos, p+1) ^ lastPredecessor; // walking right
                while (!skipWalkRight) {
                    // check if we have found the insertion position while walking to the right side
                    if (pRight == NIL) {
                        // this is the lexicographically last suffix starting with c,
                        for (unsigned int cPlus = cOrd+1; cPlus < length(lexfirstpos); ++cPlus) {
                            if (getValue(lexfirstpos, cPlus) != NIL) {
                                successor = getValue(lexfirstpos, cPlus);
                                break;
                            }
                        }
                        value(lexlastpos, cOrd) = p;
                        break;
                    } else {
                        // check if character has been found walking to the right
                        if (getValue(s, pRight-1) == c) {
                            successor = pRight-1;
                            // check if just inserted suffix is now the first suffix starting with c
                            if (getValue(lexfirstpos, cOrd) == pRight-1) value(lexfirstpos, cOrd) = p;
                            break;
                        }
                    }
                    TValue tmp = pRight;
                    pRight = getValue(lexxorpos, pRight) ^ pLast;
                    pLast = tmp;
                }
            }
            // finally insert between predecessor and successor
            _insertXorBetween(lexxorpos, p, predecessor, successor);
            lastPredecessor = predecessor;
            lastSuccessor = successor;
        }

        typename Iterator<TArray, Standard>::Type fit = begin(lexfirstpos, Standard());
        while (*fit == NIL) ++fit;
        return *fit;
    }
}
#endif
