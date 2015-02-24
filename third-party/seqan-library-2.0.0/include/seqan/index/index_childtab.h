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

#ifndef SEQAN_HEADER_INDEX_CHILDTAB_H
#define SEQAN_HEADER_INDEX_CHILDTAB_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

    struct Childtab {};

    //////////////////////////////////////////////////////////////////////////////
    // external childtab algorithm
    //////////////////////////////////////////////////////////////////////////////

    template < typename TLCPInput >
    struct Value< Pipe< TLCPInput, Childtab > > {
        typedef typename Size<TLCPInput>::Type Type;
    };


    template < typename TLCPInput, typename TDest >
    inline void _childtabProcess(TLCPInput &lcpIn, TDest &dest)
    {
        typedef typename Value<TLCPInput>::Type                TValue;
        typedef typename Size<TLCPInput>::Type                TSize;

        typedef Pair<TSize, TValue, Pack>                TPair;        // (i, lcptab[i])
        typedef std::stack<TPair>                         TStack;

        TStack stack_updown;
        TStack stack_nextl;

        #ifdef SEQAN_DEBUG_INDEX
            std::cerr << "--- CREATE CHILD TABLE ---" << std::endl;
            std::cerr << "Start stack-processing [pipelining]" << std::endl;
        #endif

        stack_updown.push(TPair(0, 0));
        stack_nextl.push(TPair(0, 0));

        dest.undefinedValue = TPair(MaxValue<TSize>::VALUE, 0);        // undefined value for unused entries
        resize(dest, length(lcpIn));
        beginRead(lcpIn);
        beginWrite(dest);

        TSize n = length(lcpIn);
        TSize undef = n + 1;
        TSize lastIndex_nextl;
        TPair lastIndex_updown = TPair(undef, 0);

        TValue lcp_i;
        for(TSize i = 0; i < n; ++i) {
            lcp_i = *lcpIn; ++lcpIn;

            //////////////////////////////////////////////////////////////////////////////
            // nextlIndex part
            lastIndex_nextl = undef;
            while (lcp_i < stack_nextl.top().i2)
                stack_nextl.pop();

            if (lcp_i == stack_nextl.top().i2) {
                lastIndex_nextl = stack_nextl.top().i1;
                push(dest, TPair(lastIndex_nextl, i + 1));                    // childtab[top].nextl

//                std::cerr << "nextl:\t\t\t" << TPair(lastIndex_nextl, i + 1) << std::endl;
                stack_nextl.pop();
            }
            stack_nextl.push(TPair(i + 1, lcp_i));


            //////////////////////////////////////////////////////////////////////////////
            // up/down part
            TPair top;
            while (lcp_i < stack_updown.top().i2) {
                lastIndex_updown = stack_updown.top();
                stack_updown.pop();

                top = stack_updown.top();
                if (lcp_i <= top.i2 && top.i2 != lastIndex_updown.i2 && top.i1 != lastIndex_nextl) {
                    push(dest, TPair(top.i1, lastIndex_updown.i1));            // childtab[top].down

//                    std::cerr << "down:\t\t" << TPair(top.i1, lastIndex_updown.i1) << std::endl;            // childtab[top].down
                }
            }
            if (lastIndex_updown.i1 != undef) {
                push(dest, TPair(i, lastIndex_updown.i1));            // childtab[top].up goes to [top - 1]

//                std::cerr << "up:\t" << TPair(i, lastIndex_updown.i1) << std::endl;
                lastIndex_updown.i1 = undef;
            }
            stack_updown.push(TPair(i + 1, lcp_i));
        }

        endRead(lcpIn);
        endWrite(dest);
    }


    //////////////////////////////////////////////////////////////////////////////
    // Enhanced class (outputs only the childtab (3rd) column of the Enhanced Suffix Array)
    template < typename TLCPInput >
    struct Pipe< TLCPInput, Childtab >
    {
        // *** SPECIALIZATION ***

        typedef typename Value<TLCPInput>::Type        TValue;
        typedef typename Size<TLCPInput>::Type        TSize;

        typedef Pair<TSize, TValue, Pack>        TCoreType;        // (i, lcptab[i])

        typedef Pool< TCoreType, MapperSpec< MapperConfigSize< filterI1<TCoreType>, TSize > > > TLinearMapper;
        typedef Pipe< TLinearMapper, Filter< filterI2<TCoreType> > > TFilter;

        TLinearMapper   mapper;
        TFilter            in;

        Pipe():
            in(mapper) {}

        Pipe(TLCPInput &_in):
            in(mapper)
        {
            process(_in);
        }

        template < typename TLcpInput_ >
        inline bool process(TLcpInput_ &_lcpIn) {

            // *** INSTANTIATION ***

            _childtabProcess(_lcpIn, mapper);
            return true;
        }

        inline typename Value<Pipe>::Type const operator*() const {
            return *in;
        }

        inline Pipe& operator++() {
            ++in;
            return *this;
        }
    };

    template < typename TInput, typename TLcpInput_ >
    inline bool operator<<(Pipe< TInput, Childtab > &me, TLcpInput_ const &in) {
         return me.process(in);
    }

    template < typename TLCPTable,
               typename TChildTable >
    void createChildtabExt(
        TChildTable &childtab,
        TLCPTable &lcp)
    {
        typedef Pipe< TLCPTable, Source<> >    TSource;
        typedef Pipe< TSource, Childtab >    TESA;

        TSource source(lcp);
        TESA    esa(source);

        childtab << esa;
    }

/*!
 * @fn IndexEsa#createChildtab
 * @headerfile <seqan/index.h>
 * @brief Creates a child table from a given lcp table.
 *
 * @signature void createChildtab(childTab, lcp[, algoTag]);
 *
 * @param[out] childTab A reference to the resulting child table.
 * @param[in]  lcp      A given lcp table.
 * @param[in]  algoTag  A tag that identifies the algorithm which is used for the creation.
 *
 * The size of <tt>childTab</tt> must be at least <tt>length(text)</tt> before calling this function.
 */
    template < typename TLCPTable,
               typename TValue,
               typename TConfig >
    inline void createChildtab(
        String<TValue, External<TConfig> > &childtab,
        TLCPTable &lcp)
    {
        createChildtabExt(childtab, lcp);
    }



    //////////////////////////////////////////////////////////////////////////////
    // internal childtab algorithm
    //////////////////////////////////////////////////////////////////////////////

    template < typename TLCPInput, typename TDest >
    inline void createChildtab(TDest &dest, TLCPInput const &lcpIn)
    {
        typedef typename Value<TLCPInput>::Type                TValue;
        typedef typename Size<TLCPInput>::Type                TSize;
        typedef typename Iterator<TLCPInput const>::Type    TIter;

        typedef Pair<TSize, TValue>                            TPair;        // (i, lcptab[i])
        typedef String<TPair>                                TStack;

        #ifdef SEQAN_DEBUG_INDEX
            std::cerr << "--- CREATE CHILD TABLE ---" << std::endl;
            std::cerr << "Start stack-processing [random access]" << std::endl;
        #endif

        TStack stack;

        appendValue(stack, TPair(0, 0));

        resize(dest, length(lcpIn));

        TSize n = length(lcpIn);
        TSize undef = n + 1;
        TPair lastIndex = TPair(undef, 0);

        TValue lcp_i;
        TIter lcpI = begin(lcpIn);

        for(TSize i = 0; i < n; ++i, ++lcpI)
        {
            lcp_i = *lcpI;

            TPair top;
            while (lcp_i < back(stack).i2)
            {
                lastIndex = back(stack);
                eraseBack(stack);

                top = back(stack);
                if (lcp_i < top.i2 && top.i2 != lastIndex.i2 /*&& top.i1 != lastIndex_nextl*/) {
                    dest[top.i1] = lastIndex.i1;            // childtab[top].down
//                    std::cerr << "down:\t" << TPair(top.i1, lastIndex.i1) << std::endl;            // childtab[top].down
                }
//                else
//                    if (lcp_i == top.i2 && top.i2 != lastIndex.i2 /*&& top.i1 != lastIndex_nextl*/)
//                        std::cerr << "down:\t" << TPair(top.i1, lastIndex.i1) << std::endl;            // childtab[top].down

            }
            if (lastIndex.i1 != undef)
            {
                dest[i] = lastIndex.i1;                        // childtab[top].up goes to [top - 1]

//                std::cerr << "up:  \t" << TPair(i+1, lastIndex.i1) << std::endl;
                lastIndex.i1 = undef;
            }
            if (lcp_i == back(stack).i2)
            {
                dest[back(stack).i1] = i + 1;                // childtab[top].nextl
//                std::cerr << "nextl:\t\t\t" << TPair(back(stack).i1, i + 1) << std::endl;
            }
            appendValue(stack, TPair(i + 1, lcp_i));
        }
    }

//}

}

#endif
