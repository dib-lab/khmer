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
// Author: Stephan Aiche <aiche@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_FIND_ABNDMALGO_H
#define SEQAN_HEADER_FIND_ABNDMALGO_H

// uncomment this for verbose output of the ABNDM ALGO
//#define SEQAN_DEBUG_ABNDM

namespace SEQAN_NAMESPACE_MAIN
{


#ifdef SEQAN_DEBUG_ABNDM
inline void _printMask(String <unsigned> const &  mask,String <char> name)
{
    unsigned len = length(mask);
    std::cout << name << ": ";
    for(unsigned int j=0;j<len;++j) {
        for(unsigned int bit_pos=0;bit_pos<BitsPerValue<unsigned int>::VALUE;++bit_pos) {
            std::cout << ((mask[j] & (1<<(bit_pos % BitsPerValue<unsigned int>::VALUE))) !=0);
        }
        std::cout << " ";
    }
    std::cout << std::endl;
}

inline void _printMask(String <unsigned> const &  mask,unsigned start, unsigned len,String <char> name)
{
    std::cout << name << ": ";
    for(unsigned int j=start;j<start + len;++j) {
        for(unsigned int bit_pos=0;bit_pos<BitsPerValue<unsigned int>::VALUE;++bit_pos) {
            std::cout << ((mask[j] & (1<<(bit_pos % BitsPerValue<unsigned int>::VALUE))) !=0);
        }
        std::cout << " ";
    }
    std::cout << std::endl;
}

#endif

//////////////////////////////////////////////////////////////////////////////
// AbndmAlgo
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class AbndmAlgoPattern
 * @extends Pattern
 * @headerfile <seqan/find.h>
 * @brief Approximate Backward Nondeterministic Dawg Matching algorithm.
 *
 * Approximate string matching using bit parallelism.
 *
 * @signature template <typename TNeedle>
 *            class Pattern<TNeedle, AbndmAlgo>;
 *
 * @tparam TNeedle The needle type.  Type @link ContainerConcept @endlink.
 *
 * @note The types of the needle and the haystack have to match.
 */


struct AbndmAlgo;

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
struct FindBeginPatternSpec< Pattern<TNeedle, AbndmAlgo> >:
    DefaultFindBeginPatternSpec<>
{
};

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
class Pattern<TNeedle, AbndmAlgo>:
    public FindBegin_<Pattern<TNeedle, AbndmAlgo> >
{
//////////////////////////////////////////////////////////////////////////////
public:
    typedef unsigned int TWord;
    Holder<TNeedle> data_host;

    String<TWord>       b_table;    // called B in Navarro
    String<TWord>       r_table;    // called R in Navarro

    TWord blockCount;
    TWord last;
    TWord needleLength;                // e.g., needleLength=33 --> blockCount=2 (iff w=32 bits)
    TWord haystackLength;
    TWord limit;                    // the maximal accepted error
    TWord cP;  // save current position of the nfa

    bool findNext;

    Pattern<TNeedle,MyersUkkonen> verifier;
//////////////////////////////////////////////////////////////////////////////

    Pattern() :
        blockCount(0), last(0), needleLength(0), haystackLength(0), limit(1), cP(0), findNext(false)
    {}

    template <typename TNeedle2>
    Pattern(TNeedle2 const & ndl) :
        blockCount(0), last(0), needleLength(0), haystackLength(0), limit(1), cP(0), findNext(false),
        verifier(ndl,-1)
    {
        setHost(*this, ndl);
    }

    template <typename TNeedle2>
    Pattern(TNeedle2 const & ndl, int _limit = -1) :
        limit(- _limit), cP(0), verifier(ndl,_limit)
    {
        setHost(*this, ndl);
    }
};

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
void _printR(Pattern<TNeedle, AbndmAlgo> & me)
{
    std::cout << "R ----------------------------- " << std::endl;
    for(unsigned i = 0;i <= me.limit;++i)
    {
        _printMask(me.r_table,me.blockCount * i, me.blockCount ," ");
    }
    std::cout << " ------------------------------ " << std::endl;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TNeedle2>
void setHost (Pattern<TNeedle, AbndmAlgo> & me, TNeedle2 const& needle)
{
SEQAN_CHECKPOINT
    typedef unsigned int TWord;
    typedef typename Value<TNeedle>::Type TValue;

    me.cP = 0;
    me.findNext = false;

    me.needleLength = length(needle);
    if (me.needleLength<1)
        me.blockCount = 1;
    else
        me.blockCount = ((me.needleLength-1) / BitsPerValue<TWord>::VALUE)+1;

    clear(me.b_table);
    resize(me.b_table, me.blockCount * ValueSize<TValue>::VALUE, 0, Exact());

    for (TWord j = 0; j < me.needleLength; ++j) {
        // Determine character position in array table
        TWord pos = convert<TWord>(getValue(needle,j));
        me.b_table[me.blockCount*pos + j / BitsPerValue<TWord>::VALUE] |= (1<<(j%BitsPerValue<TWord>::VALUE));
    }

#ifdef SEQAN_DEBUG_ABNDM
    std::cout << "Needle:   " << needle << std::endl;
    std::cout << "|Needle|: " << length(needle) << std::endl;
    std::cout << "Alphabet size: " << ValueSize<TValue>::VALUE << std::endl;

    for(unsigned i=0;i<ValueSize<TValue>::VALUE;++i) {
        if (((i<97) && (4 < i) ) || (i>122)) continue;
        std::cout << static_cast<TValue>(i) << ": ";
        for(unsigned int j=0;j<me.blockCount;++j) {
            for(int bit_pos=0;bit_pos<BitsPerValue<unsigned>::VALUE;++bit_pos) {
                std::cout << ((me.b_table[me.blockCount*i+j] & (1<<(bit_pos % BitsPerValue<unsigned>::VALUE))) !=0);
            }
            std::cout << " ";
        }
        std::cout << std::endl;
    }
#endif
    clear(me.r_table); // init is only possible if we know the the error count

    me.data_host = needle;
    setHost(me.verifier,needle);
}

template <typename TNeedle, typename TNeedle2>
void setHost (Pattern<TNeedle, AbndmAlgo> & me, TNeedle2 & needle)
{
    setHost(me, reinterpret_cast<TNeedle2 const &>(needle));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
inline void _patternInit (Pattern<TNeedle, AbndmAlgo> & me)
{
    SEQAN_CHECKPOINT
    clear(me.r_table);
    resize(me.r_table, me.blockCount * (me.limit + 1), 0, Exact());
    me.findNext = false;
    me.last = 0;
    _findBeginInit(me, needle(me));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, AbndmAlgo> >::Type &
host(Pattern<TNeedle, AbndmAlgo> & me)
{
    SEQAN_CHECKPOINT
    return value(me.data_host);
}

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, AbndmAlgo> const>::Type &
host(Pattern<TNeedle, AbndmAlgo> const & me)
{
    SEQAN_CHECKPOINT
    return value(me.data_host);
}

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn AbndmAlgoPattern#getScore
 * @headerfile <seqan/find.h>
 * @brief Score of the last found match in approximate searching.
 *
 * @signature TScoreValue getScore(pattern);
 *
 * @param[in] pattern A abndmAlgo pattern that can be used for approximate searching.
 *
 * @return TScoreValue The score of the last match found using <tt>pattern</tt>.  If no match was found then the value
 *                     is undefined.
 */


template <typename TNeedle>
int getScore(Pattern<TNeedle, AbndmAlgo > & me)
{
    SEQAN_CHECKPOINT
    return getScore(me.verifier);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFinder, typename TNeedle>
inline bool _findAbndmSmallNeedle(TFinder & finder, Pattern<TNeedle, AbndmAlgo> & me)
{
    SEQAN_CHECKPOINT
    typedef unsigned int TWord;

    typedef typename Host<TFinder>::Type    THost;
    typedef Segment<THost>                  THostSegment;
    typedef Finder<THostSegment>            THSFinder;

    TWord j,nR,oR,i,cB;
    TWord startPos = position(finder);

    // if returned from previous search we need to check the current position for additional matches
    if(me.findNext)
    {
        // reset the finder
        TWord offset = startPos - me.cP;
        finder -= offset;
        int end = me.cP + me.needleLength + me.limit;
        // adjust end if it points over the edges of host(finder)
        end = (end > static_cast<int>(length(host(finder))) ? static_cast<int>(length(host(finder))) : end);

        THostSegment s(infix(host(finder),me.cP, end));
        THSFinder f(s);
#ifdef SEQAN_DEBUG_ABNDM
        std::cout << "additional verification of current position" << std::endl;
        std::cout << "initialized on: " << s << std::endl;
        std::cout << "parameters are: " << me.cP << " " << (me.cP + me.needleLength + me.limit) << std::endl;
        std::cout << "original start pos: " << startPos << std::endl;
#endif

        while(find(f,me.verifier,- (int) me.limit)){
            TWord newP = position(finder) + position(f);
#ifdef SEQAN_DEBUG_ABNDM
            std::cout << "found on new position: " << newP << std::endl;
#endif
            if(newP > startPos){
                finder += position(f);
                return true;
            }
        }
#ifdef SEQAN_DEBUG_ABNDM
        std::cout << "additional verification done .. shifted by last=" << me.last << std::endl << std::endl;
#endif
        finder += me.last;
        me.cP += me.last;
    }

    me.cP = position(finder);
    // walk on
    while (position(finder) <= me.haystackLength - me.needleLength)
    {
        SEQAN_CHECKPOINT
            j = me.needleLength - me.limit - 1;
        me.last = j;

        // init R_0
        TWord pos = convert<TWord>(*(finder + j));

        me.r_table[0] = me.b_table[pos];
        nR = (TWord)-1;
    for(i = 1;i <= me.limit;++i) me.r_table[i] = nR;

#ifdef SEQAN_DEBUG_ABNDM
        std::cout << "reading " << *(finder + j) << std::endl;
        _printMask(me.b_table[pos],"");
        _printR(me);
#endif
        // process R_1 .. R_limit
        while((nR != 0) && (j != 0))
        {
            pos = convert<TWord>(*(finder + j - 1));
            cB = me.b_table[pos];
            oR = me.r_table[0];
            nR = (oR >> 1) & cB;
            me.r_table[0] = nR;
            for(i = 1;i <= me.limit;++i)
            {
                nR = ((me.r_table[i] >> 1) & cB) | oR | ((oR | nR) >> 1);
                oR = me.r_table[i];
                me.r_table[i] = nR;
            }
#ifdef SEQAN_DEBUG_ABNDM
            std::cout << "reading " << *(finder + j - 1) << std::endl;
            _printMask(me.b_table[pos],"");
            _printR(me);
#endif
            --j;
            if((nR & 1) != 0){
                if(j > 0){
#ifdef SEQAN_DEBUG_ABNDM
                    std::cout << "last bit is set so last is set to " << j << std::endl;
#endif
                    me.last = j;
                }
                else // verify
                {
                    // call find
                    int end = me.cP + me.needleLength + me.limit;
                    // adjust end if it points over the edges of host(finder)
                    end = (end > static_cast<int>(length(host(finder))) ? static_cast<int>(length(host(finder))) : end);

                    THostSegment s(infix(host(finder),me.cP, end));
                    THSFinder f(s);

#ifdef SEQAN_DEBUG_ABNDM
                    std::cout << "found verifiable prefix" << std::endl;
                    // init finder on infix
                    std::cout << "parameters are " << me.cP << " " << (me.cP + me.needleLength + me.limit) << std::endl;
            std::cout << "init finder on: " << s << std::endl;
#endif
                    // try to find the sequence
                    while(find(f,me.verifier,- (int) me.limit)){
                        TWord nP = position(finder) + position(f);
                        if(nP > startPos){
                            finder += position(f);
                            me.findNext = true;
                            return true;
                        }
                    }
                }
            }
        }
#ifdef SEQAN_DEBUG_ABNDM
        std::cout << "automaton runs out of active stats so window is shifted by last=" << me.last << std::endl << std::endl;
#endif

        finder += me.last;
        me.cP += me.last;
    }
    return false;
}

template <typename TFinder, typename TNeedle>
inline bool _findAbndmLargeNeedle(TFinder & finder, Pattern<TNeedle, AbndmAlgo> & me)
{
    SEQAN_CHECKPOINT
        typedef unsigned int TWord;
    TWord carryPattern = ((TWord)1 << (BitsPerValue<TWord>::VALUE - 1));
    typedef typename Host<TFinder>::Type    THost;
    typedef Segment<THost>                  THostSegment;
    typedef Finder<THostSegment>            THSFinder;

    TWord j,i;
    String<TWord> nR,oR;
    TWord startPos = position(finder);

    // if returned from previous search we need to check the current position for additional matches
    if(me.findNext)
    {
        // reset the finder
        TWord offset = startPos - me.cP;
        finder -= offset;

        int end = me.cP + me.needleLength + me.limit;
        // adjust end if it points over the edges of host(finder)
        end = (end > static_cast<int>(length(host(finder))) ? static_cast<int>(length(host(finder))) : end);

        THostSegment s(infix(host(finder),me.cP, end));
        THSFinder f(s);
#ifdef SEQAN_DEBUG_ABNDM
        std::cout << "additional verification of current position" << std::endl;
        std::cout << "initialized on: " << s << std::endl;
        std::cout << "parameters are: " << me.cP << " " << (me.cP + me.needleLength + me.limit) << std::endl;
        std::cout << "original start pos: " << startPos << std::endl;
#endif

        while(find(f,me.verifier,- (int) me.limit)){
            TWord newP = position(finder) + position(f);
#ifdef SEQAN_DEBUG_ABNDM
            std::cout << "found on new position: " << newP << std::endl;
#endif
            if(newP > startPos){
                finder += position(f);
                return true;
            }
        }
#ifdef SEQAN_DEBUG_ABNDM
        std::cout << "additional verification done .. shifted by last=" << me.last << std::endl << std::endl;
#endif
        finder += me.last;
        me.cP += me.last;
    }else me.cP = position(finder);
    // walk on
    while (position(finder) <= me.haystackLength - me.needleLength)
    {
        SEQAN_CHECKPOINT
            j = me.needleLength - me.limit - 1;
        me.last = j;

#ifdef SEQAN_DEBUG_ABNDM
    std::cout << "starting new frame ending at position " << me.cP + j << std::endl;
#endif
        // init R_0
        TWord pos = convert<TWord>(*(finder + j));
    for(int block = me.blockCount - 1;block >= 0;--block){
            me.r_table[block] = me.b_table[pos * me.blockCount + block];
    }
        //me.r_table[0] = me.b_table[pos];

    resize(nR,me.blockCount,0 - 1);
    resize(oR,me.blockCount,0);
    // nR = 0 - 1;
    //for(i = me.blockCount;i <= me.limit * me.blockCount;++i) me.r_table[i] = 0 - 1;//set all states in r_table to 1
    for(i = 1;i <= me.limit;++i){
            for(int block = me.blockCount - 1;block >= 0;--block){
                me.r_table[me.blockCount * i + block] = 0 -1 ;
            }
    }

    // track active states in nR
    TWord nRTrack = 1;

#ifdef SEQAN_DEBUG_ABNDM
        std::cout << "reading " << *(finder + j) << std::endl;
        _printMask(me.b_table,pos * me.blockCount, me.blockCount," ");
        _printR(me);
#endif
        while((nRTrack != 0) && (j != 0))
        {
            // get current letter
            pos = convert<TWord>(*(finder + j - 1));

            bool carry = 0;
            for(int block = me.blockCount -1;block >= 0;--block){
                oR[block] = me.r_table[block];
                bool nCarry=((oR[block] & 1) != 0);
                TWord toR = oR[block] >> 1;
                if(carry) toR |= carryPattern;
                carry = nCarry;
                nR[block] = toR & me.b_table[pos*me.blockCount + block];
                me.r_table[block] = nR[block];
            }


            for(i = 1;i <= me.limit;++i){
                bool rCarry = 0;
                bool nCarry = 0;
                for(int block = me.blockCount - 1;block >= 0;--block){
                    bool nRCarry = ((me.r_table[i * me.blockCount + block] & 1) != 0);
                    bool nNCarry = (((oR[block] | nR[block]) & 1) != 0);

                    TWord rPattern = (rCarry ? carryPattern : 0);
                    TWord nPattern = (nCarry ? carryPattern : 0);
                    nR[block] = (((me.r_table[i * me.blockCount + block] >> 1) | rPattern) & me.b_table[pos*me.blockCount + block]) | (((oR[block] | nR[block]) >> 1) | nPattern) | oR[block];
                    rCarry = nRCarry;
                    nCarry = nNCarry;
                    oR[block] = me.r_table[i  * me.blockCount + block];
                    me.r_table[i * me.blockCount + block] = nR[block];
                }
            }
#ifdef SEQAN_DEBUG_ABNDM
            std::cout << "reading " << *(finder + j - 1) << std::endl;
         _printMask(me.b_table,pos*me.blockCount, me.blockCount," ");
            _printR(me);
#endif
            --j;
            // track active states in nR
            nRTrack = 0;
            for(i = 0;i < me.blockCount;++i) nRTrack |= nR[i];

            if((nR[0] & 1) != 0){
                if(j > 0){
#ifdef SEQAN_DEBUG_ABNDM
                    std::cout << "last bit is set so last is set to " << j << std::endl;
#endif
                    me.last = j;
                }
                else // verify
                {
                    // call find
                    int end = me.cP + me.needleLength + me.limit;
                    // adjust end if it points over the edges of host(finder)
                    end = (end > static_cast<int>(length(host(finder))) ? static_cast<int>(length(host(finder))) : end);

                    THostSegment s(infix(host(finder),me.cP, end));
                    THSFinder f(s);

#ifdef SEQAN_DEBUG_ABNDM
                    std::cout << "found verifiable prefix" << std::endl;
                    // init finder on infix
                    std::cout << "parameters are " << me.cP << " " << (me.cP + me.needleLength + me.limit) << std::endl;
                    std::cout << "init finder on: " << s << std::endl;
#endif
                    // try to find the sequence
                    while(find(f,me.verifier,- (int) me.limit)){
                        TWord nP = position(finder) + position(f);
                        if(nP > startPos){
                            finder += position(f);
#ifdef SEQAN_DEBUG_ABNDM
                            std::cout << "found pattern at position " << position(finder) << std::endl;
#endif
                            me.findNext = true;
                            return true;
                        }
                    }
                }
            }
        }
#ifdef SEQAN_DEBUG_ABNDM
        std::cout << "automaton runs out of active stats so window is shifted by last=" << me.last << std::endl << std::endl;
#endif

        finder += me.last;
        me.cP += me.last;
    }
    return false;
}

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn AbndmAlgoPattern#scoreLimit
 * @headerfile <seqan/find.h>
 * @brief The minimal score a match must reach in approximate searching.
 *
 * @signature TScoreValue scoreLimit(pattern);
 *
 * @param[in] pattern The AbndmAlgoPattern to query.
 *
 * @return TScoreValue The score limit value.
 */
template <typename TNeedle>
inline int
scoreLimit(Pattern<TNeedle, AbndmAlgo > const & me)
{
    SEQAN_CHECKPOINT
    return - (int) me.limit;
}


//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn AbndmAlgoPattern#setSoreLimit
 * @headerfile <seqan/find.h>
 * @brief Set the minimal score a match must reach in approximate serach.
 *
 * @signature void setScoreLimit(pattern, limit);
 *
 * @param[in,out] pattern The AbndmAlgoPattern to set the limit for.
 * @param[in]     limit   The limit score value to set.
 *
 * @return TScoreValue The score limit value.
 */

template <typename TNeedle, typename TScoreValue>
inline void
setScoreLimit(Pattern<TNeedle, AbndmAlgo > & me,
              TScoreValue _limit)
{
    SEQAN_CHECKPOINT
        setScoreLimit(me.verifier,_limit);
    me.limit = (- _limit);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFinder, typename TNeedle>
inline bool find (TFinder & finder,
                  Pattern<TNeedle, AbndmAlgo > & me)
{
    SEQAN_CHECKPOINT
    if (empty(finder)) {
            _patternInit(me);
            _finderSetNonEmpty(finder);
            me.haystackLength = length(container(finder));
    }

    bool ret;
    if (me.blockCount == 1) {
        ret = _findAbndmSmallNeedle(finder, me);
    } else {
        ret = _findAbndmLargeNeedle(finder, me);
    }
    if (ret)
    {
        _setFinderEnd(finder);
    }
    return ret;

}

template <typename TFinder, typename TNeedle>
inline bool find (TFinder & finder,
                  Pattern<TNeedle, AbndmAlgo > & me,
                  int const k)
{
    SEQAN_CHECKPOINT
    setScoreLimit(me, k);
    return find(finder, me);
}

//////////////////////////////////////////////////////////////////////////////
}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_FIND_ABNDMALGO_H
