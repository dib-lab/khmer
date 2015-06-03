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
// Author: Stefan Aiche <aiche@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_FIND_WILD_SHIFTAND_H
#define SEQAN_HEADER_FIND_WILD_SHIFTAND_H

// uncomment this for detailed debug output
//#define SEQAN_WILD_SHIFTAND_DEBUG

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// ShiftAnd with Wildcards
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class WildShiftAndPattern
 * @extends Pattern
 * @headerfile <seqan/find.h>
 * @brief Exact string matching with wildcards using bit parallelism.
 *
 * The Shift-And algorithm is applicable to search small patterns in texts using small alphabets.
 *
 * @signature template <typename TNeedle>
 *            class Pattern<TNeedle, WildShiftAnd>;
 *
 * @tparam TNeedle The needle type.  Type: @link ContainerConcept @endlink.
 *
 * The supported wildcards are <tt>*</tt> (zero or more occurrence), <tt>+</tt> (one or more occurrences), <tt>?</tt>
 * (optional character), <tt>.</tt> (every character), character classes (e.g. <tt>[a-z]</tt>) and bounded repeats (e.g.
 * <tt>{n,m}</tt>).
 *
 * After the find-Method returned the Finder will point to the last position of the occurrence
 *
 * We encourage the user to intialize the Pattern with a <tt>String&lt;char&gt;</tt> (@link HostedConcept#setHost
 * @endlink or the C'tor). If you use for instance <tt>String&lt;Dna&gt;</tt> instead you won't be able to specify
 * wildcards.
 */

struct WildShiftAnd_;
typedef Tag<WildShiftAnd_> WildShiftAnd;

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
class Pattern<TNeedle, WildShiftAnd> {
//____________________________________________________________________________
public:
    typedef unsigned TWord;

    Holder<TNeedle> data_host;
    String<TWord> table;            // Look up table for each character in the alphabet (called B in "Navarro")

    String<TWord> s_table;            // marks all positions, that can remain active, after reading a specific character (called S in "Navarro")
    String<TWord> a_table;            // marks all positions of optional characters in the pattern (called A in "Navarro")
    String<TWord> i_table;            // marks all positions in the pattern, that preceed a block of optional characters (called I in "Navarro")
    String<TWord> f_table;            // marks all end-positions of blocks of optional characters (called F in "Navarro")

    String<TWord> prefSufMatch;        // Set of all the prefixes of needle that match a suffix of haystack (called D in "Navarro")
    String<TWord> df;                // additional bit mask to enable flooding of bits

    TWord needleLength;                // e.g., needleLength=33 --> blockCount=2 (iff w=32 bits)
    TWord character_count;            // number of normal characters in the needle
    TWord blockCount;                // #unsigned ints required to store needle

    bool _valid;                    // is the pattern valid or not

//____________________________________________________________________________

    Pattern()
        : _valid(false)
        {}

    template <typename TNeedle2>
    Pattern(TNeedle2 const & ndl)
        : _valid(false) {
        setHost(*this, ndl);
    }

//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

/*
// debug method to visualize bitmasks
inline void _printMask(String <unsigned> const &  mask,unsigned line,String <char> name)
{
    unsigned len = length(mask);
    std::cout << name << " " << line << "   ";
    for(unsigned int j=0;j<len;++j) {
        for(unsigned int bit_pos=0;bit_pos<BitsPerValue<unsigned int>::VALUE;++bit_pos) {
            std::cout << ((mask[j] & (1<<(bit_pos % BitsPerValue<unsigned int>::VALUE))) !=0);
        }
        std::cout << " ";
    }
    std::cout << std::endl;
}
*/

//____________________________________________________________________________
//                                VALIDATION

inline bool _isUnsigned(String<char> const & number)
{
//IOREV move to future is-module, maybe just use: return (strtol(number) >= 0) [if I understand this correctly)
    unsigned int len = length(number);
    for(unsigned int i = 0;i < len;++i){
        if(!(convert<unsigned int>(getValue(number,i)) <= 57 && convert<unsigned int>(getValue(number,i)) >= 47))
            return false;
    }
    return true;
}

// check if the pattern is valid
template <typename TNeedle2>
bool _validate(TNeedle2 const & needle)
{
    typedef unsigned TWord;

    TWord nl = length(needle);
    TWord len = nl;
    TWord i = 0;

    while (i < nl){
        if (convert<char>(getValue(needle,i)) == '['){
            // check if class ends
          while(i < nl && convert<char>(getValue(needle,i)) != ']')++i;
            if(i == nl)    return false;
        }
        else if(convert<char>(getValue(needle,i)) == '*' || convert<char>(getValue(needle,i)) == '+' || convert<char>(getValue(needle,i)) == '?'){
            // check for consecutive runs of +,*,?
            ++i;
            if(convert<char>(getValue(needle,i)) == '*' || convert<char>(getValue(needle,i)) == '+' || convert<char>(getValue(needle,i)) == '?')
                return false;
            ++i;
        }
        else if(convert<char>(getValue(needle,i)) == '{'){
SEQAN_CHECKPOINT
            String <char> number;
            TWord n,m;
            n = m = 0;
            --len;++i; // get to the first number
            while(i < nl && convert<char>(getValue(needle,i)) != '}' && convert<char>(getValue(needle,i)) != ',') {
SEQAN_CHECKPOINT
                append(number,convert<char>(getValue(needle,i)));
                --len;++i;
            }
            // class didn't ended
            if(i == nl)    return false;

            // isNumber
            if(!_isUnsigned(number)) return false;
            n = atoi(toCString(number));

            // check the second number
            if (convert<char>(getValue(needle,i)) == ','){
SEQAN_CHECKPOINT
                --len;++i;
                clear(number);
                while(i < nl && convert<char>(getValue(needle,i)) != '}') {
SEQAN_CHECKPOINT
                    append(number,convert<char>(getValue(needle,i)));
                    --len;++i;
                }

                // class didn't ended
                if(i == nl) return false;

                // isNumber
                if(!_isUnsigned(number)) return false;
                m = atoi(toCString(number));
                --len;
            }
            else --len;

            // check if optional part is greater or equal then fixed
            if(m < n && m != 0)    return false;
        }
        else if(convert<char>(getValue(needle,i)) == '\\'){
SEQAN_CHECKPOINT
            // check if there exists a next character
            if(i == nl - 1)    return false;
            else ++i;
        }
        ++i;
    }
    return true;
}

/*
  determine the original pattern length -- wo wildcard-chars
  currently supported +,*,[..],?,.,{n,m}
*/
template <typename TNeedle>
unsigned _lengthWithoutWildcards(TNeedle const & needle) {

    typedef unsigned TWord;

    TWord nl = length(needle);
    TWord len = nl;
    TWord i = 0;
    while(i < nl) {
        if(convert<char>(getValue(needle,i)) == '+'){
SEQAN_CHECKPOINT
            --len;
        }
        else if(convert<char>(getValue(needle,i)) == '*'){
SEQAN_CHECKPOINT
            --len;
        }
        else if(convert<char>(getValue(needle,i)) == '?'){
SEQAN_CHECKPOINT
            --len;
        }
        else if(convert<char>(getValue(needle,i)) == '['){
SEQAN_CHECKPOINT
            while(convert<char>(getValue(needle,i)) != ']') {
SEQAN_CHECKPOINT
                --len;++i;
            }
        }
        else if(convert<char>(getValue(needle,i)) == '{'){
SEQAN_CHECKPOINT
            String <char> number;
            TWord n,m;
            n = m = 0;
            --len;++i; // get to the first number
            while(convert<char>(getValue(needle,i)) != '}' && convert<char>(getValue(needle,i)) != ',') {
SEQAN_CHECKPOINT
                append(number,convert<char>(getValue(needle,i)));
                --len;++i;
            }
            // we also have to read the second number
            n = atoi(toCString(number));
            if (convert<char>(getValue(needle,i)) == ','){
SEQAN_CHECKPOINT
                --len;++i;
                //lets get m
                clear(number);
                while(convert<char>(getValue(needle,i)) != '}') {
SEQAN_CHECKPOINT
                    append(number,convert<char>(getValue(needle,i)));
                    --len;++i;
                }
                m = atoi(toCString(number));
                --len;
            }
            else
                --len;

            len += (m != 0 ? m : n) - 1;
        }
        else if(convert<char>(getValue(needle,i)) == '\\'){
SEQAN_CHECKPOINT
            --len;++i; // next character could be a \ too
        }
        ++i;
    }

    return len;
}

/*
    returns all character codes contained in the class
*/
template <typename TValue,typename TNeedle2>
String <unsigned> _getCharacterClass(TNeedle2 const & host,unsigned start,unsigned end){
    typedef unsigned TWord;
    String <unsigned> ret;
    unsigned pos = start;
    while (pos < end){
SEQAN_CHECKPOINT
        if(convert<char>(getValue(host,pos)) != '-')
            append(ret,convert<TWord>(convert<TValue>(getValue(host,pos))));
        else{
SEQAN_CHECKPOINT
            // could be a range
            if (pos > start && pos < end && (convert<TWord>(convert<TValue>(getValue(host,pos-1))) < convert<TWord>(convert<TValue>(getValue(host,pos+1)))) ){
                unsigned r_s = convert<TWord>(convert<TValue>(getValue(host,pos-1))) + 1;
                while(r_s < convert<TWord>(convert<TValue>(getValue(host,pos+1)))){
                    append(ret,r_s);
                    ++r_s;
                }
            }
            else // or simply a '-'
                append(ret,convert<TWord>(convert<TValue>(getValue(host,pos))));
        }
        ++pos;
    }
    return ret;
}


template <typename TNeedle, typename TNeedle2>
void setHost (Pattern<TNeedle, WildShiftAnd> & me, TNeedle2 const & needle) {
SEQAN_CHECKPOINT
    me._valid = _validate(needle);

    if(!valid(me)) return;

    typedef unsigned TWord;
    typedef typename Value<TNeedle>::Type TValue;

    me.needleLength = length(needle);
    me.character_count = _lengthWithoutWildcards(needle);

    if (me.character_count<1) me.blockCount=1;
    else me.blockCount=((me.character_count-1) / BitsPerValue<TWord>::VALUE)+1;

    clear(me.table);
    resize(me.table, me.blockCount * ValueSize<TValue>::VALUE, 0, Exact());

    clear(me.s_table);
    resize(me.s_table, me.blockCount * ValueSize<TValue>::VALUE, 0, Exact());

    clear(me.a_table);
    resize(me.a_table,me.blockCount,0,Exact());

    int i = -1;
    String <TWord> last_char; // stores the character (or characters) that were read in the last step
    TWord j=0;
    while(j < me.needleLength){
SEQAN_CHECKPOINT
        if (convert<char>(getValue(needle,j)) == '+'){
SEQAN_CHECKPOINT
            TWord len = length(last_char);
            for (unsigned int k = 0; k < len; ++k)
                me.s_table[me.blockCount*last_char[k] + i / BitsPerValue<TWord>::VALUE] |= (1<<(i%BitsPerValue<TWord>::VALUE));
        }
        else if (convert<char>(getValue(needle,j)) == '?'){
SEQAN_CHECKPOINT
            me.a_table[i / BitsPerValue<TWord>::VALUE] |= (1<<(i%BitsPerValue<TWord>::VALUE));
        }
        else if (convert<char>(getValue(needle,j)) == '*'){
SEQAN_CHECKPOINT
            TWord len = length(last_char);
            for (unsigned int k = 0; k < len; ++k)
                me.s_table[me.blockCount*last_char[k] + i / BitsPerValue<TWord>::VALUE] |= (1<<(i%BitsPerValue<TWord>::VALUE));
            me.a_table[i / BitsPerValue<TWord>::VALUE] |= (1<<(i%BitsPerValue<TWord>::VALUE));
        }
        else if(convert<char>(getValue(needle,j)) == '['){
SEQAN_CHECKPOINT
            /* find characters in class */
            TWord e = j;
            while(convert<char>(getValue(needle,e)) != ']') ++e;
            /* get character codes of class */
            last_char = _getCharacterClass<TValue>(needle,j+1,e);
            TWord len = length(last_char);

            /* add class to the mask */
            ++i;
            for (unsigned int k = 0; k < len; ++k){
SEQAN_CHECKPOINT
                me.table[me.blockCount*last_char[k] + i / BitsPerValue<TWord>::VALUE] |= (1<<(i%BitsPerValue<TWord>::VALUE));
            }
            j = e;
        }
        else if(convert<char>(getValue(needle,j)) == '.'){ // matches all characters in the current alphabet
SEQAN_CHECKPOINT
            clear(last_char);
            ++i;
            for(unsigned int l = 0;l < ValueSize<TValue>::VALUE;++l){
                append(last_char,l);
                me.table[me.blockCount*l + i / BitsPerValue<TWord>::VALUE] |= (1<<(i%BitsPerValue<TWord>::VALUE));
            }

        }
        else if(convert<char>(getValue(needle,j)) == '\\'){ // handle escape characters
SEQAN_CHECKPOINT
            /* goto next character use this for the bit mask */
            ++i;++j;
            clear(last_char);
            append(last_char, convert<TWord>(convert<TValue>(getValue(needle,j))));
            me.table[me.blockCount*last_char[0] + i / BitsPerValue<TWord>::VALUE] |= (1<<(i%BitsPerValue<TWord>::VALUE));
        }
        else if(convert<char>(getValue(needle,j)) == '{'){ // handle bounded character repeats
SEQAN_CHECKPOINT
            String <char> number;
            TWord n,m,r;
            TWord len = length(last_char);
            n = m = 0;
            ++j;
            while(convert<char>(getValue(needle,j)) != '}' && convert<char>(getValue(needle,j)) != ',') {
SEQAN_CHECKPOINT
                append(number,convert<char>(getValue(needle,j)));
                ++j;
            }
            n = atoi(toCString(number));
            if (convert<char>(getValue(needle,j)) == ','){
SEQAN_CHECKPOINT
                ++j;
                clear(number);
                while(convert<char>(getValue(needle,j)) != '}') {
SEQAN_CHECKPOINT
                    append(number,convert<char>(getValue(needle,j)));
                    ++j;
                }
                m = atoi(toCString(number));
            }
            // we already have seen one required occurrence of the character (last_char)
            n -= 1;
            r = 0;
            while(r < n){ // add n normal characters
SEQAN_CHECKPOINT
                ++i;
                for (unsigned int k = 0; k < len; ++k){
                    me.table[me.blockCount*last_char[k] + i / BitsPerValue<TWord>::VALUE] |= (1<<(i%BitsPerValue<TWord>::VALUE));
                }
                ++r;
            }
            ++r; // correct the -1 of n to get in the correct relation to m
            while (r < m){ // if there was no m specified this won't be used
SEQAN_CHECKPOINT
                // add m - n charaters and make them optional
                ++i;
                for (unsigned int k = 0; k < len; ++k){
SEQAN_CHECKPOINT
                    me.table[me.blockCount*last_char[k] + i / BitsPerValue<TWord>::VALUE] |= (1<<(i%BitsPerValue<TWord>::VALUE));
                }
                me.a_table[i / BitsPerValue<TWord>::VALUE] |= (1<<(i%BitsPerValue<TWord>::VALUE));
                ++r;
            }
        }

        else // we have a character here
        {
SEQAN_CHECKPOINT
            // determine character position in array table
            clear(last_char);
            append(last_char, convert<TWord>(convert<TValue>(getValue(needle,j))));
            ++i;
            me.table[me.blockCount*last_char[0] + i / BitsPerValue<TWord>::VALUE] |= (1<<(i%BitsPerValue<TWord>::VALUE));
        }
        ++j;
    }

    clear(me.i_table);
    resize(me.i_table,me.blockCount,0,Exact());

    clear(me.f_table);
    resize(me.f_table,me.blockCount,0,Exact());

    for (unsigned int i = 0; i < me.character_count; ++i){
SEQAN_CHECKPOINT
        if ((me.a_table[i / BitsPerValue<TWord>::VALUE] & (1 << (i % BitsPerValue<TWord>::VALUE))) != 0){
SEQAN_CHECKPOINT
            if ((me.f_table[i / BitsPerValue<TWord>::VALUE] & (1 << ((i-1) % BitsPerValue<TWord>::VALUE))) == 0){
SEQAN_CHECKPOINT

                if(i > 0)
                    me.i_table[(i-1) / BitsPerValue<TWord>::VALUE] |= 1 << ((i-1) % BitsPerValue<TWord>::VALUE);
                me.f_table[i / BitsPerValue<TWord>::VALUE] |= 1 << (i % BitsPerValue<TWord>::VALUE);
#ifdef SEQAN_WILD_SHIFTAND_DEBUG
                std::cout << "Update F and I" << std::endl;
                _printMask(me.f_table,0,"F ");
                _printMask(me.i_table,0,"I ");
                _printMask(me.a_table,0,"A ");
                std::cout << std::endl;
#endif
            }
            else{
SEQAN_CHECKPOINT
                TWord curBlock = i / BitsPerValue<TWord>::VALUE;
                for (unsigned int k = 0; k < me.blockCount; ++k){
SEQAN_CHECKPOINT
                    if(k != curBlock)
                        me.f_table[i / BitsPerValue<TWord>::VALUE] &= ~0;
                    else
                        me.f_table[i / BitsPerValue<TWord>::VALUE] &= ~(1 << ((i-1) % BitsPerValue<TWord>::VALUE));

                }
                //me.f_table[i / BitsPerValue<TWord>::VALUE] &= ~(1 << ((i-1) % BitsPerValue<TWord>::VALUE));
                me.f_table[i / BitsPerValue<TWord>::VALUE] |= 1 << (i % BitsPerValue<TWord>::VALUE);
#ifdef SEQAN_WILD_SHIFTAND_DEBUG
                std::cout << "Update F" << std::endl;
                _printMask(me.f_table,0,"F ");
                std::cout << std::endl;
#endif
            }
        }
    }

    setValue(me.data_host, needle);

#ifdef SEQAN_WILD_SHIFTAND_DEBUG
    // Debug code
    std::cout << "Alphabet size: " << ValueSize<TValue>::VALUE << std::endl;
    std::cout << "Needle length (with wildcards): " << me.needleLength << std::endl;
    std::cout << "Needle length (wo wildcards): " << me.character_count << std::endl;
    std::cout << "Block count: " << me.blockCount << std::endl;

    std::cout << "Needle:" << needle << std::endl;

    _printMask(me.f_table,0,"F ");
    _printMask(me.i_table,0,"I ");
    _printMask(me.a_table,0,"A ");
    std::cout << std::endl << std::endl;

    for(unsigned i=0;i<ValueSize<TValue>::VALUE;++i) {
        if (((i<97) && (4 < i) ) || (i>122)) continue;
        std::cout << static_cast<TValue>(i) << ": ";
        for(unsigned int j=0;j<me.blockCount;++j) {
            for(int bit_pos=0;bit_pos<BitsPerValue<unsigned>::VALUE;++bit_pos) {
                std::cout << ((me.table[me.blockCount*i+j] & (1<<(bit_pos % BitsPerValue<unsigned>::VALUE))) !=0);
            }
        }
        std::cout << std::endl;
    }
#endif

}

template <typename TNeedle, typename TNeedle2>
inline void setHost (Pattern<TNeedle, WildShiftAnd> & me, TNeedle2 & needle)
{
    setHost(me, reinterpret_cast<TNeedle2 const &>(needle));
}

//____________________________________________________________________________


template <typename TNeedle>
inline void _patternInit (Pattern<TNeedle, WildShiftAnd> & me)
{
SEQAN_CHECKPOINT
    clear(me.prefSufMatch);
    resize(me.prefSufMatch, me.blockCount, 0, Exact());

    clear(me.df);
    resize(me.df, me.blockCount, 0, Exact());

}

//____________________________________________________________________________
template <typename TNeedle>
inline bool valid(Pattern <TNeedle,WildShiftAnd> & me)
{
SEQAN_CHECKPOINT
    return me._valid;
}

//____________________________________________________________________________
template <typename TNeedle>
inline bool valid(Pattern <TNeedle,WildShiftAnd> const & me)
{
SEQAN_CHECKPOINT
    return me._valid;
}

//____________________________________________________________________________


template <typename TNeedle>
inline typename Host<Pattern<TNeedle, WildShiftAnd>const>::Type &
host(Pattern<TNeedle, WildShiftAnd> & me)
{
SEQAN_CHECKPOINT
    return value(me.data_host);
}

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, WildShiftAnd>const>::Type &
host(Pattern<TNeedle, WildShiftAnd> const & me)
{
SEQAN_CHECKPOINT
    return value(me.data_host);
}

//____________________________________________________________________________


template <typename TFinder, typename TNeedle>
inline bool _findShiftAndSmallNeedle(TFinder & finder, Pattern<TNeedle, WildShiftAnd> & me) {
SEQAN_CHECKPOINT
    typedef unsigned TWord;
    TWord compare = (1 << (me.character_count-1));
    while (!atEnd(finder)) {
SEQAN_CHECKPOINT
        TWord pos = convert<TWord>(*finder);
        /* added  | (me.prefSufMatch[0] & me.s_table[me.blockCount*pos]) at the end of the line */
        me.prefSufMatch[0] = (((me.prefSufMatch[0] << 1) | 1) & me.table[me.blockCount*pos]) | (me.prefSufMatch[0] & me.s_table[me.blockCount*pos]) ;

        /* additional bit operations */
        me.df[0] = me.prefSufMatch[0] | me.f_table[0];
        me.prefSufMatch[0] |= ((me.a_table[0] & (~(me.df[0] - me.i_table[0]))) ^ me.df[0]);
        if ((me.prefSufMatch[0] & compare) != 0) {
SEQAN_CHECKPOINT
            return true;
        }
        goNext(finder);
    }
    return false;
}

template <typename TFinder, typename TNeedle>
inline bool _findShiftAndLargeNeedle(TFinder & finder, Pattern<TNeedle, WildShiftAnd> & me) {
SEQAN_CHECKPOINT
    typedef unsigned TWord;
    const TWord all1 = ~0;
    TWord compare = (1 << ((me.character_count-1) % BitsPerValue<TWord>::VALUE));

    while (!atEnd(finder)) {
SEQAN_CHECKPOINT
        TWord pos = convert<TWord>(*finder);
        TWord carry = 1;
        TWord wc_carry = 0;
        // shift of blocks with carry
        for(TWord block=0;block<me.blockCount;++block) {
SEQAN_CHECKPOINT
            bool newCarry = ((me.prefSufMatch[block] & (1<< (BitsPerValue<TWord>::VALUE - 1)))!=0);
            me.prefSufMatch[block] = (((me.prefSufMatch[block] << 1) | carry) & me.table[me.blockCount*pos+block]) | (me.prefSufMatch[block] & me.s_table[me.blockCount*pos+block]) ;
            carry = newCarry;

            me.df[block] = me.prefSufMatch[block] | me.f_table[block];
            TWord Z = me.df[block] - me.i_table[block] - wc_carry;
            wc_carry = ((me.df[block] < Z) || (me.i_table[block]==all1 && wc_carry)) ? 1 : 0;
            me.prefSufMatch[block] |= (me.a_table[block] & (~Z ^ me.df[block]));
        }

#ifdef SEQAN_WILD_SHIFTAND_DEBUG
        std::cout << "reading " << *finder << std::endl;
        _printMask(me.prefSufMatch,position(finder),"D ");
        _printMask(me.df,position(finder),"Df");
        std::cout << std::endl;
#endif

        // check for match
        if ((me.prefSufMatch[me.blockCount-1] & compare) != 0)
            return true;

        goNext(finder);
    }
    return false;
}

template <typename TFinder, typename TNeedle>
inline bool find(TFinder & finder, Pattern<TNeedle, WildShiftAnd> & me) {
SEQAN_CHECKPOINT

    if (empty(finder)) {
        _patternInit(me);
        _finderSetNonEmpty(finder);

        if(!valid(me))
            return false;
    }
    else
        goNext(finder);

    // Fast algorithm for needles < machine word?
    if (me.blockCount == 1) {
        return _findShiftAndSmallNeedle(finder, me);
    } else {
        return _findShiftAndLargeNeedle(finder, me);
    }
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_FIND_WILD_SHIFTAND_H
