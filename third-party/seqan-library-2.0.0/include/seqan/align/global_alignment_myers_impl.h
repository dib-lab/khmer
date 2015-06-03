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
// Author: Stephan Aiche <stephan.aiche@fu-berlin.de>
// ==========================================================================

// TODO(holtgrew): Should be called _globalAlignmentScore()!

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_GLOBAL_ALIGNMENT_MYERS_IMPL_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_GLOBAL_ALIGNMENT_MYERS_IMPL_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// When using different alphabets, we will internally use the pattern alphabet
// for the comparison.  This means that the text character is converted to the
// pattern alphabet.  Here, the "pattern" is the shorter source sequence.

template <typename TAlphabetH, typename TSpecH, typename TAlphabetV, typename TSpecV>
int
_globalAlignmentScore(String<TAlphabetH, TSpecH> const & seqH,
                      String<TAlphabetV, TSpecV> const & seqV,
                      MyersBitVector const & algorithmTag)
{
    // Switch horizontal and vertical gap roles, gapsV should be the shorter one
    // to fit into less words.
    if (length(seqH) < length(seqV))
        return _globalAlignmentScore(seqV, seqH, algorithmTag);

    // Use size of unsigned int as blocksize for bit-vectors.
    const unsigned int BLOCK_SIZE = BitsPerValue<unsigned int>::VALUE;

    typedef String<TAlphabetH, TSpecH> const TSequenceH;
    typedef String<TAlphabetV, TSpecV> const TSequenceV;

    typedef typename Value<TSequenceV>::Type TPatternAlphabet;
    typedef typename Size<TSequenceH>::Type  TSourceSize;

    TSequenceH const & x = seqH;
    TSequenceV const & y = seqV;

    TSourceSize len_x = length(x);
    unsigned int pos = 0;

    // init variables
    unsigned int len_y = length(y);
    int score = (-1)*len_y;
    unsigned int patternAlphabetSize = ValueSize<TPatternAlphabet>::VALUE;
    unsigned int blockCount = (len_y + BLOCK_SIZE - 1) / BLOCK_SIZE;

    unsigned int scoreMask = 1 << ((len_y % BLOCK_SIZE) - 1);    // the mask with a bit set at the position of the last active cell

    String<unsigned> VP;
    resize(VP, blockCount, maxValue<unsigned>());
    String<unsigned> VN;
    resize(VN, blockCount, 0);
    String<unsigned> bitMask;
    resize(bitMask, patternAlphabetSize * blockCount, 0);

    // encoding the letters as bit-vectors
    for (unsigned int j = 0; j < len_y; j++)
        bitMask[blockCount * ordValue(getValue(y,j)) + j/BLOCK_SIZE] = bitMask[blockCount * ordValue(getValue(y,j)) + j/BLOCK_SIZE] | 1 << (j%BLOCK_SIZE);

    // compute score
    unsigned int X, D0, HN, HP;
    if(blockCount == 1)
    {
        while (pos < len_x) {
            X = bitMask[ordValue(static_cast<TPatternAlphabet>(getValue(x,pos)))] | VN[0];

            D0 = ((VP[0] + (X & VP[0])) ^ VP[0]) | X;
            HN = VP[0] & D0;
            HP = VN[0] | ~(VP[0] | D0);

            // customized to compute edit distance
            X = (HP << 1) | 1;
            VN[0] = X & D0;
            VP[0] = (HN << 1) | ~(X | D0);

            if (HP & scoreMask)
                score--;
            else if (HN & scoreMask)
                score++;

            ++pos;
        }
    } // end compute score - short pattern
    else
    {
        unsigned int temp, shift, currentBlock;
        unsigned int carryD0, carryHP, carryHN;

        while (pos < len_x)
        {
            // set vars
            carryD0 = carryHP = carryHN = 0;
            shift = blockCount * ordValue(static_cast<TPatternAlphabet>(getValue(x,pos)));

            // computing first the top most block
            X = bitMask[shift] | VN[0];

            temp = VP[0] + (X & VP[0]);
            carryD0 = temp < VP[0];

            D0 = (temp ^ VP[0]) | X;
            HN = VP[0] & D0;
            HP = VN[0] | ~(VP[0] | D0);

            // customized to compute edit distance
            X = (HP << 1) | 1;
            carryHP = HP >> (BLOCK_SIZE - 1);

            VN[0] = X & D0;

            temp = (HN << 1);
            carryHN = HN >> (BLOCK_SIZE - 1);

             VP[0] = temp | ~(X | D0);

            // computing the necessary blocks, carries between blocks following one another are stored
            for (currentBlock = 1; currentBlock < blockCount; currentBlock++) {
                X = bitMask[shift + currentBlock] | VN[currentBlock];

                temp = VP[currentBlock] + (X & VP[currentBlock]) + carryD0;

                carryD0 = ((carryD0) ? temp <= VP[currentBlock] : temp < VP[currentBlock]);

                D0 = (temp ^ VP[currentBlock]) | X;
                HN = VP[currentBlock] & D0;
                HP = VN[currentBlock] | ~(VP[currentBlock] | D0);

                X = (HP << 1) | carryHP;
                carryHP = HP >> (BLOCK_SIZE-1);

                VN[currentBlock] = X & D0;

                temp = (HN << 1) | carryHN;
                carryHN = HN >> (BLOCK_SIZE - 1);

                 VP[currentBlock] = temp | ~(X | D0);
            }

            // update score with the HP and HN values of the last block the last block
            if (HP & scoreMask)
                score--;
            else if (HN & scoreMask)
                score++;
            ++pos;
        }

    }  // end compute score - long pattern

    return score;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_GLOBAL_ALIGNMENT_MYERS_IMPL_H_
