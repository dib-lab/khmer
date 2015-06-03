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

#ifndef SEQAN_HEADER_RADIX_H
#define SEQAN_HEADER_RADIX_H

namespace SEQAN_NAMESPACE_MAIN
{

    // stably sort a[0..n-1] to b[0..n-1] with keys in 0..K-1 from r
    template <
        typename TSortedArray,
        typename TUnsortedArray,
        typename TCountArray,
        typename TText >
    void radixPass(
        TSortedArray &b,            // result array (sorted by 1 character)
        TUnsortedArray const &a,    // source array (unsorted)
        TText const &r,                // text to compare with
        TCountArray &c,                // temp. counter array
        unsigned K)                    // alphabet size
    {
        typedef typename Value<TCountArray>::Type TSize;

        TSize i, sum = 0, n = length(a);
        arrayFill(begin(c, Standard()), begin(c, Standard()) + K, 0);    // reset counters

        for (i = 0;  i < n;  i++)                                        // count occurrences
            c[ordValue(getValue(r,getValue(a,i)))]++;

        for (i = 0;  i < K;  i++)
        {                                                                // exclusive prefix sums
            TSize t = getValue(c,i);
            c[i] = sum;
            sum += t;
        }

        for (i = 0;  i < n;  i++)
        {
            TSize j = getValue(a,i);                                    // sort
            b[c[ordValue(getValue(r,j))]++] = j;
        }
    }

    template <
        typename TSortedArray,
        typename TUnsortedArray,
        typename TCountArray,
        typename TText,
        typename TShift >
    void radixPass(
        TSortedArray &b,            // result array (sorted by 1 character)
        TUnsortedArray const &a,    // source array (unsorted)
        TText const &r,                // text to compare with
        TCountArray &c,                // temp. counter array
        unsigned K,                    // alphabet size
        TShift shift)                // shift value
    {
        typedef typename Value<TCountArray>::Type TSize;

        TSize i, sum = 0, n = length(a), sn = length(r);
        arrayFill(begin(c, Standard()), begin(c, Standard()) + K, 0);    // reset counters

        for (i = 0;  i < n;  i++)
        {
            TSize j = getValue(a,i) + shift;                            // count occurrences
            if (j < sn) c[ordValue(getValue(r,j))]++;
            else        sum++;
        }

        for (i = 0;  i < K;  i++)
        {                                                                // exclusive prefix sums
            TSize t = getValue(c,i);
            c[i] = sum;
            sum += t;
        }

        for (i = 0, sum = 0;  i < n;  i++)
        {                                                                // sort
            TSize j = getValue(a,i);
            TSize k = j + shift;
            if (k < sn)
                b[c[ordValue(getValue(r,k))]++] = j;    // On Exception: Make sure, you have resized your sufarray
            else
                b[sum++] = j;                            // before calling createSuffixArray(..) to length(text)?
        }
    }

    // stably sort a[0..n-1] to b[0..n-1] with keys in 0..K-1 from r
    template <
        typename TSortedArray,
        typename TUnsortedArray,
        typename TCountArray,
        typename TText >
    void radixExtend(
        TSortedArray &b,            // result array (sorted by 1 character)
        TUnsortedArray const &a,    // source array (unsorted)
        TText const &r,                // text to compare with
        TCountArray &c,                // temp. counter array
        unsigned K)                    // alphabet size
    {
        typedef typename Value<TCountArray>::Type TSize;

        TSize i, sum = 0, n = length(a);
        arrayFill(begin(c, Standard()), begin(c, Standard()) + K, 0);    // reset counters

        for (i = 0;  i < n;  i++)                                        // count occurrences
            c[ordValue(getValue(r,getValue(a,i)-1))]++;

        for (i = 0;  i < K;  i++)
        {                                                                // exclusive prefix sums
            TSize t = getValue(c,i);
            c[i] = sum;
            sum += t;
        }

        for (i = 0;  i < n;  i++)
        {
            TSize j = getValue(a,i);                                    // sort
            b[c[ordValue(getValue(r,j-1))]++] = j - 1;
        }
    }

    // stably sort a[0..n-1] to b[0..n-1] with keys in 0..K-1 from r
    template <
        typename TSortedArray,
        typename TUnsortedArray,
        typename TCountArray,
        typename TText >
    void radixExtendClip(
        TSortedArray &b,            // result array (sorted by 1 character)
        TUnsortedArray const &a,    // source array (unsorted)
        TText const &r,                // text to compare with
        TCountArray &c,                // temp. counter array
        unsigned K)                    // alphabet size
    {
        typedef typename Value<TCountArray>::Type TSize;

        TSize i, sum = 0, n = length(a);
        arrayFill(begin(c, Standard()), begin(c, Standard()) + K, 0);    // reset counters

        for (i = 0;  i < n;  i++)
        {                                                                // count occurrences
            TSize j = getValue(a,i);
            if (j > 0) c[ordValue(getValue(r,j-1))]++;
        }

        for (i = 0;  i < K;  i++)
        {                                                                // exclusive prefix sums
            TSize t = getValue(c,i);
            c[i] = sum;
            sum += t;
        }

        for (i = 0;  i < n;  i++)
        {
            TSize j = getValue(a,i);                                    // sort
            if (j > 0) b[c[ordValue(getValue(r,j-1))]++] = j - 1;
        }
    }

}

#endif
