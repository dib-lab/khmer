/***************************************************************************
 *  include/stxxl/bits/algo/intksort.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002 Peter Sanders <sanders@mpi-sb.mpg.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_ALGO_INTKSORT_HEADER
#define STXXL_ALGO_INTKSORT_HEADER

#include <algorithm>
#include <cassert>
#include <stxxl/bits/common/types.h>
#include <stxxl/bits/unused.h>
#include <stxxl/bits/parallel.h>

STXXL_BEGIN_NAMESPACE

template <typename TypeKey>
static void
count(TypeKey* a, TypeKey* aEnd, int_type* bucket, int_type K,
      typename TypeKey::key_type offset, unsigned shift)
{
    // reset buckets
    std::fill(bucket, bucket + K, 0);

    // count occupancies
    for (TypeKey* p = a; p < aEnd; p++)
    {
        int_type i = (int_type)((p->key - offset) >> shift);
        /*
        if (!(i < K && i >= 0))
        {
            STXXL_ERRMSG("i: " << i);
            abort();
        }
        */
        bucket[i]++;
    }
}

static inline void
exclusive_prefix_sum(int_type* bucket, int_type K)
{
    int_type sum = 0;
    for (int_type i = 0; i < K; i++)
    {
        int_type current = bucket[i];
        bucket[i] = sum;
        sum += current;
    }
}

// distribute input a to output b using bucket for the starting indices
template <typename TypeKey>
static void
classify(TypeKey* a, TypeKey* aEnd, TypeKey* b, int_type* bucket,
         typename TypeKey::key_type offset, unsigned shift)
{
    for (TypeKey* p = a; p < aEnd; p++)
    {
        int_type i = (int_type)((p->key - offset) >> shift);
        int_type bi = bucket[i];
        b[bi] = *p;
        bucket[i] = bi + 1;
    }
}

template <class Type>
inline void
sort2(Type& a, Type& b)
{
    if (b < a)
        std::swap(a, b);
}

template <class Type>
inline void
sort3(Type& a, Type& b, Type& c)
{
    Type temp;
    if (b < a)
    {
        if (c < a)
        {                       // b , c < a
            if (b < c)
            {                   // b < c < a
                temp = a;
                a = b;
                b = c;
                c = temp;
            }
            else
            {                   // c <=b < a
                std::swap(c, a);
            }
        }
        else
        {                       // b < a <=c
            std::swap(a, b);
        }
    }
    else
    {                           // a <=b
        if (c < a)
        {                       // c < a <=b
            temp = a;
            a = c;
            c = b;
            b = temp;
        }
        else
        {                       // a <=b , c
            if (c < b)
            {                   // a <=c < b
                std::swap(b, c);
            }
        }
    }
    // Assert1 (!(b < a) && !(c < b));
}

template <class Type>
inline void
sort4(Type& a, Type& b, Type& c, Type& d)
{
    sort2(a, b);
    sort2(c, d);                // a < b ; c < d
    if (c < a)
    {                           // c minimal, a < b
        if (d < a)
        {                       // c < d < a < b
            std::swap(a, c);
            std::swap(b, d);
        }
        else
        {                       // c < a < {db}
            if (d < b)
            {                   // c < a < d < b
                Type temp = a;
                a = c;
                c = d;
                d = b;
                b = temp;
            }
            else
            {                   // c < a < b < d
                Type temp = a;
                a = c;
                c = b;
                b = temp;
            }
        }
    }
    else
    {                           // a minimal ; c < d
        if (c < b)
        {                       // c < (bd)
            if (d < b)
            {                   // c < d < b
                Type temp = b;
                b = c;
                c = d;
                d = temp;
            }
            else
            {                   // a < c < b < d
                std::swap(b, c);
            }
        }                       // else sorted
    }
    //Assert1 (!(b < a) && !(c < b) & !(d < c));
}

template <class Type>
inline void
sort5(Type& a, Type& b, Type& c, Type& d, Type& e)
{
    sort2(a, b);
    sort2(d, e);
    if (d < a)
    {
        std::swap(a, d);
        std::swap(b, e);
    }                           // a < d < e, a < b
    if (d < c)
    {
        std::swap(c, d);        // a minimal, c < {de}
        sort2(d, e);
    }
    else
    {                           // a<b, a<d<e, c<d<e
        sort2(a, c);
    }                           // a min, c < d < e
    // insert b int cde by binary search
    if (d < b)
    {                           // c < d < {be}
        if (e < b)
        {                       // c < d < e < b
            Type temp = b;
            b = c;
            c = d;
            d = e;
            e = temp;
        }
        else
        {                       // c < d < b < e
            Type temp = b;
            b = c;
            c = d;
            d = temp;
        }
    }
    else
    {                           // {cb} <=d < e
        sort2(b, c);
    }
    //Assert1 (!(b < a) && !(c < b) & !(d < c) & !(e < d));
}

template <class Type>
inline void
insertion_sort(Type* a, Type* aEnd)
{
    Type* pp;
    for (Type* p = a + 1; p < aEnd; p++)
    {
        // Invariant a..p-1 is sorted;
        Type t = *p;
        if (t < *a)
        {   // new minimum
            // move stuff to the right
            for (pp = p; pp != a; pp--)
            {
                *pp = *(pp - 1);
            }
            *pp = t;
        }
        else
        {
            // now we can use *a as a sentinel
            for (pp = p; t < *(pp - 1); pp--)
            {
                *pp = *(pp - 1);
            }
            *pp = t;
        }
    }
}

// sort each bucket
// bucket[i] is an index one off to the right from
// the end of the i-th bucket
template <class Type>
static void
cleanup(Type* b, int_type* bucket, int_type K)
{
    Type* c = b;
    for (int_type i = 0; i < K; i++)
    {
        Type* cEnd = b + bucket[i];
        switch (cEnd - c)
        {
        case 0:
            break;
        case 1:
            break;
        case 2:
            sort2(c[0], c[1]);
            break;
        case 3:
            sort3(c[0], c[1], c[2]);
            break;
        case 4:
            sort4(c[0], c[1], c[2], c[3]);
            break;
        case 5:
#if 0
            sort5(c[0], c[1], c[2], c[3], c[4]);
            break;
#endif
        case 6:
        case 7:
        case 8:
        case 9:
        case 10:
        case 11:
        case 12:
        case 13:
        case 14:
        case 15:
        case 16:
            insertion_sort(c, cEnd);
            break;
        default:
            check_sort_settings();
            potentially_parallel::
            sort(c, cEnd);
        }
        c = cEnd;
    }
}

// do a single level MDS radix sort
// using bucket[0..K-1] as a counter array
// and using (key(x) - offset) >> shift to index buckets.
// and using (key(x) - offset) >> shift to index buckets.
// the input comes from a..aEnd-1
// the output goes to b
template <typename TypeKey>
void
l1sort(TypeKey* a,
       TypeKey* aEnd,
       TypeKey* b, int_type* bucket, int_type K,
       typename TypeKey::key_type offset, int shift)
{
    count(a, aEnd, bucket, K, offset, shift);
    exclusive_prefix_sum(bucket, K);
    classify(a, aEnd, b, bucket, offset, shift);
    cleanup(b, bucket, K);
}

template <typename Type, typename TypeKey, typename KeyExtractor>
void classify_block(Type* begin, Type* end, TypeKey*& out,
                    int_type* bucket, typename KeyExtractor::key_type offset, unsigned shift, KeyExtractor keyobj)
{
    assert(shift < (sizeof(typename KeyExtractor::key_type) * 8 + 1));
    for (Type* p = begin; p < end; p++, out++)  // count & create references
    {
        out->ptr = p;
        typename KeyExtractor::key_type key = keyobj(*p);
        int_type ibucket = (int_type)((key - offset) >> shift);
        out->key = key;
        bucket[ibucket]++;
    }
}
template <typename Type, typename TypeKey, typename KeyExtractor>
void classify_block(Type* begin, Type* end, TypeKey*& out, int_type* bucket, typename Type::key_type offset, unsigned shift,
                    const int_type K, KeyExtractor keyobj)
{
    assert(shift < (sizeof(typename Type::key_type) * 8 + 1));
    for (Type* p = begin; p < end; p++, out++)  // count & create references
    {
        out->ptr = p;
        typename Type::key_type key = keyobj(*p);
        int_type ibucket = (key - offset) >> shift;
        /*
        if (!(ibucket < K && ibucket >= 0))
        {
            STXXL_ERRMSG("ibucket: " << ibucket << " K:" << K);
            abort();
        }
        */
        out->key = key;
        bucket[ibucket]++;
    }
    STXXL_UNUSED(K);
}

STXXL_END_NAMESPACE

#endif // !STXXL_ALGO_INTKSORT_HEADER
