/***************************************************************************
 *  include/stxxl/bits/algo/adaptor.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_ALGO_ADAPTOR_HEADER
#define STXXL_ALGO_ADAPTOR_HEADER

#include <stxxl/bits/mng/bid.h>
#include <stxxl/bits/mng/adaptor.h>

STXXL_BEGIN_NAMESPACE

template <unsigned BlockSize, typename RunType, class PosType = int_type>
struct runs2bid_array_adaptor : public two2one_dim_array_adapter_base<RunType*, BID<BlockSize>, PosType>
{
    typedef runs2bid_array_adaptor<BlockSize, RunType, PosType> self_type;
    typedef BID<BlockSize> data_type;

    enum {
        block_size = BlockSize
    };

    unsigned_type dim_size;

    typedef two2one_dim_array_adapter_base<RunType*, BID<BlockSize>, PosType> parent_type;
    using parent_type::array;
    using parent_type::pos;

    runs2bid_array_adaptor(RunType** a, PosType p, unsigned_type d)
        : two2one_dim_array_adapter_base<RunType*, BID<BlockSize>, PosType>(a, p), dim_size(d)
    { }
    runs2bid_array_adaptor(const self_type& a)
        : two2one_dim_array_adapter_base<RunType*, BID<BlockSize>, PosType>(a), dim_size(a.dim_size)
    { }

    const self_type& operator = (const self_type& a)
    {
        array = a.array;
        pos = a.pos;
        dim_size = a.dim_size;
        return *this;
    }

    data_type& operator * ()
    {
        CHECK_RUN_BOUNDS(pos);
        return (BID<BlockSize>&)((*(array[(pos) % dim_size]))[(pos) / dim_size].bid);
    }

    const data_type* operator -> () const
    {
        CHECK_RUN_BOUNDS(pos);
        return &((*(array[(pos) % dim_size])[(pos) / dim_size].bid));
    }

    data_type& operator [] (PosType n) const
    {
        n += pos;
        CHECK_RUN_BOUNDS(n);
        return (BID<BlockSize>&)((*(array[(n) % dim_size]))[(n) / dim_size].bid);
    }
};

BLOCK_ADAPTOR_OPERATORS(runs2bid_array_adaptor)

template <unsigned BlockSize, typename RunType, class PosType = int_type>
struct runs2bid_array_adaptor2
    : public two2one_dim_array_adapter_base<RunType*, BID<BlockSize>, PosType>
{
    typedef runs2bid_array_adaptor2<BlockSize, RunType, PosType> self_type;
    typedef BID<BlockSize> data_type;

    typedef two2one_dim_array_adapter_base<RunType*, BID<BlockSize>, PosType> base_type;

    using base_type::pos;
    using base_type::array;

    enum {
        block_size = BlockSize
    };

    PosType w, h, K;

    runs2bid_array_adaptor2(RunType** a, PosType p, int_type _w, int_type _h)
        : two2one_dim_array_adapter_base<RunType*, BID<BlockSize>, PosType>(a, p),
          w(_w), h(_h), K(_w * _h)
    { }

    runs2bid_array_adaptor2(const self_type& a)
        : two2one_dim_array_adapter_base<RunType*, BID<BlockSize>, PosType>(a),
          w(a.w), h(a.h), K(a.K)
    { }

    const self_type& operator = (const self_type& a)
    {
        array = a.array;
        pos = a.pos;
        w = a.w;
        h = a.h;
        K = a.K;
        return *this;
    }

    data_type& operator * ()
    {
        PosType i = pos - K;
        if (i < 0)
            return (BID<BlockSize>&)((*(array[(pos) % w]))[(pos) / w].bid);

        PosType _w = w;
        _w--;
        return (BID<BlockSize>&)((*(array[(i) % _w]))[h + (i / _w)].bid);
    }

    const data_type* operator -> () const
    {
        PosType i = pos - K;
        if (i < 0)
            return &((*(array[(pos) % w])[(pos) / w].bid));

        PosType _w = w;
        _w--;
        return &((*(array[(i) % _w])[h + (i / _w)].bid));
    }

    data_type& operator [] (PosType n) const
    {
        n += pos;
        PosType i = n - K;
        if (i < 0)
            return (BID<BlockSize>&)((*(array[(n) % w]))[(n) / w].bid);

        PosType _w = w;
        _w--;
        return (BID<BlockSize>&)((*(array[(i) % _w]))[h + (i / _w)].bid);
    }
};

BLOCK_ADAPTOR_OPERATORS(runs2bid_array_adaptor2)

template <typename trigger_iterator_type>
struct trigger_entry_iterator
{
    typedef trigger_entry_iterator<trigger_iterator_type> self_type;
    typedef typename std::iterator_traits<trigger_iterator_type>::value_type::bid_type bid_type;

    // STL typedefs
    typedef bid_type value_type;
    typedef std::random_access_iterator_tag iterator_category;
    typedef int_type difference_type;
    typedef value_type* pointer;
    typedef value_type& reference;

    trigger_iterator_type value;

    trigger_entry_iterator(const self_type& a) : value(a.value) { }
    trigger_entry_iterator(trigger_iterator_type v) : value(v) { }

    bid_type& operator * ()
    {
        return value->bid;
    }
    bid_type* operator -> () const
    {
        return &(value->bid);
    }
    const bid_type& operator [] (int_type n) const
    {
        return (value + n)->bid;
    }
    bid_type& operator [] (int_type n)
    {
        return (value + n)->bid;
    }

    self_type& operator ++ ()
    {
        value++;
        return *this;
    }
    self_type operator ++ (int)
    {
        self_type tmp = *this;
        value++;
        return tmp;
    }
    self_type& operator -- ()
    {
        value--;
        return *this;
    }
    self_type operator -- (int)
    {
        self_type tmp = *this;
        value--;
        return tmp;
    }
    bool operator == (const self_type& a) const
    {
        return value == a.value;
    }
    bool operator != (const self_type& a) const
    {
        return value != a.value;
    }
    self_type operator += (int_type n)
    {
        value += n;
        return *this;
    }
    self_type operator -= (int_type n)
    {
        value -= n;
        return *this;
    }
    int_type operator - (const self_type& a) const
    {
        return value - a.value;
    }
    int_type operator + (const self_type& a) const
    {
        return value + a.value;
    }
};

template <typename Iterator>
inline
trigger_entry_iterator<Iterator>
make_bid_iterator(Iterator iter)
{
    return trigger_entry_iterator<Iterator>(iter);
}

STXXL_END_NAMESPACE

#endif // !STXXL_ALGO_ADAPTOR_HEADER
