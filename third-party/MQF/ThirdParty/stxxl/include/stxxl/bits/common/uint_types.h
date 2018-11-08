/***************************************************************************
 *  include/stxxl/bits/common/uint_types.h
 *
 *  Class representing a 40-bit or 48-bit unsigned integer encoded in five or
 *  six bytes.
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_COMMON_UINT_TYPES_HEADER
#define STXXL_COMMON_UINT_TYPES_HEADER

#include <stxxl/bits/config.h>
#include <stxxl/bits/namespace.h>
#include <stxxl/bits/common/utils.h>
#include <stxxl/bits/common/types.h>

#include <limits>
#include <ostream>
#include <cassert>

STXXL_BEGIN_NAMESPACE

/*!
 * Construct an 40-bit or 48-bit unsigned integer stored in five or six bytes.
 *
 * The purpose of this class is to provide integers with smaller data storage
 * footprints when more than 32-bit, but less than 64-bit indexes are
 * needed. This is commonly the case for storing file offsets and indexes. Here
 * smaller types currently suffice for files < 1 TiB or < 16 TiB.
 *
 * The class combines a 32-bit integer with a HighType (either 8-bit or 16-bit)
 * to get a larger type. Only unsigned values are supported, which fits the
 * general application of file offsets.
 *
 * Calculation in uint_pair are generally done by transforming everything to
 * 64-bit data type, so that 64-bit register arithmetic can be used. The
 * exception here is \b increment and \b decrement, which is done directly on
 * the lower/higher part. Not all arithmetic operations are supported, patches
 * welcome if you really need the operations.
 */
#if STXXL_MSVC
#pragma pack(push, 1)
#endif
template <typename HighType>
class uint_pair
{
public:
    //! lower part type, always 32-bit
    typedef uint32 low_type;
    //! higher part type, currently either 8-bit or 16-bit
    typedef HighType high_type;

private:
    //! member containing lower significant integer value
    low_type low;
    //! member containing higher significant integer value
    high_type high;

    //! return highest value storable in lower part, also used as a mask.
    static unsigned_type low_max()
    {
        return std::numeric_limits<low_type>::max();
    }

    //! number of bits in the lower integer part, used a bit shift value.
    static const size_t low_bits = 8 * sizeof(low_type);

    //! return highest value storable in higher part, also used as a mask.
    static unsigned_type high_max()
    {
        return std::numeric_limits<high_type>::max();
    }

    //! number of bits in the higher integer part, used a bit shift value.
    static const size_t high_bits = 8 * sizeof(high_type);

public:
    //! number of binary digits (bits) in uint_pair
    static const size_t digits = low_bits + high_bits;

    //! number of bytes in uint_pair
    static const size_t bytes = sizeof(low_type) + sizeof(high_type);

    //! empty constructor, does not even initialize to zero!
    inline uint_pair()
    {
        // compile-time assertions about size of low_type
        STXXL_STATIC_ASSERT(8 * sizeof(low_type) == 32);
        // compile-time assertions about size of our data structure, this tests
        // packing of structures by the compiler
        STXXL_STATIC_ASSERT(sizeof(uint_pair) == bytes);
        STXXL_STATIC_ASSERT(sizeof(uint_pair) == digits / 8);
        STXXL_STATIC_ASSERT(digits / 8 == bytes);
    }

    //! construct unit pair from lower and higher parts.
    inline uint_pair(const low_type& l, const high_type& h)
        : low(l), high(h)
    { }

    //! copy constructor
    inline uint_pair(const uint_pair& a)
        : low(a.low), high(a.high)
    { }

    //! const from a simple 32-bit unsigned integer
    inline uint_pair(const uint32& a)
        : low(a), high(0)
    { }

    //! const from a simple 32-bit signed integer
    inline uint_pair(const int32& a)
        : low(a), high(0)
    {
        if (a >= 0)
            low = a;
        else
            low = a, high = (high_type)high_max();
    }

    //! construct from an uint64 (unsigned long long)
    inline uint_pair(const uint64& a)
        : low((low_type)(a & low_max())),
          high((high_type)((a >> low_bits) & high_max()))
    {
        // check for overflow
        assert((a >> (low_bits + high_bits)) == 0);
    }

    //! return the number as an uint64 (unsigned long long)
    inline uint64 ull() const
    {
        return ((uint64)high) << low_bits | (uint64)low;
    }

    //! implicit cast to an unsigned long long
    inline operator uint64 () const
    {
        return ull();
    }

    //! return the number as a uint64
    inline uint64 u64() const
    {
        return ((uint64)high) << low_bits | (uint64)low;
    }

    //! prefix increment operator (directly manipulates the integer parts)
    inline uint_pair& operator ++ ()
    {
        if (UNLIKELY(low == low_max()))
            ++high, low = 0;
        else
            ++low;
        return *this;
    }

    //! prefix decrement operator (directly manipulates the integer parts)
    inline uint_pair& operator -- ()
    {
        if (UNLIKELY(low == 0))
            --high, low = (low_type)low_max();
        else
            --low;
        return *this;
    }

    //! addition operator (uses 64-bit arithmetic)
    inline uint_pair& operator += (const uint_pair& b)
    {
        uint64 add = (uint64)low + b.low;
        low = (low_type)(add & low_max());
        high = (high_type)(high + b.high + ((add >> low_bits) & high_max()));
        return *this;
    }

    //! equality checking operator
    inline bool operator == (const uint_pair& b) const
    {
        return (low == b.low) && (high == b.high);
    }

    //! inequality checking operator
    inline bool operator != (const uint_pair& b) const
    {
        return (low != b.low) || (high != b.high);
    }

    //! less-than comparison operator
    inline bool operator < (const uint_pair& b) const
    {
        return (high < b.high) || (high == b.high && low < b.low);
    }

    //! less-or-equal comparison operator
    inline bool operator <= (const uint_pair& b) const
    {
        return (high < b.high) || (high == b.high && low <= b.low);
    }

    //! greater comparison operator
    inline bool operator > (const uint_pair& b) const
    {
        return (high > b.high) || (high == b.high && low > b.low);
    }

    //! greater-or-equal comparison operator
    inline bool operator >= (const uint_pair& b) const
    {
        return (high > b.high) || (high == b.high && low >= b.low);
    }

    //! make a uint_pair outputtable via iostreams, using unsigned long long.
    friend std::ostream& operator << (std::ostream& os, const uint_pair& a)
    {
        return os << a.ull();
    }

    //! return an uint_pair instance containing the smallest value possible
    static uint_pair min()
    {
        return uint_pair(std::numeric_limits<low_type>::min(),
                         std::numeric_limits<high_type>::min());
    }

    //! return an uint_pair instance containing the largest value possible
    static uint_pair max()
    {
        return uint_pair(std::numeric_limits<low_type>::max(),
                         std::numeric_limits<high_type>::max());
    }
}
#if STXXL_MSVC
;
#pragma pack(pop)
#else
__attribute__ ((packed));
#endif

//! \addtogroup support
//! \{

//! Construct a 40-bit unsigned integer stored in five bytes.
typedef uint_pair<uint8> uint40;

//! Construct a 48-bit unsigned integer stored in six bytes.
typedef uint_pair<uint16> uint48;

//! \}

STXXL_END_NAMESPACE

namespace std {

//! template class providing some numeric_limits fields for uint_pair types.
template <typename HighType>
class numeric_limits<stxxl::uint_pair<HighType> >
{
public:
    //! yes we have information about uint_pair
    static const bool is_specialized = true;

    //! return an uint_pair instance containing the smallest value possible
    static stxxl::uint_pair<HighType> min()
    { return stxxl::uint_pair<HighType>::min(); }

    //! return an uint_pair instance containing the largest value possible
    static stxxl::uint_pair<HighType> max()
    { return stxxl::uint_pair<HighType>::max(); }

    //! return an uint_pair instance containing the smallest value possible
    static stxxl::uint_pair<HighType> lowest()
    { return min(); }

    //! unit_pair types are unsigned
    static const bool is_signed = false;

    //! uint_pair types are integers
    static const bool is_integer = true;

    //! unit_pair types contain exact integers
    static const bool is_exact = true;

    //! unit_pair radix is binary
    static const int radix = 2;

    //! number of binary digits (bits) in uint_pair
    static const int digits = stxxl::uint_pair<HighType>::digits;

    //! epsilon is zero
    static const stxxl::uint_pair<HighType> epsilon()
    { return stxxl::uint_pair<HighType>(0, 0); }

    //! rounding error is zero
    static const stxxl::uint_pair<HighType> round_error()
    { return stxxl::uint_pair<HighType>(0, 0); }

    //! no exponent
    static const int min_exponent = 0;

    //! no exponent
    static const int min_exponent10 = 0;

    //! no exponent
    static const int max_exponent = 0;

    //! no exponent
    static const int max_exponent10 = 0;

    //! no infinity
    static const bool has_infinity = false;
};

} // namespace std

#endif // !STXXL_COMMON_UINT_TYPES_HEADER
