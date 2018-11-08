/***************************************************************************
 *  include/stxxl/bits/common/utils.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002-2006 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2007-2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2008 Johannes Singler <singler@ira.uka.de>
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_COMMON_UTILS_HEADER
#define STXXL_COMMON_UTILS_HEADER

#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <sstream>
#include <limits>

#include <stxxl/bits/config.h>
#include <stxxl/bits/namespace.h>
#include <stxxl/bits/common/types.h>
#include <stxxl/bits/compat/type_traits.h>
#include <stxxl/bits/msvc_compatibility.h>

STXXL_BEGIN_NAMESPACE

////////////////////////////////////////////////////////////////////////////

#if defined(__GNUC__) && ((__GNUC__ > 4) || ((__GNUC__ == 4) && (__GNUC_MINOR__ >= 7)))
#  define STXXL_ATTRIBUTE_UNUSED __attribute__ ((unused))
#else
#  define STXXL_ATTRIBUTE_UNUSED
#endif

////////////////////////////////////////////////////////////////////////////

#if defined(__GXX_EXPERIMENTAL_CXX0X__)
#define STXXL_STATIC_ASSERT(x) static_assert(x, #x)
#else
#define STXXL_STATIC_ASSERT(x) { typedef int static_assert_dummy_type[(x) ? 1 : -1] STXXL_ATTRIBUTE_UNUSED; }
#endif

////////////////////////////////////////////////////////////////////////////

//! Split a string by given separator string. Returns a vector of strings with
//! at least min_fields and at most limit_fields
static inline std::vector<std::string>
split(const std::string& str, const std::string& sep,
      unsigned int min_fields = 0,
      unsigned int limit_fields = std::numeric_limits<unsigned int>::max())
{
    std::vector<std::string> result;
    if (str.empty()) {
        result.resize(min_fields);
        return result;
    }

    std::string::size_type CurPos(0), LastPos(0);
    while (1)
    {
        if (result.size() + 1 == limit_fields)
            break;

        CurPos = str.find(sep, LastPos);
        if (CurPos == std::string::npos)
            break;

        result.push_back(
            str.substr(LastPos,
                       std::string::size_type(CurPos - LastPos))
            );

        LastPos = CurPos + sep.size();
    }

    std::string sub = str.substr(LastPos);
    result.push_back(sub);

    if (result.size() < min_fields)
        result.resize(min_fields);

    return result;
}

////////////////////////////////////////////////////////////////////////////

//! Format any ostream-able type into a string
template <typename Type>
std::string to_str(const Type& t)
{
    std::ostringstream oss;
    oss << t;
    return oss.str();
}

////////////////////////////////////////////////////////////////////////////

//! Parse a string like "343KB" or "44 GiB" into the corresponding size in
//! bytes. Returns the number of bytes and sets ok = true if the string could
//! be parsed correctly. If no units indicator is given, use def_unit in
//! k/m/g/t/p (powers of ten) or in K/M/G/T/P (power of two).
bool parse_SI_IEC_size(const std::string& str, uint64& size, char def_unit = 0);

//! Format a byte size using SI (K, M, G, T) suffixes (powers of ten). Returns
//! "123 M" or similar.
std::string format_SI_size(uint64 number);

//! Format a byte size using IEC (Ki, Mi, Gi, Ti) suffixes (powers of
//! two). Returns "123 Ki" or similar.
std::string format_IEC_size(uint64 number);

////////////////////////////////////////////////////////////////////////////

inline stxxl::int64 atoi64(const char* s)
{
#if STXXL_MSVC
    return _atoi64(s);
#else
    return atoll(s);
#endif
}

////////////////////////////////////////////////////////////////////////////

inline stxxl::uint64 atouint64(const char* s)
{
#if STXXL_MSVC
    return _strtoui64(s, NULL, 10);
#else
    return strtoull(s, NULL, 10);
#endif
}

////////////////////////////////////////////////////////////////////////////

template <typename Type>
inline const Type&
STXXL_MIN(const Type& a, const Type& b)
{
    return std::min<Type>(a, b);
}

template <typename Type>
inline const Type&
STXXL_MAX(const Type& a, const Type& b)
{
    return std::max<Type>(a, b);
}

////////////////////////////////////////////////////////////////////////////

//! calculate the log2 floor of an integral type using math.h
template <typename Integral>
inline Integral log2_ceil(Integral i)
{
    return Integral(ceil(log2(i)));
}

//! calculate the log2 ceiling of an integral type using math.h
template <typename Integral>
inline Integral log2_floor(Integral i)
{
    return Integral(log2(i));
}

////////////////////////////////////////////////////////////////////////////

//! calculate the log2 floor of an integer type (by repeated bit shifts)
template <typename IntegerType>
unsigned int ilog2_floor(IntegerType i)
{
    unsigned int p = 0;
    while (i >= 256) i >>= 8, p += 8;
    while (i >>= 1) ++p;
    return p;
}

//! calculate the log2 ceiling of an integer type (by repeated bit shifts)
template <typename IntegerType>
unsigned int ilog2_ceil(const IntegerType& i)
{
    if (i <= 1) return 0;
    return ilog2_floor(i - 1) + 1;
}

////////////////////////////////////////////////////////////////////////////

template <typename Integral, typename Integral2>
inline
typename compat::remove_const<Integral>::type
div_ceil(Integral n, Integral2 d)
{
#if 0  // ambiguous overload for std::div(unsigned_anything, unsigned_anything)
    typedef __typeof__ (std::div (n, d)) div_type;
    div_type result = std::div(n, d);
    return result.quot + (result.rem != 0);
#else
    return n / d + ((n % d) != 0);
#endif
}

////////////////////////////////////////////////////////////////////////////

#ifdef __GNUC__
#define HAVE_BUILTIN_EXPECT
#endif

#ifdef HAVE_BUILTIN_EXPECT
 #define LIKELY(c)   __builtin_expect((c), 1)
#else
 #define LIKELY(c)   c
#endif

#ifdef HAVE_BUILTIN_EXPECT
 #define UNLIKELY(c)   __builtin_expect((c), 0)
#else
 #define UNLIKELY(c)   c
#endif

////////////////////////////////////////////////////////////////////////////

inline size_t longhash1(uint64 key_)
{
    key_ += ~(key_ << 32);
    key_ ^= (key_ >> 22);
    key_ += ~(key_ << 13);
    key_ ^= (key_ >> 8);
    key_ += (key_ << 3);
    key_ ^= (key_ >> 15);
    key_ += ~(key_ << 27);
    key_ ^= (key_ >> 31);
    return (size_t)key_;
}

////////////////////////////////////////////////////////////////////////////

template <class Type>
inline void swap_1D_arrays(Type* a, Type* b, unsigned_type size)
{
    for (unsigned_type i = 0; i < size; ++i)
        std::swap(a[i], b[i]);
}

////////////////////////////////////////////////////////////////////////////

template <typename Integral>
static inline Integral round_up_to_power_of_two(Integral n)
{
    --n;
    for (int k = 1; k != 8 * sizeof(n); k <<= 1)
        n |= n >> k;
    ++n;
    return n;
}

//! round n up to next larger multiple of 2^power. example: (48,4) = 64, (48,3) = 48.
template <typename Integral>
inline Integral round_up_to_power_of_two(Integral n, unsigned_type power)
{
    Integral pot = Integral(1) << power, // = 0..0 1 0^power
        mask = pot - 1;                  // = 0..0 0 1^power
    if (n & mask)                        // n not divisible by pot
        return (n & ~mask) + pot;
    else
        return n;
}

////////////////////////////////////////////////////////////////////////////

template <class Container>
inline typename Container::value_type pop(Container& c)
{
    typename Container::value_type r = c.top();
    c.pop();
    return r;
}

template <class Container>
inline typename Container::value_type pop_front(Container& c)
{
    typename Container::value_type r = c.front();
    c.pop_front();
    return r;
}

template <class Container>
inline typename Container::value_type pop_back(Container& c)
{
    typename Container::value_type r = c.back();
    c.pop_back();
    return r;
}

template <class Container>
inline typename Container::value_type pop_begin(Container& c)
{
    typename Container::value_type r = *c.begin();
    c.erase(c.begin());
    return r;
}

////////////////////////////////////////////////////////////////////////////

STXXL_END_NAMESPACE

#endif // !STXXL_COMMON_UTILS_HEADER
// vim: et:ts=4:sw=4
