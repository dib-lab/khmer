/***************************************************************************
 *  include/stxxl/bits/common/rand.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002, 2003, 2005 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2007 Andreas Beckmann <beckmann@mpi-inf.mpg.de>
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_COMMON_RAND_HEADER
#define STXXL_COMMON_RAND_HEADER

#include <cmath>

#include <stxxl/bits/config.h>
#include <stxxl/bits/common/types.h>
#include <stxxl/bits/common/seed.h>
#include <stxxl/bits/namespace.h>

#if STXXL_STD_RANDOM
 #include <random>
#elif STXXL_BOOST_RANDOM
 #include <boost/random.hpp>
#endif

// Recommended seeding procedure:
// by default, the global seed is initialized from a high resolution timer and the process id
// 1. stxxl::set_seed(seed); // optionally, do this if you wan't to us a specific seed to replay a certain program run
// 2. seed = stxxl::get_next_seed(); // store/print/... this value can be used for step 1 to replay the program with a specific seed
// 3. stxxl::srandom_number32(); // seed the global state of stxxl::random_number32
// 4. create all the other prngs used.

STXXL_BEGIN_NAMESPACE

extern unsigned ran32State;

//! \addtogroup support
//! \{

//! Fast uniform [0, 2^32) pseudo-random generator with period 2^32, random
//! bits: 32.
//! \warning Uses a global state and is not reentrant or thread-safe!
struct random_number32
{
    typedef unsigned value_type;

    //! Returns a random number from [0, 2^32)
    inline value_type operator () () const
    {
        return (ran32State = 1664525 * ran32State + 1013904223);
    }

    //! Returns a random number from [0, N)
    inline value_type operator () (const value_type& N) const
    {
        return operator () () % N;
    }
};

//! Set a seed value for \c random_number32.
inline void srandom_number32(unsigned seed = get_next_seed())
{
    ran32State = seed;
}

//! Fast uniform [0, 2^32) pseudo-random generator with period 2^32, random
//! bits: 32.
//! Reentrant variant of random_number32 that keeps it's private state.
struct random_number32_r
{
    typedef unsigned value_type;
    mutable unsigned state;

    random_number32_r(unsigned seed = get_next_seed())
    {
        state = seed;
    }

    //! Change the current seed
    void set_seed(unsigned seed)
    {
        state = seed;
    }

    //! Returns a random number from [0, 2^32)
    inline value_type operator () () const
    {
        return (state = 1664525 * state + 1013904223);
    }
};

//! Fast uniform [0, 255] pseudo-random generator with period 2^8, random bits:
//! 8 (one byte).
class random_number8_r
{
    random_number32_r m_rnd32;
    uint32 m_value;
    unsigned int m_pos;

public:
    typedef uint8 value_type;

    random_number8_r(unsigned seed = get_next_seed())
        : m_rnd32(seed), m_pos(4)
    { }

    //! Returns a random byte from [0, 255]
    inline value_type operator () ()
    {
        if (++m_pos >= 4) {
            m_value = m_rnd32();
            m_pos = 0;
        }
        return ((uint8*)&m_value)[m_pos];
    }
};

//! Fast uniform [0.0, 1.0) pseudo-random generator
//! \warning Uses a global state and is not reentrant or thread-safe!
struct random_uniform_fast
{
    typedef double value_type;
    random_number32 rnd32;

    random_uniform_fast(unsigned /*seed*/ = get_next_seed())
    { }

    //! Returns a random number from [0.0, 1.0)
    inline value_type operator () () const
    {
        return (double(rnd32()) * (0.5 / 0x80000000));
    }
};

#if STXXL_MSVC
#pragma warning(push)
#pragma warning(disable:4512) // assignment operator could not be generated
#endif

//! Slow and precise uniform [0.0, 1.0) pseudo-random generator
//! period: at least 2^48, random bits: at least 31
//!
//! \warning Seed is not the same as in the fast generator \c random_uniform_fast
struct random_uniform_slow
{
    typedef double value_type;
#if STXXL_STD_RANDOM
    typedef std::default_random_engine gen_type;
    mutable gen_type gen;
    typedef std::uniform_real_distribution<> uni_type;
    mutable uni_type uni;

    random_uniform_slow(unsigned seed = get_next_seed())
        : gen(seed), uni(0.0, 1.0)
    { }
#elif STXXL_BOOST_RANDOM
    typedef boost::minstd_rand base_generator_type;
    base_generator_type generator;
    boost::uniform_real<> uni_dist;
    mutable boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni;

    random_uniform_slow(unsigned seed = get_next_seed())
        : uni(generator, uni_dist)
    {
        uni.engine().seed(seed);
    }
#else
    mutable unsigned short state48[3];
/*
 * embedded erand48.c
 *
 * Copyright (c) 1993 Martin Birgmeier
 * All rights reserved.
 *
 * You may redistribute unmodified or modified versions of this source
 * code provided that the above copyright notice and this and the
 * following conditions are retained.
 *
 * This software is provided ``as is'', and comes with no warranties
 * of any kind. I shall in no event be liable for anything that happens
 * to anyone/anything when using this software.
 */
    static void
    _dorand48(unsigned short xseed[3])
    {
        unsigned long accu;
        unsigned short temp[2];

        static const unsigned short _mult[3] = { 0xe66d, 0xdeec, 0x0005 };
        static const unsigned short _add = 0x000b;

        accu = (unsigned long)_mult[0] * (unsigned long)xseed[0]
               + (unsigned long)_add;
        temp[0] = (unsigned short)accu;         /* lower 16 bits */
        accu >>= sizeof(unsigned short) * 8;
        accu += (unsigned long)_mult[0] * (unsigned long)xseed[1]
                + (unsigned long)_mult[1] * (unsigned long)xseed[0];
        temp[1] = (unsigned short)accu;         /* middle 16 bits */
        accu >>= sizeof(unsigned short) * 8;
        accu += _mult[0] * xseed[2] + _mult[1] * xseed[1] + _mult[2] * xseed[0];
        xseed[0] = temp[0];
        xseed[1] = temp[1];
        xseed[2] = (unsigned short)accu;
    }

    static double
    _erand48(unsigned short xseed[3])
    {
        _dorand48(xseed);
        return ldexp((double)xseed[0], -48)
               + ldexp((double)xseed[1], -32)
               + ldexp((double)xseed[2], -16);
    }
/* end erand48.c */

    random_uniform_slow(unsigned seed = get_next_seed())
    {
        state48[0] = (unsigned short)(seed & 0xffff);
        state48[1] = (unsigned short)(seed >> 16);
        state48[2] = 42;
        _dorand48(state48);
    }
#endif

    //! Returns a random number from [0.0, 1.0)
    inline value_type operator () () const
    {
#if STXXL_STD_RANDOM
        return uni(gen);
#elif STXXL_BOOST_RANDOM
        return uni();
#else
        return _erand48(state48);
#endif
    }
};

//! Uniform [0, N) pseudo-random generator
template <class UniformRGen = random_uniform_fast>
struct random_number
{
    typedef unsigned value_type;
    UniformRGen uniform;

    random_number(unsigned seed = get_next_seed())
        : uniform(seed)
    { }

    //! Returns a random number from [0, N)
    inline value_type operator () (value_type N) const
    {
        return static_cast<value_type>(uniform() * double(N));
    }
};

//! Slow and precise uniform [0, 2^64) pseudo-random generator
struct random_number64
{
    typedef stxxl::uint64 value_type;
    random_uniform_slow uniform;

    random_number64(unsigned seed = get_next_seed())
        : uniform(seed)
    { }

    //! Returns a random number from [0, 2^64)
    inline value_type operator () () const
    {
        return static_cast<value_type>(uniform() * (18446744073709551616.));
    }

    //! Returns a random number from [0, N)
    inline value_type operator () (value_type N) const
    {
        return static_cast<value_type>(uniform() * double(N));
    }
};

#if STXXL_MSVC
#pragma warning(pop) // assignment operator could not be generated
#endif

//! \}

STXXL_END_NAMESPACE

#endif // !STXXL_COMMON_RAND_HEADER
