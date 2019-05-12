/***************************************************************************
 *  include/stxxl/bits/mng/block_alloc_interleaved.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002, 2003 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2007-2009, 2011 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_MNG_BLOCK_ALLOC_INTERLEAVED_HEADER
#define STXXL_MNG_BLOCK_ALLOC_INTERLEAVED_HEADER

#include <vector>

#include <stxxl/bits/mng/block_manager.h>
#include <stxxl/bits/common/rand.h>

STXXL_BEGIN_NAMESPACE

#define CHECK_RUN_BOUNDS(pos)

struct interleaved_striping
{
protected:
    int_type nruns;
    unsigned_type begindisk;
    unsigned_type diff;

    interleaved_striping(int_type nruns, unsigned_type begindisk, unsigned_type diff)
        : nruns(nruns), begindisk(begindisk), diff(diff)
    { }

public:
    interleaved_striping(int_type _nruns, const striping& strategy)
        : nruns(_nruns), begindisk(strategy.begin), diff(strategy.diff)
    { }

    unsigned_type operator () (unsigned_type i) const
    {
        return begindisk + (i / nruns) % diff;
    }
};

struct interleaved_FR : public interleaved_striping
{
    typedef random_number<random_uniform_fast> rnd_type;
    rnd_type rnd;

    interleaved_FR(int_type _nruns, const FR& strategy)
        : interleaved_striping(_nruns, strategy.begin, strategy.diff)
    { }

    unsigned_type operator () (unsigned_type /*i*/) const
    {
        return begindisk + rnd(rnd_type::value_type(diff));
    }
};

struct interleaved_SR : public interleaved_striping
{
    typedef random_number<random_uniform_fast> rnd_type;
    std::vector<int> offsets;

    interleaved_SR(int_type _nruns, const SR& strategy)
        : interleaved_striping(_nruns, strategy.begin, strategy.diff)
    {
        rnd_type rnd;
        for (int_type i = 0; i < nruns; i++)
            offsets.push_back(rnd(rnd_type::value_type(diff)));
    }

    unsigned_type operator () (unsigned_type i) const
    {
        return begindisk + (i / nruns + offsets[i % nruns]) % diff;
    }
};

struct interleaved_RC : public interleaved_striping
{
    std::vector<std::vector<unsigned_type> > perms;

    interleaved_RC(int_type _nruns, const RC& strategy)
        : interleaved_striping(_nruns, strategy.begin, strategy.diff),
          perms(nruns, std::vector<unsigned_type>(diff))
    {
        for (int_type i = 0; i < nruns; i++)
        {
            for (unsigned_type j = 0; j < diff; j++)
                perms[i][j] = j;

            random_number<random_uniform_fast> rnd;
            std::random_shuffle(perms[i].begin(), perms[i].end(), rnd _STXXL_FORCE_SEQUENTIAL);
        }
    }

    unsigned_type operator () (unsigned_type i) const
    {
        return begindisk + perms[i % nruns][(i / nruns) % diff];
    }
};

struct first_disk_only : public interleaved_striping
{
    first_disk_only(int_type _nruns, const single_disk& strategy)
        : interleaved_striping(_nruns, strategy.disk, 1)
    { }

    unsigned_type operator () (unsigned_type) const
    {
        return begindisk;
    }
};

template <typename scheme>
struct interleaved_alloc_traits
{ };

template <>
struct interleaved_alloc_traits<striping>
{
    typedef interleaved_striping strategy;
};

template <>
struct interleaved_alloc_traits<FR>
{
    typedef interleaved_FR strategy;
};

template <>
struct interleaved_alloc_traits<SR>
{
    typedef interleaved_SR strategy;
};

template <>
struct interleaved_alloc_traits<RC>
{
    typedef interleaved_RC strategy;
};

template <>
struct interleaved_alloc_traits<RC_disk>
{
    // FIXME! HACK!
    typedef interleaved_RC strategy;
};

template <>
struct interleaved_alloc_traits<RC_flash>
{
    // FIXME! HACK!
    typedef interleaved_RC strategy;
};

template <>
struct interleaved_alloc_traits<single_disk>
{
    typedef first_disk_only strategy;
};

STXXL_END_NAMESPACE

#endif // !STXXL_MNG_BLOCK_ALLOC_INTERLEAVED_HEADER
// vim: et:ts=4:sw=4
