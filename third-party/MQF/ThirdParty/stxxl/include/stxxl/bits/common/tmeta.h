/***************************************************************************
 *  include/stxxl/bits/common/tmeta.h
 *
 *  Template Metaprogramming Tools
 *  (from the Generative Programming book Krysztof Czarnecki, Ulrich Eisenecker)
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2003 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2008 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_COMMON_TMETA_HEADER
#define STXXL_COMMON_TMETA_HEADER

#include <stxxl/bits/namespace.h>
#include <stxxl/bits/common/types.h>

STXXL_BEGIN_NAMESPACE

//! IF template metaprogramming statement.
//!
//! If \c Flag is \c true then \c IF<>::result is of type Type1
//! otherwise of \c IF<>::result is of type Type2
template <bool Flag, class Type1, class Type2>
struct IF
{
    typedef Type1 result;
};

template <class Type1, class Type2>
struct IF<false, Type1, Type2>
{
    typedef Type2 result;
};

//! If \c Flag is \c true then \c IF<>::result is Num1
//! otherwise of \c IF<>::result is Num2
template <bool Flag, unsigned Num1, unsigned Num2>
struct IF_N
{
    enum
    {
        result = Num1
    };
};

template <unsigned Num1, unsigned Num2>
struct IF_N<false, Num1, Num2>
{
    enum
    {
        result = Num2
    };
};

const int DEFAULT = ~(~0u >> 1); // initialize with the smallest int

struct NilCase { };

template <int tag_, class Type_, class Next_ = NilCase>
struct CASE
{
    enum { tag = tag_ };
    typedef Type_ Type;
    typedef Next_ Next;
};

template <int tag, class Case>
class SWITCH
{
    typedef typename Case::Next NextCase;
    enum
    {
        caseTag = Case::tag,
        found = (caseTag == tag || caseTag == DEFAULT)
    };

public:
    typedef typename IF<found,
                        typename Case::Type,
                        typename SWITCH<tag, NextCase>::result
                        >::result result;
};

template <int tag>
class SWITCH<tag, NilCase>
{
public:
    typedef NilCase result;
};

//! \internal, use LOG2 instead
template <unsigned_type Input>
class LOG2_floor
{
public:
    enum
    {
        value = LOG2_floor<Input / 2>::value + 1
    };
};

template <>
class LOG2_floor<1>
{
public:
    enum
    {
        value = 0
    };
};

template <>
class LOG2_floor<0>
{
public:
    enum
    {
        value = 0
    };
};

template <unsigned_type Input>
class LOG2
{
public:
    enum
    {
        floor = LOG2_floor<Input>::value,
        ceil = LOG2_floor<Input - 1>::value + 1
    };
};

template <>
class LOG2<1>
{
public:
    enum
    {
        floor = 0,
        ceil = 0
    };
};

template <>
class LOG2<0>
{
public:
    enum
    {
        floor = 0,
        ceil = 0
    };
};

STXXL_END_NAMESPACE

#endif // !STXXL_COMMON_TMETA_HEADER
