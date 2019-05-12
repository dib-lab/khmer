/***************************************************************************
 *  include/stxxl/bits/stream/unique.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2003-2005 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_STREAM_UNIQUE_HEADER
#define STXXL_STREAM_UNIQUE_HEADER

#include <stxxl/bits/namespace.h>

STXXL_BEGIN_NAMESPACE

//! Stream package subnamespace.
namespace stream {

////////////////////////////////////////////////////////////////////////
//     UNIQUE                                                         //
////////////////////////////////////////////////////////////////////////

struct dummy_cmp_unique { };

//! Equivalent to std::unique algorithms.
//!
//! Removes consecutive duplicates from the stream.
//! Uses BinaryPredicate to compare elements of the stream
template <class Input, class BinaryPredicate = dummy_cmp_unique>
class unique
{
    Input& input;
    BinaryPredicate binary_pred;
    typename Input::value_type current;

public:
    //! Standard stream typedef.
    typedef typename Input::value_type value_type;

    unique(Input& input_, BinaryPredicate binary_pred_)
        : input(input_), binary_pred(binary_pred_)
    {
        if (!input.empty())
            current = *input;
    }

    //! Standard stream method.
    unique& operator ++ ()
    {
        value_type old_value = current;
        ++input;
        while (!input.empty() && (binary_pred(current = *input, old_value)))
            ++input;
        return *this;
    }

    //! Standard stream method.
    const value_type& operator * () const
    {
        return current;
    }

    //! Standard stream method.
    const value_type* operator -> () const
    {
        return &current;
    }

    //! Standard stream method.
    bool empty() const
    {
        return input.empty();
    }
};

//! Equivalent to std::unique algorithms.
//!
//! Removes consecutive duplicates from the stream.
template <class Input>
class unique<Input, dummy_cmp_unique>
{
    Input& input;
    typename Input::value_type current;

public:
    //! Standard stream typedef.
    typedef typename Input::value_type value_type;

    unique(Input& input_) : input(input_)
    {
        if (!input.empty())
            current = *input;
    }

    //! Standard stream method.
    unique& operator ++ ()
    {
        value_type old_value = current;
        ++input;
        while (!input.empty() && ((current = *input) == old_value))
            ++input;
        return *this;
    }

    //! Standard stream method.
    const value_type& operator * () const
    {
        return current;
    }

    //! Standard stream method.
    const value_type* operator -> () const
    {
        return &current;
    }

    //! Standard stream method.
    bool empty() const
    {
        return input.empty();
    }
};

//! \}

} // namespace stream

STXXL_END_NAMESPACE

#endif // !STXXL_STREAM_UNIQUE_HEADER
// vim: et:ts=4:sw=4
