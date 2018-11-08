/***************************************************************************
 *  include/stxxl/bits/stream/choose.h
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

#ifndef STXXL_STREAM_CHOOSE_HEADER
#define STXXL_STREAM_CHOOSE_HEADER

#include <stxxl/bits/namespace.h>

STXXL_BEGIN_NAMESPACE

//! Stream package subnamespace.
namespace stream {

////////////////////////////////////////////////////////////////////////
//     CHOOSE                                                         //
////////////////////////////////////////////////////////////////////////

template <class Input, int Which>
class choose
{ };

//! Creates stream from a tuple stream taking the first component of each tuple.
//!
//! \tparam Input type of the input tuple stream
//!
//! \remark Tuple stream is a stream which \c value_type is \c stxxl::tuple .
template <class Input>
class choose<Input, 1>
{
    Input& in;

    typedef typename Input::value_type tuple_type;

public:
    //! Standard stream typedef.
    typedef typename tuple_type::first_type value_type;

    //! Construction.
    choose(Input& in_) : in(in_)
    { }

    //! Standard stream method.
    const value_type& operator * () const
    {
        return (*in).first;
    }

    const value_type* operator -> () const
    {
        return &(*in).first;
    }

    //! Standard stream method.
    choose& operator ++ ()
    {
        ++in;
        return *this;
    }

    //! Standard stream method.
    bool empty() const
    {
        return in.empty();
    }
};

//! Creates stream from a tuple stream taking the second component of each tuple.
//!
//! \tparam Input type of the input tuple stream
//!
//! \remark Tuple stream is a stream which \c value_type is \c stxxl::tuple .
template <class Input>
class choose<Input, 2>
{
    Input& in;

    typedef typename Input::value_type tuple_type;

public:
    //! Standard stream typedef.
    typedef typename tuple_type::second_type value_type;

    //! Construction.
    choose(Input& in_) : in(in_)
    { }

    //! Standard stream method.
    const value_type& operator * () const
    {
        return (*in).second;
    }

    const value_type* operator -> () const
    {
        return &(*in).second;
    }

    //! Standard stream method.
    choose& operator ++ ()
    {
        ++in;
        return *this;
    }

    //! Standard stream method.
    bool empty() const
    {
        return in.empty();
    }
};

//! Creates stream from a tuple stream taking the third component of each tuple.
//!
//! \tparam Input type of the input tuple stream
//!
//! \remark Tuple stream is a stream which \c value_type is \c stxxl::tuple .
template <class Input>
class choose<Input, 3>
{
    Input& in;

    typedef typename Input::value_type tuple_type;

public:
    //! Standard stream typedef.
    typedef typename tuple_type::third_type value_type;

    //! Construction.
    choose(Input& in_) : in(in_)
    { }

    //! Standard stream method.
    const value_type& operator * () const
    {
        return (*in).third;
    }

    const value_type* operator -> () const
    {
        return &(*in).third;
    }

    //! Standard stream method.
    choose& operator ++ ()
    {
        ++in;
        return *this;
    }

    //! Standard stream method.
    bool empty() const
    {
        return in.empty();
    }
};

//! Creates stream from a tuple stream taking the fourth component of each tuple.
//!
//! \tparam Input type of the input tuple stream
//!
//! \remark Tuple stream is a stream which \c value_type is \c stxxl::tuple .
template <class Input>
class choose<Input, 4>
{
    Input& in;

    typedef typename Input::value_type tuple_type;

public:
    //! Standard stream typedef.
    typedef typename tuple_type::fourth_type value_type;

    //! Construction.
    choose(Input& in_) : in(in_)
    { }

    //! Standard stream method.
    const value_type& operator * () const
    {
        return (*in).fourth;
    }

    const value_type* operator -> () const
    {
        return &(*in).fourth;
    }

    //! Standard stream method.
    choose& operator ++ ()
    {
        ++in;
        return *this;
    }

    //! Standard stream method.
    bool empty() const
    {
        return in.empty();
    }
};

//! Creates stream from a tuple stream taking the fifth component of each tuple.
//!
//! \tparam Input type of the input tuple stream
//!
//! \remark Tuple stream is a stream which \c value_type is \c stxxl::tuple .
template <class Input>
class choose<Input, 5>
{
    Input& in;

    typedef typename Input::value_type tuple_type;

public:
    //! Standard stream typedef.
    typedef typename tuple_type::fifth_type value_type;

    //! Construction.
    choose(Input& in_) : in(in_)
    { }

    //! Standard stream method.
    const value_type& operator * () const
    {
        return (*in).fifth;
    }

    const value_type* operator -> () const
    {
        return &(*in).fifth;
    }

    //! Standard stream method.
    choose& operator ++ ()
    {
        ++in;
        return *this;
    }

    //! Standard stream method.
    bool empty() const
    {
        return in.empty();
    }
};

//! Creates stream from a tuple stream taking the sixth component of each tuple.
//!
//! \tparam Input type of the input tuple stream
//!
//! \remark Tuple stream is a stream which \c value_type is \c stxxl::tuple .
template <class Input>
class choose<Input, 6>
{
    Input& in;

    typedef typename Input::value_type tuple_type;

public:
    //! Standard stream typedef.
    typedef typename tuple_type::sixth_type value_type;

    //! Construction.
    choose(Input& in_) : in(in_)
    { }

    //! Standard stream method.
    const value_type& operator * () const
    {
        return (*in).sixth;
    }

    const value_type* operator -> () const
    {
        return &(*in).sixth;
    }

    //! Standard stream method.
    choose& operator ++ ()
    {
        ++in;
        return *this;
    }

    //! Standard stream method.
    bool empty() const
    {
        return in.empty();
    }
};

//! \}

} // namespace stream

STXXL_END_NAMESPACE

#include <stxxl/bits/stream/unique.h>

#endif // !STXXL_STREAM_CHOOSE_HEADER
// vim: et:ts=4:sw=4
