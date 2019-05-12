/***************************************************************************
 *  include/stxxl/bits/common/exceptions.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2006 Roman Dementiev <dementiev@ira.uka.de>
 *  Copyright (C) 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_COMMON_EXCEPTIONS_HEADER
#define STXXL_COMMON_EXCEPTIONS_HEADER

#include <iostream>
#include <string>
#include <stdexcept>

#include <stxxl/bits/namespace.h>

STXXL_BEGIN_NAMESPACE

class io_error : public std::ios_base::failure
{
public:
    io_error() throw ()
        : std::ios_base::failure("")
    { }

    io_error(const std::string& message) throw ()
        : std::ios_base::failure(message)
    { }
};

class resource_error : public std::runtime_error
{
public:
    resource_error() throw ()
        : std::runtime_error("")
    { }

    resource_error(const std::string& message) throw ()
        : std::runtime_error(message)
    { }
};

class bad_ext_alloc : public std::runtime_error
{
public:
    bad_ext_alloc() throw ()
        : std::runtime_error("")
    { }

    bad_ext_alloc(const std::string& message) throw ()
        : std::runtime_error(message)
    { }
};

class bad_parameter : public std::runtime_error
{
public:
    bad_parameter() throw ()
        : std::runtime_error("")
    { }

    bad_parameter(const std::string& message) throw ()
        : std::runtime_error(message)
    { }
};

class unreachable : public std::runtime_error
{
public:
    unreachable() throw ()
        : std::runtime_error("")
    { }

    unreachable(const std::string& message) throw ()
        : std::runtime_error(message)
    { }
};

STXXL_END_NAMESPACE

#endif // !STXXL_COMMON_EXCEPTIONS_HEADER
// vim: et:ts=4:sw=4
