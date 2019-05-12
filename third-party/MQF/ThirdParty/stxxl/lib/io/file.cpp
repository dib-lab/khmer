/***************************************************************************
 *  lib/io/file.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/bits/io/file.h>
#include "ufs_platform.h"

STXXL_BEGIN_NAMESPACE

int file::unlink(const char* path)
{
    return ::unlink(path);
}

STXXL_END_NAMESPACE
