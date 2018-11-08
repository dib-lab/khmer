/***************************************************************************
 *  include/stxxl/bits/common/settings.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2007 Johannes Singler <singler@ira.uka.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_COMMON_SETTINGS_HEADER
#define STXXL_COMMON_SETTINGS_HEADER

/*!
 * @file stxxl/bits/common/settings.h
 * Provides a static class to store runtime tuning parameters.
 */

#include <stxxl/bits/namespace.h>

STXXL_BEGIN_NAMESPACE

template <typename MustBeInt = int>
class settings
{
public:
    static bool native_merge;
};

template <typename MustBeInt>
bool settings<MustBeInt>::native_merge = false;

typedef settings<> SETTINGS;

STXXL_END_NAMESPACE

#endif // !STXXL_COMMON_SETTINGS_HEADER
