/***************************************************************************
 *  include/stxxl/bits/parallel/tags.h
 *
 *  Tags for compile-time options.
 *  Extracted from MCSTL http://algo2.iti.uni-karlsruhe.de/singler/mcstl/
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2007 Johannes Singler <singler@ira.uka.de>
 *  Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_PARALLEL_TAGS_HEADER
#define STXXL_PARALLEL_TAGS_HEADER

#include <stxxl/bits/namespace.h>

STXXL_BEGIN_NAMESPACE

namespace parallel {

/** Makes the parametrized class actually do work, i. e. actives it. */
class active_tag
{ };
/** Makes the parametrized class do nothing, i. e. deactives it. */
class inactive_tag
{ };

} // namespace parallel

STXXL_END_NAMESPACE

#endif // !STXXL_PARALLEL_TAGS_HEADER
