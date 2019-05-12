/***************************************************************************
 *  include/stxxl/bits/io/io.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2008 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_IO_IO_HEADER
#define STXXL_IO_IO_HEADER

#include <stxxl/request>
#include <stxxl/bits/io/file.h>
#include <stxxl/bits/io/syscall_file.h>
#include <stxxl/bits/io/mmap_file.h>
#include <stxxl/bits/io/simdisk_file.h>
#include <stxxl/bits/io/wincall_file.h>
#include <stxxl/bits/io/boostfd_file.h>
#include <stxxl/bits/io/mem_file.h>
#include <stxxl/bits/io/fileperblock_file.h>
#include <stxxl/bits/io/wbtl_file.h>
#include <stxxl/bits/io/linuxaio_file.h>
#include <stxxl/bits/io/create_file.h>
#include <stxxl/bits/io/disk_queues.h>
#include <stxxl/bits/io/iostats.h>
#include <stxxl/bits/namespace.h>

//! \c STXXL library namespace
STXXL_BEGIN_NAMESPACE

    STXXL_END_NAMESPACE

#endif // !STXXL_IO_IO_HEADER
