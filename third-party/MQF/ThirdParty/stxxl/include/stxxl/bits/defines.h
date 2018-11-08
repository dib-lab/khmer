/***************************************************************************
 *  include/stxxl/bits/defines.h
 *
 *  Document all defines that may change the behavior of stxxl.
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2008-2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_DEFINES_HEADER
#define STXXL_DEFINES_HEADER

//#define STXXL_HAVE_MMAP_FILE 0/1
//#define STXXL_HAVE_SIMDISK_FILE 0/1
//#define STXXL_HAVE_BOOSTFD_FILE 0/1
//#define STXXL_HAVE_WINCALL_FILE 0/1
//#define STXXL_HAVE_WBTL_FILE 0/1
//#define STXXL_HAVE_LINUXAIO_FILE 0/1
// default: 0/1 (platform and type dependent)
// used in: io/*_file.h, io/*_file.cpp, mng/mng.cpp
// affects: library
// effect:  enables/disables some file implementations

//#define STXXL_CHECK_BLOCK_ALIGNING
// default: not defined
// used in: io/*_file.cpp
// effect:  call request::check_alignment() from request::request(...)

//#define STXXL_CHECK_FOR_PENDING_REQUESTS_ON_SUBMISSION 0/1
// default: 1
// used in: io/*_queue*.cpp
// affects: library
// effect:  check (and warn) for multiple concurrently pending I/O requests
//          for the same block, usually causing coherency problems on
//          out-of-order execution

//#define STXXL_DO_NOT_COUNT_WAIT_TIME
// default: not defined
// used in: io/iostats.{h,cpp}
// effect:  makes calls to wait time counting functions no-ops

//#define STXXL_WAIT_LOG_ENABLED
// default: not defined
// used in: common/log.cpp, io/iostats.cpp
// effect:  writes wait timing information to the file given via environment
//          variable STXXLWAITLOGFILE, does nothing if this is not defined

//#define STXXL_PRINT_TIMESTAMP_ALWAYS
// default: not defined
// used in: common/verbose.cpp
// affects: library
// effect:  prefix all MSG/ERRMSG/VERBOSE with elapsed time since program start

//#define STXXL_SORT_OPTIMAL_PREFETCHING 0/1
// default: 1
// used in: algo/*sort.h, stream/sort_stream.h
// effect if defined to 0: does not reorder prefetch requests to a disk
//          optimal schedule (Hutchinson, Sanders, Vitter: Duality between
//          prefetching and queued writing on parallel disks, 2005)

//#define STXXL_CHECK_ORDER_IN_SORTS 0/1
// default: 0
// used in: algo/*sort.h, stream/sort_stream.h, containers/priority_queue.h
// effect if set to 1: perform additional checking of sorted results

//#define STXXL_NO_WARN_RECURSIVE_SORT
// default: not defined
// used in: algo/sort_base.h
// affects: programs
// effect if defined: does not print error messages about possibly inefficient
//          recursive merging

//#define STXXL_HACK_SINGLE_IO_THREAD
// default: not defined
// used in: io/disk_queues.h
// affects: programs
// effect if defined: uses only a single I/O thread instead of one per disk
//          used e.g. by EcoSort which puts input file, output file and
//          scratch on a single disk (RAID0)

//#define STXXL_MNG_COUNT_ALLOCATION 0/1
// default: 1
// used in: mng/block_manager.h
// effect if defined: counts current, total and maximum allocation of bytes in
// block manager. The numbers are exported via block_manager's get_
// functions. This can be used to determine the maximum disk space required by
// an application.

//#define STXXL_NO_DEPRECATED 0/1
// default: 0
// used in deprecated.h
// turns off deprecated warnings for some forced template instantiations

#endif // !STXXL_DEFINES_HEADER
