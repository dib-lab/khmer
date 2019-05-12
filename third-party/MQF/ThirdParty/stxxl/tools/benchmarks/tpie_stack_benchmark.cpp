/***************************************************************************
 *  tools/benchmarks/tpie_stack_benchmark.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2006 Roman Dementiev <dementiev@ira.uka.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

//! \example containers/tpie_stack_benchmark.cpp
//! This is a benchmark mentioned in the paper
//! R. Dementiev, L. Kettner, P. Sanders "STXXL: standard template library for XXL data sets"
//! Software: Practice and Experience
//! Volume 38, Issue 6, Pages 589-637, May 2008
//! DOI: 10.1002/spe.844

#include "app_config.h"

#include <portability.h>
#include <versions.h>

// Get the AMI_stack definition.
#include <ami_stack.h>

// Utilities for ASCII output.
#include <ami_scan_utils.h>

#include <stxxl/bits/common/utils.h>
#include <stxxl/bits/verbose.h>
#include <stxxl/timer>

#define MEM_2_RESERVE    (768 * 1024 * 1024)

#define BLOCK_SIZE       (2 * 1024 * 1024)

#ifndef DISKS
 #define DISKS 1
#endif

template <unsigned RECORD_SIZE>
struct my_record_
{
    char data[RECORD_SIZE];
    my_record_() { }
};

template <class my_record>
void run_stack(stxxl::int64 volume)
{
    typedef AMI_stack<my_record> stack_type;

    MM_manager.set_memory_limit(BLOCK_SIZE * DISKS * 8);

    STXXL_MSG("Record size: " << sizeof(my_record) << " bytes");

    stack_type Stack;

    stxxl::int64 ops = volume / sizeof(my_record);

    stxxl::int64 i;

    my_record cur;

    stxxl::timer Timer;
    Timer.start();

    for (i = 0; i < ops; ++i)
    {
        Stack.push(cur);
    }

    Timer.stop();

    STXXL_MSG("Records in Stack: " << Stack.stream_len());
    if (i != Stack.stream_len())
    {
        STXXL_MSG("Size does not match");
        abort();
    }

    STXXL_MSG("Insertions elapsed time: " << (Timer.mseconds() / 1000.) <<
              " seconds : " << (double(volume) / (1024. * 1024. * Timer.mseconds() / 1000.)) <<
              " MiB/s");

    ////////////////////////////////////////////////
    Timer.reset();
    Timer.start();

    my_record* out;

    for (i = 0; i < ops; ++i)
    {
        Stack.pop(&out);
    }

    Timer.stop();

    STXXL_MSG("Records in Stack: " << Stack.stream_len());
    if (Stack.stream_len() != 0)
    {
        STXXL_MSG("Stack must be empty");
        abort();
    }

    STXXL_MSG("Deletions elapsed time: " << (Timer.mseconds() / 1000.) <<
              " seconds : " << (double(volume) / (1024. * 1024. * Timer.mseconds() / 1000.)) <<
              " MiB/s");
}

int main(int argc, char* argv[])
{
    using std::cout;
    using std::endl;
#ifdef BTE_COLLECTION_IMP_MMAP
    cout << "BTE_COLLECTION_IMP_MMAP is defined" << endl;
#endif
#ifdef BTE_COLLECTION_IMP_UFS
    cout << "BTE_COLLECTION_IMP_UFS is defined" << endl;
#endif
#ifdef BTE_STREAM_IMP_UFS
    cout << "BTE_STREAM_IMP_UFS is defined" << endl;
    cout << "BTE_STREAM_UFS_BLOCK_FACTOR is " << BTE_STREAM_UFS_BLOCK_FACTOR << endl;
    cout << "Actual block size is " << (TPIE_OS_BLOCKSIZE() * BTE_STREAM_UFS_BLOCK_FACTOR / 1024) << " KiB" << endl;
#endif
#ifdef BTE_STREAM_IMP_MMAP
    cout << "BTE_STREAM_IMP_MMAP is defined" << endl;
    cout << "BTE_STREAM_MMAP_BLOCK_FACTOR is " << BTE_STREAM_MMAP_BLOCK_FACTOR << endl;
    cout << "Actual block size is " << (TPIE_OS_BLOCKSIZE() * BTE_STREAM_MMAP_BLOCK_FACTOR / 1024) << " KiB" << endl;
#endif
#ifdef BTE_STREAM_IMP_STDIO
    cout << "BTE_STREAM_IMP_STDIO is defined" << endl;
#endif
    cout << "TPIE_OS_BLOCKSIZE() is " << TPIE_OS_BLOCKSIZE() << endl;

    if (argc < 3)
    {
        STXXL_MSG("Usage: " << argv[0] << " version #volume");
        STXXL_MSG("\t version = 1: TPIE stack with 4 byte records");
        STXXL_MSG("\t version = 2: TPIE stack with 32 byte records");
        return -1;
    }

    int version = atoi(argv[1]);
    stxxl::uint64 volume = stxxl::atouint64(argv[2]);

    STXXL_MSG("Allocating array with size " << MEM_2_RESERVE
                                            << " bytes to prevent file buffering.");
    //int * array = new int[MEM_2_RESERVE/sizeof(int)];
    int* array = (int*)malloc(MEM_2_RESERVE);
    std::fill(array, array + (MEM_2_RESERVE / sizeof(int)), 0);

    STXXL_MSG("Running version: " << version);
    STXXL_MSG("Data volume    : " << volume << " bytes");

    switch (version)
    {
    case 1:
        run_stack<my_record_<4> >(volume);
        break;
    case 2:
        run_stack<my_record_<32> >(volume);
        break;
    default:
        STXXL_MSG("Unsupported version " << version);
    }

    //delete [] array;
    free(array);
}
