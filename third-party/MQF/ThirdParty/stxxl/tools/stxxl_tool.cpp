/***************************************************************************
 *  tools/stxxl_tool.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2007, 2009-2011 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/io>
#include <stxxl/mng>
#include <stxxl/version.h>
#include <stxxl/bits/common/utils.h>
#include <stxxl/bits/common/cmdline.h>
#include <stxxl/bits/parallel.h>

int stxxl_info(int, char**)
{
    stxxl::config::get_instance();
    stxxl::block_manager::get_instance();
    stxxl::stats::get_instance();
    stxxl::disk_queues::get_instance();

#if STXXL_PARALLEL
    STXXL_MSG("STXXL_PARALLEL, max threads = " << omp_get_max_threads());
#endif
    STXXL_MSG("sizeof(unsigned int)   = " << sizeof(unsigned int));
    STXXL_MSG("sizeof(unsigned_type)  = " << sizeof(stxxl::unsigned_type));
    STXXL_MSG("sizeof(uint64)         = " << sizeof(stxxl::uint64));
    STXXL_MSG("sizeof(long)           = " << sizeof(long));
    STXXL_MSG("sizeof(size_t)         = " << sizeof(size_t));
    STXXL_MSG("sizeof(off_t)          = " << sizeof(off_t));
    STXXL_MSG("sizeof(void*)          = " << sizeof(void*));

#if defined(STXXL_HAVE_LINUXAIO_FILE)
    STXXL_MSG("STXXL_HAVE_LINUXAIO_FILE = " << STXXL_HAVE_LINUXAIO_FILE);
#endif

    return 0;
}

extern int create_files(int argc, char* argv[]);
extern int benchmark_disks(int argc, char* argv[]);
extern int benchmark_files(int argc, char* argv[]);
extern int benchmark_sort(int argc, char* argv[]);
extern int benchmark_disks_random(int argc, char* argv[]);
extern int benchmark_pqueue(int argc, char* argv[]);
extern int do_mlock(int argc, char* argv[]);
extern int do_mallinfo(int argc, char* argv[]);

struct SubTool
{
    const char* name;
    int (* func)(int argc, char* argv[]);
    bool shortline;
    const char* description;
};

struct SubTool subtools[] = {
    {
        "info", &stxxl_info, false,
        "Print out information about the build system and which optional "
        "modules where compiled into STXXL."
    },
    {
        "create_files", &create_files, false,
        "Precreate large files to keep file system allocation time out to measurements."
    },
    {
        "benchmark_disks", &benchmark_disks, false,
        "This program will benchmark the disks configured by the standard "
        ".stxxl disk configuration files mechanism."
    },
    {
        "benchmark_files", &benchmark_files, false,
        "Benchmark different file access methods, e.g. syscall or mmap_files."
    },
    {
        "benchmark_sort", &benchmark_sort, false,
        "Run benchmark tests of different sorting methods in STXXL"
    },
    {
        "benchmark_disks_random", &benchmark_disks_random, false,
        "Benchmark random block access time to .stxxl configured disks."
    },
    {
        "benchmark_pqueue", &benchmark_pqueue, false,
        "Benchmark priority queue implementation using sequence of operations."
    },
    {
        "mlock", &do_mlock, true,
        "Lock physical memory."
    },
    {
        "mallinfo", &do_mallinfo, true,
        "Show mallinfo statistics."
    },
    { NULL, NULL, false, NULL }
};

int main_usage(const char* arg0)
{
    STXXL_MSG(stxxl::get_version_string_long());

    std::cout << "Usage: " << arg0 << " <subtool> ..." << std::endl
              << "Available subtools: " << std::endl;

    int shortlen = 0;

    for (unsigned int i = 0; subtools[i].name; ++i)
    {
        if (!subtools[i].shortline) continue;
        shortlen = std::max(shortlen, (int)strlen(subtools[i].name));
    }

    for (unsigned int i = 0; subtools[i].name; ++i)
    {
        if (subtools[i].shortline) continue;
        std::cout << "  " << subtools[i].name << std::endl;
        stxxl::cmdline_parser::output_wrap(std::cout, subtools[i].description, 80, 6, 6);
        std::cout << std::endl;
    }

    for (unsigned int i = 0; subtools[i].name; ++i)
    {
        if (!subtools[i].shortline) continue;
        std::cout << "  " << std::left << std::setw(shortlen + 2)
                  << subtools[i].name << subtools[i].description << std::endl;
    }
    std::cout << std::endl;

    return 0;
}

int main(int argc, char** argv)
{
    char progsub[256];

    if (stxxl::check_library_version() != 0)
        STXXL_ERRMSG("version mismatch between headers and library");

    if (argc > 1)
    {
        for (unsigned int i = 0; subtools[i].name; ++i)
        {
            if (strcmp(subtools[i].name, argv[1]) == 0)
            {
                // replace argv[1] with call string of subtool.
                snprintf(progsub, sizeof(progsub), "%s %s", argv[0], argv[1]);
                argv[1] = progsub;
                return subtools[i].func(argc - 1, argv + 1);
            }
        }
        std::cout << "Unknown subtool '" << argv[1] << "'" << std::endl;
    }

    return main_usage(argv[0]);
}
