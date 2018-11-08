/***************************************************************************
 *  lib/io/simdisk_file.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002-2003 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2008 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/bits/io/simdisk_file.h>

#if STXXL_HAVE_SIMDISK_FILE

#include <stxxl/bits/io/iostats.h>
#include <stxxl/bits/common/error_handling.h>
#include "ufs_platform.h"
#include <sys/mman.h>

STXXL_BEGIN_NAMESPACE

const double simdisk_geometry::s_average_speed = (15 * 1024 * 1024);

void simdisk_geometry::add_zone(int& first_cyl, int last_cyl,
                                int sec_per_track, int& first_sect)
{
    double rate =
        nsurfaces * sec_per_track * bytes_per_sector /
        ((nsurfaces - 1) * head_switch_time +
         cyl_switch_time +
         nsurfaces * revolution_time);
    int sectors =
        (last_cyl - first_cyl +
         1) * nsurfaces * sec_per_track;
    zones.insert(Zone(first_sect, sectors, rate));
    first_sect += sectors;
    first_cyl = last_cyl + 1;
}

// returns delay in s
double simdisk_geometry::get_delay(file::offset_type offset, file::size_type size)
{
#if 0
    int first_sect = offset / bytes_per_sector;
    int last_sect = (offset + size) / bytes_per_sector;
    int sectors = size / bytes_per_sector;
    double delay =
        cmd_ovh + seek_time + rot_latency +
        double(bytes_per_sector) /
        double(interface_speed);

    std::set<Zone, ZoneCmp>::iterator zone = zones.lower_bound(first_sect);
    //std::cout << __FUNCTION__ << " " << (*zone).first_sector << std::endl;
    while (1)
    {
        int from_this_zone =
            last_sect - ((*zone).first_sector +
                         (*zone).sectors);
        if (from_this_zone <= 0)
        {
            delay += sectors * bytes_per_sector /
                     ((*zone).sustained_data_rate);
            break;
        }
        else
        {
            delay += from_this_zone *
                     bytes_per_sector /
                     ((*zone).sustained_data_rate);
            zone++;
            stxxl_nassert(zone == zones.end());
            sectors -= from_this_zone;
        }
    }

    return delay;
#else
    STXXL_UNUSED(offset);
    return double(size) / s_average_speed;
#endif
}

IC35L080AVVA07::IC35L080AVVA07()
{
    std::cout << "Creating IBM 120GXP IC35L080AVVA07" <<
        std::endl;

    nsurfaces = 4;
    bytes_per_sector = 512;
    cmd_ovh = 0.0002;                           // in s
    seek_time = 0.0082;                         // in s
    rot_latency = 0.00417;                      // in s
    head_switch_time = 0.0015;                  // in s
    cyl_switch_time = 0.002;                    // in s
    revolution_time = 0.0083;                   // in s
    interface_speed = 100000000;                // in byte/s

    int first_sect = 0;
    int last_cyl = 0;
    add_zone(last_cyl, 1938, 928, first_sect);
    add_zone(last_cyl, 3756, 921, first_sect);
    add_zone(last_cyl, 5564, 896, first_sect);
    add_zone(last_cyl, 7687, 896, first_sect);
    add_zone(last_cyl, 9526, 888, first_sect);
    add_zone(last_cyl, 11334, 883, first_sect);
    add_zone(last_cyl, 13331, 864, first_sect);
    add_zone(last_cyl, 15128, 850, first_sect);
    add_zone(last_cyl, 16925, 840, first_sect);
    add_zone(last_cyl, 18922, 822, first_sect);
    add_zone(last_cyl, 20709, 806, first_sect);
    add_zone(last_cyl, 22601, 792, first_sect);
    add_zone(last_cyl, 24138, 787, first_sect);
    add_zone(last_cyl, 26024, 768, first_sect);
    add_zone(last_cyl, 27652, 752, first_sect);
    add_zone(last_cyl, 29501, 740, first_sect);
    add_zone(last_cyl, 31234, 725, first_sect);
    add_zone(last_cyl, 33009, 698, first_sect);
    add_zone(last_cyl, 34784, 691, first_sect);
    add_zone(last_cyl, 36609, 672, first_sect);
    add_zone(last_cyl, 38374, 648, first_sect);
    add_zone(last_cyl, 40139, 630, first_sect);
    add_zone(last_cyl, 41904, 614, first_sect);
    add_zone(last_cyl, 43519, 595, first_sect);
    add_zone(last_cyl, 45250, 576, first_sect);
    add_zone(last_cyl, 47004, 552, first_sect);
    add_zone(last_cyl, 48758, 533, first_sect);
    add_zone(last_cyl, 50491, 512, first_sect);
    add_zone(last_cyl, 52256, 493, first_sect);
    add_zone(last_cyl, 54010, 471, first_sect);
    add_zone(last_cyl, 55571, 448, first_sect);

#if 0
    set<Zone, ZoneCmp>::iterator it = zones.begin();
    int i = 0;
    for ( ; it != zones.end(); it++, i++)
    {
        //const int block_size = 128*3*1024* 4;  // one cylinder

        std::cout << "Zone " << i << " first sector: " << (*it).first_sector;
        std::cout << " sectors: " << (*it).sectors << " sustained rate: ";
        std::cout << (*it).sustained_data_rate / 1024 / 1024 << " MiB/s" << std::endl;
    }

    std::cout << "Last sector     : " << first_sect << std::endl;
    std::cout << "Approx. capacity: " << (first_sect / 1024 / 1024) * bytes_per_sector << " MiB" << std::endl;
#endif

    std::cout << "Transfer 16 MiB from zone 0 : " <<
        get_delay(0, 16 * 1024 * 1024) << " s" << std::endl;
    std::cout << "Transfer 16 MiB from zone 30: " <<
        get_delay(file::offset_type(158204036) * file::offset_type(bytes_per_sector), 16 * 1024 * 1024) << " s" << std::endl;
}

////////////////////////////////////////////////////////////////////////////

void sim_disk_file::serve(void* buffer, offset_type offset, size_type bytes,
                          request::request_type type)
{
    scoped_mutex_lock fd_lock(fd_mutex);

    double op_start = timestamp();

    stats::scoped_read_write_timer read_write_timer(bytes, type == request::WRITE);

    void* mem = mmap(NULL, bytes, PROT_READ | PROT_WRITE, MAP_SHARED, file_des, offset);
    if (mem == MAP_FAILED)
    {
        STXXL_THROW_ERRNO
            (io_error,
            " mmap() failed." <<
            " Page size: " << sysconf(_SC_PAGESIZE) <<
            " offset modulo page size " << (offset % sysconf(_SC_PAGESIZE)));
    }
    else if (mem == 0)
    {
        STXXL_THROW_ERRNO(io_error, "mmap() returned NULL");
    }
    else
    {
        if (type == request::READ)
        {
            memcpy(buffer, mem, bytes);
        }
        else
        {
            memcpy(mem, buffer, bytes);
        }
        STXXL_THROW_ERRNO_NE_0(munmap(mem, bytes), io_error,
                               "munmap() failed");
    }

    double delay = get_delay(offset, bytes);

    delay = delay - timestamp() + op_start;

    assert(delay > 0.0);

    int seconds_to_wait = static_cast<int>(floor(delay));
    if (seconds_to_wait)
        sleep(seconds_to_wait);

    usleep((useconds_t)((delay - seconds_to_wait) * 1000000.));
}

const char* sim_disk_file::io_type() const
{
    return "simdisk";
}

////////////////////////////////////////////////////////////////////////////

void sim_disk_file::set_size(offset_type newsize)
{
    scoped_mutex_lock fd_lock(fd_mutex);
    if (newsize > _size())
    {
        STXXL_THROW_ERRNO_LT_0(::lseek(file_des, newsize - 1, SEEK_SET), io_error,
                               "lseek() fd=" << file_des << " pos=" << newsize - 1);
        STXXL_THROW_ERRNO_LT_0(::write(file_des, "", 1), io_error,
                               "write() fd=" << file_des << " size=1");
    }
}

STXXL_END_NAMESPACE

#endif  // #if STXXL_HAVE_SIMDISK_FILE
// vim: et:ts=4:sw=4
