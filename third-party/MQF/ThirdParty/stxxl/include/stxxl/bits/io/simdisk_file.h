/***************************************************************************
 *  include/stxxl/bits/io/simdisk_file.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002-2003 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2008 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2009 Johannes Singler <singler@ira.uka.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_IO_SIMDISK_FILE_HEADER
#define STXXL_IO_SIMDISK_FILE_HEADER

#include <stxxl/bits/config.h>

#ifndef STXXL_HAVE_SIMDISK_FILE
// use mmap call
#define STXXL_HAVE_SIMDISK_FILE STXXL_HAVE_MMAP_FILE
#endif

#if STXXL_HAVE_SIMDISK_FILE

#include <set>
#include <cmath>

#include <stxxl/bits/io/ufs_file_base.h>
#include <stxxl/bits/io/disk_queued_file.h>

STXXL_BEGIN_NAMESPACE

//! \weakgroup fileimpl
//! \{

class simdisk_geometry : private noncopyable
{
    struct Zone
    {
        // manufactured data
#if 0
        int last_cyl;
        int sect_per_track;
#endif
        // derived data
        int first_sector;
        int sectors;
        double sustained_data_rate;  // in MiB/s
        inline Zone(int _first_sector) : first_sector(_first_sector)
        { }                          // constructor for zone search

        inline Zone(
#if 0
            int _last_cyl,
            int _sect_per_track,
#endif
            int _first_sector,
            int _sectors,
            double _rate)
            :
#if 0
              last_cyl(_last_cyl),
              sect_per_track(_sect_per_track),
#endif
              first_sector(_first_sector),
              sectors(_sectors),
              sustained_data_rate(_rate)
        { }
    };
    struct ZoneCmp
    {
        inline bool operator () (const Zone& a, const Zone& b) const
        {
            return a.first_sector < b.first_sector;
        }
    };

protected:
    int nsurfaces;
    int bytes_per_sector;
    double cmd_ovh;                     // in s
    double seek_time;                   // in s
    double rot_latency;                 // in s
    double head_switch_time;            // in s
    double cyl_switch_time;             // in s
    double revolution_time;             // in s
    double interface_speed;             // in byte/s
    std::set<Zone, ZoneCmp> zones;

    void add_zone(int& first_cyl, int last_cyl,
                  int sec_per_track, int& first_sect);

public:
    inline simdisk_geometry()
    { }
    double get_delay(file::offset_type offset, file::size_type size);                // returns delay in s

    inline ~simdisk_geometry()
    { }

    static const double s_average_speed;
};

class IC35L080AVVA07 : public simdisk_geometry              // IBM series 120GXP
{
public:
    IC35L080AVVA07();
};

//! Implementation of disk emulation.
//! \remark It is emulation of IBM IC35L080AVVA07 disk's timings
class sim_disk_file : public ufs_file_base, public disk_queued_file, public IC35L080AVVA07
{
public:
    //! Constructs file object.
    //! \param filename path of file
    //! \attention filename must be resided at memory disk partition
    //! \param mode open mode, see \c stxxl::file::open_modes
    //! \param queue_id disk queue identifier
    //! \param allocator_id linked disk_allocator
    //! \param device_id physical device identifier
    inline sim_disk_file(
        const std::string& filename,
        int mode,
        int queue_id = DEFAULT_QUEUE,
        int allocator_id = NO_ALLOCATOR,
        unsigned int device_id = DEFAULT_DEVICE_ID)
        : file(device_id),
          ufs_file_base(filename, mode),
          disk_queued_file(queue_id, allocator_id)
    {
        std::cout << "Please, make sure that '" << filename <<
            "' is resided on swap memory partition!" <<
            std::endl;
    }
    void serve(void* buffer, offset_type offset, size_type bytes,
               request::request_type type);
    void set_size(offset_type newsize);
    const char * io_type() const;
};

//! \}

STXXL_END_NAMESPACE

#endif // #if STXXL_HAVE_SIMDISK_FILE

#endif // !STXXL_IO_SIMDISK_FILE_HEADER
