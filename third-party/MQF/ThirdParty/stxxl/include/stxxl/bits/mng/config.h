/***************************************************************************
 *  include/stxxl/bits/mng/config.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002-2005 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2008, 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_MNG_CONFIG_HEADER
#define STXXL_MNG_CONFIG_HEADER

#include <vector>
#include <string>
#include <cstdlib>
#include <cassert>

#include <stxxl/version.h>
#include <stxxl/bits/singleton.h>
#include <stxxl/bits/common/log.h>

STXXL_BEGIN_NAMESPACE

//! \addtogroup mnglayer
//! \{

//! Encapsulate the configuration of one "disk". The disk is actually a file
//! I/O object which block_manager uses to read/write blocks.
class disk_config
{
public:
    //! \name Basic Disk Configuration Parameters
    //! \{

    //! the file path used by the io implementation
    std::string path;

    //! file size to initially allocate
    uint64 size;

    //! io implementation to access file
    std::string io_impl;

    //! \}

public:
    //! default constructor
    disk_config();

    //! initializing constructor, also parses fileio parameter
    disk_config(const std::string& path, uint64 size, const std::string& fileio);

    //! initializing constructor, parse full line as in config files
    disk_config(const std::string& line);

    //! parse a disk=\<path>,\<size>,\<fileio> options line into disk_config,
    //! throws std::runtime_error on parse errors.
    void parse_line(const std::string& line);

    //! parse the "io_impl" parameter into the optional parameter fields.
    void parse_fileio();

    //! return formatted fileio name and optional configuration parameters
    std::string fileio_string() const;

public:
    //! \name Optional Disk / File I/O Implementation Parameters
    //! \{

    //! autogrow file if more disk space is needed, automatically set if size == 0.
    bool autogrow;

    //! delete file on program exit (default for autoconfigurated files)
    bool delete_on_exit;

    //! tristate variable: direct=0 -> force direct OFF, direct=1 -> try direct
    //! ON, if fails print warning and open without direct, direct=2 -> force
    //! direct ON, fail if unavailable.
    enum direct_type { DIRECT_OFF = 0, DIRECT_TRY = 1, DIRECT_ON = 2 } direct;

    //! marks flash drives (configuration entries with flash= instead of disk=)
    bool flash;

    //! select request queue for disk. Use different queues for files on
    //! different disks. queue=-1 -> default queue (one for each disk).
    int queue;

    //! the selected physical device id (e.g. for calculating prefetching
    //! sequences). If -1 then the device id is chosen automatically.
    unsigned int device_id;

    //! turned on by syscall fileio when the path points to a raw block device
    bool raw_device;

    //! unlink file immediately after opening (available on most Unix)
    bool unlink_on_open;

    //! desired queue length for linuxaio_file and linuxaio_queue
    int queue_length;

    //! \}
};

//! Access point to disks properties. Since 1.4.0: no config files are read
//! automatically!
//! \remarks is a singleton
class config : public singleton<config>
{
    friend class singleton<config>;

    //! typedef of list of configured disks
    typedef std::vector<disk_config> disk_list_type;

    //! list of configured disks
    disk_list_type disks_list;

    //! In disks_list, flash devices come after all regular disks
    unsigned first_flash;

    //! Finished initializing config
    bool is_initialized;

    //! Constructor: this must be inlined to print the header version
    //! string.
    inline config()
        : is_initialized(false)
    {
        logger::get_instance();
        STXXL_MSG(get_version_string_long());
        print_library_version_mismatch();
    }

    //! deletes autogrow files
    ~config();

    //! Search several places for a config file.
    void find_config();

    //! If disk list is empty, then search different locations for a disk
    //! configuration file, or load a default config if everything fails.
    void initialize();

public:
    //! \name Initialization Functions
    //! \{

    //! Check that initialize() was called.
    //! \note This function need not be called by the user, block_manager will
    //! always call it.
    void check_initialized()
    {
        if (!is_initialized) initialize();
    }

    //! Load disk configuration file.
    void load_config_file(const std::string& config_path);

    //! Load default configuration.
    void load_default_config();

    //! Add a disk to the configuration list.
    //!
    //! \warning This function should only be used during initialization, as it
    //! has no effect after construction of block_manager.
    inline config & add_disk(const disk_config& cfg)
    {
        disks_list.push_back(cfg);
        return *this;
    }

    //! \}

protected:
    //! \name Automatic Disk Enumeration Functions
    //! \{

    //! static counter for automatic physical device enumeration
    unsigned int m_max_device_id;

public:
    //! Returns automatic physical device id counter
    unsigned int get_max_device_id();

    //! Returns next automatic physical device id counter
    unsigned int get_next_device_id();

    //! Update the automatic physical device id counter
    void update_max_device_id(unsigned int devid);

    //! \}

public:
    //! \name Query Functions
    //! \{

    //! Returns number of disks available to user.
    //! \return number of disks
    inline size_t disks_number()
    {
        check_initialized();
        return disks_list.size();
    }

    //! Returns contiguous range of regular disks w/o flash devices in the array of all disks.
    //! \return range [begin, end) of regular disk indices
    inline std::pair<unsigned, unsigned> regular_disk_range() const
    {
        assert(is_initialized);
        return std::pair<unsigned, unsigned>(0, first_flash);
    }

    //! Returns contiguous range of flash devices in the array of all disks.
    //! \return range [begin, end) of flash device indices
    inline std::pair<unsigned, unsigned> flash_range() const
    {
        assert(is_initialized);
        return std::pair<unsigned, unsigned>(first_flash, (unsigned)disks_list.size());
    }

    //! Returns mutable disk_config structure for additional disk parameters
    inline disk_config & disk(size_t disk)
    {
        check_initialized();
        return disks_list[disk];
    }

    //! Returns constant disk_config structure for additional disk parameters
    inline const disk_config & disk(size_t disk) const
    {
        assert(is_initialized);
        return disks_list[disk];
    }

    //! Returns path of disks.
    //! \param disk disk's identifier
    //! \return string that contains the disk's path name
    inline const std::string & disk_path(size_t disk) const
    {
        assert(is_initialized);
        return disks_list[disk].path;
    }

    //! Returns disk size.
    //! \param disk disk's identifier
    //! \return disk size in bytes
    inline stxxl::uint64 disk_size(size_t disk) const
    {
        assert(is_initialized);
        return disks_list[disk].size;
    }

    //! Returns name of I/O implementation of particular disk.
    //! \param disk disk's identifier
    inline const std::string & disk_io_impl(size_t disk) const
    {
        assert(is_initialized);
        return disks_list[disk].io_impl;
    }

    //! Returns the total size over all disks
    uint64 total_size() const;

    //! \}
};

//! \}

STXXL_END_NAMESPACE

#endif // !STXXL_MNG_CONFIG_HEADER
// vim: et:ts=4:sw=4
