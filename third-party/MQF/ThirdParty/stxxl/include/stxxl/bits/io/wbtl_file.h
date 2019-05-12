/***************************************************************************
 *  include/stxxl/bits/io/wbtl_file.h
 *
 *  a write-buffered-translation-layer pseudo file
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2008-2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2009 Johannes Singler <singler@ira.uka.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_IO_WBTL_FILE_HEADER
#define STXXL_IO_WBTL_FILE_HEADER

#ifndef STXXL_HAVE_WBTL_FILE
#define STXXL_HAVE_WBTL_FILE 1
#endif

#if STXXL_HAVE_WBTL_FILE

#include <map>

#include <stxxl/bits/io/disk_queued_file.h>

STXXL_BEGIN_NAMESPACE

//! \addtogroup fileimpl
//! \{

//! Implementation of file based on buffered writes and block remapping via a
//! translation layer.
class wbtl_file : public disk_queued_file
{
    typedef std::pair<offset_type, offset_type> place;
    typedef std::map<offset_type, offset_type> sortseq;
    typedef std::map<offset_type, place> place_map;

    // the physical disk used as backend
    file* storage;
    offset_type sz;
    size_type write_block_size;

    mutex mapping_mutex;
    // logical to physical address translation
    sortseq address_mapping;
    // physical to (logical address, size) translation
    place_map reverse_mapping;
    // list of free (physical) regions
    sortseq free_space;
    offset_type free_bytes;

    // the write buffers:
    // write_buffer[curbuf] is the current write buffer
    // write_buffer[1-curbuf] is the previous write buffer
    // buffer_address if the start offset on the backend file
    // curpos is the next writing position in write_buffer[curbuf]
    mutex buffer_mutex;
    char* write_buffer[2];
    offset_type buffer_address[2];
    int curbuf;
    size_type curpos;
    request_ptr backend_request;

    struct FirstFit : public std::binary_function<place, offset_type, bool>
    {
        bool operator () (
            const place& entry,
            const offset_type size) const
        {
            return (entry.second >= size);
        }
    };

public:
    //! Constructs file object.
    //! param backend_file file object used as storage backend, will be deleted in ~wbtl_file()
    wbtl_file(
        file* backend_file,
        size_type write_buffer_size,
        int write_buffers = 2,
        int queue_id = DEFAULT_QUEUE,
        int allocator_id = NO_ALLOCATOR);
    ~wbtl_file();
    offset_type size();
    void set_size(offset_type newsize);
    void lock();
    void serve(void* buffer, offset_type offset, size_type bytes,
               request::request_type type);
    void discard(offset_type offset, offset_type size);
    const char * io_type() const;

private:
    void _add_free_region(offset_type offset, offset_type size);

protected:
    void sread(void* buffer, offset_type offset, size_type bytes);
    void swrite(void* buffer, offset_type offset, size_type bytes);
    offset_type get_next_write_block();
    void check_corruption(offset_type region_pos, offset_type region_size,
                          sortseq::iterator pred, sortseq::iterator succ);
};

//! \}

STXXL_END_NAMESPACE

#endif // #if STXXL_HAVE_WBTL_FILE

#endif // !STXXL_IO_WBTL_FILE_HEADER
