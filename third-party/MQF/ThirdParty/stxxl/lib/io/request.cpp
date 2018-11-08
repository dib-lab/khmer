/***************************************************************************
 *  lib/io/request.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2008 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2009 Johannes Singler <singler@ira.uka.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <ostream>

#include <stxxl/bits/io/request.h>
#include <stxxl/bits/io/file.h>

STXXL_BEGIN_NAMESPACE

request::request(
    const completion_handler& on_compl,
    file* file,
    void* buffer,
    offset_type offset,
    size_type bytes,
    request_type type)
    : m_on_complete(on_compl),
      m_file(file),
      m_buffer(buffer),
      m_offset(offset),
      m_bytes(bytes),
      m_type(type)
{
    STXXL_VERBOSE3_THIS("request::(...), ref_cnt=" << get_reference_count());
    m_file->add_request_ref();
}

request::~request() noexcept(false)
{
    STXXL_VERBOSE3_THIS("request::~request(), ref_cnt=" << get_reference_count());
}

void request::check_alignment() const
{
    if (m_offset % STXXL_BLOCK_ALIGN != 0)
        STXXL_ERRMSG("Offset is not aligned: modulo " <<
                     STXXL_BLOCK_ALIGN << " = " << m_offset % STXXL_BLOCK_ALIGN);

    if (m_bytes % STXXL_BLOCK_ALIGN != 0)
        STXXL_ERRMSG("Size is not a multiple of " <<
                     STXXL_BLOCK_ALIGN << ", = " << m_bytes % STXXL_BLOCK_ALIGN);

    if (unsigned_type(m_buffer) % STXXL_BLOCK_ALIGN != 0)
        STXXL_ERRMSG("Buffer is not aligned: modulo " <<
                     STXXL_BLOCK_ALIGN << " = " << unsigned_type(m_buffer) % STXXL_BLOCK_ALIGN <<
                     " (" << m_buffer << ")");
}

void request::check_nref_failed(bool after)
{
    STXXL_ERRMSG("WARNING: serious error, reference to the request is lost " <<
                 (after ? "after" : "before") << " serve()" <<
                 " nref=" << get_reference_count() <<
                 " this=" << this <<
                 " offset=" << m_offset <<
                 " buffer=" << m_buffer <<
                 " bytes=" << m_bytes <<
                 " type=" << ((m_type == READ) ? "READ" : "WRITE") <<
                 " file=" << m_file <<
                 " iotype=" << m_file->io_type()
                 );
}

const char* request::io_type() const
{
    return m_file->io_type();
}

std::ostream& request::print(std::ostream& out) const
{
    out << "File object address: " << static_cast<void*>(m_file);
    out << " Buffer address: " << static_cast<void*>(m_buffer);
    out << " File offset: " << m_offset;
    out << " Transfer size: " << m_bytes << " bytes";
    out << " Type of transfer: " << ((m_type == READ) ? "READ" : "WRITE");
    return out;
}

STXXL_END_NAMESPACE
// vim: et:ts=4:sw=4
