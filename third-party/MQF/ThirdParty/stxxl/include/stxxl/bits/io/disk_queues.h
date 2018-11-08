/***************************************************************************
 *  include/stxxl/bits/io/disk_queues.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2008-2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2009 Johannes Singler <singler@ira.uka.de>
 *  Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_IO_DISK_QUEUES_HEADER
#define STXXL_IO_DISK_QUEUES_HEADER

#include <map>

#include <stxxl/bits/namespace.h>
#include <stxxl/bits/singleton.h>
#include <stxxl/bits/io/iostats.h>
#include <stxxl/bits/io/request.h>
#include <stxxl/bits/io/request_queue_impl_qwqr.h>
#include <stxxl/bits/io/linuxaio_queue.h>
#include <stxxl/bits/io/linuxaio_request.h>
#include <stxxl/bits/io/serving_request.h>

STXXL_BEGIN_NAMESPACE

//! \addtogroup reqlayer
//! \{

//! Encapsulates disk queues.
//! \remark is a singleton
class disk_queues : public singleton<disk_queues>
{
    friend class singleton<disk_queues>;

    typedef stxxl::int64 DISKID;
    typedef std::map<DISKID, request_queue*> request_queue_map;

protected:
    request_queue_map queues;
    disk_queues()
    {
        stxxl::stats::get_instance(); // initialize stats before ourselves
    }

public:
    void add_request(request_ptr& req, DISKID disk)
    {
#ifdef STXXL_HACK_SINGLE_IO_THREAD
        disk = 42;
#endif
        request_queue_map::iterator qi = queues.find(disk);
        request_queue* q;
        if (qi == queues.end())
        {
            // create new request queue
#if STXXL_HAVE_LINUXAIO_FILE
            if (dynamic_cast<linuxaio_request*>(req.get()))
                q = queues[disk] = new linuxaio_queue(
                        dynamic_cast<linuxaio_file*>(req->get_file())->get_desired_queue_length()
                        );
            else
#endif
            q = queues[disk] = new request_queue_impl_qwqr();
        }
        else
            q = qi->second;

        q->add_request(req);
    }

    //! Cancel a request.
    //! The specified request is canceled unless already being processed.
    //! However, cancelation cannot be guaranteed.
    //! Cancelled requests must still be waited for in order to ensure correct
    //! operation.
    //! \param req request to cancel
    //! \param disk disk number for disk that \c req was scheduled on
    //! \return \c true iff the request was canceled successfully
    bool cancel_request(request_ptr& req, DISKID disk)
    {
#ifdef STXXL_HACK_SINGLE_IO_THREAD
        disk = 42;
#endif
        if (queues.find(disk) != queues.end())
            return queues[disk]->cancel_request(req);
        else
            return false;
    }

    request_queue * get_queue(DISKID disk)
    {
        if (queues.find(disk) != queues.end())
            return queues[disk];
        else
            return NULL;
    }

    ~disk_queues()
    {
        // deallocate all queues
        for (request_queue_map::iterator i = queues.begin(); i != queues.end(); i++)
            delete (*i).second;
    }

    //! Changes requests priorities.
    //! \param op one of:
    //!                 - READ, read requests are served before write requests within a disk queue
    //!                 - WRITE, write requests are served before read requests within a disk queue
    //!                 - NONE, read and write requests are served by turns, alternately
    void set_priority_op(request_queue::priority_op op)
    {
        for (request_queue_map::iterator i = queues.begin(); i != queues.end(); i++)
            i->second->set_priority_op(op);
    }
};

//! \}

STXXL_END_NAMESPACE

#endif // !STXXL_IO_DISK_QUEUES_HEADER
// vim: et:ts=4:sw=4
