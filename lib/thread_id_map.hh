//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt.
// Contact: khmer-project@idyll.org
//

#ifndef THREAD_ID_MAP_HH
#   define THREAD_ID_MAP_HH

#include <exception>
#include <map>

// TODO? Just use 'pthread_t' everywhere.
// Linux uses standard PIDs for thread IDs.
#if defined (__linux__)
#   include <sys/types.h>
// MacOS X uses Mach kernel ports for thread IDs.
#elif defined (__APPLE__) && defined (__MACH__)
#   include <mach/mach.h>
#   include <pthread.h>
// Else, hope that some POSIX threads implementation is available.
// If so, then try to use 'pthread_t' instances as thread IDs.
#else
#   include <pthread.h>
#endif

#include "khmer.hh"

namespace khmer
{


struct InvalidNumberOfThreadsRequested : public std:: exception {
};

struct TooManyThreads : public std:: exception {
};


struct ThreadIDMap {

    ThreadIDMap( uint32_t number_of_threads );
    ~ThreadIDMap( );

    uint32_t const get_thread_id( );

private:

    uint32_t				_number_of_threads;
    uint32_t				_thread_counter;
#if defined (__linux__)
    std:: map< pid_t, uint32_t >	_thread_id_map;
#elif defined (__APPLE__) && defined (__MACH__)
    std:: map< mach_port_t, uint32_t >	_thread_id_map;
#else
    std:: map< pthread_t, uint32_t >	_thread_id_map;
#endif
    uint32_t				_tid_map_spin_lock;

}; // struct ThreadIDMap


} // namespace khmer

#endif // THREAD_ID_MAP_HH

// vim: set ft=cpp sts=4 sw=4 tw=79:
