//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt.
// Contact: khmer-project@idyll.org
//

#if defined (__linux__)
#   include <unistd.h>
#   include <sys/syscall.h>
#endif

#include "thread_id_map.hh"

namespace khmer
{


ThreadIDMap::
ThreadIDMap( uint32_t number_of_threads )
    :   _number_of_threads( number_of_threads ),
        _thread_counter( 0 ),
        _tid_map_spin_lock( 0 )
{
    if (0 == number_of_threads) {
        throw InvalidNumberOfThreadsRequested( );
    }
}


ThreadIDMap::
~ThreadIDMap( )
{
    _thread_id_map.clear( );
}

} // namespace khmer

// vim: set ft=cpp sts=4 sw=4 tw=79:
