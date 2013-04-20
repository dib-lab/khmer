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
    if (0 == number_of_threads) throw InvalidNumberOfThreadsRequested( );
}


ThreadIDMap::
~ThreadIDMap( )
{
    _thread_id_map.clear( );
}


uint32_t const
ThreadIDMap::
get_thread_id( )
{
#if defined (__linux__)
    // Note: No error handling because this call always succeeds, allegedly.
    pid_t native_thread_id = syscall( SYS_gettid ); 
    std:: map< pid_t, uint32_t > :: iterator match;
#elif defined (__APPLE__) && defined (__MACH__)
    // TODO? Error handling.
    mach_port_t native_thread_id = pthread_mach_thread_np( pthread_self( ) );
    std:: map< mach_port_t, uint32_t > :: iterator match;
#else
    // TODO? Error handling.
    pthread_t native_thread_id = pthread_self( );
    std:: map< pthread_t, uint32_t > :: iterator match;
#endif

    match = _thread_id_map.find( native_thread_id );
    if (match == _thread_id_map.end( ))
    {
	uint32_t thread_id;

	while (!__sync_bool_compare_and_swap( &_tid_map_spin_lock, 0, 1 ));

	thread_id = _thread_counter++;

	try
	{
	    if (_number_of_threads < _thread_counter)
		throw TooManyThreads( );
	    _thread_id_map[ native_thread_id ] = thread_id;
	}
	catch (...)
	{
	    _thread_counter--;
	    __sync_bool_compare_and_swap( &_tid_map_spin_lock, 1, 0 );
	    throw;
	}

	__sync_bool_compare_and_swap( &_tid_map_spin_lock, 1, 0 );

	return thread_id;
    }

    return (*match).second;
} // get_thread_id


} // namespace khmer

// vim: set ft=cpp sts=4 sw=4 tw=79:
