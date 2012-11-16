#include <cstdlib>
#include <cstdio>

#include "khmer_config.hh"

namespace khmer
{
  

// Create default config at initialization.
static Config	    the_config	      = Config( );


Config &
get_active_config( void )
{
    return the_config;
}


void
set_active_config( Config &config )
{
    the_config = config;
}


Config::
Config( void )
{
#ifdef KHMER_EXTRA_SANITY_CHECKS
    _has_extra_sanity_checks	= true;
#else
    _has_extra_sanity_checks	= false;
#endif
    _number_of_threads		= 1;
    _reads_input_buffer_size	= 512U*1024*1024;
    // TODO? Precalculate hash count thresholds.
}


bool const
Config::
has_extra_sanity_checks( void )
{ return _has_extra_sanity_checks; }


uint32_t const
Config::
get_number_of_threads( void ) 
{ return _number_of_threads; } 


void
Config::
set_number_of_threads( uint32_t const number_of_threads )
{
    if (0 == number_of_threads)
	// TODO: Throw exception,
	//	 or interpret as number of cores on current machine.
	;

    // TODO: Move this check to hash table instantiation.
    // TODO: Scale according to size of counter type.
    if (number_of_threads > 127)
	// TODO: Throw exception.
	;

    // Require at least 64 KiB buffer sizes per thread.
    if (_reads_input_buffer_size < (65536U * (uint64_t)number_of_threads))
	// TODO: Throw exception.
	;

    _number_of_threads = number_of_threads;
}


uint64_t const
Config::
get_reads_input_buffer_size( void )
{ return _reads_input_buffer_size; }


void
Config::
set_reads_input_buffer_size( uint64_t const reads_input_buffer_size )
{
    // Require at least 64 KiB per thread.
    if (reads_input_buffer_size < (65536U * (uint64_t)_number_of_threads))
	// TODO: Throw exception.
	;

    _reads_input_buffer_size = reads_input_buffer_size;
}

} // namespace khmer

// vim: set ft=cpp sts=4 sw=4 tw=79:
