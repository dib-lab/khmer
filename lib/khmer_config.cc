//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt. 
// Contact: khmer-project@idyll.org
//

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
Config( )
:
#ifdef KHMER_EXTRA_SANITY_CHECKS
    _has_extra_sanity_checks( true ),
#else
    _has_extra_sanity_checks( false ),
#endif
    _number_of_threads( 1 ),
    _reads_input_buffer_size( 512U*1024*1024 ),
    _ibmgr_trace_level( 255U ),
    _rparser_trace_level( 255U )
{ }


bool const
Config::
has_extra_sanity_checks( void ) const
{ return _has_extra_sanity_checks; }


uint32_t const
Config::
get_number_of_threads( void ) const
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
get_reads_input_buffer_size( void ) const
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


uint8_t const
Config::
get_input_buffer_trace_level( void ) const
{ return _ibmgr_trace_level; }


void
Config::
set_input_buffer_trace_level( uint8_t trace_level )
{ _ibmgr_trace_level = trace_level; }


uint8_t const
Config::
get_reads_parser_trace_level( void ) const
{ return _rparser_trace_level; }


void
Config::
set_reads_parser_trace_level( uint8_t trace_level )
{ _rparser_trace_level = trace_level; }


} // namespace khmer

// vim: set ft=cpp sts=4 sw=4 tw=79:
