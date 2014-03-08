//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2014. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt. 
// Contact: khmer-project@idyll.org
//

#ifndef KHMER_CONFIG_HH
#   define KHMER_CONFIG_HH

#   include "khmer.hh"

namespace khmer
{

// Special class for holding configuration values 
// and providing accessors to them.
// There is no need for there to be a singleton.
// Objects of this class can be produced by a factory, 
// allowing multiple extant configurations, should they ever be needed.
struct Config
{

    Config( );

    bool const has_extra_sanity_checks( void ) const;

    uint32_t const get_number_of_threads( void ) const;
    void set_number_of_threads( uint32_t const );

    uint64_t const get_reads_input_buffer_size( void ) const;
    void set_reads_input_buffer_size( uint64_t const );

    uint8_t const get_input_buffer_trace_level( void ) const;
    void set_input_buffer_trace_level( uint8_t const );
    uint8_t const get_reads_parser_trace_level( void ) const;
    void set_reads_parser_trace_level( uint8_t const );

private:
    
    bool	_has_extra_sanity_checks;

    uint32_t	_number_of_threads;

    uint64_t	_reads_input_buffer_size;

    uint8_t	_ibmgr_trace_level;
    uint8_t	_rparser_trace_level;

};

// Get and set default configuration instance.
Config & get_active_config( void );
void set_active_config( Config & );

} // namespace khmer

#endif // KHMER_CONFIG_HH
// vim: set ft=cpp sts=4 sw=4 tw=79:
