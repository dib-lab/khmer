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
  set_active_config( Config & config )
  {
    the_config = config;
  }


  inline
  unsigned int const
  Config::
  detected_number_of_threads( void ) const
  { return 1; }


  Config::
  Config( void )
  {
    char      buf[ 1024 ];
    
    _config_data[ "is_threaded" ]	      = string( "true" );
    sprintf( buf, "%u", detected_number_of_threads( ) );
    _config_data[ "number_of_threads" ]	      = string( buf );
    _config_data[ "reads_parser_threading" ]  = string( "true" );
#ifdef KHMER_EXTRA_SANITY_CHECKS
    _config_data[ "has_extra_sanity_checks" ] = string( "true" );
#else
    _config_data[ "has_extra_sanity_checks" ] = string( "false" );
#endif
    sprintf( buf, "%ull", 512U*1024*1024 );
    _config_data[ "reads_file_chunk_size" ]   = string( buf );
  }


  bool const
  Config::
  is_threaded( void ) 
  { return !_config_data[ "is_threaded" ].compare( "true" ); }


  bool const
  Config::
  has_extra_sanity_checks( void )
  { return !_config_data[ "has_extra_sanity_checks" ].compare( "true" ); }


  unsigned int const
  Config::
  get_number_of_threads( void ) 
  { 
    unsigned int  number_of_threads;

    sscanf(
      (const char *)_config_data[ "number_of_threads" ].c_str( ), 
      "%u", &number_of_threads
    );
    return number_of_threads;
  }


  void
  Config::
  set_number_of_threads( const unsigned int number_of_threads )
  {
    // TODO: Throw exception if 0 threads are requested.
    // TODO: Sanity check to make sure the value is not so large 
    //	     that it would produce unreasonable hash count slop thresholds.

    char	  buf[ 1024 ];

    sprintf( buf, "%u", number_of_threads );
    _config_data[ "number_of_threads" ] = buf;
  }
  

  bool const
  Config::
  get_reads_parser_threading( void )
  { return !_config_data[ "reads_parser_threading" ].compare( "true" ); }


  void
  Config::
  set_reads_parser_threading( const bool threading )
  {
    _config_data[ "reads_parser_threading" ] = 
      threading ? string( "true" ) : string( "false" );
  }


  unsigned long long const
  Config::
  get_reads_file_chunk_size( void )
  {
    unsigned long long  reads_file_chunk_size;

    sscanf(
      (const char *)_config_data[ "reads_file_chunk_size" ].c_str( ), 
      "%llu", &reads_file_chunk_size
    );
    return reads_file_chunk_size;
  }


  void
  Config::
  set_reads_file_chunk_size( unsigned long long reads_file_chunk_size )
  {
    // TODO: Throw exception if chunk size is 0.

    char	  buf[ 1024 ];

    sprintf( buf, "%llu", reads_file_chunk_size );
    _config_data[ "reads_file_chunk_size" ] = buf;
  }


  Byte const
  Config::
  get_hash_count_threshold( void )
  {
    unsigned int    number_of_threads	    = get_number_of_threads( );
    if (1 == number_of_threads) return (Byte)MAX_COUNT;
    return (Byte)(MAX_COUNT - get_number_of_threads( ) + 1);
  }


  BoundedCounterType const
  Config::
  get_hash_bigcount_threshold( void )
  {
    unsigned int    number_of_threads	    = get_number_of_threads( );
    if (1 == number_of_threads) return (BoundedCounterType)MAX_BIGCOUNT;
    return (BoundedCounterType)(MAX_BIGCOUNT - get_number_of_threads( ) + 1);
  }

}

// vim: set sts=2 sw=2:
