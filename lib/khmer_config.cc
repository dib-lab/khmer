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
  const
  unsigned int
  Config::
  detected_number_of_threads( void ) const
  {
#ifdef KHMER_THREADED
    // NOTE: Assume OpenMP for now.
    //return omp_get_num_threads( );
    return 1;
#else
    return 1;
#endif
  }


  Config::
  Config( void )
  {
    char      buf[ 1024 ];
    
#ifdef KHMER_THREADED
    _config_data[ "is_threaded" ]	      = string( "true" );
    sprintf( buf, "%u", detected_number_of_threads( ) );
    _config_data[ "number_of_threads" ]	      = string( buf );
#else
    _config_data[ "is_threaded" ]	      = string( "false" );
    sprintf( buf, "%u", 1 );
    _config_data[ "number_of_threads" ]	      = string( buf );
#endif
#ifdef KHMER_EXTRA_SANITY_CHECKS
    _config_data[ "has_extra_sanity_checks" ] = string( "true" );
#else
    _config_data[ "has_extra_sanity_checks" ] = string( "false" );
#endif
    
  }


  const
  bool
  Config::
  is_threaded( void ) 
  { return !_config_data[ "is_threaded" ].compare( "true" ); }


  const
  bool
  Config::
  has_extra_sanity_checks( void )
  { return !_config_data[ "has_extra_sanity_checks" ].compare( "true" ); }


  const
  unsigned int
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


#ifdef KHMER_THREADED
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
#endif


  const
  Byte
  Config::
  get_hash_count_threshold( void )
  {
#ifdef KHMER_THREADED
    unsigned int    number_of_threads	    = get_number_of_threads( );
    if (1 == number_of_threads) return (Byte)MAX_COUNT;
    return (Byte)(MAX_COUNT - get_number_of_threads( ) + 1);
#else
    return (Byte)MAX_COUNT;
#endif
  }


  const
  BoundedCounterType
  Config::
  get_hash_bigcount_threshold( void )
  {
#ifdef KHMER_THREADED
    unsigned int    number_of_threads	    = get_number_of_threads( );
    if (1 == number_of_threads) return (BoundedCounterType)MAX_BIGCOUNT;
    return (BoundedCounterType)(MAX_BIGCOUNT - get_number_of_threads( ) + 1);
#else
    return (BoundedCounterType)MAX_BIGCOUNT;
#endif
  }

}

// vim: set sts=2 sw=2:
