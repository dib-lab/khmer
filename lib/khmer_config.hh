#ifndef KHMER_CONFIG_HH
# define KHMER_CONFIG_HH

# include <string>
# include <map>

using namespace std;

# include "khmer.hh"

namespace khmer
{

  // Special class for holding ocnfiguration values 
  // and providing accessors to them.
  // There is no need for there to be a singleton.
  // Objects of this class can be produced by a factory, 
  // allowing multiple extant configurations, should they ever be needed.
  // If we were content with a singleton, then it would be eaiser just to 
  // do things in a khmer:: config namespace.
  class Config
  {

  private:

    map< string, string >	      _config_data;
    
    const unsigned int detected_number_of_threads( void ) const;

  public:
    
    Config( void );
    
    const bool is_threaded( void );
    const bool has_extra_sanity_checks( void );
    
    const unsigned int get_number_of_threads( void );
    void set_number_of_threads( const unsigned int );
    
    // Control multi-threaded parsing of reads files.
    const bool get_reads_parser_threading( void );
    void set_reads_parser_threading( const bool );
    const unsigned long long get_reads_file_chunk_size( void );
    void set_reads_file_chunk_size( const unsigned long long );

    // Calculate the saturation thresholds for hash counts.
    // This changes depending on whether multi-threading is enabled.
    const Byte get_hash_count_threshold( void );
    const BoundedCounterType get_hash_bigcount_threshold( void );

  };

  Config & get_active_config( void );
  void set_active_config( Config & );

}

#endif
// vim: set sts=2 sw=2:
