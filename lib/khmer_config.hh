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

    Config( void );

    bool const has_extra_sanity_checks( void );

    uint32_t const get_number_of_threads( void );
    void set_number_of_threads( uint32_t const );

    uint64_t const get_reads_input_buffer_size( void );
    void set_reads_input_buffer_size( uint64_t const );

private:
    
    bool	_has_extra_sanity_checks;
    uint32_t	_number_of_threads;
    uint64_t	_reads_input_buffer_size;

};

// Get and set default configuration instance.
Config & get_active_config( void );
void set_active_config( Config & );

} // namespace khmer

#endif // KHMER_CONFIG_HH
// vim: set ft=cpp sts=4 sw=4 tw=79:
