//
// This file is part of khmer, https://github.com/dib-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2015. It is licensed under
// the three-clause BSD license; see LICENSE.
// Contact: khmer-project@idyll.org
//

#include <fcntl.h>

#include "trace_logger.hh"
#include "khmer_exception.hh"

namespace khmer
{


#ifdef WITH_INTERNAL_TRACING
TraceLogger::
TraceLogger( uint8_t const level, FILE * stream_handle )
    : _level( level ), _shared_stream( true ), _stream_handle( stream_handle )
{
    if( !(NULL != stream_handle) ) {
        throw khmer_exception();
    }
}
#endif


TraceLogger::
TraceLogger( uint8_t const level, char const * const file_name_format, ... )
#ifdef WITH_INTERNAL_TRACING
    : _level( level ), _shared_stream( false )
{
    char	tfile_name[ FILENAME_MAX + 1 ];
    va_list	varargs;

    va_start( varargs, file_name_format );
    vsnprintf( tfile_name, FILENAME_MAX, file_name_format, varargs );
    va_end( varargs );

    _stream_handle = fopen( tfile_name, "w" );
    if (NULL == _stream_handle) {
        std::ostringstream msg;
        const char *err_str = strerror(errno);
        msg << "Invalid trace filename " << tfile_name << std::endl;
        msg << err_str << std::endl;
        throw InvalidStream(msg.str());
    }

}
#else	// WITH_INTERNAL_TRACING
{ }
#endif	// !WITH_INTERNAL_TRACING


TraceLogger::
~TraceLogger( )
#ifdef WITH_INTERNAL_TRACING
{

    if ((!_shared_stream) && (NULL != _stream_handle)) {
        fclose( _stream_handle );
        _stream_handle = NULL;
    }

}
#else	// WITH_INTERNAL_TRACING
{ }
#endif	// !WITH_INTENRAL_TRACING

} // namespace khmer


// vim: set ft=cpp sts=4 sw=4 tw=80:
