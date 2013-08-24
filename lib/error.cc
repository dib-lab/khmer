#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdarg>


void error( int status, int errnum, char const * format, ... )
{
    va_list posargs;

    fflush( stdout );
    if (errnum)
        fprintf( stderr, "\n%s: ", strerror( errnum ) );
    va_start( posargs, format );
    vfprintf( stderr, format, posargs );
    va_end( posargs );
    fprintf( stderr, "\n" );
    fflush( stderr );
    exit( status );
}


// vim: set ft=cpp sts=4 sw=4 et tw=79:
