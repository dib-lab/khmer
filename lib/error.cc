//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2014. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt. 
// Contact: khmer-project@idyll.org
//

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
