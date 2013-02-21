# Tests for the ReadParser and Read classes.


import khmer
from khmer import ReadParser
import khmer_tst_utils as utils


def test_DEFAULT_ARGUMENTS( ):
    
    read_names = [ ]
    # Note: Using a data file where read names are just integers on [0,99).
    rparser = ReadParser( utils.get_test_data( "random-20-a.fa" ) )

    for read in rparser:
        read_names.append( int( read.name ) )
    
    # "Derandomize".
    read_names.sort( )

    # Each read number should match the corresponding name.
    for m, n in enumerate( read_names ):
        assert m == n


# TODO: Write more tests.


# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
