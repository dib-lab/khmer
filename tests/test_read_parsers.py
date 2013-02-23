# Tests for the ReadParser and Read classes.


import khmer
from khmer import ReadParser
import khmer_tst_utils as utils


def test_read_properties( ):
    
    # Note: Using a data file with only one read.
    rparser = ReadParser( utils.get_test_data( "single-read.fq" ) )

    # Check the properties of all one reads in data set.
    for read in rparser:
        assert read.name == "895:1:1:1246:14654 1:N:0:NNNNN"
        assert read.sequence == "CAGGCGCCCACCACCGTGCCCTCCAACCTGATGGT"
        assert read.annotations == ""
        assert read.accuracy == """][aaX__aa[`ZUZ[NONNFNNNNNO_____^RQ_"""


def test_with_default_arguments( ):
    
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


def test_gzip_decompression( ):
    
    reads_count = 0
    rparser = ReadParser( utils.get_test_data( "100-reads.fq.gz" ) )
    for read in rparser:
        reads_count += 1

    assert 100 == reads_count


def test_bzip2_decompression( ):
    
    reads_count = 0
    rparser = ReadParser( utils.get_test_data( "100-reads.fq.bz2" ) )
    for read in rparser:
        reads_count += 1

    assert 100 == reads_count


def test_with_multiple_threads( ):
    
    import operator
    import threading

    reads_count_1thr = 0
    rparser = ReadParser( utils.get_test_data( "test-reads.fq.bz2" ) )
    for read in rparser: reads_count_1thr += 1

    def count_reads( rparser, counters, tnum ):
        counters[ tnum ] = reduce( operator.add, ( 1 for read in rparser ) )

    N_THREADS = 4
    config = khmer.get_config( )
    bufsz = config.get_reads_input_buffer_size( )
    config.set_reads_input_buffer_size( N_THREADS * 64 * 1024 )
    threads = [ ]
    reads_counts_per_thread = [ 0 ] * N_THREADS
    # Note: This file, when used in conjunction with a 64 KiB per-thread
    #       prefetch buffer, also tests the paired read mating logic with the
    #       Casava >= 1.8 read name format.
    rparser = ReadParser( utils.get_test_data( "test-reads.fq.bz2" ), N_THREADS )
    for tnum in xrange( N_THREADS ):
        t = \
        threading.Thread(
            target = count_reads,
            args = [ rparser, reads_counts_per_thread, tnum ]
        )
        threads.append( t )
        t.start( )
    for t in threads:
        t.join( )
    config.set_reads_input_buffer_size( bufsz )

    assert reads_count_1thr == sum( reads_counts_per_thread )


# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
