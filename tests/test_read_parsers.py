# Tests for the ReadParser and Read classes.


import khmer
from khmer import ReadParser
import khmer_tst_utils as utils
from nose.plugins.attrib import attr

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


# CTB note: one reason this test can fail is if you've completely broken
# multithreaded parsing ;).
def test_with_multiple_threads( ):
    import operator
    import threading

    reads_count_1thr = 0
    rparser = ReadParser( utils.get_test_data( "test-reads.fq.bz2" ) )
    for read in rparser: reads_count_1thr += 1

    def count_reads( rparser, counters, tnum ):
        counters[ tnum ] = sum([ 1 for read in rparser ])

    N_THREADS = 4
    config = khmer.get_config( )
    bufsz = config.get_reads_input_buffer_size( )
    config.set_reads_input_buffer_size( N_THREADS * 64 * 1024 )
    threads = [ ]
    reads_counts_per_thread = [ 0 ] * N_THREADS
    rparser = ReadParser( utils.get_test_data( "test-reads.fq.bz2" ), N_THREADS )
    for tnum in xrange( N_THREADS ):
        t = \
        threading.Thread(
            target = count_reads,
            args = [ rparser, reads_counts_per_thread, tnum ]
        )
        print 'starting'
        threads.append( t )
        t.start( )
    for t in threads:
        print 'joining'
        t.join( )
    config.set_reads_input_buffer_size( bufsz )

    #print reads_counts_per_thread

    assert reads_count_1thr == sum( reads_counts_per_thread ), \
           reads_counts_per_thread


def test_old_illumina_pair_mating( ):
    import threading

    config = khmer.get_config( )
    bufsz = config.get_reads_input_buffer_size( )
    config.set_reads_input_buffer_size( 65600 * 2 )
    # Note: This file, when used in conjunction with a 65600 byte per-thread
    #       prefetch buffer, tests the paired read mating logic with the
    #       old Illumina read name format.
    rparser = ReadParser( utils.get_test_data( "test-reads.fa" ), 2 )

    def thread_1_runtime( rparser ):
        for read in rparser: pass

    def thread_2_runtime( rparser ):
        for readnum, read in enumerate( rparser ):
            if 0 == readnum:
                assert "850:2:1:1198:16820/1" == read.name, (read.name,)

    t1 = threading.Thread( target = thread_1_runtime, args = [ rparser ] )
    t2 = threading.Thread( target = thread_2_runtime, args = [ rparser ] )

    t1.start( )
    t2.start( )

    t1.join( )
    t2.join( )

    config.set_reads_input_buffer_size( bufsz )


def test_casava_1_8_pair_mating( ):
    import threading

    config = khmer.get_config( )
    bufsz = config.get_reads_input_buffer_size( )
    config.set_reads_input_buffer_size( 128 * 1024 )
    # Note: This file, when used in conjunction with a 64 KiB per-thread
    #       prefetch buffer, tests the paired read mating logic with the
    #       Casava >= 1.8 read name format.
    rparser = ReadParser( utils.get_test_data( "test-reads.fq.bz2" ), 2 )

    def thread_1_runtime( rparser ):
        for read in rparser: pass

    def thread_2_runtime( rparser ):
        for readnum, read in enumerate( rparser ):
            if 0 == readnum:
                assert "895:1:1:1761:13189 2:N:0:NNNNN" == read.name

    t1 = threading.Thread( target = thread_1_runtime, args = [ rparser ] )
    t2 = threading.Thread( target = thread_2_runtime, args = [ rparser ] )

    t1.start( )
    t2.start( )

    t1.join( )
    t2.join( )

    config.set_reads_input_buffer_size( bufsz )


def test_iterator_identities( ):
    
    rparser = \
    ReadParser( utils.get_test_data( "test-abund-read-paired.fa" ) )
    assert rparser is rparser.__iter__( )
    assert rparser is rparser.iter_reads( )


@attr('known_failing')
def test_read_pair_iterator_in_error_mode( ):
    assert 0
    
    rparser = \
    ReadParser( utils.get_test_data( "test-abund-read-paired.fa" ) )

    # If walks like an iterator and quacks like an iterator...
    rpi = rparser.iter_read_pairs( )
    assert "__iter__" in dir( rpi )
    assert "next" in dir( rpi )

    # Are the alleged pairs actually pairs?
    read_pairs_1 = [ ]
    for read_1, read_2 in rpi:
        read_pairs_1.append( [ read_1, read_2 ] )
        assert read_1.name[ : 19 ] == read_2.name[ : 19 ]

    # Reload parser.
    # Note: No 'rewind' or 'reset' capability at the time of this writing.
    rparser = \
    ReadParser( utils.get_test_data( "test-abund-read-paired.fa" ) )

    # Ensure that error mode is the default mode.
    read_pairs_2 = [ ]
    for read_1, read_2 \
    in rparser.iter_read_pairs( ReadParser.PAIR_MODE_ERROR_ON_UNPAIRED ):
        read_pairs_2.append( [ read_1, read_2 ] )
    matches = \
    map(
        lambda rp1, rp2: rp1[ 0 ].name == rp2[ 0 ].name, 
        read_pairs_1, read_pairs_2
    )
    assert all( matches )  # Assert ALL the matches. :-]


def test_read_pair_iterator_in_error_mode_xfail( ):

    rparser = \
    ReadParser( utils.get_test_data( "test-abund-read-impaired.fa" ) )

    failed = True
    try:
        for rpair in rparser.iter_read_pairs( ): pass
        failed = False
    except IOError as exc: pass
    assert failed


@attr('known_failing')
def test_read_pair_iterator_in_ignore_mode( ):
    assert 0

    rparser = \
    ReadParser( utils.get_test_data( "test-abund-read-impaired.fa" ) )
    
    read_pairs = [ ]
    for read_1, read_2 \
    in rparser.iter_read_pairs( ReadParser.PAIR_MODE_IGNORE_UNPAIRED ):
        read_pairs.append( [ read_1, read_2 ] )
        assert read_1.name[ : 19 ] == read_2.name[ : 19 ]
    assert 2 == len( read_pairs )


# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
