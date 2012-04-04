"""
    Tests various aspects wrapper for C++ API configuration interface.
"""


import khmer
# NOTE: Currently the wrapper only supports a config singleton.
#	In the future, manipulation of multiple configs may be allowed.
#	The following alias is a hedge against the future.
from khmer	    import get_config		as get_active_config


def test_EXISTENCE_is_threaded( ):
    """
	Verify that 'is_threaded' exists.
	An exception should be thrown if a config object cannot be obtained.
    """
    config = get_active_config( )
    assert "is_threaded" in dir( config )


def test_RETURN_TYPE_is_threaded( ):
    """
	Verify that 'is_threaded" returns a boolean.
	Note: This test is performed here because it determines the nature of some of the later tests.
    """
    config = get_active_config( )
    assert bool is type( config.is_threaded( ) ) 


def check_attribute_exists( config, attr_name ):
    """
	Helper function for testing attribute existence.
    """
    assert True == hasattr( config, attr_name ), attr_name

def test_EXISTENCE_OTHERS( ):
    """
	Verify that all of the various attributes exist.
    """
    config = get_active_config( )
    for attr_name in \
	[
	    "has_extra_sanity_checks", "get_number_of_threads", 
	    "get_reads_parser_threading",
	    "get_reads_file_chunk_size", "set_reads_file_chunk_size",
	    "get_hash_count_threshold", "get_hash_bigcount_threshold",
	]:
	yield check_attribute_exists, config, attr_name
    if config.is_threaded( ):
	for attr_name in \
	    [
		"set_number_of_threads", "set_reads_parser_threading",
	    ]:
	    yield check_attribute_exists, config, attr_name


def test_1_ARGS_set_number_of_threads( ):
    """
	Verify that the number of threads cannot be set to a negative number.
    """
    config = get_active_config( )
    if config.is_threaded( ):
	try: config.set_number_of_threads( -1 );
	except: pass
	else: assert False, "config.set_number_of_threads( -1 )"


def test_2_ARGS_set_number_of_threads( ):
    """
	Verify that the number of threads cannot be set to zero.
    """
    config = get_active_config( )
    if config.is_threaded( ):
	try: config.set_number_of_threads( 0 );
	except: pass
	else: assert False, "config.set_number_of_threads( 0 )"


def test_USE_set_number_of_threads( ):
    """
	Verify that the number of threads set is what is reported.
    """
    config = get_active_config( )
    if config.is_threaded( ):
	config.set_number_of_threads( 8 )
	assert 8 == config.get_number_of_threads( )


def test_USE_set_reads_parser_threading( ):
    """
	Verify that the reads parser threading can be toggled.
    """
    config = get_active_config( )
    if config.is_threaded( ):
	config.set_reads_parser_threading( True )
	assert True == config.get_reads_parser_threading( )


def test_USE_set_reads_file_chunk_size( ):
    """
	Verify that the reads file chunk size is what is reported.
    """
    config = get_active_config( )
    if config.is_threaded( ):
	config.set_reads_file_chunk_size( 123456789L )
	assert 123456789L == config.get_reads_file_chunk_size( )


def check_hash_count_threshold( config, i, max_count, count_func_name ):
    """
	Verify that sloppiness thresholds for counting hashes change appropriately with number of threads.
    """
    number_of_threads = 2**i
    config.set_number_of_threads( number_of_threads )
    if 0 == i:	assert max_count == getattr( config, count_func_name )( )
    else:	assert (max_count - number_of_threads + 1) == getattr( config, count_func_name )( )

def test_EFFECTS_set_number_of_threads( ):
    """
	When multi-threading is enabled, verify that the effects of changing the number of threads are proper.
    """
    config = get_active_config( )
    if config.is_threaded( ):
	config.set_number_of_threads( 1 )
	max_count	= config.get_hash_count_threshold( )
	max_bigcount	= config.get_hash_bigcount_threshold( )
	# Descend from greater number of threads to 1 thread as extra assurance that same answer is produced again after various ops.
	for i in reversed( xrange( 4 ) ):
	    yield check_hash_count_threshold, config, i, max_count, "get_hash_count_threshold"
	    yield check_hash_count_threshold, config, i, max_bigcount, "get_hash_bigcount_threshold"


# vim: set sts=4 sw=4:
