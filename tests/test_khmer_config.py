#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name
"""
    Tests various aspects wrapper for C++ API configuration interface.
"""


# NOTE: Currently the wrapper only supports a config singleton.
#       In the future, manipulation of multiple configs may be allowed.
#       The following alias is a hedge against the future.
from khmer import get_config as get_active_config


def test_EXISTENCE_has_extra_sanity_checks():
    """
        Verify that 'has_extra_sanity_checks' exists.
        An exception should be thrown if a config object cannot be obtained.
    """
    config = get_active_config()
    assert "has_extra_sanity_checks" in dir(config)
    assert isinstance(config.has_extra_sanity_checks(), bool)


def check_attribute_exists(config, attr_name):
    """
        Helper function for testing attribute existence.
    """
    assert True == hasattr(config, attr_name), attr_name


def test_EXISTENCE_OTHERS():
    """
        Verify that all of the various attributes exist.
    """
    config = get_active_config()
    for attr_name in \
        ["set_number_of_threads", "get_number_of_threads",
         "get_reads_input_buffer_size", "set_reads_input_buffer_size", ]:
        yield check_attribute_exists, config, attr_name


# def test_1_ARGS_set_number_of_threads( ):
#    """
#       Verify that the number of threads cannot be set to a negative number.
#    """
#    config = get_active_config( )
#    if config.is_threaded( ):
#       try: config.set_number_of_threads( -1 );
#       except: pass
#       else: assert False, "config.set_number_of_threads( -1 )"


# def test_2_ARGS_set_number_of_threads( ):
#    """
#       Verify that the number of threads cannot be set to zero.
#    """
#    config = get_active_config( )
#    if config.is_threaded( ):
#       try: config.set_number_of_threads( 0 );
#       except: pass
#       else: assert False, "config.set_number_of_threads( 0 )"


def test_USE_set_number_of_threads():
    """
        Verify that the number of threads set is what is reported.
    """
    config = get_active_config()
    tnum = config.get_number_of_threads()
    config.set_number_of_threads(8)
    assert 8 == config.get_number_of_threads()
    config.set_number_of_threads(tnum)
    assert tnum == config.get_number_of_threads()


def test_WRONGETYPE_set_number_of_threads():
    config = get_active_config()
    try:
        config.set_number_of_threads("a")
        assert 0, "Shouldn't be able to set the number of threads to a letter"
    except TypeError, err:
        print str(err)


def test_WRONGTYPE_set_reads_input_buffer():
    config = get_active_config()
    try:
        config.set_reads_input_buffer_size("a")
        assert 0, ("Shouldn't be able to set the reads input buffer size to a "
                   "letter")
    except TypeError, err:
        print str(err)


def test_USE_set_reads_input_buffer_size():
    """
        Verify that the reads file chunk size is what is reported.
    """
    config = get_active_config()
    bufsz = config.get_reads_input_buffer_size()
    config.set_reads_input_buffer_size(123456789L)
    assert 123456789L == config.get_reads_input_buffer_size()
    config.set_reads_input_buffer_size(bufsz)
    assert bufsz == config.get_reads_input_buffer_size()


# vim: set ft=python sts=4 sw=4 tw=79:
