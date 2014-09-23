import screed
import os
from screed.DBConstants import fileExtension

def test_nodb():
    """
    Tests if screed throws an appropriate exception if it is
    asked to open a non-existant screed database
    """
    try:
        db = screed.ScreedDB('foo')
        assert 1 == 0 # Previous line should throw an error
    except ValueError:
        pass

def test_wrongdb():
    """
    Tests if screed throws an appropriate exception if it is
    asked to open a file that isn't a screed database
    """
    try:
        blah = 'blah_screed'
        blah_file = open(blah, 'wb')
        blah_file.close()
        
        db = screed.ScreedDB(blah)
        os.unlink(blah)
        assert 1 == 0
    except TypeError:
        os.unlink(blah)
        pass
