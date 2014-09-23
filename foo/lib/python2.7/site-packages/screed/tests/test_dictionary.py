import os
import screed
from screed.DBConstants import fileExtension

class Test_dict_methods(object):
    """
    Make sure that screed returns sensible results for standard dictionary
    queries.
    """
    def setup(self):
        self._testfa = os.path.join(os.path.dirname(__file__), 'test.fa')
        screed.read_fasta_sequences(self._testfa)
        self.db = screed.ScreedDB(self._testfa)

    def teardown(self):
        os.unlink(self._testfa + fileExtension)

    def test_iter_stuff(self):
        db = self.db
        keys = db.keys()
        ikeys = list(db.iterkeys())
        assert sorted(keys) == sorted(ikeys)

        values = db.values()
        ivalues = list(db.itervalues())
        assert sorted(values) == sorted(ivalues)

        items = db.items()
        iitems = list(db.iteritems())
        assert sorted(items) == sorted(iitems)

    def test_contains(self):
        for k in self.db:
            assert self.db.has_key(k)

        assert db.get('FOO') == None

        assert not self.db.has_key('FOO')
            
    def test_contains(self):
        for k in self.db:
            assert k in self.db

        assert not 'FOO' in self.db

    def test_get(self):
        for k in self.db:
            record = self.db.get(k)
            assert record.name == k

            record = self.db[k]
            assert record.name == k

        try:
            self.db['FOO']
            assert False, "the previous line should raise a KeyError"
        except KeyError:
            pass

    def test_missing(self):
        """
        Make sure that unsupported dict attributes are actually missing.
        """
        db = self.db

        try:
            db.clear()
            assert 0
        except AttributeError:
            pass

        try:
            db.update({})
            assert 0
        except AttributeError:
            pass

        try:
            db.setdefault(None)
            assert 0
        except AttributeError:
            pass

        try:
            db.pop()
            assert 0
        except AttributeError:
            pass

        try:
            db.popitem()
            assert 0
        except AttributeError:
            pass
