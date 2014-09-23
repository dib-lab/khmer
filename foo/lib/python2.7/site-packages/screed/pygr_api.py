# Copyright (c) 2008-2010, Michigan State University

"""
A simple wrapper implementing a pygr-compatible SequenceDB based on screed.

There are two implementions:
 - ScreedSequenceDB
 - ScreedSequenceDB_ByIndex

ScreedSequenceDB uses the sequence name as the sequence ID, which
mimics the behavior of pygr's SequenceFileDB and is good for
small-to-medium sized collections of sequences.
ScreedSequenceDB_ByIndex uses the sequence's index (0...size of
database) as a sequence ID, rather than the sequence name; this is
much faster for databases with many, many sequences.

Unlike the normal seqdb, screed will load the entire sequence record
into memory on request, so it's not good for large sequences.

All screed records are guaranteed to have an 'index', a 'name', and a
'sequence' attribute; anything else is specific to the database writer
you use.  The raw screed record (which contains any other information)
is available under seqObj.record.

Note: the underlying screed database must already have been built with
fadbm or fqdbm.

CTB 3/20/09
"""

import UserDict

from screed import ScreedDB

from pygr.sequence import SequenceBase
from pygr.seqdb import SequenceDB
from pygr.sequtil import DNA_SEQTYPE

###

class ScreedSequence(SequenceBase):
    """Sequence implementation based on screed; stores screed record info.

    Attributes:
      - 'id' and 'db' are the standard pygr-ish name/database attrs.
      - 'record' is the screed 'record' object, containing name, etc.
      - 'name' is the record name, which can be the same as 'id' but
        can also be different (see ScreedSequenceDB_ByIndex).
      - 'seq' is the sequence.

    """
    def __init__(self, db, id):
        self.id = id
        SequenceBase.__init__(self)
        info = db.seqInfoDict[id]
        
        self.record = info.record
        self.name = info.record.name
        self.seq = info.record.sequence

class ScreedSequenceDB(SequenceDB):
    """SequenceDB implementation based on screed; retrieve seqs by name."""
    itemClass = ScreedSequence

    def __init__(self, filepath):
        self.filepath = filepath
        self.seqInfoDict = _ScreedSeqInfoDict_ByName(filepath)
        SequenceDB.__init__(self)

    def _set_seqtype(self):
        self._seqtype = DNA_SEQTYPE

    def __repr__(self):
        return "<%s '%s'>" % (self.__class__.__name__, self.filepath)

    # override inherited __reduce__/__getstate__/__setstate__ from SequenceDB.
    def __reduce__(self):
        return (ScreedSequenceDB, (self.filepath,))
    
class ScreedSequenceDB_ByIndex(SequenceDB):
    """SequenceDB implementation based on screed; retrieve seqs by index."""
    itemClass = ScreedSequence
    
    def __init__(self, filepath):
        self.filepath = filepath
        self.seqInfoDict = _ScreedSeqInfoDict_ByIndex(filepath)
        SequenceDB.__init__(self)
        
    def _set_seqtype(self):
        self._seqtype = DNA_SEQTYPE

    def __repr__(self):
        return "<%s '%s'>" % (self.__class__.__name__, self.filepath)

    # override inherited __reduce__/__getstate__/__setstate__ from SequenceDB.
    def __reduce__(self):
        return (ScreedSequenceDB_ByIndex, (self.filepath,))
    
class _ScreedSequenceInfo(object):
    """Objects to put in seqInfoDict values, for holding screed record info."""
    def __init__(self, id, record):
        self.id = id
        self.record = record
        self.length = len(record.sequence)

class _ScreedSeqInfoDict_ByName(object, UserDict.DictMixin):
    """seqInfoDict implementation that uses names to retrieve records."""
    def __init__(self, filepath):
        self.sdb = ScreedDB(filepath)

    def __getitem__(self, k):
        v = self.sdb[k]
        return _ScreedSequenceInfo(k, v)

    def keys(self):
        return self.sdb.keys()

    def itervalues(self):
        i = 0
        max_index = len(self.sdb)
        while i < max_index:
            v = self.sdb.loadRecordByIndex(i)
            yield _ScreedSequenceInfo(v.name, v)
            i += 1

    def iteritems(self):
        for v in self.itervalues():
            yield v.record.name, v
        

class _ScreedSeqInfoDict_ByIndex(object, UserDict.DictMixin):
    """seqInfoDict implementation that uses indices to retrieve records."""
    def __init__(self, filepath):
        self.sdb = ScreedDB(filepath)

    def __getitem__(self, k):
        n = int(k) 
        v = self.sdb.loadRecordByIndex(n)
        return _ScreedSequenceInfo(k, v)

    def keys(self):
        return xrange(0, len(self.sdb))

    def iterkeys(self):
        i = 0
        max_index = len(self.sdb)
        while i < max_index:
            yield i
            i += 1

###

if __name__ == '__main__':
    import sys
    filename = sys.argv[1]

    db = ScreedSequenceDB(filename)
    for k in db:
        print k, repr(db[k]), db[k].name

    db = ScreedSequenceDB_ByIndex(filename)
    for k in db:
        print k, repr(db[k]), db[k].name
