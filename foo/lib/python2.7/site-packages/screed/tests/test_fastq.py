import screed
from screed.DBConstants import fileExtension
import os
from cStringIO import StringIO

def test_new_record():
    # test for a bug where the record dict was not reset after each
    # sequence load, leading to all records being identical if you
    # kept a handle on the returned dictionary.
    
    s = StringIO("@1\nACTG\n+\nAAAA\n@2\nACGG\n+\nAAAA\n")

    records = list(iter(screed.fastq.fastq_iter(s)))
    assert records[0]['name'] == '1'
    assert records[1]['name'] == '2'

def test_parse_description_true():
    # test for a bug where the record dict was not reset after each
    # sequence load, leading to all records being identical if you
    # kept a handle on the returned dictionary.
    
    s = StringIO("@1 FOO\nACTG\n+\nAAAA\n@2\nACGG\n+\nAAAA\n")

    records = list(iter(screed.fastq.fastq_iter(s, parse_description=True)))
    assert records[0]['name'] == '1'
    assert records[1]['name'] == '2'

    # also is default behavior
    s = StringIO("@1 FOO\nACTG\n+\nAAAA\n@2\nACGG\n+\nAAAA\n")

    records = list(iter(screed.fastq.fastq_iter(s)))
    assert records[0]['name'] == '1'
    assert records[1]['name'] == '2'

def test_parse_description_false():
    # test for a bug where the record dict was not reset after each
    # sequence load, leading to all records being identical if you
    # kept a handle on the returned dictionary.
    
    s = StringIO("@1 FOO\nACTG\n+\nAAAA\n@2\nACGG\n+\nAAAA\n")

    records = list(iter(screed.fastq.fastq_iter(s, parse_description=False)))
    assert records[0]['name'] == '1 FOO'
    assert records[1]['name'] == '2'

class Test_fastq(object):
    def setup(self):
        self._testfq = os.path.join(os.path.dirname(__file__), 'test.fastq')
        screed.read_fastq_sequences(self._testfq)
        self.db = screed.ScreedDB(self._testfq)

    def teardown(self):
        os.unlink(self._testfq + fileExtension)

    def test_length(self):
        assert len(self.db) == 125

    def test_keys(self):
        for key in self.db:
            assert key == self.db[key].name

    def test_id_retrieval(self):
        for key in self.db:
            record = self.db[key]
            intRcrd = self.db.loadRecordByIndex(record.id)
            assert record == intRcrd

    def test_contains_front(self):
        first = self.db[self.db.keys()[0]]
        assert first.id == 0
        assert first.name == 'HWI-EAS_4_PE-FC20GCB:2:1:492:573/2'
        assert first.sequence == 'ACAGCAAAATTGTGATTGAGGATGAAGAACTGCTGT'
        assert first.accuracy == 'AA7AAA3+AAAAAA.AAA.;7;AA;;;;*;<1;<<<'

    def test_contains_middle(self):
        middle = self.db[self.db.keys()[62]]
        assert middle.id == 62
        assert middle.name == 'HWI-EAS_4_PE-FC20GCB:2:1:245:483/2'
        assert middle.sequence == 'TGTCGAGCAAAGCAAAACAGGCGTAAAAATTGCCAT'
        assert middle.accuracy == 'AAAAAAAAAAAAAAAAAAAAA>AAAAAAAA?9>6><'

    def test_contains_end(self):
        end = self.db[self.db.keys()[124]]
        assert end.id == 124
        assert end.name == 'HWI-EAS_4_PE-FC20GCB:2:1:350:588/2'
        assert end.sequence == 'GGTACAAAATAGATGCTGGACTCTCCGAATCCTATA'
        assert end.accuracy == ';?5AAAAAAAAAA?A??;?AA;AAA>AAAA?4?844'

    def test_contains(self):
        for k in self.db:
            assert self.db.has_key(k)

        assert self.db.get('FOO') == None

        assert not 'FOO' in self.db

    def test_iterv(self):
        entries = []
        for entry in self.db:
            entries.append(self.db[entry])

        ivalues = list(self.db.itervalues())
        assert sorted(entries) == sorted(ivalues)

    def test_iteri(self):
        for id, entry in self.db.iteritems():
            assert id == self.db[entry.name].id
            assert entry == self.db[entry.name]

def test_writer():
    fp = StringIO()
    w = screed.fastq.FASTQ_Writer("", fp)

    class FakeRecord(object):
        pass
    
    read = FakeRecord()
    read.name = 'foo'
    read.description = 'bar'
    read.sequence = 'ATCG'
    read.accuracy = '####'

    w.write(read)

    assert fp.getvalue() == '@foo bar\nATCG\n+\n####\n'
    
def test_writer_2():
    fp = StringIO()
    w = screed.fastq.FASTQ_Writer("", fp)

    class FakeRecord(object):
        pass
    
    read = FakeRecord()
    read.name = 'foo'
    read.description = 'bar'
    read.sequence = 'ATCG'
    read.accuracy = '####'

    read_iter = [read]

    w.consume(read_iter)

    assert fp.getvalue() == '@foo bar\nATCG\n+\n####\n'
