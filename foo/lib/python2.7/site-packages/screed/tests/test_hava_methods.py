import screed
import screed.seqparse
from screed.DBConstants import fileExtension
import os

testha = os.path.join(os.path.dirname(__file__), 'test.hava')

class test_hava(object):
    def setup(self):
        screed.seqparse.read_hava_sequences(testha)
        self._db = screed.ScreedDB(testha)

    def teardown(self):
        b = 7
        #os.unlink(testha + fileExtension)

    def test_contains(self):
        assert 'test_006' in self._db

    def test_beginning_key_retrieval(self):
        result = self._db['test_000']
        assert result.hava == 'test_000'
        assert result.quarzk == 'ACGGTGACGGTCACCGTCGACGGCCCAAGCCCATCGAACG'\
               'TACCACCCCCACCTATCGTCACGCTGGTGGAGAGCCAATG'
        assert result.muchalo == 'AFPPCLHBCCILGMMOCHKNNDBKCCPNHAMKJOCCDJA'\
               'OEPNMHFHCBAJOKEMMMBHCPHIOAEPFFCAOJPGIMKGK'
        assert result.fakours == '218583165871861127719451483455294521865'\
               '68176931571171542294878855181415261425688'
        assert result.selimizicka == 'bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb'\
               'bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb'
        assert result.marshoon == 'C7AF246AC7AAEABE5A557FCBC6FD5F5263BCDE'\
               '4E745BEF1GG7DD1AB511GBC63A4GF1F4E1A154B35D'

    def test_middle_key_retrieval(self):
        result = self._db['test_0063']
        assert result.hava == 'test_0063'
        assert result.quarzk == 'CAACACGATCAAGTTTGGTAAGAATTCCGCCTTAAGCTTT'\
               'CTAGAACGATAGTTGCCCCCAATCTGGTTCGAAATCTCTT'
        assert result.muchalo == 'GMDAPLMOOFANDHHMLBPIKGHIAFFFOABFMNNJNIJ'\
               'ILEEFEPOCAJLNDLIFBPMGKOFJIEFAHNJPIOFAJMLM'
        assert result.fakours == '392363971393898522756138876485334274384'\
               '39122136418369146118333919885587613673488'
        assert result.selimizicka == 'bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb'\
               'bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb'
        assert result.marshoon == 'FC25E2CFC2BAFA7A2AA4757F3GFFFEE37G7752'\
               'FCDBAEADBA1AC7374FB5C15552E6E2GG6GFF62C6GE'

    def test_end_key_retrieval(self):
        result = self._db['test_00124']
        assert result.hava == 'test_00124'
        assert result.quarzk == 'ATCGCAACCGTTTCCCCTATCTGGCAATTGAATCCGCGTC'\
               'CTAAAACGAAAGCTTATCCCTGGCGAGGCACGCTAGGCCT'
        assert result.muchalo == 'CIHNCECANFNLKGCHNOEHJDHADHPAEMMNKGMMMPD'\
               'OBMOCKNBCMCPHEBEOINHMBMMGCHEMOIOAPEFPDDJP'
        assert result.fakours == '327364511483537131695325595876269716778'\
               '14946924334424648676283848861393812686731'
        assert result.selimizicka == 'bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb'\
               'bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb'
        assert result.marshoon == '4FE5FDD76CC5DE4DC2F25AA2GFBD7BEG326C6D'\
               '7AB5B71GA67BAFD63AE1A562CDC1C2D157G6EF17CD'

