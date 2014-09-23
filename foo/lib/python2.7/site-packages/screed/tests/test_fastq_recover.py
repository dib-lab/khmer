import test_fastq
import os
import screed
from screed.DBConstants import fileExtension

class test_fq_recover(test_fastq.Test_fastq):
    def setup(self):
        thisdir = os.path.dirname(__file__)
        self._fileName = os.path.join(thisdir, 'fastqRecovery')
        self._testfq = os.path.join(thisdir, 'test.fastq')
        screed.read_fastq_sequences(self._testfq)
        screed.ToFastq(self._testfq, self._fileName)
        screed.read_fastq_sequences(self._fileName)
        self.db = screed.ScreedDB(self._fileName)

    def teardown(self):
        os.unlink(self._fileName)
        os.unlink(self._fileName + fileExtension)
        os.unlink(self._testfq + fileExtension)
