import sys
import os
import shutil
from cStringIO import StringIO
import traceback

import khmer_tst_utils as utils
import khmer
import khmer.file
import screed


def scriptpath(script):
    return script


def test_norm_by_median():
    infile = utils.get_test_data('R.pe.qc.fq.gz')
    hashfile = utils.get_test_data('normC20k20.kh')
    script = scriptpath('normalize-by-median.py')
    args = ['--loadtable', hashfile, infile]
    (status, out, err) = utils.runscript(script, args)
    assert status == 0, (out, err)
    print (out, err)
    assert os.path.exists(hashfile), hashfile
