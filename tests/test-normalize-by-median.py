#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import khmer
from screed.fasta import fasta_iter
from nose.plugins.attrib import attr

import khmer_tst_utils as utils

def t_normalize_by_median_fpr():
    MIN_TABLESIZE_PARAM = 1
    infile = utils.get_temp_filename('test-fpr.fq')
    in_dir = os.path.dirname(infile)
    shutil.copyfile(utils.get_test_data('test-fastq-reads.fq'), infile)
    script = scriptpath('normalize-by-median.py')
    args = ['-f', '-k 17', '-x ' + str(MIN_TABLESIZE_PARAM), infile]
    (status, out, err) = utils.runscript(script, args, in_dir, fail_ok=True)
    assert os.path.exists(infile)
