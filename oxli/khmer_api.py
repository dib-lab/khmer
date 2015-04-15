#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#

from screed.screedRecord import _screed_record_dict
import os
from khmer.utils import write_record, write_record_pair

def diginorm(input_stream, ct, coverage):
    n = 0
    discard = 0
    for _, is_pair, read1, read2 in input_stream:
        if is_pair:
            med1, _, _ = ct.get_median_count(read1.sequence)
            med2, _, _ = ct.get_median_count(read2.sequence)

            if med1 < coverage or med2 < coverage:
                ct.consume(read1.sequence)
                ct.consume(read2.sequence)
                yield n, True, read1, read2
                n += 2
        else:
            med, _, _ = ct.get_median_count(read1.sequence)
            if med < coverage:
                ct.consume(read1.sequence)
                yield n, False, read1, None
                n += 1