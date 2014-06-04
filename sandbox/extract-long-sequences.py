#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import argparse
import screed
import sys

min_length = int(sys.argv[1])

for filename in sys.argv[2:]:
    for record in screed.open(filename):
        if len(record['sequence']) >= min_length:
            #if fasta
                print '>%s\n%s' % (record['name'], record['sequence'],)
            #elif fastq
                if hasattr(record, 'accuracy'):
                    outfp.write(
                        '@{name}\n{seq}\n'
                        '+\n{acc}\n'.format(name=record.name,
                                            seq=record.sequence,
                                            acc=record.accuracy))
                else:
                    outfp.write(
                        '>{name}\n{seq}\n'.format(name=record.name,
                                                    seq=record.sequence))
