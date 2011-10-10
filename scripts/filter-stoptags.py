"""
Trim sequences at k-mers in the given stoptags file.  Output sequences
will be placed in 'infile.stopfilt'.

% python scripts/filter-stoptags.py <stoptags> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""

import sys, screed.fasta, os
import khmer
import argparse
from khmer.thread_utils import ThreadedSequenceProcessor, verbose_loader

# @CTB K should be loaded from file...
DEFAULT_K = 32

###

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-k', default=DEFAULT_K, type=int, help='k-mer size',
                        dest='ksize')
    parser.add_argument('stoptags_file')
    parser.add_argument('input_filenames', nargs='+')

    args = parser.parse_args()
    K = args.ksize

    stoptags = args.stoptags_file
    infiles = args.input_filenames

    print 'loading stop tags, with K', K
    ht = khmer.new_hashbits(K, 1, 1)
    ht.load_stop_tags(stoptags)

    def process_fn(record):
        name = record['name']
        seq = record['sequence']
        if 'N' in seq:
            return None, None

        trim_seq, trim_at = ht.trim_on_stoptags(seq)

        if trim_at >= K:
            return name, trim_seq

        return None, None

    ### the filtering loop
    for infile in infiles:
       print 'filtering', infile
       outfile = os.path.basename(infile) + '.stopfilt'

       outfp = open(outfile, 'w')

       tsp = ThreadedSequenceProcessor(process_fn)
       tsp.start(verbose_loader(infile), outfp)

       print 'output in', outfile

if __name__ == '__main__':
    main()
