#! /usr/bin/env python
"""
Test script for threaded reading/writing.  Should output files identical
in content (but potentially shuffled) to those given to it.  -x turns
off threading, -p turns on paired-end required.

% python sandbox/thread-rw.py [ -p ] [ -x ] [ --threads N ] <input files>

Use '-h' for parameter help.
"""
import sys
import os
import threading
import time
import argparse
import screed

import khmer

from khmer.threading_args import add_threading_args
from khmer import thread_utils
from khmer.thread_utils import ThreadedProcessor, PairThreadedProcessor

###

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-x', action="store_false", dest="do_threading")
    parser.add_argument('-p', action="store_true", dest="paired")
    parser.add_argument('input_filenames', nargs='+')
    add_threading_args(parser)

    args = parser.parse_args()

    infiles = args.input_filenames
    n_threads = max(int(args.n_threads), 2)

    if args.do_threading:
        print 'n_threads:', n_threads
    else:
        print 'NO THREADS FOR YOU!'

    config = khmer.get_config()
    config.set_reads_input_buffer_size(n_threads * 64 * 1024)

    ### the filtering function - return sequence, rest is taken care of.
    def filter_fn(name, seq):
        time.sleep(0.00001)
        return seq
    
    def pair_filter_fn(name, seq, name2, seq2):
        time.sleep(0.00001)
        return seq, seq2
        
    ### the filtering loop
    for n, filename in enumerate(infiles):
        print 'filtering', filename
        outfile = os.path.basename(filename) + '.foo'
        
        print 'output in', outfile
        outfp = open(outfile, 'w')

        # create and start a threaded writer
        if args.paired:
            tw = PairThreadedProcessor(outfp).start()
            ffn = pair_filter_fn
        else:
            tw = ThreadedProcessor(outfp).start()
            ffn = filter_fn

        threads = []
        if args.do_threading:
            # create multithreaded readparser
            rparser = khmer.ReadParser(filename, n_threads - 1)
            threads = thread_utils.start_threads(n_threads - 1,
                                                 target=tw.process_fn,
                                                 args=(rparser, ffn))

        else:
            # no threads except for the writer thread
            rparser = khmer.ReadParser(filename, 1)
            tw.process_fn(rparser, ffn)
                        
        # wait for threads to finish & flush out any remaining records.
        tw.join(threads)

if __name__ == '__main__':
    main()
