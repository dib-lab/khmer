#! /usr/bin/env python
"""
Eliminate reads with median k-mer abundance higher than
DESIRED_COVERAGE.  Output sequences will be placed in 'infile.keep'.

% python scripts/normalize-by-median.py [ -C <cutoff> ] <data1> <data2> ...

Use '-h' for parameter help.
"""

import sys
import screed
import os
import argparse

import khmer
from khmer.counting_args import build_construct_args, DEFAULT_MIN_HASHSIZE
from khmer.threading_args import add_threading_args
from khmer import thread_utils
from khmer.thread_utils import ThreadedProcessor, PairThreadedProcessor, \
                               FilterReporter

DEFAULT_DESIRED_COVERAGE = 10

def normalize_by_median(input_filename, outfp, ht, args, n_threads,
                        report_fp=None):
    def single_filter_fn(name, seq):
        if len(seq) < K:
            return None

        seq = seq.replace('N', 'A')
        med, _, _ = ht.get_median_count(seq)

        if med < DESIRED_COVERAGE:
            ht.consume(seq)
            return seq

    def pair_filter_fn(n1, s1, n2, s2):
        if single_filter_fn(n1, s1) or single_filter_fn(n2, s2):
            return s1, s2
        else:
            return None, None

    DESIRED_COVERAGE = args.cutoff
    K = ht.ksize()

    if args.paired:
        processor_class = PairThreadedProcessor
        filter_fn = pair_filter_fn
    else:
        processor_class = ThreadedProcessor
        filter_fn = single_filter_fn

    reporter = FilterReporter()

    tw = processor_class(outfp, reporter=reporter).start()
    rparser = khmer.ReadParser(input_filename, n_threads - 1)
    threads = thread_utils.start_threads(n_threads - 1,
                                         target=tw.process_fn,
                                         args=(rparser, filter_fn))

    if not tw.join(threads):
        print >>sys.stderr, "** failure.  See message above."
        raise IOError

    return reporter.n_read, reporter.n_read - reporter.n_saved

def handle_error(error, output_name, input_name, ht):
    print >>sys.stderr, '** ERROR:', error
    print >>sys.stderr, '** Failed on {}: '.format(input_name)
    hashname = os.path.basename(input_name) + '.ht.failed'
    print >>sys.stderr, '** ...dumping hashtable to {}'.format(hashname)
    ht.save(hashname)
    try:
        os.remove(output_name)
    except:
        print >>sys.stderr, '** ERROR: problem removing erroneous .keep file'

def main():
    parser = build_construct_args()
    add_threading_args(parser)
    parser.add_argument('-C', '--cutoff', type=int, dest='cutoff',
                        default=DEFAULT_DESIRED_COVERAGE)
    parser.add_argument('-p', '--paired', action='store_true')
    parser.add_argument('-s', '--savehash', dest='savehash', default='')
    parser.add_argument('-l', '--loadhash', dest='loadhash',
                        default='')
    parser.add_argument('-R', '--report-to-file', dest='report_file',
                        type=argparse.FileType('w'))
    parser.add_argument('-f', '--force-processing', dest='force',
                        help='continue on next file if read errors are \
                         encountered', action='store_true')
    parser.add_argument('-d', '--dump-frequency', dest='dump_frequency',
                        type=int, help='dump hashtable every d files',
                        default=-1)
    parser.add_argument('input_filenames', nargs='+')

    args = parser.parse_args()

    if not args.quiet:
        if args.min_hashsize == DEFAULT_MIN_HASHSIZE:
            print >>sys.stderr, \
                "** WARNING: hashsize is default!  " \
                "You absodefly want to increase this!\n** " \
                "Please read the docs!"

        print >>sys.stderr, '\nPARAMETERS:'
        print >>sys.stderr, \
            ' - kmer size =    {:d} \t\t(-k)'.format(args.ksize)
        print >>sys.stderr, \
            ' - n hashes =     {:d} \t\t(-N)'.format(args.n_hashes)
        print >>sys.stderr, \
            ' - min hashsize = {:-5.2g} \t(-x)'.format(args.min_hashsize)
        print >>sys.stderr, ' - paired = {} \t\t(-p)'.format(args.paired)
        print >>sys.stderr, ''
        print >>sys.stderr, \
            'Estimated memory usage is {:.2g} bytes \
            (n_hashes x min_hashsize)'.format(args.n_hashes*args.min_hashsize)
        print >>sys.stderr, '-' * 8

    K = args.ksize
    HT_SIZE = args.min_hashsize
    N_HT = args.n_hashes

    report_fp = args.report_file
    filenames = args.input_filenames
    n_threads = max(int(args.n_threads), 2) # min one reader, one writer.

    config = khmer.get_config()
    config.set_reads_input_buffer_size(n_threads * 64 * 1024)

    force=args.force
    dump_frequency = args.dump_frequency
    
    # list to save error files along with throwing exceptions
    if force == True:
        corrupt_files = []
    
    if args.loadhash:
        print 'loading hashtable from', args.loadhash
        ht = khmer.load_counting_hash(args.loadhash)
    else:
        print 'making hashtable'
        ht = khmer.new_counting_hash(K, HT_SIZE, N_HT)

    total = 0
    discarded = 0
    for n, input_filename in enumerate(filenames):
        output_name = os.path.basename(input_filename) + '.keep'
        outfp = open(output_name, 'w')
        
        total_acc = 0
        discarded_acc = 0
        
        try:
            total_acc, discarded_acc = normalize_by_median(input_filename, 
                                                           outfp, ht, args,
                                                           n_threads,
                                                           report_fp)
        except IOError as e:
            handle_error(e, output_name, input_filename, ht)
            if not force:
                print >>sys.stderr, '** Exiting!'
                sys.exit(-1)
            else:
                print >>sys.stderr, '*** Skipping error file, moving on...'
                corrupt_files.append(input_filename)
        else:
            if total_acc == 0 and discarded_acc == 0:
                print 'SKIPPED empty file', input_filename
            else:
                total += total_acc
                discarded += discarded_acc
                print 'DONE with {inp}; kept {kept} of {total} or {perc:2}%'\
                    .format(inp=input_filename,
                    kept=total-discarded, total=total,
                    perc=int(100. - discarded / float(total) * 100.))
                print 'output in', output_name

        if dump_frequency > 0 and n > 0 and n % dump_frequency == 0:
            print 'Backup: Saving hashfile through', input_filename
            if args.savehash:
                hashname = args.savehash
                print '...saving to', hashname
            else:
                hashname = 'backup.ht'
                print 'Nothing given for savehash, saving to', hashname
            ht.save(hashname)
            
    if args.savehash:
        print 'Saving hashfile through', input_filename
        print '...saving to', args.savehash
        ht.save(args.savehash)

    # Change 0.2 only if you really grok it.  HINT: You don't.
    fp_rate = khmer.calc_expected_collisions(ht)
    print 'fp rate estimated to be {:1.3f}'.format(fp_rate)
    
    if force and len(corrupt_files) > 0:
        print >>sys.stderr, "** WARNING: Finished with errors!"
        print >>sys.stderr, "** IOErrors occurred in the following files:"
        print >>sys.stderr, "\t", " ".join(corrupt_files)
    
    if fp_rate > 0.20:
        print >>sys.stderr, "**"
        print >>sys.stderr, "** ERROR: the counting hash is too small for"
        print >>sys.stderr, "** this data set.  Increase hashsize/num ht."
        print >>sys.stderr, "**"
        print >>sys.stderr, "** Do not use these results!!"
        sys.exit(-1)

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
