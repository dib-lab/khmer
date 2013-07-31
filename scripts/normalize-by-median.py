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
import khmer
from itertools import izip
from khmer.counting_args import build_construct_args, DEFAULT_MIN_HASHSIZE
import argparse

DEFAULT_DESIRED_COVERAGE = 5

# Iterate a collection in arbitrary batches
# from: http://stackoverflow.com/questions/4628290/pairs-from-single-list


def batchwise(t, size):
    it = iter(t)
    return izip(*[it] * size)

# Returns true if the pair of records are properly pairs
def validpair(r0, r1):
    return r0.name[-1] == "1" and \
        r1.name[-1] == "2" and \
        r0.name[0:-1] == r1.name[0:-1]

def normalize_by_median(input_filename, outfp, ht, K, DESIRED_COVERAGE, args):

    # In paired mode we read two records at a time
    batch_size = 1
    if args.paired:
        batch_size = 2

    n = -1
    total = 0
    discarded = 0
    for n, batch in enumerate(batchwise(screed.open(
            input_filename), batch_size)):
        if n > 0 and n % 100000 == 0:
            print '... kept', total - discarded, 'of', total, ', or', \
                int(100. - discarded / float(total) * 100.), '%'
            print '... in file', input_filename

            if report_fp:
                print>>report_fp, total, total - discarded, \
                    1. - (discarded / float(total))
                report_fp.flush()

        total += batch_size

        # If in paired mode, check that the reads are properly interleaved
        if args.paired:
            if not validpair(batch[0], batch[1]):
                raise IOError('Error: Improperly interleaved pairs %s %s' \
                               % (batch[0].name, batch[1].name))

        # Emit the batch of reads if any read passes the filter
        # and all reads are longer than K
        passed_filter = False
        passed_length = True
        for record in batch:
            if len(record.sequence) < K:
                passed_length = False
                continue

            seq = record.sequence.replace('N', 'A')
            med, _, _ = ht.get_median_count(seq)

            if med < DESIRED_COVERAGE:
                ht.consume(seq)
                passed_filter = True

        # Emit records if any passed
        if passed_length and passed_filter:
            for record in batch:
                if hasattr(record, 'accuracy'):
                    outfp.write('@%s\n%s\n+\n%s\n' % (record.name,
                                                      record.sequence,
                                                      record.accuracy))
                else:
                    outfp.write(
                        '>%s\n%s\n' % (record.name, record.sequence))
        else:
            discarded += batch_size
    
    return total, discarded

def main():
    parser = build_construct_args()
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
                        type=int, help='dump hashtable every F files',
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
        print >>sys.stderr, ' - kmer size =    %d \t\t(-k)' % args.ksize
        print >>sys.stderr, ' - n hashes =     %d \t\t(-N)' % args.n_hashes
        print >>sys.stderr, \
            ' - min hashsize = %-5.2g \t(-x)' % args.min_hashsize
        print >>sys.stderr, ' - paired =	      %s \t\t(-p)' % args.paired
        print >>sys.stderr, ''
        print >>sys.stderr, \
            'Estimated memory usage is %.2g bytes (n_hashes x min_hashsize)' \
            % (args.n_hashes * args.min_hashsize)
        print >>sys.stderr, '-' * 8

    K = args.ksize
    HT_SIZE = args.min_hashsize
    N_HT = args.n_hashes
    DESIRED_COVERAGE = args.cutoff
    report_fp = args.report_file
    filenames = args.input_filenames
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

        try:
            t, d = normalize_by_median(input_filename, outfp, ht, K, DESIRED_COVERAGE, args)
            total += t
            discarded += d

            if t == 0 and d == 0:
                print 'SKIPPED empty file', input_filename
            else:
                print \
                    'DONE with', input_filename, '; kept', total - discarded, \
                    'of', total, 'or', \
                    int(100. - discarded / float(total) * 100.), '%'
                print 'output in', output_name

        except IOError as e:
            print >>sys.stderr, '** ERROR:', e
            print >>sys.stderr, '** Failed on {}: '.format(input_filename)
            hashname = os.path.basename(input_filename) + '.ht.failed'
            print >>sys.stderr, '** ...dumping hashtable to {}'.format(hashname)
            
            ht.save(hashname)
            os.remove(output_name)
            
            if not force:
                print >>sys.stderr, '** Exiting!'
                sys.exit(-1)
            else:
                print >>sys.stderr, '*** Skipping error file, moving on...'
                corrupt_files.append(input_filename)
                pass

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
    print 'fp rate estimated to be %1.3f' % fp_rate
    
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
