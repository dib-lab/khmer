#! /usr/bin/env python
"""
XXX

Eliminate reads with minimum k-mer abundance higher than
DESIRED_COVERAGE.  Output sequences will be placed in 'infile.keep'.

% python scripts/normalize-by-min.py [ -C <cutoff> ] <data1> <data2> ...

Use '-h' for parameter help.
"""

import sys
import screed
import os
import khmer
from khmer.counting_args import build_construct_args, DEFAULT_MIN_HASHSIZE

DEFAULT_MINIMUM_COVERAGE = 5


def main():
    parser = build_construct_args()
    parser.add_argument("-t", "--trusted-cutoff", dest="trusted_cutoff", type=int, default=3)
    parser.add_argument("--bits-theta", help="Tuning parameter controlling trade off of speed vs alignment sensitivity", default=1.0, type=float, dest="bits_theta")
    parser.add_argument('-C', '--cutoff', type=int, dest='cutoff',
                        default=DEFAULT_MINIMUM_COVERAGE)
    parser.add_argument('-s', '--savehash', dest='savehash', default='')
    parser.add_argument('-l', '--loadhash', dest='loadhash',
                        default='')
    parser.add_argument('input_filenames', nargs='+')

    args = parser.parse_args()

    if not args.quiet:
        if args.min_hashsize == DEFAULT_MIN_HASHSIZE:
            print >>sys.stderr, "** WARNING: hashsize is default!  " \
                "You absodefly want to increase this!\n** " \
                "Please read the docs!"

        print >>sys.stderr, '\nPARAMETERS:'
        print >>sys.stderr, ' - kmer size =    %d \t\t(-k)' % args.ksize
        print >>sys.stderr, ' - n hashes =     %d \t\t(-N)' % args.n_hashes
        print >>sys.stderr, ' - min hashsize = %-5.2g \t(-x)' % \
            args.min_hashsize
        print >>sys.stderr, ''
        print >>sys.stderr, 'Estimated memory usage is %.2g bytes ' \
            '(n_hashes x min_hashsize)' % (
            args.n_hashes * args.min_hashsize)
        print >>sys.stderr, '-' * 8

    K = args.ksize
    HT_SIZE = args.min_hashsize
    N_HT = args.n_hashes
    DESIRED_COVERAGE = args.cutoff

    filenames = args.input_filenames

    if args.loadhash:
        print 'loading hashtable from', args.loadhash
        ht = khmer.load_counting_hash(args.loadhash)
    else:
        print 'making hashtable'
        ht = khmer.new_counting_hash(K, HT_SIZE, N_HT)

    aligner = khmer.new_readaligner(ht, args.trusted_cutoff, args.bits_theta)
            
    total = 0
    discarded = 0
    for input_filename in filenames:
        output_name = os.path.basename(input_filename) + '.keepalign'
        outfp = open(output_name, 'w')

        for n, record in enumerate(screed.open(input_filename)):
            if n > 0 and n % 10000 == 0:
                print '... kept', total - discarded, 'of', total, ', or', \
                    int(100. - discarded / float(total) * 100.), '%'
                print '... in file', input_filename

            total += 1

            if len(record.sequence) < K:
                continue

            seq = record.sequence.replace('N', 'A')

            ##
            score, graph_alignment, read_alignment, truncated = aligner.align(record.sequence)
            print >>sys.stderr, "{0}\t{1:0.2f}\t{2}".format(record.name, score, truncated);

            keep = False
            if truncated:
                keep = True
            else:
                graph_seq = graph_alignment.replace('-', '')
                mincount = ht.get_min_count(graph_seq)

                if mincount < DESIRED_COVERAGE:
                    keep = True
                    seq = graph_seq
                else:
                    assert not keep
                    
            if keep:
                ht.consume(seq)
                outfp.write('>%s\n%s\n' % (record.name, seq))
            else:
                discarded += 1

        print 'DONE with', input_filename, '; kept', total - discarded, 'of',\
            total, 'or', int(100. - discarded / float(total) * 100.), '%'
        print 'output in', output_name

    if args.savehash:
        print 'Saving hashfile through', input_filename
        print '...saving to', args.savehash
        ht.save(args.savehash)

    # Change 0.2 only if you really grok it.  HINT: You don't.
    fp_rate = khmer.calc_expected_collisions(ht)
    print 'fp rate estimated to be %1.3f' % fp_rate

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
