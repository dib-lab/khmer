#! /usr/bin/env python
import os.path
import khmer
import argparse
import screed

K=20
N=4
X=1e8

C1=20
C2=100
PCNT=50

def output_single(r):
    if hasattr(r, 'accuracy'):
        return "@%s\n%s\n+\n%s\n" % (r.name, r.sequence, r.accuracy)
    else:
        return ">%s\n%s\n" % (r.name, r.sequence)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("readfiles", nargs="+")

    args = parser.parse_args()

    kh = khmer.new_counting_hash(K, X, N)
    elim = khmer.new_hashbits(K, X, N)
    outp_index = 1

    for filename in args.readfiles:
        outputfilename = os.path.basename(filename) + '.bink.%d' % outp_index
        fp = open(outputfilename, 'w')

        n = 0
        n_low = 0
        n_med = 0
        n_high = 0
        n_elim = 0
        recent_elim = [0] * 100
        
        for n, record in enumerate(screed.open(filename)):
            if n % 25000 == 0:
                print '...', n
                
            seq = record.sequence.upper()
            if 'N' in seq:
                seq = seq.replace('N', 'G')

            if K <= len(seq):
                do_elim, _, _ = elim.get_median_count(seq)
                if do_elim:
                    n_elim += 1
                    continue
                
                a, _, _ = kh.get_median_count(seq)

                if a < C1:              # low coverage? keep.
                    n_low += 1
                    kh.consume(seq)
                    do_output = True
                elif a >= C2:           # high coverage? discard.
                    n_high += 1
                    elim.consume(seq)
                    recent_elim.pop(0)
                    recent_elim.append(1)
                    if sum(recent_elim) >= PCNT:
                        print 'PING'
                        recent_elim = [0]*100
                        outp_index += 1
                        outputfilename = os.path.basename(filename) + '.bink.%d' % outp_index
                        fp = open(outputfilename, 'w')
                    do_output = False
                else:                   # medium coverage? keep, sort of :)
                    n_med += 1
                    kh.consume_high_abund_kmers(seq, C1)
                    do_output = True

            if do_output:
                fp.write(output_single(record))

        print n, n_low, n_med, n_high, n_elim

if __name__ == '__main__':
    main()
