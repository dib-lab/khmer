#! /usr/bin/env python
import argparse
import khmer
import screed

# @CTB ideas: load in the high-abund k-mers only, as in streaming, past C=20.
# @CTB ideas: count the number of tags; requires merging in running stuff.
# @CTB ideas: extend to include internal/isolated errors

N_HT = 4
HASHSIZE = 1e8
K=20
C=20
MAX_SEQ_LEN = 65535


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('filenames', nargs='+')
    
    args = parser.parse_args()

    ht = khmer.new_counting_hash(K, HASHSIZE, N_HT)

    positions = [0]*MAX_SEQ_LEN
    lengths = []
    n_consumed = 0
    bp_consumed = 0
    total = 0

    for filename in args.filenames:
        print 'opening', filename
        for n, record in enumerate(screed.open(filename)):
            total += 1
            if n % 10000 == 0:
                print '...', n
                if n > 300000:
                    break
            seq = record.sequence.replace('N', 'A')
            med, _, _ = ht.get_median_count(seq)

            if med < C:
                n_consumed += 1
                bp_consumed += len(seq)
                ht.consume(seq)
            else:
                pos = ht.find_first_low_abund_kmer(seq, 2)
                lengths.append(len(seq))
                if pos != len(seq):
                    positions[pos] += 1

    fp = open('out.hist', 'w')
    for n, i in enumerate(positions):
        print >>fp, n, i
    fp.close()

    print 'total sequences:', total
    print 'n consumed:', n_consumed
    print 'bp consumed:', bp_consumed, bp_consumed / float(C)

if __name__ == '__main__':
    main()
