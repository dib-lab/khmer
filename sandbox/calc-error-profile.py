#! /usr/bin/env python
import argparse
import khmer
import screed

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
                posns = ht.find_low_abund_kmers(seq, 2)
                lengths.append(len(seq))

                for p in posns:
                    positions[p] += 1

    lengths.sort()
    max_length = lengths[-1]

    length_count = [0]*max_length
    for j in range(max_length):
        length_count[j] = sum([ 1 for i in lengths if i >= j ])


    last_zero = len(positions) - 1
    while positions[last_zero] == 0:
        last_zero -= 1

    fp = open('out.hist', 'w')
    for n, i in enumerate(positions[:last_zero + 1]):
        print >>fp, n, i, float(i) / float(length_count[n])
    fp.close()

    print 'total sequences:', total
    print 'n consumed:', n_consumed
    print 'bp consumed:', bp_consumed, bp_consumed / float(C)

if __name__ == '__main__':
    main()
