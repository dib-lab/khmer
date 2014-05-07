#! /usr/bin/env python2
import sys
import khmer
import argparse
import screed

K = 20
N = 4
X = 1e8


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("readfiles", nargs="+")
    args = parser.parse_args()
    outfp = open('out', 'w')

    kh = khmer.new_counting_hash(K, X, N)

    for filename in args.readfiles:
        for n, record in enumerate(screed.open(filename)):
            if n % 1000 == 0:
                print '...', n
                p = kh.do_subset_partition_with_abundance(10, 255)

                if p.partition_sizes()[0]:
                    p_sizes = dict(p.partition_sizes()[0])
                    p_cov = dict(p.partition_average_coverages(kh))

                    p_sizes_items = p_sizes.items()
                    p_sizes_items.sort(key=lambda x: x[1], reverse=True)

                    for k, size in p_sizes_items[:5]:
                        print '\tp %d: size %d, avg cov %d' % (k, size, p_cov[k])

                    largest_p = p_sizes_items[0][0]
                    if p_cov[largest_p] > 100:
                        pass
# break

            seq = record.sequence.upper()

            a, _, _ = kh.get_median_count(seq)

            if a < 20:              # low coverage? keep.
                kh.consume_and_tag(seq)
            elif a > 100:
                pass
            else:                   # medium coverage? keep, sort of :)
                kh.consume_high_abund_kmers(seq, 20)

            outfp.write('>%s\n%s\n' % (record.name, record.sequence))

if __name__ == '__main__':
    main()
