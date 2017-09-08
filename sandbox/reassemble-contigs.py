#! /usr/bin/env python
import argparse
import screed
import khmer

K = 31


def main():
    p = argparse.ArgumentParser()
    p.add_argument('contig_files', nargs='+')
    args = p.parse_args()

    ng = khmer.Nodegraph(K, 1e8, 4)
    starts = []

    for filename in args.contig_files:
        for n, record in enumerate(screed.open(filename)):
            if n and n % 10000 == 0:
                print('...', n)
            ng.consume(record.sequence)
            starts.append(record.sequence[:K])

    hdn = khmer.HashSet(K)
    for filename in args.contig_files:
        for n, record in enumerate(screed.open(filename)):
            if n and n % 10000 == 0:
                print('...', n)
            hdn += ng.find_high_degree_nodes(record.sequence)

    lh = khmer._GraphLabels(ng)
    for filename in args.contig_files:
        for n, record in enumerate(screed.open(filename)):
            if n and n % 10000 == 0:
                print('...', n)
                lh.label_across_high_degree_nodes(record.sequence, hdn, n)


    counter = 0
    for k in starts:
        contigs = lh.assemble_labeled_path(k)
        if not contigs:
            print('nada...')
        for c in contigs:
            print('>%d\n%s' % (counter, c))
            counter += 1


if __name__ == '__main__':
    main()
